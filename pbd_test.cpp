#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>

// 简单的3D向量类（如果没有glm库）
struct vec3 {
    float x, y, z;
    
    vec3(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}
    
    vec3 operator+(const vec3& other) const { return vec3(x + other.x, y + other.y, z + other.z); }
    vec3 operator-(const vec3& other) const { return vec3(x - other.x, y - other.y, z - other.z); }
    vec3 operator*(float scalar) const { return vec3(x * scalar, y * scalar, z * scalar); }
    vec3& operator+=(const vec3& other) { x += other.x; y += other.y; z += other.z; return *this; }
    vec3& operator*=(float scalar) { x *= scalar; y *= scalar; z *= scalar; return *this; }
    
    float length() const { return std::sqrt(x*x + y*y + z*z); }
    vec3 normalized() const { 
        float len = length(); 
        return len > 0 ? vec3(x/len, y/len, z/len) : vec3(0,0,0); 
    }
};

vec3 mix(const vec3& a, const vec3& b, float t) {
    return a * (1.0f - t) + b * t;
}

struct Particle {
    vec3 position;
    vec3 predicted_position;
    vec3 velocity;
    float mass;
    float inv_mass;
    bool is_static;
    
    Particle(vec3 pos, float m) : 
        position(pos), predicted_position(pos), velocity(0,0,0), 
        mass(m), inv_mass(m > 0 ? 1.0f/m : 0.0f), is_static(false) {}
};

// 基础距离约束（作为应力约束的基础）
struct DistanceConstraint {
    int particle1, particle2;
    float rest_length;
    float stiffness;
    
    DistanceConstraint(int p1, int p2, float length, float k = 1.0f) :
        particle1(p1), particle2(p2), rest_length(length), stiffness(k) {}
};

// 应力约束 - 考虑材料的弹性和塑性
struct StressConstraint {
    int particle1, particle2;
    float rest_length;
    float young_modulus;        // 杨氏模量
    float yield_strength;       // 屈服强度
    float ultimate_strength;    // 极限强度
    float cross_section_area;   // 横截面积
    float current_stress;       // 当前应力
    bool is_broken;            // 是否断裂
    
    StressConstraint(int p1, int p2, float length, float E, float yield, float ultimate, float area) :
        particle1(p1), particle2(p2), rest_length(length), 
        young_modulus(E), yield_strength(yield), ultimate_strength(ultimate),
        cross_section_area(area), current_stress(0.0f), is_broken(false) {}
};

class PBDStressSystem {
private:
    std::vector<Particle> particles;
    std::vector<StressConstraint> stress_constraints;
    float damping;
    int solver_iterations;
    
public:
    PBDStressSystem(float damp = 0.99f, int iterations = 5) :
        damping(damp), solver_iterations(iterations) {}
    
    // 添加粒子
    int addParticle(vec3 position, float mass) {
        particles.emplace_back(position, mass);
        return particles.size() - 1;
    }
    
    // 创建粒子链
    void createParticleChain(vec3 start, vec3 end, int num_particles, 
                           float mass, float E, float yield, float ultimate, float area) {
        std::vector<int> chain_particles;
        
        // 创建链上的粒子
        for (int i = 0; i < num_particles; i++) {
            float t = static_cast<float>(i) / (num_particles - 1);
            vec3 pos = mix(start, end, t);
            int particle_id = addParticle(pos, mass);
            chain_particles.push_back(particle_id);
            
            // 第一个粒子固定
            if (i == 0) {
                particles[particle_id].is_static = true;
                particles[particle_id].inv_mass = 0.0f;
            }
        }
        
        // 创建相邻粒子间的应力约束
        for (int i = 0; i < chain_particles.size() - 1; i++) {
            int p1 = chain_particles[i];
            int p2 = chain_particles[i + 1];
            float length = (particles[p2].position - particles[p1].position).length();
            
            stress_constraints.emplace_back(p1, p2, length, E, yield, ultimate, area);
        }
    }
    
    // 计算应力约束
    void solveStressConstraint(StressConstraint& constraint) {
        if (constraint.is_broken) return;
        
        Particle& p1 = particles[constraint.particle1];
        Particle& p2 = particles[constraint.particle2];
        
        vec3 delta = p2.predicted_position - p1.predicted_position;
        float current_length = delta.length();
        
        if (current_length < 1e-6f) return;
        
        // 计算应变
        float strain = (current_length - constraint.rest_length) / constraint.rest_length;
        
        // 计算应力 (σ = E * ε 在弹性范围内)
        float stress = constraint.young_modulus * strain;
        constraint.current_stress = std::abs(stress);
        
        // 检查是否超过极限强度（断裂）
        if (constraint.current_stress > constraint.ultimate_strength) {
            constraint.is_broken = true;
            return;
        }
        
        // 计算修正向量
        vec3 correction_direction = delta * (1.0f / current_length);
        float correction_magnitude;
        
        // 弹性范围内的线性修正
        if (constraint.current_stress <= constraint.yield_strength) {
            correction_magnitude = constraint.rest_length - current_length;
        }
        // 塑性范围内的非线性修正
        else {
            float plastic_factor = 1.0f - (constraint.current_stress - constraint.yield_strength) / 
                                  (constraint.ultimate_strength - constraint.yield_strength);
            correction_magnitude = (constraint.rest_length - current_length) * plastic_factor;
        }
        
        vec3 correction = correction_direction * (correction_magnitude * 0.5f);
        
        // 应用位置修正
        float total_inv_mass = p1.inv_mass + p2.inv_mass;
        if (total_inv_mass > 0) {
            vec3 correction1 = correction * (-p1.inv_mass / total_inv_mass);
            vec3 correction2 = correction * (p2.inv_mass / total_inv_mass);
            
            if (!p1.is_static) p1.predicted_position += correction1;
            if (!p2.is_static) p2.predicted_position += correction2;
        }
    }
    
    // 时间步进
    void timeStep(float dt, vec3 gravity = vec3(0, -9.81f, 0)) {
        // 1. 预测位置
        for (auto& particle : particles) {
            if (!particle.is_static) {
                particle.velocity += gravity * dt;
                particle.predicted_position = particle.position + particle.velocity * dt;
            }
        }
        
        // 2. 求解约束
        for (int iter = 0; iter < solver_iterations; iter++) {
            for (auto& constraint : stress_constraints) {
                solveStressConstraint(constraint);
            }
        }
        
        // 3. 更新位置和速度
        for (auto& particle : particles) {
            if (!particle.is_static) {
                particle.velocity = (particle.predicted_position - particle.position) * (1.0f / dt);
                particle.velocity *= damping;
                particle.position = particle.predicted_position;
            }
        }
    }
    
    // 获取系统状态信息
    void getSystemInfo() {
        int active_constraints = 0;
        float max_stress = 0.0f;
        
        for (const auto& constraint : stress_constraints) {
            if (!constraint.is_broken) {
                active_constraints++;
                max_stress = std::max(max_stress, constraint.current_stress);
            }
        }
        
        std::cout << "Active constraints: " << active_constraints << "/" << stress_constraints.size() << std::endl;
        std::cout << "Max stress: " << max_stress << std::endl;
    }
    
    // 添加外力
    void applyForce(int particle_id, vec3 force, float dt) {
        if (particle_id < particles.size() && !particles[particle_id].is_static) {
            particles[particle_id].velocity += force * (particles[particle_id].inv_mass * dt);
        }
    }
    
    // 获取粒子位置
    std::vector<vec3> getParticlePositions() const {
        std::vector<vec3> positions;
        for (const auto& particle : particles) {
            positions.push_back(particle.position);
        }
        return positions;
    }
    
    // 打印链的状态
    void printChainStatus() {
        std::cout << "\n=== 粒子链状态 ===" << std::endl;
        for (size_t i = 0; i < particles.size(); i++) {
            const auto& p = particles[i];
            std::cout << "粒子 " << i << ": 位置(" 
                     << std::fixed << std::setprecision(3)
                     << p.position.x << ", " 
                     << p.position.y << ", " 
                     << p.position.z << ")";
            if (p.is_static) std::cout << " [固定]";
            std::cout << std::endl;
        }
        
        std::cout << "\n=== 约束状态 ===" << std::endl;
        for (size_t i = 0; i < stress_constraints.size(); i++) {
            const auto& c = stress_constraints[i];
            vec3 delta = particles[c.particle2].position - particles[c.particle1].position;
            float current_length = delta.length();
            float strain = (current_length - c.rest_length) / c.rest_length * 100;
            
            std::cout << "约束 " << i << " (粒子" << c.particle1 << "-" << c.particle2 << "): ";
            if (c.is_broken) {
                std::cout << "[已断裂]";
            } else {
                std::cout << "应力=" << std::scientific << std::setprecision(2) << c.current_stress 
                         << " Pa, 应变=" << std::fixed << std::setprecision(1) << strain << "%";
            }
            std::cout << std::endl;
        }
    }
    
    // 获取统计信息
    void printStatistics() {
        int broken_constraints = 0;
        float max_stress = 0.0f;
        float total_length = 0.0f;
        
        for (const auto& constraint : stress_constraints) {
            if (constraint.is_broken) {
                broken_constraints++;
            } else {
                max_stress = std::max(max_stress, constraint.current_stress);
            }
            
            vec3 delta = particles[constraint.particle2].position - particles[constraint.particle1].position;
            total_length += delta.length();
        }
        
        std::cout << "\n=== 系统统计 ===" << std::endl;
        std::cout << "总粒子数: " << particles.size() << std::endl;
        std::cout << "总约束数: " << stress_constraints.size() << std::endl;
        std::cout << "断裂约束: " << broken_constraints << std::endl;
        std::cout << "完整约束: " << stress_constraints.size() - broken_constraints << std::endl;
        std::cout << "最大应力: " << std::scientific << std::setprecision(2) << max_stress << " Pa" << std::endl;
        std::cout << "总链长: " << std::fixed << std::setprecision(3) << total_length << " m" << std::endl;
    }
    
    // 检查是否有断裂
    bool hasAnyBrokenConstraints() const {
        for (const auto& constraint : stress_constraints) {
            if (constraint.is_broken) return true;
        }
        return false;
    }
    
    // 获取断裂约束数量
    int getBrokenConstraintCount() const {
        int count = 0;
        for (const auto& constraint : stress_constraints) {
            if (constraint.is_broken) count++;
        }
        return count;
    }
};

// 测试场景1：基础重力下垂测试
void testBasicGravity() {
    std::cout << "\n" << "="*50 << std::endl;
    std::cout << "测试场景1：基础重力下垂测试" << std::endl;
    std::cout << "="*50 << std::endl;
    
    PBDStressSystem system(0.99f, 5);
    
    // 使用较弱的材料参数便于观察变形
    float young_modulus = 1e6f;      // 1 MPa (较软)
    float yield_strength = 1e5f;     // 100 kPa  
    float ultimate_strength = 2e5f;  // 200 kPa
    float cross_section = 1e-4f;     // 1 cm²
    
    // 创建垂直链条
    system.createParticleChain(
        vec3(0, 2, 0),    // 起点(顶部)
        vec3(0, 0, 0),    // 终点(底部) 
        5,                // 5个粒子
        0.1f,             // 质量
        young_modulus, yield_strength, ultimate_strength, cross_section
    );
    
    std::cout << "初始状态:" << std::endl;
    system.printChainStatus();
    
    // 仿真50步
    float dt = 1.0f / 60.0f;
    for (int i = 0; i < 50; i++) {
        system.timeStep(dt);
    }
    
    std::cout << "\n重力作用50步后:" << std::endl;
    system.printChainStatus();
    system.printStatistics();
}

// 测试场景2：外力拉伸测试
void testExternalForce() {
    std::cout << "\n" << "="*50 << std::endl;
    std::cout << "测试场景2：外力拉伸测试" << std::endl;
    std::cout << "="*50 << std::endl;
    
    PBDStressSystem system(0.95f, 8);
    
    // 中等强度材料
    float young_modulus = 10e6f;     // 10 MPa
    float yield_strength = 5e5f;     // 500 kPa
    float ultimate_strength = 1e6f;  // 1 MPa
    float cross_section = 1e-4f;
    
    // 创建水平链条
    system.createParticleChain(
        vec3(0, 1, 0),    // 起点(左)
        vec3(2, 1, 0),    // 终点(右)
        6,                // 6个粒子
        0.2f,             // 质量
        young_modulus, yield_strength, ultimate_strength, cross_section
    );
    
    std::cout << "开始拉伸测试..." << std::endl;
    
    float dt = 1.0f / 60.0f;
    
    // 阶段1：轻微拉伸
    std::cout << "\n阶段1：轻微拉伸 (前30步)" << std::endl;
    for (int i = 0; i < 30; i++) {
        system.applyForce(5, vec3(2.0f, 0, 0), dt);  // 向右拉最后一个粒子
        system.timeStep(dt);
    }
    system.printStatistics();
    
    // 阶段2：强力拉伸
    std::cout << "\n阶段2：强力拉伸 (30-60步)" << std::endl;
    for (int i = 30; i < 60; i++) {
        system.applyForce(5, vec3(10.0f, 0, 0), dt);  // 更大的拉力
        system.timeStep(dt);
    }
    system.printStatistics();
    
    // 阶段3：极限拉伸
    std::cout << "\n阶段3：极限拉伸 (60-100步)" << std::endl;
    for (int i = 60; i < 100; i++) {
        system.applyForce(5, vec3(50.0f, 0, 0), dt);  // 极大拉力
        system.timeStep(dt);
    }
    
    std::cout << "\n最终状态:" << std::endl;
    system.printChainStatus();
    system.printStatistics();
}

// 测试场景3：材料断裂测试
void testMaterialFailure() {
    std::cout << "\n" << "="*50 << std::endl;
    std::cout << "测试场景3：材料断裂测试" << std::endl;
    std::cout << "="*50 << std::endl;
    
    PBDStressSystem system(0.9f, 10);
    
    // 脆性材料（容易断裂）
    float young_modulus = 50e6f;     // 50 MPa (较硬)
    float yield_strength = 1e5f;     // 100 kPa
    float ultimate_strength = 1.5e5f;// 150 kPa (低断裂强度)
    float cross_section = 5e-5f;     // 较小截面积
    
    // 创建长链条
    system.createParticleChain(
        vec3(0, 3, 0),    // 起点
        vec3(0, 0, 0),    // 终点
        8,                // 8个粒子
        0.3f,             // 较大质量
        young_modulus, yield_strength, ultimate_strength, cross_section
    );
    
    std::cout << "测试断裂过程..." << std::endl;
    
    float dt = 1.0f / 60.0f;
    
    // 逐步增加外力直到断裂
    for (int step = 0; step < 200; step++) {
        // 在链末端施加向下的力
        float force_magnitude = 1.0f + step * 0.5f;  // 逐渐增大的力
        system.applyForce(7, vec3(0, -force_magnitude, 0), dt);
        
        system.timeStep(dt);
        
        // 每20步输出一次状态
        if (step % 20 == 0 || step < 10) {
            std::cout << "\n步骤 " << step << " (力=" << std::fixed << std::setprecision(1) 
                     << force_magnitude << "N):" << std::endl;
            system.printStatistics();
        }
        
        // 检查是否有断裂发生
        if (system.hasAnyBrokenConstraints() && step > 50) {  // 忽略前50步的可能误报
            std::cout << "\n!!! 检测到断裂，停止测试 !!!" << std::endl;
            break;
        }
    }
    
    std::cout << "\n最终链条状态:" << std::endl;
    system.printChainStatus();
}

// 主函数
int main() {
    std::cout << "PBD应力约束系统测试程序" << std::endl;
    std::cout << "========================" << std::endl;
    
    try {
        // 运行所有测试
        testBasicGravity();
        
        std::cout << "\n按回车继续下一个测试...";
        std::cin.get();
        
        testExternalForce();
        
        std::cout << "\n按回车继续下一个测试...";
        std::cin.get();
        
        testMaterialFailure();
        
        std::cout << "\n" << "="*50 << std::endl;
        std::cout << "所有测试完成！" << std::endl;
        std::cout << "="*50 << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "\n按回车退出...";
    std::cin.get();
    return 0;
}