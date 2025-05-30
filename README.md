# pbd_stress_constraint_test

PBD算法的stress应力约束测试

约束求解公式

> 应变 ε = (当前长度 - 原长) / 原长
> 应力 σ = E × ε (弹性范围内)
> 根据应力状态决定修正力度

优化测试设置

## 链式连接的实现

**创建粒子链：**

- `createParticleChain()` 在起点和终点间创建指定数量的粒子
- 相邻粒子间自动建立应力约束
- 第一个粒子设为静态（固定点）

**约束传递：**

- 每个约束连接两个相邻粒子
- 力通过约束在整条链上传递
- 支持链的任意位置断裂

## 使用建议

**材料参数调整：**

- 调大`young_modulus`使链更硬
- 调小`yield_strength`使更容易塑性变形
- 调大`ultimate_strength`使更难断裂

**性能优化：**

- 增加`solver_iterations`提高稳定性但降低性能
- 适当的`damping`防止震荡
- 合理设置时间步长`dt`

## 编译和运行说明

**编译命令：**

```bash
g++ -std=c++11 -O2 -o pbd_test pbd_stress_constraint.cpp
```

或者使用更详细的编译选项：

```bash
g++ -std=c++11 -Wall -Wextra -O2 -o pbd_test pbd_stress_constraint.cpp
```

**运行程序：**

```bash
./pbd_test
```