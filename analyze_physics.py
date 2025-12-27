#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os

print("物理意义分析...")

# 理论临界温度
THEORETICAL_TC = {
    1: 0.0,      # 1D无相变
    2: 2.269,    # 2D Onsager精确解
    3: 4.511,    # 3D数值解
    4: 6.68      # 4D平均场论
}

# 理论β值
THEORETICAL_BETA = {
    1: 0.0,      # 1D无相变
    2: 0.125,    # 2D精确解
    3: 0.326,    # 3D数值解
    4: 0.5       # 4D平均场论
}

plt.figure(figsize=(12, 8))

# 读取并绘制每个维度的数据
colors = ['blue', 'red', 'green', 'purple']
markers = ['o', 's', '^', 'D']

for dim in range(1, 5):
    filename = f'm_vs_t_{dim}d.dat'
    if os.path.exists(filename):
        t, m = np.loadtxt(filename, unpack=True)
        
        # 绘制数据点
        plt.plot(t, m, marker=markers[dim-1], color=colors[dim-1], 
                linestyle='-', label=f'{dim}D Simulated', markersize=3, alpha=0.7)
        
        # 标记理论临界温度
        if dim > 1:  # 1D没有相变
            plt.axvline(x=THEORETICAL_TC[dim], color=colors[dim-1], 
                       linestyle='--', alpha=0.5, 
                       label=f'{dim}D Tc≈{THEORETICAL_TC[dim]:.3f}')

# 美化图表
plt.xlabel('Temperature T', fontsize=12)
plt.ylabel('|M|/N', fontsize=12)
plt.title('Ising Model: Critical Temperature Increases with Dimension\n' + 
          'Physical Explanation: Higher dimensions = more neighbors = stronger cooperation',
          fontsize=14)
plt.legend(fontsize=10, loc='best')
plt.grid(True, alpha=0.3)

# 添加理论解释文本
explanation = (
    "Critical Temperature (Tc) increases with dimension:\n"
    "• 1D: Tc = 0 (no phase transition)\n"
    "• 2D: Tc ≈ 2.269 (Onsager exact solution)\n"
    "• 3D: Tc ≈ 4.511 (numerical calculation)\n"
    "• 4D: Tc ≈ 6.68 (mean field theory)\n\n"
    "Why? Higher dimensions → more neighbors → stronger cooperation"
)
plt.text(0.02, 0.98, explanation, transform=plt.gca().transAxes,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
         fontsize=9)

plt.tight_layout()
plt.savefig('physics_analysis.pdf', dpi=300, bbox_inches='tight')
plt.close()

print("物理分析图表已保存: physics_analysis.pdf")

# 创建临界温度对比图
plt.figure(figsize=(10, 6))
dimensions = list(range(1, 5))
tc_values = [THEORETICAL_TC[d] for d in dimensions]
beta_values = [THEORETICAL_BETA[d] for d in dimensions]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# 临界温度图
ax1.bar(dimensions, tc_values, color='lightblue', edgecolor='black')
ax1.set_xlabel('Dimension')
ax1.set_ylabel('Critical Temperature Tc')
ax1.set_title('Critical Temperature vs Dimension')
ax1.set_xticks(dimensions)
ax1.grid(True, alpha=0.3, axis='y')

for i, v in enumerate(tc_values):
    if v > 0:  # 1D的Tc=0不显示
        ax1.text(i+1, v + 0.1, f'{v:.3f}', ha='center', va='bottom')

# 临界指数图
ax2.bar(dimensions, beta_values, color='lightcoral', edgecolor='black')
ax2.set_xlabel('Dimension')
ax2.set_ylabel('Critical Exponent β')
ax2.set_title('Critical Exponent β vs Dimension')
ax2.set_xticks(dimensions)
ax2.grid(True, alpha=0.3, axis='y')

for i, v in enumerate(beta_values):
    ax2.text(i+1, v + 0.02, f'{v:.3f}', ha='center', va='bottom')

plt.tight_layout()
plt.savefig('critical_parameters.pdf', dpi=300, bbox_inches='tight')
plt.close()

print("临界参数对比图已保存: critical_parameters.pdf")

print("\n=== 物理解释 ===")
print("1. 临界温度随维度增加是完全正确的！")
print("2. 物理原因：")
print("   - 高维度有更多邻居（2D:4个, 3D:6个, 4D:8个）")
print("   - 更多邻居=更强的协同效应")
print("   - 需要更高温度来破坏长程有序")
print("3. 临界指数β也随维度增加，趋近于平均场论的0.5")
print("4. 这是统计物理中的基本规律！")