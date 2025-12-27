#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os

print("Starting simple multi-dimensional Ising model comparison...")

# 理论β值
theoretical_beta = {1: 0.0, 2: 0.125, 3: 0.326, 4: 0.5}

# 1. 磁化强度对比图
plt.figure(figsize=(8, 6))
colors = ['blue', 'red', 'green', 'purple']
markers = ['o', 's', '^', 'D']

for dim in range(1, 5):
    filename = f'm_vs_t_{dim}d.dat'
    if os.path.exists(filename):
        t, m = np.loadtxt(filename, unpack=True)
        plt.plot(t, m, marker=markers[dim-1], color=colors[dim-1], 
                linestyle='-', label=f'{dim}D', markersize=4)

plt.xlabel('Temperature T')
plt.ylabel('|M|/N')
plt.title('Ising Model Magnetization Comparison')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('m_vs_t_comparison.pdf')
plt.close()
print('Magnetization comparison saved: m_vs_t_comparison.pdf')

# 2. β值对比图
dimensions = [1, 2, 3, 4]
theoretical_betas = [theoretical_beta[d] for d in dimensions]

plt.figure(figsize=(8, 6))
x = np.arange(len(dimensions))
plt.bar(x, theoretical_betas, color='lightblue', edgecolor='black')
plt.xlabel('Dimension')
plt.ylabel('Critical Exponent β')
plt.title('Theoretical Critical Exponent β Values')
plt.xticks(x, [f'{d}D' for d in dimensions])
plt.grid(True, alpha=0.3, axis='y')

# 添加数值标签
for i, v in enumerate(theoretical_betas):
    plt.text(i, v + 0.01, f'{v:.3f}', ha='center', va='bottom')

plt.savefig('beta_theoretical.pdf')
plt.close()
print('Theoretical β values saved: beta_theoretical.pdf')

print("\nSimple plots completed!")
print("Generated files:")
print("  - m_vs_t_comparison.pdf")
print("  - beta_theoretical.pdf")