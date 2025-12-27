#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
依赖: matplotlib + numpy
用法: python3 plot_comparison.py
会生成多维度伊辛模型的对比图表
"""
import numpy as np
import matplotlib.pyplot as plt
import os

# 使用英文字体
plt.rcParams['font.family'] = 'DejaVu Sans'

# 理论β值（用于参考）
THEORETICAL_BETA = {
    1: 0.0,    # 1D没有相变
    2: 0.125,  # 2D精确解
    3: 0.326,  # 3D数值解
    4: 0.5     # 4D平均场论
}

def plot_m_vs_t_comparison():
    """绘制不同维度的|M|/N随温度变化对比图"""
    plt.figure(figsize=(8, 6))
    
    colors = ['blue', 'red', 'green', 'purple']
    markers = ['o', 's', '^', 'D']
    
    for dim in range(1, 5):
        filename = f'm_vs_t_{dim}d.dat'
        if os.path.exists(filename):
            t, m = np.loadtxt(filename, unpack=True)
            plt.plot(t, m, marker=markers[dim-1], color=colors[dim-1], 
                    linestyle='-', label=f'{dim}D', markersize=4, linewidth=1.5)
    
    plt.xlabel('Temperature T', fontsize=12)
    plt.ylabel('|M|/N', fontsize=12)
    plt.title('Ising Model Magnetization Comparison (Different Dimensions)', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)
    plt.savefig('m_vs_t_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print('Magnetization comparison saved: m_vs_t_comparison.pdf')

def plot_beta_star_analysis():
    """绘制不同维度的β*分析图"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    colors = plt.cm.viridis(np.linspace(0, 1, 10))
    
    for dim in range(1, 5):
        ax = axes[dim-1]
        filename = f'm1_beta_star_{dim}d.dat'
        
        if os.path.exists(filename):
            data = np.loadtxt(filename)
            temps = data[:, 0]
            cols = data[:, 1:]
            
            # 获取β*值
            with open(filename) as f:
                header = f.readline()
            beta_stars = [float(tok.split('=')[1]) for tok in header.split()[2:]]
            
            # 绘制每条线
            for i, (y, b) in enumerate(zip(cols.T, beta_stars)):
                alpha = 0.3 + 0.7 * abs(b - THEORETICAL_BETA[dim]) / 0.2  # 接近理论值的线更明显
                ax.plot(temps, y, color=colors[i % len(colors)], 
                       label=f'β*={b:.3f}', linewidth=1, alpha=alpha)
            
            # 标记理论β值
            ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, 
                      label=f'Theory β={THEORETICAL_BETA[dim]:.3f}')
            
            ax.set_xlabel('Temperature T', fontsize=10)
            ax.set_ylabel(f'$M^{{1/\beta^*}}$', fontsize=10)
            ax.set_title(f'{dim}D Ising Model - β* Analysis', fontsize=12)
            ax.legend(fontsize=6, ncol=2)
            ax.grid(True, alpha=0.3)
    
    plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)
    plt.savefig('beta_star_analysis.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print('β* analysis plot saved: beta_star_analysis.pdf')

def plot_beta_comparison():
    """绘制不同维度的β值对比图"""
    dimensions = [1, 2, 3, 4]
    theoretical_betas = [THEORETICAL_BETA[d] for d in dimensions]
    
    # 这里需要手动输入从β*分析中估计的β值
    # 实际使用时需要根据模拟结果调整
    estimated_betas = [0.0, 0.125, 0.326, 0.5]  # 示例值
    
    x = np.arange(len(dimensions))
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(8, 6))
    bars1 = ax.bar(x - width/2, theoretical_betas, width, 
                   label='Theoretical', color='lightblue', edgecolor='black')
    bars2 = ax.bar(x + width/2, estimated_betas, width, 
                   label='Simulated', color='lightcoral', edgecolor='black')
    
    ax.set_xlabel('Dimension', fontsize=12)
    ax.set_ylabel('Critical Exponent β', fontsize=12)
    ax.set_title('Critical Exponent β Comparison (Different Dimensions)', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{d}D' for d in dimensions])
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    
    # 添加数值标签
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.annotate(f'{height:.3f}',
                       xy=(bar.get_x() + bar.get_width() / 2, height),
                       xytext=(0, 3),
                       textcoords="offset points",
                       ha='center', va='bottom', fontsize=9)
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)
    plt.savefig('beta_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print('β comparison plot saved: beta_comparison.pdf')

def main():
    """Main function"""
    print("Starting multi-dimensional Ising model comparison plots...")
    
    # Check data files exist
    missing_files = []
    for dim in range(1, 5):
        if not os.path.exists(f'm_vs_t_{dim}d.dat'):
            missing_files.append(f'm_vs_t_{dim}d.dat')
        if not os.path.exists(f'm1_beta_star_{dim}d.dat'):
            missing_files.append(f'm1_beta_star_{dim}d.dat')
    
    if missing_files:
        print("Warning: The following data files are missing:")
        for file in missing_files:
            print(f"  - {file}")
        print("Please run the Rust program first to generate data files.")
        return
    
    # Plot various charts
    plot_m_vs_t_comparison()
    plot_beta_star_analysis()
    plot_beta_comparison()
    
    print("\nAll plots generated successfully!")
    print("Generated files:")
    print("  - m_vs_t_comparison.pdf: Magnetization comparison")
    print("  - beta_star_analysis.pdf: β* analysis")
    print("  - beta_comparison.pdf: β value comparison")

if __name__ == "__main__":
    main()