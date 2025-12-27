#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
依赖: matplotlib + numpy
用法: python plot_ising.py
会在当前目录生成 m_vs_t.pdf 和 m1_beta_star.pdf
"""
import numpy as np
import matplotlib.pyplot as plt

# 1. 读取 |M|/N 随 T 数据
t, m = np.loadtxt('m_vs_t.dat', unpack=True)

plt.figure(figsize=(4,3))
plt.plot(t, m, 'o-', ms=4, lw=1)
plt.xlabel(r'$T$')
plt.ylabel(r'$|M|/N$')
plt.tight_layout()
plt.savefig('m_vs_t.pdf')
plt.close()

# 2. 读取 M^(1/β*) 多列数据
data = np.loadtxt('m1_beta_star.dat')
temps = data[:, 0]          # 第一列是温度
cols  = data[:, 1:]         # 其余每列对应一个 β*

# 取 β* 值作为图例标签
with open('m1_beta_star.dat') as f:
    header = f.readline()
beta_stars = [float(tok.split('=')[1]) for tok in header.split()[2:]]

plt.figure(figsize=(5,4))
for y, b in zip(cols.T, beta_stars):
    plt.plot(temps, y, label=f'β*={b:.3f}', lw=1)
plt.xlabel(r'$T$')
plt.ylabel(r'$M^{1/\beta^*}$')
plt.legend(fontsize=8, ncol=2)
plt.tight_layout()
plt.savefig('m1_beta_star.pdf')
plt.close()

print('绘图完成：m_vs_t.pdf  和  m1_beta_star.pdf')