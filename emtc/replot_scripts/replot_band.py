# -*- coding: utf-8 -*-
"""
@author: XY Ding
mail to: dxy_vasp@163.com
python3: replot_band.py
"""

import numpy as np 
from matplotlib import pyplot as plt
import os 

def band_plot(fname, colorb):
    ax = plt.subplot(111)
    data = np.loadtxt(fname, skiprows=1, dtype=float)
    for b in range(1, len(data[0, :])):
        ax.plot(data[:, 0], data[:, b], color=colorb)
    return ax

cls = ['blue', 'red']
ffs = os.listdir()
bd_f = []
if 'BAND_OUTPUT.dat' in ffs:
    bd_f.append('BAND_OUTPUT.dat')
elif 'BAND_OUTPUT_UP.dat' in ffs:
    bd_f.append('BAND_OUTPUT_UP.dat')
    bd_f.append('BAND_OUTPUT_DW.dat')

klabels = open('KLINES.dat', 'r').readlines()[0].split(':')[1].split('-')
klabel_ticks =[float(item) for item in open('KLINES.dat', 'r').readlines()[1].split(':')[1].split()]
print(klabels)
print(klabel_ticks)
kdata = np.loadtxt('KLINES.dat', skiprows=2, dtype=float)
fig = plt.figure(figsize=(8, 6))
for i in range(len(bd_f)):
    ax = band_plot(bd_f[i], cls[i])
ax.plot(kdata[:, 0], kdata[:, 1], linewidth=1.0, linestyle="--", color='gray')
ax.axhline(0, xmin=0, xmax=np.max(kdata[:, 0]), linewidth=1.0, linestyle="--", color='gray')
ax.set_ylim([-10, 5])
ax.set_xlim([0, np.max(kdata[:, 0])])
ax.set_xticks(klabel_ticks)
ax.set_xticklabels(klabels)

if len(bd_f) == 2:
    ax.plot([],[], c=cls[0], label='SPIN UP')
    ax.plot([],[], c=cls[0], label='SPIN DW')
leg = ax.legend()
plt.ylabel("Energy(eV)", fontsize=20, fontname='Arial')
labels = ax.get_xticklabels() + ax.get_yticklabels() 
[label.set_fontname('Arial') for label in labels]
[label.set_fontsize(20) for label in labels]
labelsl = leg.get_texts()
[label.set_fontname('Arial') for label in labelsl]
[label.set_fontsize(14) for label in labelsl]
plt.savefig('band.png', dpi=300, bbox_inches = 'tight', pad_inches=0.2)
plt.close()