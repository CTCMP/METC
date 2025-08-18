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
bd_f = 'PHON_BAND.dat'

klabels = open('QLINES.dat', 'r').readlines()[0].split(':')[1].split('-')
klabel_ticks =[float(item) for item in open('QLINES.dat', 'r').readlines()[1].split(':')[1].split()]
print(klabels)
print(klabel_ticks)
kdata = np.loadtxt('QLINES.dat', skiprows=2, dtype=float)
fig = plt.figure(figsize=(8, 6))
ax = band_plot(bd_f, cls[0])
ax.plot(kdata[:, 0], kdata[:, 1], linewidth=1.0, linestyle="--", color='gray')
ax.set_ylim([-5, np.max(kdata[:, 1])])
ax.set_xlim([0, np.max(kdata[:, 0])])
ax.set_xticks(klabel_ticks)
ax.set_xticklabels(klabels)

plt.ylabel("Frequency (meV)", fontsize=20, fontname='Arial')
labels = ax.get_xticklabels() + ax.get_yticklabels() 
[label.set_fontname('Arial') for label in labels]
[label.set_fontsize(20) for label in labels]
plt.savefig('phonon.png', dpi=300, bbox_inches = 'tight', pad_inches=0.1)
plt.close()