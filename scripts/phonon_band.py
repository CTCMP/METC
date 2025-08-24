# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: band.py
"""
import numpy as np
from matplotlib import pyplot as plt
import os, sys
import os.path as osp
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(osp.join(curPath, '../'))
import configfile.inputpara as parameters
sys.path.append(os.getcwd())
import AutoinputQE as inpp
# import matplotlib as mpl
from matplotlib import font_manager
fontpath = osp.join(curPath, "../font/arial.ttf")
# fontpath = "/work/wangr/data/dxy/scripts/program/QE/font/arial.ttf"
# ------------------ font setup ----------------------#
font_properties = font_manager.FontProperties(fname=fontpath)
styles = ['normal', 'italic', 'oblique']
weights = ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
font = {'fontproperties': font_properties,  # 'family': 'Liberation Sans'
		'style': styles[0],
		'color': 'black',
		'weight': weights[1],
		'size': 12, 
        }
sys.dont_write_bytecode = True       ## not generate __pycache__ files

#------------------- Data manipulating ----------------------
#-------------------------  QE  -----------------------------
def band_data_pro(filepath, filename):
    if osp.exists(osp.join(filepath, filename)): 
        file1 = open(osp.join(filepath, filename), 'r').readlines()
        file2 = open(osp.join(filepath, filename + ".gp"), 'r').readlines()
    else:
        file1 = open(osp.join(filepath, parameters.PHONON_mat['flfrq']), 'r').readlines()
        file2 = open(osp.join(filepath, parameters.PHONON_mat['flfrq'] + ".gp"), 'r').readlines()
    nbands = int(file1[0].split(',')[0].split('=')[1])
    nks = int(file1[0].split(',')[1].split("=")[1].split("/")[0])
    band_data = np.loadtxt(file2, skiprows=0, dtype=float)
    band_k = band_data[:nks, 0]
    kk = band_data[:, 0]
    energy = band_data[:, 1:len(band_data[0, :])]/8.065
    hsp_labels = [str(items.split("!")[1]).strip() for items in inpp.kpoints_band.split("\n")[1:-1]]
    for i in range(len(hsp_labels)):
        # if str(hsp_labels[i]).lower() == "\Gamma".lower():
        #   hsp_labels[i] = u"Γ"
        if '\\' in hsp_labels[i].lower():
            hsp_labels[i] = "$" + hsp_labels[i] + "$"
    hsp_nums = len(hsp_labels)-1
    hsp_data = band_k[0:nks:(int(nks/(hsp_nums)))]
    hsp = {
        'klabels' : hsp_labels,
        'hspdata' : hsp_data,
    }
    # print(hsp['klabels'])
    return nbands, nks, hsp, band_k, kk, energy 

def plot_band(filepath, filename):
    nbands, nks, hsp, band_k, kk, energy  = band_data_pro(filepath, filename)
    fig = plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)
    for kk in range(0, len(energy[0, :])):
        ax.plot(band_k, energy[:, kk], color='#B8860B', lw=2.0)
    max_value = np.max(energy)
    for i in range(0,len(hsp['klabels'])):
        ax.vlines(hsp['hspdata'][i], ymin=0, ymax=max_value*1.5,linewidth=0.8,linestyle="--", color='gray')
    ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
    ax.set_xticks(hsp['hspdata'])
    ax.set_xticklabels(hsp['klabels'])
    font_size = font['size']
    ax.set_ylabel("Frequency (meV)", fontdict=font)
    ax.set_ylim(-0.5, max_value*1.03)
    ax.set_xlim(0, np.max(band_k))
    labels = ax.get_xticklabels() + ax.get_yticklabels() 
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(font_size) for label in labels]
    plt.savefig("phonon.png", dpi=400, bbox_inches = 'tight', pad_inches=0.1)

if __name__ == '__main__':
    # energy = [0, 3]
    plot_band(os.getcwd(), "pwscf.freq")

