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
sys.path.append(curPath)
sys.path.append(os.getcwd())
import AutoinputQE as inpp
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
		'size': 24, 
        }
font_size = 24
sys.dont_write_bytecode = True       ## not generate __pycache__ files

#------------------- Data manipulating ----------------------
#-------------------------  QE  -----------------------------
def get_fermi_level(filepath, filename="scf.log"):
    try:
        data = open(osp.join(filepath, filename), 'r').readlines()
        fermi_all = []
        for i in range(len(data)):
            if "the Fermi energy" in data[i]:
                fermi_all.append(float(data[i].split()[4]))
        # print(fermi)
        return fermi_all[-1]
    except:
        return 0.0
def get_electrons(filepath, filename='scf.log'):
    data = open(osp.join(filepath, filename), 'r').readlines()
    for i in range(len(data)):
        if 'number of electrons' in data[i]:
            return int(data[i].split('=')[1])
        
def band_data_pro(filepath, filename):
    file1 = open(osp.join(filepath, filename), 'r').readlines()
    file2 = open(osp.join(filepath, filename + ".gnu"), 'r').readlines()
    nbands = int(file1[0].split(',')[0].split('=')[1])
    nks = int(file1[0].split(',')[1].split("=")[1].split("/")[0])
    band_data = np.loadtxt(file2, dtype=float)
    band_k = band_data[:nks, 0]
    kk = band_data[:, 0]
    fermi = get_fermi_level(filepath)
    energy = band_data[:, 1:len(band_data[0, :])] - fermi
    hsp_labels = [str(items.split("!")[1]).strip() for items in inpp.kpoints_band.split("\n")[1:-1]]
    for i in range(len(hsp_labels)):
        # if str(hsp_labels[i]).lower() == "\Gamma".lower():
        #     hsp_labels[i] = u"Γ"
        if '\\' in hsp_labels[i].lower():
            hsp_labels[i] = "$" + hsp_labels[i] + "$"
    hsp_nums = len(hsp_labels)-1
    hsp_data = band_k[0:nks:(int(nks/(hsp_nums)))]
    if hsp_data[-1] != band_k[-1]:
        hsp_data = np.append(hsp_data, band_k[-1])
    hsp = {
        'klabels' : hsp_labels,
        'hspdata' : hsp_data,
    }
    return nbands, nks, hsp, band_k, kk, energy 

def plot_band(filepath, energyrange=[-3, 3]):
    fig = plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)
    linewidth = 2
    if inpp.SOC or inpp.noncolin:
        logical = True 
    else:
        logical = False
    if inpp.ISPIN == 2 and not logical:
        nbands, nks, hsp, band_k, kk, energy  = band_data_pro(filepath, "band_up.dat")
        for js in range(nbands-2):
            ax.plot(band_k, energy[js*nks:(js+1)*nks, 0], color="blue", lw=linewidth)
        ax.plot(band_k, energy[(nbands-1)*nks:((nbands-1)+1)*nks, 0], color="blue", lw=linewidth, label="Spin-up")
        nbands1, nks1, hsp1, band_k1, kk1, energy1  = band_data_pro(filepath, "band_dw.dat")
        for jd in range(nbands1-2):
            ax.plot(band_k1, energy1[jd*nks1:(jd+1)*nks1, 0], color="green", lw=linewidth)
        ax.plot(band_k1, energy1[(nbands1-1)*nks1:((nbands1-1)+1)*nks1, 0], color="green", lw=linewidth, label="Spin-down")
        leg = ax.legend(loc='upper right')
        labelss = leg.get_texts()
        [label.set_fontproperties(font_properties) for label in labelss]
        [label.set_fontsize(font_size-6) for label in labelss]
    else:
        nbands, nks, hsp, band_k, kk, energy  = band_data_pro(filepath, "band.dat")
        for j in range(nbands-1):
            ax.plot(band_k, energy[j*nks:(j+1)*nks, 0], color="blue", lw=linewidth)
        ax.plot(band_k, energy[(nbands-1)*nks:((nbands-1)+1)*nks, 0], color="red", lw=linewidth)
    print("High Symmerty Points: ")
    print(hsp['hspdata'], len(hsp['hspdata']))
    print(hsp['klabels'], len(hsp['klabels']))
    for i in range(0,len(hsp['klabels'])):
        ax.vlines(hsp['hspdata'][i], ymin=energyrange[0], ymax=energyrange[-1],linewidth=1,linestyle="--", color='gray')
    ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1], linewidth=1.0, linestyle="--", color='gray')
    ax.set_xticks(hsp['hspdata'])
    ax.set_xticklabels(hsp['klabels'])
    ax.set_ylabel("Energy (eV)", fontdict=font)
    ax.set_ylim(energyrange[0], energyrange[1])
    ax.set_xlim(0, np.max(band_k))
    labels = ax.get_xticklabels() + ax.get_yticklabels() 
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(font_size) for label in labels]
    plt.savefig("band.png", dpi=400, bbox_inches = 'tight', pad_inches=0.2)
    return band_k, energy.reshape(nbands, len(band_k)), hsp

def plot_band_mark(filepath, energyrange=[-3, 3], index_band=[]):
    fig = plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)
    linewidth = 2
    electrons = get_electrons()
    if inpp.SOC or inpp.noncolin:
        logical = True 
    else:
        logical = False
    if inpp.ISPIN == 2 and not logical:
        nbands, nks, hsp, band_k, kk, energy  = band_data_pro(filepath, "band_up.dat")
        for js in range(nbands-2):
            ax.plot(band_k, energy[js*nks:(js+1)*nks, 0], color="blue", lw=linewidth)
        ax.plot(band_k, energy[(nbands-1)*nks:((nbands-1)+1)*nks, 0], color="blue", lw=linewidth, label="Spin-up")
        nbands1, nks1, hsp1, band_k1, kk1, energy1  = band_data_pro(filepath, "band_dw.dat")
        for jd in range(nbands1-2):
            ax.plot(band_k1, energy1[jd*nks1:(jd+1)*nks1, 0], color="green", lw=linewidth)
        ax.plot(band_k1, energy1[(nbands1-1)*nks1:((nbands1-1)+1)*nks1, 0], color="green", lw=linewidth, label="Spin-down")
        leg = ax.legend(loc='upper right')
        labelss = leg.get_texts()
        [label.set_fontproperties(font_properties) for label in labelss]
        [label.set_fontsize(font_size-6) for label in labelss]
    else:
        nbands, nks, hsp, band_k, kk, energy  = band_data_pro(filepath, "band.dat")
        for j in range(nbands-1):
            ax.plot(band_k, energy[j*nks:(j+1)*nks, 0], color="blue", lw=linewidth)
        ax.plot(band_k, energy[(nbands-1)*nks:((nbands-1)+1)*nks, 0], color="red", lw=linewidth)
    mark_band = energy.reshape(nbands, len(band_k))
    if inpp.ISPIN == 2 and not logical:
        if index_band == []:
            index_band.append(int(electrons/2))
            ax.plot(band_k, energy[(nbands-1)*nks:((nbands-1)+1)*nks, 0], color="blue", lw=linewidth, label="Spin-up")
            ax.text(band_k[:], mark_band[index_band[0], kmark],
                    jj+1,
                    fontdict=font,
                    color=self.lcolor[index_c])
    else:
        pass
    print("High Symmerty Points: ")
    print(hsp['hspdata'], len(hsp['hspdata']))
    print(hsp['klabels'], len(hsp['klabels']))
    for i in range(0,len(hsp['klabels'])):
        ax.vlines(hsp['hspdata'][i], ymin=energyrange[0], ymax=energyrange[-1],linewidth=1,linestyle="--", color='gray')
    ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1], linewidth=1.0, linestyle="--", color='gray')
    ax.set_xticks(hsp['hspdata'])
    ax.set_xticklabels(hsp['klabels'])
    ax.set_ylabel("Energy (eV)", fontdict=font)
    ax.set_ylim(energyrange[0], energyrange[1])
    ax.set_xlim(0, np.max(band_k))
    labels = ax.get_xticklabels() + ax.get_yticklabels() 
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(font_size) for label in labels]
    plt.savefig("band.png", dpi=400, bbox_inches = 'tight', pad_inches=0.2)
    return band_k, energy.reshape(nbands, len(band_k)), hsp

if __name__ == '__main__':
    energy = [-6, 6]
    plot_band(os.getcwd(), energy)

