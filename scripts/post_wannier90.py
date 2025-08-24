# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: band.py
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import font_manager
import os, sys
import os.path as osp
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
import configfile.inputpara as parameters
sys.path.append(os.getcwd())
import AutoinputQE as inpp
sys.dont_write_bytecode = True       ## not generate __pycache__ files

# fontpath='/work/wangr/dxy/scripts/font/times.ttf'
fontpath=osp.join(osp.join(curPath, "../font"), 'arial.ttf')
#fontpath="D:\\JianGuoYun\\others\\keda\\old\\fermi\\font\\times.ttf"
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
#------------------- Data manipulating ----------------------
#-------------------------  QE  -----------------------------
def get_fermi_level(filepath, filename="scf.log"):
    try:
        data = open(osp.join(filepath, filename), 'r').readlines()
        fermi_all = []
        for i in range(len(data)):
            if "Fermi" in data[i]:
                fermi_all.append(float(data[i].split()[4]))
        # print(fermi)
        return fermi_all[-1]
    except:
        return 0.0
def band_data_pro(filepath, filename):
    if osp.exists(osp.join(curPath, "../configfile")): 
        file1 = open(osp.join(filepath, parameters.band_plot['filband']), 'r').readlines()
        file2 = open(osp.join(filepath, parameters.band_plot['filband'] + ".gnu"), 'r').readlines()
    else:
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
    hsp = {
        'klabels' : hsp_labels,
        'hspdata' : hsp_data,
    }
    # print(hsp['klabels'])
    return nbands, nks, hsp, band_k, kk, energy 

def w90_data_pro(filepath, filename='pwscf_band.dat'):
    file = osp.join(filepath, filename)
    data = np.loadtxt(file, skiprows=0, dtype=float)
    data_wannier90 = open(file, 'r').readlines()
    k_num_wannier90 = 0
    fermi = get_fermi_level(filepath)
    for index_wannier90 in range(0, len(data_wannier90)):
        if len(str(data_wannier90[index_wannier90]).strip()) == 0:
            k_num_wannier90 = index_wannier90
            break
    kpt = np.array(data[:k_num_wannier90, 0])
    energy = data[:, 1] - fermi
    return kpt, energy.reshape(int(len(data[:, 1])/k_num_wannier90), k_num_wannier90)

def plot_nihe(filepath, filename, energyrange=[-5, 5], label_fit="Wannier90-fitting"):
    nbands, nks, hsp, band_k, kk, energy  = band_data_pro(filepath, filename)
    kpt_w90, energy_w90 = w90_data_pro(filepath)
    kpt_w90_rescale = np.array(kpt_w90)/(np.max(kpt_w90)/hsp['hspdata'][-1])
    fig = plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)
    for j in range(nbands-2):
        ax.plot(band_k, energy[j*nks:(j+1)*nks, 0], color="blue", lw=1.5)
    ax.plot(band_k, energy[(nbands-1)*nks:((nbands-1)+1)*nks, 0], color="blue", lw=1.5, label="PBE")
    ########### wannier90
    for jj in range(len(energy_w90[:, 0])-1):
        ax.plot(kpt_w90_rescale, energy_w90[jj, :], '--', color="red", lw=1.5)
    ax.plot(kpt_w90_rescale, energy_w90[len(energy_w90[:, 0])-1, :], '--', color="red", label=label_fit, lw=1.5)
    for i in range(0,len(hsp['klabels'])):
        ax.vlines(hsp['hspdata'][i], ymin=energyrange[0], ymax=energyrange[-1],linewidth=0.8,linestyle="--", color='gray')
    ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
    ax.set_xticks(hsp['hspdata'])
    ax.set_xticklabels(hsp['klabels'])
    leg = ax.legend(loc="upper right")
    font_size = 24
    ax.set_ylabel("Energy (eV)", fontdict=font)
    ax.set_ylim(energyrange[0], energyrange[1])
    ax.set_xlim(0, np.max(band_k))
    labels = ax.get_xticklabels() + ax.get_yticklabels() 
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(font_size) for label in labels]
    labels1 = leg.get_texts() 
    [label.set_fontproperties(font_properties) for label in labels1]
    [label.set_fontsize(font_size-6) for label in labels1]
    plt.savefig("nihe.png", dpi=400, bbox_inches = 'tight', pad_inches=0.2)
    plt.close()

class Pw90:
    def __init__(self):
        self.file_path = os.getcwd()
        self.lcolor = ['#B8860B', '#00FFFF', '#808000', '#00008B', '#006400', '#FF00FF', '#D2691E', '#DC143C', \
                        '#B8860B', '#00FFFF', '#00FF00', "#8B4513", "#FFB6C1", "#6A5ACD", "#1E90FF", "#40E0D0", "#008080", \
                        "#3CB371", "#00FF7F", "#6B8E23", "#DAA520", "#DC143C", "#FF4500"] 
        self.defaultcolor = 'black'
        self.font_size=20
        self.lw = 2.0
    def w90_data(self, filepath, filename='pwscf_band.dat'):
        file = osp.join(filepath, filename)
        data = np.loadtxt(file, skiprows=0, dtype=float)
        data_wannier90 = open(file, 'r').readlines()
        k_num_wannier90 = 0
        fermi = get_fermi_level(filepath)
        for index_wannier90 in range(0, len(data_wannier90)):
            if len(str(data_wannier90[index_wannier90]).strip()) == 0:
                k_num_wannier90 = index_wannier90
                break
        kpt = np.array(data[:k_num_wannier90, 0])
        energy = data[:, 1] - fermi
        return k_num_wannier90, kpt, energy.reshape(int(len(data[:, 1])/k_num_wannier90), k_num_wannier90)

    def plot_nihe(self, filepath, filename, energyrange=[-5, 5], label_fit="Wannier90-fitting"):
        nbands, nks, hsp, band_k, kk, energy  = band_data_pro(filepath, filename)
        knum, kpt_w90, energy_w90 = self.w90_data(filepath)
        kpt_w90_rescale = np.array(kpt_w90)/(np.max(kpt_w90)/hsp['hspdata'][-1])
        fig = plt.figure(figsize=(8, 6))
        ax = plt.subplot(111)
        for j in range(nbands-2):
            ax.plot(band_k, energy[j*nks:(j+1)*nks, 0], color="blue", lw=1.5)
        ax.plot(band_k, energy[(nbands-1)*nks:((nbands-1)+1)*nks, 0], color="blue", lw=1.5, label="PBE")
        ########### wannier90
        for jj in range(len(energy_w90[:, 0])-1):
            ax.plot(kpt_w90_rescale, energy_w90[jj, :], '--', color="red", lw=1.5)
        ax.plot(kpt_w90_rescale, energy_w90[len(energy_w90[:, 0])-1, :], '--', color="red", label=label_fit, lw=1.5)
        for i in range(0,len(hsp['klabels'])):
            ax.vlines(hsp['hspdata'][i], ymin=energyrange[0], ymax=energyrange[-1],linewidth=0.8,linestyle="--", color='gray')
        ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
        ax.set_xticks(hsp['hspdata'])
        ax.set_xticklabels(hsp['klabels'])
        leg = ax.legend(loc="upper right")
        ax.set_ylabel("Energy (eV)", fontdict=font)
        ax.set_ylim(energyrange[0], energyrange[1])
        ax.set_xlim(0, np.max(band_k))
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size) for label in labels]
        labels1 = leg.get_texts() 
        [label.set_fontproperties(font_properties) for label in labels1]
        [label.set_fontsize(self.font_size-6) for label in labels1]
        plt.savefig("nihe.png", dpi=400, bbox_inches = 'tight', pad_inches=0.2)
        plt.close()

    def plot_occupiedBand(self, filepath, filename, energyrange=[-5, 5], label_fit="Wannier90-fitting"):
        nbands, nks, hsp, band_k, kk, energy  = band_data_pro(filepath, filename)
        knum, kpt_w90, energy_w90 = self.w90_data(filepath)
        kpt_w90_rescale = np.array(kpt_w90)/(np.max(kpt_w90)/hsp['hspdata'][-1])
        fig = plt.figure(figsize=(8, 6))
        ax = plt.subplot(111)
        for j in range(nbands-2):
            ax.plot(band_k, energy[j*nks:(j+1)*nks, 0], color="blue", lw=1.5)
        ax.plot(band_k, energy[(nbands-1)*nks:((nbands-1)+1)*nks, 0], color="blue", lw=1.5, label="PBE")
        ########### wannier90
        mark = [0, int(knum/20)]
        for jj in range(len(energy_w90[:, 0])-1):
            kmark = int(knum/2) - mark[int(jj%2)]
            index_c = jj
            if jj > len(self.lcolor)-1:
                index_c = 0
            # print(index_c, len(self.lcolor))
            ax.plot(kpt_w90_rescale, energy_w90[jj, :], '--', color=self.lcolor[index_c], lw=1.5)

            if energy_w90[jj, kmark] < energyrange[0] or energy_w90[jj, kmark] > energyrange[1]:
                continue
            else:
                ax.text(kpt_w90_rescale[kmark], energy_w90[jj, kmark],
                        jj+1,
                        fontdict=font,
                        color=self.lcolor[index_c])
        ax.plot(kpt_w90_rescale, energy_w90[len(energy_w90[:, 0])-1, :], '--', color=self.lcolor[0], label=label_fit, lw=1.5)
        if energy_w90[len(energy_w90[:, 0])-1, int(knum/2)] < energyrange[0] or energy_w90[len(energy_w90[:, 0])-1, int(knum/2)] > energyrange[1]:
            print("")
        else:
            ax.text(kpt_w90_rescale[int(knum/2)], energy_w90[len(energy_w90[:, 0])-1, int(knum/2)],
                    len(energy_w90[:, 0]),
                    fontdict=font,
                    color=self.lcolor[0])
        for i in range(0,len(hsp['klabels'])):
            ax.vlines(hsp['hspdata'][i], ymin=energyrange[0], ymax=energyrange[-1],linewidth=0.8,linestyle="--", color='gray')
        ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
        ax.set_xticks(hsp['hspdata'])
        ax.set_xticklabels(hsp['klabels'])
        leg = ax.legend(loc="upper right")
        ax.set_ylabel("Energy (eV)", fontdict=font)
        ax.set_ylim(energyrange[0], energyrange[1])
        ax.set_xlim(0, np.max(band_k))
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size) for label in labels]
        labels1 = leg.get_texts() 
        [label.set_fontproperties(font_properties) for label in labels1]
        [label.set_fontsize(self.font_size-6) for label in labels1]
        plt.savefig("Occupied_band.png", dpi=400, bbox_inches = 'tight', pad_inches=0.2)
        plt.close()