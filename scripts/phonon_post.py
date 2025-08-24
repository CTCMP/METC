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
import module.funcs as bas
from matplotlib import font_manager
sys.path.append(os.getcwd())
import AutoinputQE as inpp

curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
import configfile.inputpara as parameters

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
sys.dont_write_bytecode = True       ## not generate __pycache__ files

class phononPost:
    def __init__(self):
        self.filepath = os.getcwd()
        self.latt, self.element, self.num_element, self.tot_ele, self.cor = bas.get_poscar(self.filepath, "POSCAR")
        # self.linecolor = ['#B8860B', '#00FFFF', '#808000', '#00008B', '#006400', '#FF00FF', '#D2691E', '#DC143C', \
        #                 '#B8860B', '#00FFFF', '#00FF00', "#8B4513", "#FFB6C1", "#6A5ACD", "#1E90FF", "#40E0D0", "#008080", \
        #                 "#3CB371", "#00FF7F", "#6B8E23", "DAA520", "#DC143C", "#FF4500"] 
        self.linecolor = ['black', '#00FFFF', '#808000', '#00008B', '#006400', '#FF00FF', '#D2691E', '#DC143C', \
                        '#B8860B', '#00FFFF', '#00FF00', "#8B4513", "#FFB6C1", "#6A5ACD", "#1E90FF", "#40E0D0", "#008080", \
                        "#3CB371", "#00FF7F", "#6B8E23", "DAA520", "#DC143C", "#FF4500"] 
        self.linewidth = 1.5
        self.font_size = 22
        # self.yminvalue = -0.2
        self.yminvalue = 0
        self.cm2meV = 8.065
        # self.cm2THz = 33.35641
        # self.cm2meV = 33.35641

    def band_data_pro(self, filepath, filename):
        if osp.exists(osp.join(filepath, filename)): 
            print(osp.join(filepath, filename + ".gp"))
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
        # energy = band_data[:, 1:len(band_data[0, :])]/33.35641   # cm-1 to THz
        energy = band_data[:, 1:len(band_data[0, :])]/self.cm2meV    # cm-1 to meV
        hsp_labels = [str(items.split("!")[1]).strip() for items in inpp.kpoints_band.split("\n")[1:-1]]
        print(hsp_labels)
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

    def plot_band(self):
        filepath = self.filepath
        filename = "pwscf.freq"
        nbands, nks, hsp, band_k, kk, energy  = self.band_data_pro(filepath, filename)
        fig = plt.figure(figsize=(5, 4))
        ax = plt.subplot(111)
        for kk in range(0, len(energy[0, :])):
            ax.plot(band_k, energy[:, kk], color=self.linecolor[0], lw=self.linewidth)
        max_value = np.max(energy)
        for i in range(0,len(hsp['klabels'])):
            ax.vlines(hsp['hspdata'][i], ymin=0, ymax=max_value*1.5,linewidth=0.8,linestyle="--", color='gray')
        ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
        ax.set_xticks(hsp['hspdata'])
        ax.set_xticklabels(hsp['klabels'])
        ax.set_ylabel("Frequency (meV)", fontdict=font)
        if self.yminvalue !=0:
            ax.set_ylim(self.yminvalue, max_value*1.03)
        else:
            ax.set_ylim(np.min(energy), max_value*1.03)
        ax.set_xlim(0, np.max(band_k))
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-2) for label in labels]
        plt.savefig(osp.join(self.filepath, "phonon.png"), dpi=300, bbox_inches = 'tight', pad_inches=0.1)
        plt.close()
        return band_k, energy, hsp

    def phonon_dos_data_processing(self):
        file = osp.join(self.filepath, "phonon.dos")
        data = np.loadtxt(file, skiprows=1, dtype=float)
        phondos = []
        count = 0
        for i in range(len(self.element)):
            tmp = np.zeros(shape=(len(data[:, 0])))
            for j in range(self.num_element[i]):
                tmp = tmp + data[:, 2 + count]
                count = count + 1
            phondos.append(tmp)
        return data[:, 0]/self.cm2meV, phondos, data[:, 1] 
        # return data[:, 0]/8.065, phondos, data[:, 1]
    
    def phonon_dos_plotting(self):
        x, phondos, tot = self.phonon_dos_data_processing()
        fig = plt.figure(figsize=(5, 3))
        ax = plt.subplot(111)
        for pp in range(len(phondos)):
            ax.plot(x, phondos[pp], label=self.element[pp], color=self.linecolor[pp], lw=self.linewidth)
        ax.plot(x, tot, label='tot', color=self.linecolor[len(phondos)+1], lw=self.linewidth)
        ax.set_xlabel("Frequency (meV)", fontdict=font)
        ax.set_ylabel("DOS", fontdict=font)
        ax.set_xlim(self.yminvalue, np.max(x))

        leg = ax.legend(frameon=False, loc='upper right').get_texts()
        [label.set_fontproperties(font_properties) for label in leg]
        [label.set_fontsize(self.font_size-6) for label in leg]

        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-6) for label in labels]
        plt.savefig(osp.join(os.getcwd(), "phdos.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
        plt.close()

    def phonon_band_dos_plotting(self):
        import matplotlib.gridspec as gridspec
        gs = gridspec.GridSpec(1, 6)
        nbands, nks, hsp, band_k, kk, energy  = self.band_data_pro(self.filepath, "pwscf.freq")
        x, phondos, tot = self.phonon_dos_data_processing()
        fig = plt.figure(figsize=(10, 7))
        ############### phonon band 
        self.linewidth = 2
        ax = plt.subplot(gs[0:4])
        for kk in range(0, len(energy[0, :])):
            ax.plot(band_k, energy[:, kk], color=self.linecolor[0], lw=self.linewidth)
        max_value = np.max(energy)
        for i in range(0,len(hsp['klabels'])):
            ax.vlines(hsp['hspdata'][i], ymin=0, ymax=max_value*1.5, linewidth=0.8, linestyle="-", color='gray')
        ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="-", color='gray')
        ax.set_xticks(hsp['hspdata'])
        ax.set_xticklabels(hsp['klabels'])
        ax.set_ylabel("Frequency (meV)", fontdict=font)
        # ax.set_ylabel("Frequency (meV)", fontdict=font)
        egg = np.max(x)*1.02
        if self.yminvalue !=0:
            ax.set_ylim(self.yminvalue, egg)
        else:
            ax.set_ylim(np.min(x), egg)
        ax.set_xlim(0, np.max(band_k))
        ax.tick_params(bottom=False, top=False, left=False, right=False)
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-2) for label in labels]
        ############### phonon dos
        ax1 = plt.subplot(gs[4:6])
        for pp in range(len(phondos)):
            ax1.plot(phondos[pp], x, label=self.element[pp], color=self.linecolor[pp], lw=self.linewidth)
        ax1.plot(tot, x, label='tot', color=self.linecolor[len(phondos)+1], lw=self.linewidth)
        ax1.set_yticks([])
        ax1.set_xticks([])
        ax1.set_xlabel("DOS", fontdict=font)
        # ax1.set_ylim(self.yminvalue, np.max(x)*1.05)
        if self.yminvalue !=0:
            ax.set_ylim(self.yminvalue, egg)
        else:
            ax.set_ylim(np.min(x), egg)
        leg = ax1.legend(frameon=False, loc='lower right').get_texts()
        ax1.tick_params(bottom=False, top=False, left=False, right=False)
        [label.set_fontproperties(font_properties) for label in leg]
        [label.set_fontsize(self.font_size-4) for label in leg]

        labels = ax1.get_xticklabels() + ax1.get_yticklabels()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-4) for label in labels]
        plt.savefig(osp.join(os.getcwd(), "phon_banddos.png"), dpi=300, bbox_inches='tight', pad_inches=0.05)
        plt.close()

def manipulate_phonon_band():
    php = phononPost()
    kpts, phonon, hsp = php.plot_band()
    return kpts, phonon, hsp

def manipulate_phonon_dos():
    php = phononPost()
    php.phonon_dos_plotting()

def manipulate_phonon_band_dos():
    php = phononPost()
    php.phonon_band_dos_plotting()

if __name__ == '__main__':
    pass

