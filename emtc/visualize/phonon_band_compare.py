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
import configfile.inputpara as parameters
sys.path.append(os.getcwd())
import AutoinputQE as inpp
import matplotlib as mpl
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
		'size': 20, 
        }
sys.dont_write_bytecode = True       ## not generate __pycache__ files

class phononcompare:
    def __init__(self, filepathlist, labellist):
        self.filepath = os.getcwd()
        self.filepathlist = filepathlist 
        self.filename = "pwscf.freq"
        self.labellist = labellist
        self.linecolor = ['red', 'blue', 'orange', 'black', '#B8860B', '#00FFFF', '#808000', '#00008B', '#006400', '#FF00FF', '#D2691E', '#DC143C', \
                        '#B8860B', '#00FFFF', '#00FF00', "#8B4513", "#FFB6C1", "#6A5ACD", "#1E90FF", "#40E0D0", "#008080", \
                        "#3CB371", "#00FF7F", "#6B8E23", "DAA520", "#DC143C", "#FF4500"] 
        self.font_size = 20
        self.cm2meV = 8.065

    def band_data_pro(self, filepath, filename):
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
        energy = band_data[:, 1:len(band_data[0, :])]/self.cm2meV
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

    def plot_band_compare(self):
        fig = plt.figure(figsize=(5, 4))
        ax = plt.subplot(111)
        max_eng = []
        ori_path = os.getcwd()
        for ff in range(len(self.filepathlist)):
            print("Filepath: ", self.filepathlist[ff])
            nbands, nks, hsp, band_k, kk, energy  = self.band_data_pro(self.filepathlist[ff], self.filename)
            for kk in range(0, len(energy[0, :])):
                ax.plot(band_k, energy[:, kk], color=self.linecolor[ff], lw=1.0)
            max_eng.append(np.max(energy))
        max_value = max(max_eng)
        for i in range(0,len(hsp['klabels'])):
            ax.vlines(hsp['hspdata'][i], ymin=0, ymax=max_value*1.5,linewidth=0.8,linestyle="--", color='gray')
        ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
        ax.set_xticks(hsp['hspdata'])
        ax.set_xticklabels(hsp['klabels'])
        ax.set_ylabel("Frequency (meV)", fontdict=font)
        ax.set_ylim(-5, max_value*1.05)
        # ax.set_ylim(-1, 20)
        ax.set_xlim(0, np.max(band_k))
        for i in range(len(self.labellist)):
            ax.plot([],[], c=self.linecolor[i], label=self.labellist[i])
        leg = ax.legend(frameon=True, loc='upper right').get_texts()
        [label.set_fontproperties(font_properties) for label in leg]
        [label.set_fontsize(self.font_size-2) for label in leg]

        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-2) for label in labels]
        plt.savefig(osp.join(self.filepath, "phon_compare.png"), dpi=300, bbox_inches = 'tight', pad_inches=0.1)
        plt.close()

def manipulate():
    filepath = input("please input filepath for comparing: \n").split()
    labellist = input("please input label for comparing: \n").split()
    filepathlist = [osp.join(os.getcwd(), item) for item in filepath]
    phcom = phononcompare(filepathlist, labellist)
    phcom.plot_band_compare()

if __name__ == '__main__':
    manipulate()