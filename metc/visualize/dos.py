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
sys.path.append(os.getcwd())
import AutoinputQE as inpp
import re 
from matplotlib import font_manager
curPath = os.path.abspath(os.path.dirname(__file__))
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
		'size': 16, 
        }
sys.dont_write_bytecode = True       ## not generate __pycache__ files

def get_poscar(filepath="", filename='POSCAR'):
    filepath = osp.join(filepath, filename)
    poscar1 = open(filepath, 'r')
    poscar_lines = poscar1.readlines()
    pos = []
    for flag_lines in range(2, 5):
        for i in range(0, 3):
            pos.append(float(poscar_lines[flag_lines].split()[i]))
    latt = np.array(pos).reshape(3, 3)
    ############# get element ###################
    poscar = open(filepath, 'r').read().strip('\n').splitlines()
    elements = poscar[5].lstrip().split()
    element = []
    num_element = []
    numbers = poscar[6].lstrip().split()
    for flag_element in range(0, len(elements)):
        ele = elements[flag_element]
        element.append(ele)
    for flag_element_num in range(0, len(numbers)):
        ele = int(numbers[flag_element_num])
        num_element.append(ele)
    tot_ele = sum(num_element)
    ############# get coordinate ###################
    # cor = np.loadtxt(filepath, skiprows=8, dtype=np.str_, encoding="UTF-8")
    cor_tmp = np.loadtxt(filepath, skiprows=8)
    if tot_ele == 1:
        tmp = [cor_tmp]
        cor = np.array(tmp)
    else:
        cor = cor_tmp[:tot_ele, 0:3]
    return np.round(latt, 10), element, num_element, tot_ele, np.round(cor, 10)

class densityofstate:
    def __init__(self, xrange, energyrange):
        self.filepath = os.getcwd()
        self.latt, self.element, self.num_element, self.tot_ele, self.cor = get_poscar(os.getcwd(), 'POSCAR')
        self.xrange = xrange
        self.linecolor = ['#B8860B', '#00FFFF', '#808000', '#00008B', '#006400', '#FF00FF', '#D2691E', '#DC143C', \
                        '#B8860B', '#00FFFF', '#00FF00', "#8B4513", "#FFB6C1", "#6A5ACD", "#1E90FF", "#40E0D0", "#008080", \
                        "#3CB371", "#00FF7F", "#6B8E23", "DAA520", "#DC143C", "#FF4500"] 
        self.energyrange = energyrange
        self.font_size = 16
        self.linewidth = 1

    def get_fermi_level(self, filename="scf.log"):
        try:
            data = open(osp.join(self.filepath, filename), 'r').readlines()
            fermi_all = []
            for i in range(len(data)):
                if "Fermi" in data[i]:
                    fermi_all.append(float(data[i].split()[4]))
            return fermi_all[-1]
        except:
            return 0.0

    def write_data(self, energy, data, tot):
        orbit = ['s', 'p', 'd']
        for ii in range(len(self.element)):
            conf = open(osp.join(os.getcwd(), 'PDOS_' + self.element[ii] + '.dat'), 'w')
            conf.write("Energy".rjust(10, ' ') + "  " + 's'.rjust(10, ' ') + "  " + 'p'.rjust(10, ' ') + "  " + 'd'.rjust(10, ' ') + "  "  + 'tot'.rjust(10, ' ') + '\n')
            for kk in range(len(data[ii][0])):
                conf.write(str(np.round(energy[kk], 8)).rjust(10, " ") )
                for orb in range(len(data[ii])):
                    conf.write("  " + str(np.round(data[ii][orb][kk], 8)).rjust(10, " ") + "  ")
                conf.write(str(np.round(tot[ii][kk], 8)).rjust(10, " ") + '\n')
            conf.close()

    def data_process_nospin(self):
        filelist = os.listdir()
        flist = []
        for ff in range(len(filelist)):
            if 'pwscf.pdos_atm' in filelist[ff]:
                flist.append(filelist[ff])
        pdos = []
        file_element = []
        for fe in range(len(self.element)):
            tmp = []
            for fee in range(len(flist)):
                if '('+ self.element[fe] +')' in flist[fee]:
                    tmp.append(flist[fee])
            file_element.append(tmp)
        index = []
        for ft in range(len(file_element)):
            index_tmp = 0
            for itf in range(len(file_element[ft])):
                orbit = re.findall('\(\D+\)', file_element[ft][itf])[1]
                if 'd' in orbit:
                    index_tmp = index_tmp + 1
            if index_tmp >= 1:
                index.append(1)
            else:
                index.append(0)
        pdos = []
        pdos_ele = []
        tmpdata = np.loadtxt(osp.join(os.getcwd(), file_element[0][0]))
        len_data = len(tmpdata[:, 0])
        for ii in range(len(file_element)):
            if index[ii] >=1:
                tmp_pdos = np.zeros(shape=(3, len_data))
                tmp_pdos_ele = np.zeros(shape=(len_data))
                for item in range(len(file_element[ii])):
                    orbital = re.findall('\(\D+\)', file_element[ii][item])[1]
                    data = np.loadtxt(osp.join(file_element[ii][item]), skiprows=1, dtype=float)
                    if 's' in orbital:
                        tmp_pdos[0, :] = tmp_pdos[0, :] + data[:, 2] 
                    elif 'p' in orbital:
                        tmp_pdos[1, :] = tmp_pdos[1, :] + data[:, 2] + data[:, 3] + data[:, 4]
                    elif 'd' in orbital:
                        tmp_pdos[2, :] = tmp_pdos[2, :] + data[:, 2] + data[:, 3] + data[:, 4] + data[:, 5] + data[:, 6]
                tmp_pdos_ele = tmp_pdos[0, :] + tmp_pdos[1, :] + tmp_pdos[2, :]
                pdos.append(tmp_pdos)
                pdos_ele.append(tmp_pdos_ele)
            else:
                tmp_pdos1 = np.zeros(shape=(2, len_data))
                tmp_pdos_ele1 = np.zeros(shape=(len_data))
                for item in range(len(file_element[ii])):
                    orbital1 = re.findall('\(\D+\)', file_element[ii][item])[1]
                    data1 = np.loadtxt(osp.join(file_element[ii][item]), skiprows=1, dtype=float)
                    if 's' in orbital1:
                        tmp_pdos1[0, :] = tmp_pdos1[0, :] + data1[:, 2] 
                    elif 'p' in orbital1:
                        tmp_pdos1[1, :] = tmp_pdos1[1, :] + data1[:, 2] + data1[:, 3] + data1[:, 4] 
                tmp_pdos_ele1 = tmp_pdos1[0, :] + tmp_pdos1[1, :]
                pdos.append(tmp_pdos1)
                pdos_ele.append(tmp_pdos_ele1)
        self.plot_pdos_spd(tmpdata[:, 0] - self.get_fermi_level(), pdos)

        self.plot_pdos_element(tmpdata[:, 0] - self.get_fermi_level(), pdos_ele)
        self.write_data(tmpdata[:, 0] - self.get_fermi_level(), pdos, pdos_ele)

    def plot_pdos_spd(self, eng, pdos):
        fig = plt.figure(figsize=(5, 3))
        ax = plt.subplot(111)
        orbit = ['s', 'p', 'd']
        count = 0
        for ii in range(len(pdos)):
            for jj in range(len(pdos[ii])):
                ax.plot(eng, pdos[ii][jj], label=self.element[ii] + '_' + orbit[jj], color=self.linecolor[count], lw=self.linewidth)
                count = count + 1
        tdos = np.loadtxt(osp.join(os.getcwd(), "pwscf.pdos_tot"), skiprows=1, dtype=float)
        ax.plot(eng, tdos[:, 1], label='tot', color='black', lw=self.linewidth)
        ax.vlines(0, ymin=0, ymax=100,linewidth=0.8,linestyle="--", color='gray')
        # ax.axhline(0, xmin=-10, xmax=10, linewidth=1.0,linestyle="--", color='gray')
        ax.set_xlabel("Energy (eV)", fontdict=font)
        ax.set_yticks([])
        ax.set_ylabel('DOS (states/eV)', fontdict=font)
        ax.set_ylim(self.energyrange)
        ax.set_xlim(self.xrange)
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-4) for label in labels]

        leg = ax.legend(frameon=True, loc='upper right').get_texts()
        [label.set_fontproperties(font_properties) for label in leg]
        [label.set_fontsize(self.font_size-6) for label in leg]
        plt.savefig(osp.join(os.getcwd(), "pdos_spd.png"), dpi=300, bbox_inches = 'tight', pad_inches=0.1)
        plt.close()

    def plot_pdos_element(self, eng, pdos):
        fig = plt.figure(figsize=(5, 3))
        ax = plt.subplot(111)
        orbit = ['s', 'p', 'd']
        count = 0
        tot = np.zeros(len(eng))
        for ii in range(len(pdos)):
            ax.plot(eng, pdos[ii], label=self.element[ii], color=self.linecolor[ii], lw=self.linewidth)
            tot = tot + pdos[ii]
        ax.plot(eng, tot, label='tot', color='black', lw=self.linewidth)
        ax.vlines(0, ymin=0, ymax=100,linewidth=0.8,linestyle="--", color='gray')
        # ax.axhline(0, xmin=-10, xmax=10, linewidth=1.0,linestyle="--", color='gray')
        ax.set_xlabel("Energy (eV)", fontdict=font)
        ax.set_yticks([])
        ax.set_ylabel('DOS (states/eV)', fontdict=font)
        ax.set_ylim(self.energyrange)
        ax.set_xlim(self.xrange)
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-4) for label in labels]

        leg = ax.legend(frameon=True, loc='upper right').get_texts()
        [label.set_fontproperties(font_properties) for label in leg]
        [label.set_fontsize(self.font_size-6) for label in leg]
        plt.savefig(osp.join(os.getcwd(), "pdos_element.png"), dpi=300, bbox_inches = 'tight', pad_inches=0.1)
        plt.close()

def manipulate_post():
    ds = densityofstate([-5, 5], [0, 10])
    ds.data_process_nospin()
def manipulate():
    print("===================================================================")
    print("           Please input x and y range for dos plot                 ")
    print("------------------->>")
    tmp_x = input("Input x: ").split()
    x_range = [int(item) for item in tmp_x]
    tmp_y = input("Input y: ").split()
    y_range = [int(item) for item in tmp_y]
    ds = densityofstate(x_range, y_range)
    ds.data_process_nospin()

if __name__ == '__main__':
    manipulate_post()

