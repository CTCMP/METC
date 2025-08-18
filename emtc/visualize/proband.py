# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: band.py
"""
import numpy as np
import pandas as pd 
import matplotlib as mpl 
from matplotlib import pyplot as plt
import os, sys
import os.path as osp
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
sys.path.append(os.getcwd())
import AutoinputQE as inpp
import module.funcs as bas 
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

color_proband = ['red', 'green', 'blue', 'black', 'magenta', 'orange', 'cyan'] 

color_proband_down = [color_proband[len(color_proband)-itt-1] for itt in range(len(color_proband))]
#------------------- Data manipulating ----------------------
#-------------------------  QE  -----------------------------
def get_fermi_level(filepath, filename="scf.log"):
    data = open(osp.join(filepath, filename), 'r').readlines()
    fermi_all = []
    for i in range(len(data)):
        if "the Fermi energy" in data[i]:
            fermi_all.append(float(data[i].split()[4]))
    # print(fermi)
    # return fermi_all[-1]
    try:
        return fermi_all[-1]
    except:
        return 0

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

orbitals = {
    's' : ['s'],    
    'p' : ['pz', 'px', 'py'],
    'd' : ['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy'],
    'f' : ['f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7'],
}



class proBandStructure():
    def __init__(self, eng):
        self.filepath = os.getcwd()
        self.ff, self.bandff= self.get_ffs()
        self.nbands, self.nks, self.hsp, self.band_k, self.kk, self.energy = band_data_pro(os.getcwd(), self.bandff[0]) 
        self.latt, self.element, self.num_element, self.tot_ele, self.cor = bas.get_poscar(os.getcwd(), "POSCAR")
        self.num_lines = self.get_lineNum() + sum(self.num_element) + len(self.element)
        self.nkpts, self.nbnds = self.get_info()
        self.k_nbd_lines = self.nkpts * self.nbands
        self.efermi = bas.get_fermi_level(self.filepath, 'scf.log')
        self.eng = eng
        self.font_size = 20
        self.scalar = 60
        self.alpha = [1.0]
        for i in range(len(self.element)):
            self.alpha.append(0.7)

    def get_lineNum(self):
        data = open(osp.join(self.filepath, 'proband.projwfc_up'), 'r').readlines()
        for i in range(len(data)):
            if self.element[0] in data[i]:
                return i
    def get_info(self):
        pdb = pd.read_csv(osp.join(self.filepath, 'proband.projwfc_up'), sep='\s+', \
                          header=None, nrows=self.num_lines+1, engine='python')
        info = pdb.loc[self.num_lines-1].to_numpy()
        # print(info)
        return int(info[1]), int(info[2])
    
    def get_ffs(self):
        ffs = os.listdir(self.filepath)
        if 'proband.projwfc_down' in ffs:
            return ['proband.projwfc_up', 'proband.projwfc_down'], ['band_up.dat', 'band_dw.dat']
        else:
            return ['proband.projwfc_up'], ['band.dat']
        
    def write_tocsv(self, band, pbanddata, orbs, filename):
        df3 = pd.DataFrame(pbanddata, columns=orbs)
        df3['tot'] = df3.apply(lambda x: x.sum(), axis=1)
        df4 = pd.concat([band, df3], axis=1, sort=False, ignore_index=False)
        df4.to_csv(osp.join(self.filepath, "prodata", filename), sep='\t', index=None, float_format='%.6f')

    def data_process(self, filename, spinlabel):
        data = pd.read_csv(osp.join(self.filepath, filename), sep='\s+', \
                         header=None, skiprows=self.num_lines+1, \
                         names=[1, 2, 3, 4, 5, 6, 7], engine='python')
        df2 = data.drop(columns=[1, 2], axis=1)
        orb_info = df2.dropna(thresh=2, inplace=False)
        orbs_tmp = orb_info[3] + '-' + orb_info[4]
        orbs_ATOMS = orbs_tmp.to_numpy()
        orbs = []
        p = 0
        d = 0
        f = 0
        print('orbitals: ', orbs_ATOMS)
        for i in range(len(orbs_ATOMS)):
            if 's' in orbs_ATOMS[i].split('-')[1].lower():
                orbs.append(orbs_ATOMS[i][:-1] + orbitals['s'][0])
            elif 'p' in orbs_ATOMS[i].split('-')[1].lower():
                orbs.append(orbs_ATOMS[i][:-1] + orbitals['p'][p])
                p = p + 1
            elif 'd' in orbs_ATOMS[i].split('-')[1].lower():
                orbs.append(orbs_ATOMS[i][:-1] + orbitals['d'][d])
                d = d + 1
            elif 'f' in orbs_ATOMS[i].split('-')[1].lower():
                orbs.append(orbs_ATOMS[i][:-1] + orbitals['f'][f])
                f = f + 1
            else:
                print("Something Wrong with orbitals !")
            if p >=3:
                p=0
            if d >=5:
                d=0
            if f >= 7:
                f=0
        count_orbs = np.zeros(shape=(len(self.element)), dtype=int)
        for j in range(len(self.element)):
            for k in range(len(orbs)):
                if self.element[j] == orbs[k].split('-')[0]:
                    count_orbs[j] = count_orbs[j] + 1
        orbs_len = np.zeros(shape=(len(self.element)), dtype=int)
        for g in range(len(self.element)):
            orbs_len[g] = int(count_orbs[g]/self.num_element[g])
        orb_dict = {}
        if len(orbs_len) == 1:
            orb_dict[self.element[0]] = orbs[:orbs_len[0]]
        else:
            tmp = 0
            for e in range(0, len(orbs_len)):
                orb_dict[self.element[e]] = orbs[tmp:tmp+orbs_len[e]]
                tmp = orbs_len[e]*self.num_element[e] + tmp
        df2.dropna(thresh=1, inplace=True) # 删除空列
        df2.dropna(axis=1, inplace=True)   # 删除空列
        df2.drop(df2[df2[3].isin(self.element)].index, inplace=True) # 删除元素所在的行
        npdata1 = df2.to_numpy(dtype=float)
        npdata2 = npdata1.reshape(-1, self.k_nbd_lines)
        nums_orbs = len(orbs_ATOMS)
        npdata3 = np.zeros(shape=(nums_orbs, self.k_nbd_lines))
        for i in range(nums_orbs):
            # ttmp = np.array(npdata2[i, :].reshape(-1, self.nbnds)).flatten(order='F')    # 将每个轨道的数据分为 nks, nbnds的形状，之后按行展平
            ttmp = npdata2[i, :].reshape(-1, self.nbnds).flatten(order='F')    # 将每个轨道的数据分为 nks, nbnds的形状，之后按行展平
            npdata3[i, :] = ttmp
        npdata = npdata3.flatten(order='C').reshape(-1, self.k_nbd_lines).T

        if not osp.exists(osp.join(os.getcwd(), 'prodata')):
            os.makedirs(osp.join(os.getcwd(), 'prodata'))
        if spinlabel=='_up':
            banddata_tmp = pd.read_csv(osp.join(self.filepath, "band_up.dat.gnu"), delim_whitespace=True, \
                                    names=['Kpts', "Energy"])
        elif spinlabel=='_down':
            banddata_tmp = pd.read_csv(osp.join(self.filepath, "band_dw.dat.gnu"), delim_whitespace=True, \
                                    names=['Kpts', "Energy"])
        else:
            banddata_tmp = pd.read_csv(osp.join(self.filepath, "band.dat.gnu"), delim_whitespace=True, \
                            names=['Kpts', "Energy"])
        banddata1 = banddata_tmp.to_numpy()
        banddata1[:, 1] = banddata1[:, 1] - self.efermi
        banddata = pd.DataFrame(banddata1, columns=['Kpts', 'Energy'])
        dict_prodata = {}
        count_pro = 0
        for ii in range(len(orbs_len)):
            tmp_data = np.zeros(shape=(self.k_nbd_lines, orbs_len[ii]))
            for jj in range(self.num_element[ii]):
                tmp_data = tmp_data + npdata[:, count_pro:count_pro+orbs_len[ii]]
                count_pro = count_pro + orbs_len[ii]
            dict_prodata[self.element[ii]] = tmp_data
            self.write_tocsv(banddata, tmp_data, orb_dict[self.element[ii]], self.element[ii] + "_pband" + spinlabel + ".dat")

    def plot_single_fig(self, ax, filename='band.dat', clor='gray'):
        linewidth = 0.5
        nbands, nks, hsp, band_k, kk, energy  = band_data_pro(self.filepath, filename)
        for j in range(nbands):
            ax.plot(band_k, energy[j*nks:(j+1)*nks, 0], color='gray', lw=0.5, alpha=0.4)
        for i in range(0,len(hsp['klabels'])):
            ax.vlines(hsp['hspdata'][i], ymin=self.eng[0], ymax=self.eng[-1],linewidth=0.8,linestyle="--", color='gray')
        ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1], linewidth=1.0, linestyle="--", color='gray')
        ax.set_xticks(hsp['hspdata'])
        ax.set_xticklabels(hsp['klabels'])
        ax.set_ylabel("Energy (eV)", fontdict=font)
        return ax 
    
    def plot_fig(self, filename, spinlabel):
        print("-------- begin plotting projected bands for every atoms and every orbitals -------")
        if spinlabel:
            data = pd.read_csv(osp.join(self.filepath, "prodata", filename + '_up.dat'), sep='\s+', \
                            engine='python')
            data_dw = pd.read_csv(osp.join(self.filepath, "prodata", filename + '_down.dat'), sep='\s+', \
                            engine='python')
            orblabels = list(data)
            banddata = data.to_numpy()
            banddata_dw = data_dw.to_numpy()
        else:
            data = pd.read_csv(osp.join(self.filepath, "prodata", filename + '.dat'), sep='\s+', \
                            engine='python')
            orblabels = list(data)
            banddata = data.to_numpy()
        for i in range(2, len(banddata[0, :])-1):
            plt.figure(figsize=(4, 5))
            ax = plt.subplot(111)
            if spinlabel:
                ax = self.plot_single_fig(ax, "band_up.dat", 'red')
                ax = self.plot_single_fig(ax, "band_dw.dat", 'blue')
                ax.scatter(banddata[:, 0], banddata[:, 1], color='red', s=banddata[:, i]*self.scalar, marker='.', edgecolors='none', label="SPIN UP")
                ax.scatter(banddata_dw[:, 0], banddata_dw[:, 1], color='blue', s=banddata_dw[:, i]*self.scalar, marker='.', edgecolors='none', label="SPIN DOWN")
                leg = ax.legend(frameon=True, loc='upper right').get_texts()
                [label.set_fontproperties(font_properties) for label in leg]
                [label.set_fontsize(self.font_size-8) for label in leg]
            else:
                ax = self.plot_single_fig(ax, "band.dat")
                ax.scatter(banddata[:, 0], banddata[:, 1], color='red', s=banddata[:, i]*self.scalar, marker='.', edgecolors='none')
            ax.set_xlim(0, np.max(banddata[:, 0]))
            ax.set_ylim(self.eng)
            ax.set_title(orblabels[i])
            labels = ax.get_xticklabels() + ax.get_yticklabels()
            [label.set_fontproperties(font_properties) for label in labels]
            [label.set_fontsize(font_size) for label in labels]
            plt.savefig(osp.join(self.filepath, 'prodata', str(orblabels[i]) + '.png'), dpi=150, bbox_inches = 'tight')
            plt.close()

    def plot_fig_element(self, filenamelist, spinlabel):       # 元素投影能带图
        print("-------- begin plotting projected band for elements -------")
        plt.figure(figsize=(4, 5))
        ax = plt.subplot(111)
        for j in range(len(filenamelist)):
            if spinlabel:
                ax = self.plot_single_fig(ax, 'band_up.dat', 'red')
                ax = self.plot_single_fig(ax, 'band_dw.dat', 'blue')
                data2 = pd.read_csv(osp.join(self.filepath, "prodata", filenamelist[j]+'_up.dat'), sep='\s+', \
                                engine='python')
                data = data2.to_numpy()
                data3 = pd.read_csv(osp.join(self.filepath, "prodata", filenamelist[j]+'_down.dat'), sep='\s+', \
                                engine='python')
                data_dw = data3.to_numpy()
                ax.scatter(data[:, 0], data[:, 1], s=data[:, -1]*self.scalar, marker='o', edgecolors=color_proband[j], \
                        color='none', label=self.element[j] + " : UP", linewidths=0.5, alpha=self.alpha[j])
                ax.scatter(data_dw[:, 0], data_dw[:, 1], s=data_dw[:, -1]*self.scalar, marker='o', edgecolors=color_proband_down[j], \
                        color='none', label=self.element[j] + " : DW", linewidths=0.5, alpha=self.alpha[j])
            else:
                ax = self.plot_single_fig(ax, 'band.dat')
                data2 = pd.read_csv(osp.join(self.filepath, "prodata", filenamelist[j]+'.dat'), sep='\s+', \
                                engine='python')  
                data = data2.to_numpy()
                ax.scatter(data[:, 0], data[:, 1], s=data[:, -1]*self.scalar, marker='o', edgecolors=color_proband[j], \
                        color='none', label=self.element[j], linewidths=0.5, alpha=self.alpha[j])
        ax.set_xlim(0, np.max(data[:, 0]))
        ax.set_ylim(self.eng)
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-2) for label in labels]

        leg = ax.legend(frameon=True, loc='upper right').get_texts()
        [label.set_fontproperties(font_properties) for label in leg]
        [label.set_fontsize(self.font_size-8) for label in leg]
        plt.savefig(osp.join(self.filepath, 'prodata', 'Proband_elements.png'), dpi=150, bbox_inches = 'tight')
        plt.close()

    def manipulate_data(self):
        spinlabel = ['_up', '_down']
        print(self.ff)
        for ii in range(len(self.ff)):
            if len(self.ff) == 2:
                print("-------- begin processing %s -------" %(self.ff[ii]))
                self.data_process(self.ff[ii], spinlabel[ii])
            else:
                print("-------- begin processing %s -------" %(self.ff[ii]))
                self.data_process(self.ff[ii], '')
        sp = False
        if len(self.ff) == 2:
            sp = True 
            for jj in range(len(self.element)):
                self.plot_fig(self.element[jj]+"_pband", sp)
            self.plot_fig_element([self.element[item]+"_pband" for item in range(len(self.element))], sp)
        else:
            for jj in range(len(self.element)):
                self.plot_fig(self.element[jj]+"_pband", sp)
            self.plot_fig_element([self.element[item]+"_pband" for item in range(len(self.element))], sp)

def manipulate(eng):
    # eng = [-20, 8]
    pbs = proBandStructure(eng)
    pbs.manipulate_data()

if __name__ == '__main__':
    manipulate()

