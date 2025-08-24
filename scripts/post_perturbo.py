# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: band.py
"""
import os, sys
import os.path as osp
import numpy as np 
from matplotlib import pyplot as plt 
from matplotlib import font_manager
import module.funcs as bas
import AutoinputQE as inpp
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
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

class perturboPost:
    def __init__(self, eng):
        self.file_path = os.getcwd()
        self.scatterpoints_size = 8
        self.pointcolor = ['#B8860B', '#00FFFF', '#808000', '#00008B', '#006400', '#FF00FF', '#D2691E', '#DC143C', \
                        '#B8860B', '#00FFFF', '#00FF00', "#8B4513", "#FFB6C1", "#6A5ACD", "#1E90FF", "#40E0D0", "#008080", \
                        "#3CB371", "#00FF7F", "#6B8E23", "DAA520", "#DC143C", "#FF4500"] 
        self.font_size=20
        self.mev2freq = 1/(4.1356676969)
        self.lw = 2.0
        self.cmTomeV = 1/8.06554429
        if not osp.exists(osp.join(self.file_path, '../../bands')):
            self.file_path_bands = osp.join(self.file_path, '../../w90')
        else:
            self.file_path_bands = osp.join(self.file_path, '../../bands')
        self.efermi = bas.get_fermi_level(self.file_path_bands, 'scf.log')
        self.filepath_phonon = '../../phonon'
        self.eng = eng
        self.fig_sz = (5, 4)
        self.fig_sz_large = (10, 4)
        self.log = True
        self.bwith = 1.5
        self.font_diff = 6

    def perturbo_bands(self, filename):
        data = np.loadtxt(osp.join(self.file_path, filename), skiprows=0, dtype=float)
        kpt = set(data[:, 0])
        nk = len(set(kpt))
        kpt = data[:nk, 0]
        energy = data[:, -1]-self.efermi
        return kpt, energy.reshape(-1, nk)

    def perturbo_phonon(self, filename):
        data = np.loadtxt(osp.join(self.file_path, filename), skiprows=0, dtype=float)
        kpt = set(data[:, 0])
        nk = len(list(kpt))
        kpt = data[:nk, 0]
        energy = data[:, -1]
        return kpt, energy.reshape(-1, nk)

    def plot_bands(self, name):
        fig = plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        kpt, bands = self.perturbo_bands('pwscf.bands')
        for i in range(len(bands)):
            ax.plot(kpt, bands[i], color='blue', lw=1.5)
        ax.set_xlim(0, np.max(kpt))
        ax.set_ylim(self.eng)
        ax.set_ylabel('Energy (eV)')
        ax.set_xlabel('High Sym Kpts')
        ax.set_xticks([])
        ax.set_xticklabels([])
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size) for label in labels]
        plt.savefig(osp.join(self.file_path, name+".png"), dpi=400, bbox_inches = 'tight', pad_inches=0.2)
        plt.close()

    def plot_phonon_single(self, name):
        fig = plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        kpt, bands = self.perturbo_phonon('pwscf.phdisp')
        for i in range(len(bands)):
            ax.plot(kpt, bands[i], color='blue', lw=1.5)
        ax.set_xlim(0, np.max(kpt))
        ax.set_ylim(self.eng)
        ax.set_ylabel('Energy (meV)')
        ax.set_xlabel('High Sym Kpts')
        ax.set_xticks([])
        ax.set_xticklabels([])
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size) for label in labels]
        plt.savefig(osp.join(self.file_path, name+".png"), dpi=400, bbox_inches = 'tight', pad_inches=0.2)
        plt.close()

    def plot_bands_compare(self, filename, label_fit="MLWFs"):
        nbands, nks, hsp, band_k, kk, energy  = bas.band_data_pro_band(self.file_path_bands, filename)
        kpt_perturbo, energy_perturbo = self.perturbo_bands('pwscf.bands')
        kpt_mlwf_rescale = np.array(kpt_perturbo)/(np.max(kpt_perturbo)/hsp['hspdata'][-1])
        fig = plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        for j in range(nbands-2):
            ax.plot(band_k, energy[j*nks:(j+1)*nks, 0], color="blue", lw=1.5)
        ax.plot(band_k, energy[(nbands-1)*nks:((nbands-1)+1)*nks, 0], color="blue", lw=1.5, label="PBE")
        ########### Perturbo bands
        for jj in range(len(energy_perturbo)-1):
            ax.plot(kpt_mlwf_rescale, energy_perturbo[jj, :], '--', color="red", lw=1.5)
        ax.plot(kpt_mlwf_rescale, energy_perturbo[-1, :], '--', color="red", label=label_fit, lw=1.5)
        ###########################
        for i in range(0,len(hsp['klabels'])):
            ax.vlines(hsp['hspdata'][i], ymin=self.eng[0], ymax=self.eng[-1],linewidth=0.8,linestyle="--", color='gray')
        ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
        ax.set_xticks(hsp['hspdata'])
        ax.set_xticklabels(hsp['klabels'])
        leg = ax.legend(loc="upper right")
        ax.set_ylabel("Energy (eV)", fontdict=font)
        ax.set_ylim(self.eng)
        ax.set_xlim(0, np.max(band_k))
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size) for label in labels]
        labels1 = leg.get_texts() 
        [label.set_fontproperties(font_properties) for label in labels1]
        [label.set_fontsize(self.font_size-8) for label in labels1]
        plt.savefig(osp.join(self.file_path, "band_nihe.png"), dpi=400, bbox_inches = 'tight', pad_inches=0.2)
        plt.close()

    def plot_phonon(self):
        import scripts.phonon_band as pb
        nbands, nks, hsp, band_k, kk, energy  = pb.band_data_pro(self.filepath_phonon, 'pwscf.freq')
        fig = plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        kpt, bands = self.perturbo_phonon('pwscf.phdisp')
        kpt_mlwf_rescale = np.array(kpt)/(np.max(kpt)/hsp['hspdata'][-1])
        for kk in range(0, len(energy[0, :])-1):
            ax.plot(band_k, energy[:, kk], color='black', lw=2.0)
        ax.plot(band_k, energy[:, -1], color='black', lw=2.0, label="PBE")
        max_value = np.max(energy)
        ####### EPW plot
        for bp in range(len(bands)-1):
            ax.plot(kpt_mlwf_rescale, bands[bp, :], color='red', ls='--')
        ax.plot(kpt_mlwf_rescale, bands[-1, :], color='red', ls='--', label='Perturbo')
        ####### EPW plot
        for i in range(0,len(hsp['klabels'])):
            ax.vlines(hsp['hspdata'][i], ymin=0, ymax=max_value*1.5,linewidth=0.8,linestyle="--", color='gray')
        ax.set_title("Phonon", fontdict=font)
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
        leg = ax.legend(frameon=True, loc='upper right')
        labels = leg.get_texts()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-6) for label in labels]
        plt.savefig(osp.join(self.file_path, "phonon_perturbo.png"), dpi=400, bbox_inches = 'tight', pad_inches=0.1)
    
    def write_it(self, label, E, data):
        wt_modes = open(osp.join(os.getcwd(), 'rt_data', str(label) + '_IT'+'.dat'), 'w')
        wt_modes.write("E(ibnd)(eV)".ljust(12, ' ')+ '    ')
        for k in range(len(data)):
            wt_modes.write(("it: " + str(k+1)).rjust(12, ' ') + '    ')
        wt_modes.write('\n')
        for i in range(len(E)):
            wt_modes.write(str(np.round(E[i], 8)).rjust(12, ' '))
            for j in range(len(data)):
                wt_modes.write("    " + str(np.round(data[j, i], 8)).rjust(12, ' '))
            wt_modes.write('\n')
        wt_modes.close()
    def write_modes_resolved_data(self, iband_len, label, E, data):
        tmp_E = []
        tmp_data = []
        len_data = int(len(E[0])/iband_len)
        for i in range(len(E)):
            tmp_xE = np.zeros(shape=(iband_len, len_data))
            tmp_ySR = np.zeros(shape=(iband_len, len_data))
            for j in range(len_data):
                for mode in range(iband_len):
                    tmp_xE[mode, j] = E[i, j*iband_len + mode]
                    tmp_ySR[mode, j] = data[i, j*iband_len + mode]
            tmp_E.append(tmp_xE) 
            tmp_data.append(tmp_ySR)
        for ii in range(iband_len):
            wt_modes_resolved = open(osp.join(os.getcwd(), 'rt_data', str(label)+"_modes_iband_"+str(ii+1)+".dat"), 'w')
            wt_modes_resolved.write("E(ibnd)(eV)".ljust(12, ' ')+ '    ')
            for k in range(len(data)):
                wt_modes_resolved.write(("it: " + str(k+1)).rjust(12, ' ') + '    ')
            wt_modes_resolved.write('\n')
            for jj in range(len(tmp_E[0][0])):
                wt_modes_resolved.write(str(np.round(tmp_E[0][ii][jj], 8)).rjust(12, ' '))
                for jjj in range(len(tmp_E)):
                    wt_modes_resolved.write("    " + str(np.round(tmp_data[jjj][ii][jj], 8)).rjust(12, ' '))
                wt_modes_resolved.write('\n')
            wt_modes_resolved.close()
    def relaxation_time(self):
        data = np.loadtxt(osp.join(os.getcwd(), 'relaxation_time.dat'), skiprows=1, dtype=float)
        ik = set(data[:, 0])
        ik_len = len(ik)
        iband = set(data[:, 2])
        iband_len = len(iband)
        modes_ik_E = data[:, 3].reshape(ik_len, -1)
        modes_ik_rt = data[:, 4].reshape(ik_len, -1)
        modes_ik_sr = data[:, 5].reshape(ik_len, -1)
        if not osp.exists(osp.join(os.getcwd(), 'rt_data')):
            os.makedirs(osp.join(os.getcwd(), 'rt_data'))
        self.write_it("RT", modes_ik_E[0], modes_ik_rt)
        self.write_it("SR", modes_ik_E[0], modes_ik_sr)
        self.write_modes_resolved_data(iband_len, 'RT', modes_ik_E, modes_ik_rt)
        self.write_modes_resolved_data(iband_len, 'SR', modes_ik_E, modes_ik_sr)
    def plot_relaxationTime_RT_SR(self):
        plt.figure(figsize=self.fig_sz_large)
        ax = plt.subplot(121)
        ax.spines['top'].set_linewidth(self.bwith)
        ax.spines['right'].set_linewidth(self.bwith)
        ax.spines['left'].set_linewidth(self.bwith)
        ax.spines['bottom'].set_linewidth(self.bwith)
        # RT
        data_RT = np.loadtxt(osp.join(os.getcwd(), 'rt_data', 'RT_IT.dat'), skiprows=1, dtype=float)
        for i in range(1, len(data_RT[0, :])):
            ax.scatter(data_RT[:, 0], data_RT[:, i], s=5, color='black')
        if self.log:
            ax.set_yscale('log')
        else:
            ax.set_ylim(np.min(data_RT[:, 1:]), np.max(data_RT[:, 1:]))
        ax.set_ylabel("Relaxation time (fs)", fontdict=font)
        ax.set_xlabel("Energy (eV)", fontdict=font)
        ax.set_xlim(np.min(data_RT[:, 0]), np.max(data_RT[:, 0]))
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-self.font_diff) for label in labels]
        # leg = ax.legend(frameon=True, loc='upper right')
        # labels = leg.get_texts()
        # [label.set_fontproperties(font_properties) for label in labels]
        # [label.set_fontsize(self.font_size-6) for label in labels]
        # SR
        ax1 = plt.subplot(122)
        ax1.spines['top'].set_linewidth(self.bwith)
        ax1.spines['right'].set_linewidth(self.bwith)
        ax1.spines['left'].set_linewidth(self.bwith)
        ax1.spines['bottom'].set_linewidth(self.bwith)
        data_SR = np.loadtxt(osp.join(os.getcwd(), 'rt_data', 'SR_IT.dat'), skiprows=1, dtype=float)
        for j in range(1, len(data_SR[0, :])):
            ax1.scatter(data_SR[:, 0], data_SR[:, j], s=5, color='blue')
        if self.log:
            ax1.set_yscale('log')
        else:
            ax1.set_ylim(np.min(data_SR[:, 1:]), np.max(data_SR[:, 1:]))
        ax1.set_ylabel("Scattering Rate (THz)", fontdict=font)
        ax1.set_xlabel("Energy (eV)", fontdict=font)
        ax1.set_xlim(np.min(data_SR[:, 0]), np.max(data_SR[:, 0]))
        labels = ax1.get_xticklabels() + ax1.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-self.font_diff) for label in labels]
        # leg = ax.legend(frameon=True, loc='upper right')
        # labels = leg.get_texts()
        # [label.set_fontproperties(font_properties) for label in labels]
        # [label.set_fontsize(self.font_size-6) for label in labels]
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.22, hspace=0.2)
        plt.savefig(osp.join(os.getcwd(), 'rt_data', "RT_SR.png"), dpi=300, bbox_inches='tight')
        plt.close()

    def write_velocity(self, iband_len, label, E, data):
        tmp_E = []
        tmp_data = []
        len_data = int(len(E[0])/iband_len)
        for i in range(len(E)):
            tmp_xE = np.zeros(shape=(iband_len, len_data))
            tmp_ySR = np.zeros(shape=(iband_len, len_data))
            for j in range(len_data):
                for mode in range(iband_len):
                    tmp_xE[mode, j] = E[i, j*iband_len + mode]
                    tmp_ySR[mode, j] = data[i, j*iband_len + mode]
            tmp_E.append(tmp_xE) 
            tmp_data.append(tmp_ySR)
        for ii in range(iband_len):
            wt_modes_resolved = open(osp.join(os.getcwd(), 'rt_data', str(label)+"_modes_iband_"+str(ii+1)+".dat"), 'w')
            wt_modes_resolved.write("E(ibnd)(eV)".ljust(12, ' ')+ '    ')
            for k in range(len(data)):
                wt_modes_resolved.write(("it: " + str(k+1)).rjust(12, ' ') + '    ')
            wt_modes_resolved.write('\n')
            for jj in range(len(tmp_E[0][0])):
                wt_modes_resolved.write(str(np.round(tmp_E[0][ii][jj], 8)).rjust(12, ' '))
                for jjj in range(len(tmp_E)):
                    wt_modes_resolved.write("    " + str(np.round(tmp_data[jjj][ii][jj], 8)).rjust(12, ' '))
                wt_modes_resolved.write('\n')
            wt_modes_resolved.close()
    def velocity_vel(self):
        data = np.loadtxt(osp.join(os.getcwd(), 'pwscf.vel'), skiprows=4, dtype=float)
        ik = set(data[:, 0])
        ik_len = len(ik)
        iband = set(data[:, 1])
        iband_len = len(iband)
        ikk = data[:, 0].reshape(iband_len, -1)
        kx = data[:, 3].reshape(iband_len, -1)
        ky = data[:, 4].reshape(iband_len, -1)
        kz = data[:, 5].reshape(iband_len, -1)
        vel = data[:, -1].reshape(iband_len, -1)
        modes_ik_sr = data[:, 5].reshape(iband_len, -1)
        if not osp.exists(osp.join(os.getcwd(), 'rt_data')):
            os.makedirs(osp.join(os.getcwd(), 'rt_data'))
        self.write_velocity("vel", ikk, modes_ik_sr)

    def plot_mode_RT_ST(self):
        data = np.loadtxt(osp.join(os.getcwd(), 'relaxation_time.dat'), skiprows=1, dtype=float)
        iband = set(data[:, 2])
        iband_len = len(iband)
        plt.figure(figsize=self.fig_sz_large)
        # RT
        ax = plt.subplot(121)
        ax.spines['top'].set_linewidth(self.bwith)
        ax.spines['right'].set_linewidth(self.bwith)
        ax.spines['left'].set_linewidth(self.bwith)
        ax.spines['bottom'].set_linewidth(self.bwith)
        min_x = []
        max_x = []
        min_y = []
        max_y = []
        for i in range(1, iband_len+1):
            tmp_data = np.loadtxt(osp.join(os.getcwd(), 'rt_data', 'RT_modes_iband_'+str(i)+'.dat'), skiprows=1, dtype=float)
            min_x.append(np.min(tmp_data[:, 0]))
            max_x.append(np.max(tmp_data[:, 0]))
            min_y.append(np.min(tmp_data[:, 1:]))
            max_y.append(np.max(tmp_data[:, 1:]))
            for j in range(1, len(tmp_data[0, :])):
                ax.scatter(tmp_data[:, 0], tmp_data[:, j], s=5, color=self.pointcolor[i-1])
        if self.log:
            ax.set_yscale('log')
        else:
            ax.set_ylim(min(min_y), max(max_y))
        ax.set_xlim(min(min_x), max(max_x))
        ax.set_ylabel("Relaxation time (fs)", fontdict=font)
        ax.set_xlabel("Energy (eV)", fontdict=font)
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-self.font_diff) for label in labels]
        for kk0 in range(1, iband_len+1):
            ax.scatter([], [], color=self.pointcolor[kk0-1], s=10, label="Iband: "+str(kk0))
        leg = ax.legend(frameon=True, loc='upper right')
        labels = leg.get_texts()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-6) for label in labels]
        # SR
        ax1 = plt.subplot(122)
        ax1.spines['top'].set_linewidth(self.bwith)
        ax1.spines['right'].set_linewidth(self.bwith)
        ax1.spines['left'].set_linewidth(self.bwith)
        ax1.spines['bottom'].set_linewidth(self.bwith)
        min_x1 = []
        max_x1 = []
        min_y1 = []
        max_y1 = []
        for ii in range(1, iband_len+1):
            tmp_data1 = np.loadtxt(osp.join(os.getcwd(), 'rt_data', 'SR_modes_iband_'+str(ii)+'.dat'), skiprows=1, dtype=float)
            min_x1.append(np.min(tmp_data1[:, 0]))
            max_x1.append(np.max(tmp_data1[:, 0]))
            min_y1.append(np.min(tmp_data1[:, 1:]))
            max_y1.append(np.max(tmp_data1[:, 1:]))
            for jj in range(1, len(tmp_data[0, :])):
                ax1.scatter(tmp_data1[:, 0], tmp_data1[:, jj], s=5, color=self.pointcolor[ii-1])
        if self.log:
            ax1.set_yscale('log')
        else:
            ax1.set_ylim(min(min_y1), max(max_y1))
        ax1.set_xlim(min(min_x1), max(max_x1))
        ax1.set_ylabel("Scattering Rate (THz)", fontdict=font)
        ax1.set_xlabel("Energy (eV)", fontdict=font)
        labels = ax1.get_xticklabels() + ax1.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-self.font_diff) for label in labels]
        for kk in range(1, iband_len+1):
            ax1.scatter([], [], color=self.pointcolor[kk-1], s=10, label="Iband: "+str(kk))
        leg1 = ax1.legend(frameon=True, loc='upper right')
        labels = leg1.get_texts()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-6) for label in labels]   
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.22, hspace=0.2)
        plt.savefig(osp.join(os.getcwd(), 'rt_data', "RT_SR_modes.png"), dpi=300, bbox_inches='tight')
        plt.close()
    def plot_iband_SR(self):
        data = np.loadtxt(osp.join(os.getcwd(), 'relaxation_time.dat'), skiprows=1, dtype=float)
        iband = set(data[:, 2])
        iband_len = len(iband)
        plt.figure(figsize=self.fig_sz)
        # SR
        ax1 = plt.subplot(111)
        ax1.spines['top'].set_linewidth(self.bwith)
        ax1.spines['right'].set_linewidth(self.bwith)
        ax1.spines['left'].set_linewidth(self.bwith)
        ax1.spines['bottom'].set_linewidth(self.bwith)
        min_x1 = []
        max_x1 = []
        min_y1 = []
        max_y1 = []
        for ii in range(1, iband_len+1):
            tmp_data1 = np.loadtxt(osp.join(os.getcwd(), 'rt_data', 'SR_modes_iband_'+str(ii)+'.dat'), skiprows=1, dtype=float)
            min_x1.append(np.min(tmp_data1[:, 0]))
            max_x1.append(np.max(tmp_data1[:, 0]))
            min_y1.append(np.min(tmp_data1[:, 1:]))
            max_y1.append(np.max(tmp_data1[:, 1:]))
            for jj in range(1, len(tmp_data1[0, :])):
                ax1.scatter(tmp_data1[:, 0], tmp_data1[:, jj], s=5, color=self.pointcolor[ii-1])
        if self.log:
            ax1.set_yscale('log')
        else:
            ax1.set_ylim(min(min_y1), max(max_y1))
        ax1.set_xlim(min(min_x1), max(max_x1))
        ax1.set_ylabel("Scattering Rate (THz)", fontdict=font)
        ax1.set_xlabel("Energy (eV)", fontdict=font)
        labels = ax1.get_xticklabels() + ax1.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-self.font_diff) for label in labels]
        for kk in range(1, iband_len+1):
            ax1.scatter([], [], color=self.pointcolor[kk-1], s=10, label="Iband: "+str(kk))
        leg1 = ax1.legend(frameon=True, loc='lower right')
        labels = leg1.get_texts()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-6) for label in labels]   
        plt.savefig(osp.join(os.getcwd(), 'rt_data', "SR_modes.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    def get_data_SR(self):
        data = open(osp.join(os.getcwd(), "pwscf.imsigma_mode"), 'r').readlines()
        wr_dat = open(osp.join(os.getcwd(), "pwscf.imsigma_mode_data"), 'w')
        [nk, nbnd, nT, nmodes] = [int(item.split()[0]) for item in data[3].split(':')[1:]]
        wr_dat.write(data[3])
        for i in range(len(data)):
            if not "#" in data[i]:
                wr_dat.write(data[i])
        wr_dat.close()
        return nk, nbnd, nT, nmodes
    def write_ITemp(self, label, E, data):
        wt_modes = open(osp.join(os.getcwd(), 'sr_data', str(label) + '_IT'+'.dat'), 'w')
        wt_modes.write("E(ibnd)(eV)".ljust(12, ' ')+ '    ')
        for k in range(len(data)):
            wt_modes.write(("it: " + str(k+1)).rjust(12, ' ') + '    ')
        wt_modes.write('\n')
        for i in range(len(E)):
            wt_modes.write(str(np.round(E[i], 8)).rjust(12, ' '))
            for j in range(len(data)):
                wt_modes.write("    " + str(np.round(data[j, i], 8)).rjust(12, ' '))
            wt_modes.write('\n')
        wt_modes.close()
    def write_modes_resolved_data_SR(self, nmodes, label, E, data):
        tmp_E = []
        tmp_data = []
        len_data = int(len(E[0])/nmodes)
        for i in range(len(E)):
            tmp_xE = np.zeros(shape=(nmodes, len_data))
            tmp_ySR = np.zeros(shape=(nmodes, len_data))
            for j in range(len_data):
                for mode in range(nmodes):
                    tmp_xE[mode, j] = E[i, j*nmodes + mode]
                    tmp_ySR[mode, j] = data[i, j*nmodes + mode]
            tmp_E.append(tmp_xE) 
            tmp_data.append(tmp_ySR)
        for ii in range(nmodes):
            wt_modes_resolved = open(osp.join(os.getcwd(), 'sr_data', str(label)+"_modes_"+str(ii+1)+".dat"), 'w')
            wt_modes_resolved.write("E(ibnd)(eV)".ljust(12, ' ')+ '    ')
            for k in range(len(data)):
                wt_modes_resolved.write(("it: " + str(k+1)).rjust(12, ' ') + '    ')
            wt_modes_resolved.write('\n')
            for jj in range(len(tmp_E[0][0])):
                wt_modes_resolved.write(str(np.round(tmp_E[0][ii][jj], 8)).rjust(12, ' '))
                for jjj in range(len(tmp_E)):
                    wt_modes_resolved.write("    " + str(np.round(tmp_data[jjj][ii][jj], 8)).rjust(12, ' '))
                wt_modes_resolved.write('\n')
            wt_modes_resolved.close()
    def scattering_rate(self):
        nk, nbnd, nT, nmodes = self.get_data_SR()
        print("NO.k: ", nk, " NO.bands: ", nbnd, " NO.T: ", nT, ' NO.modes: ', nmodes)
        hbar = 0.65821195 # eV
        data= np.loadtxt(osp.join(os.getcwd(), "pwscf.imsigma_mode_data"), skiprows=1, dtype=float)
        modes_ik_E = data[:, 3].reshape(nT, -1)
        modes_ik_sr = ((2/hbar*10) *(data[:, 5])).reshape(nT, -1)
        if not osp.exists(osp.join(os.getcwd(), 'sr_data')):
            os.makedirs(osp.join(os.getcwd(), 'sr_data'))
        self.write_ITemp("SR", modes_ik_E[0], modes_ik_sr)
        self.write_modes_resolved_data_SR(nmodes, 'SR', modes_ik_E, modes_ik_sr)
    def plot_mode_SR(self, fermi, eng):
        data = np.loadtxt(osp.join(os.getcwd(), 'pwscf.imsigma_mode_data'), skiprows=1, dtype=float)
        iband = set(data[:, 4])
        iband_len = len(iband)
        plt.figure(figsize=self.fig_sz)
        # SR
        ax1 = plt.subplot(111)
        ax1.spines['top'].set_linewidth(self.bwith)
        ax1.spines['right'].set_linewidth(self.bwith)
        ax1.spines['left'].set_linewidth(self.bwith)
        ax1.spines['bottom'].set_linewidth(self.bwith)
        min_x1 = []
        max_x1 = []
        min_y1 = []
        max_y1 = []
        for ii in range(1, iband_len+1):
            tmp_data1 = np.loadtxt(osp.join(os.getcwd(), 'sr_data', 'SR_modes_'+str(ii)+'.dat'), skiprows=1, dtype=float)
            min_x1.append(np.min(tmp_data1[:, 0]))
            max_x1.append(np.max(tmp_data1[:, 0]))
            min_y1.append(np.min(tmp_data1[:, 1:]))
            max_y1.append(np.max(tmp_data1[:, 1:]))
            for jj in range(1, len(tmp_data1[0, :])):
                ax1.scatter(tmp_data1[:, 0]-fermi, tmp_data1[:, jj], s=5, color=self.pointcolor[ii-1])
        if self.log:
            ax1.set_yscale('log')
        else:
            ax1.set_ylim(min(min_y1), max(max_y1))
        ax1.set_xlim(eng)
        ax1.set_ylabel("Scattering Rate (THz)", fontdict=font)
        ax1.set_xlabel("Energy (eV)", fontdict=font)
        labels = ax1.get_xticklabels() + ax1.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-self.font_diff) for label in labels]
        for kk in range(1, iband_len+1):
            ax1.scatter([], [], color=self.pointcolor[kk-1], s=10, label="modes: "+str(kk))
        leg1 = ax1.legend(frameon=True, loc='lower right')
        labels = leg1.get_texts()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-6) for label in labels]   
        plt.savefig(osp.join(os.getcwd(), 'sr_data', "SR_modes.png"), dpi=300, bbox_inches='tight')
        plt.close()

    def write_ITemp_mfp(self, label, E, data):
        wt_modes = open(osp.join(os.getcwd(), 'mfp_data', str(label) + '_IT'+'.dat'), 'w')
        wt_modes.write("E(ibnd)(eV)".ljust(12, ' ')+ '    ')
        for k in range(len(data)):
            wt_modes.write(("it: " + str(k+1)).rjust(12, ' ') + '    ')
        wt_modes.write('\n')
        for i in range(len(E)):
            wt_modes.write(str(np.round(E[i], 8)).rjust(12, ' '))
            for j in range(len(data)):
                wt_modes.write("    " + str(np.round(data[j, i], 8)).rjust(12, ' '))
            wt_modes.write('\n')
        wt_modes.close()
    def write_modes_resolved_data_MFP(self, nmodes, label, E, data):
        tmp_E = []
        tmp_data = []
        len_data = int(len(E[0])/nmodes)
        for i in range(len(E)):
            tmp_xE = np.zeros(shape=(nmodes, len_data))
            tmp_ySR = np.zeros(shape=(nmodes, len_data))
            for j in range(len_data):
                for mode in range(nmodes):
                    tmp_xE[mode, j] = E[i, j*nmodes + mode]
                    tmp_ySR[mode, j] = data[i, j*nmodes + mode]
            tmp_E.append(tmp_xE) 
            tmp_data.append(tmp_ySR)
        for ii in range(nmodes):
            wt_modes_resolved = open(osp.join(os.getcwd(), 'mfp_data', str(label)+"_modes_"+str(ii+1)+".dat"), 'w')
            wt_modes_resolved.write("E(ibnd)(eV)".ljust(12, ' ')+ '    ')
            for k in range(len(data)):
                wt_modes_resolved.write(("it: " + str(k+1)).rjust(12, ' ') + '    ')
            wt_modes_resolved.write('\n')
            for jj in range(len(tmp_E[0][0])):
                wt_modes_resolved.write(str(np.round(tmp_E[0][ii][jj], 8)).rjust(12, ' '))
                for jjj in range(len(tmp_E)):
                    wt_modes_resolved.write("    " + str(np.round(tmp_data[jjj][ii][jj], 8)).rjust(12, ' '))
                wt_modes_resolved.write('\n')
            wt_modes_resolved.close()
    def get_data_mfp(self):
        data = open(osp.join(os.getcwd(), "pwscf.mfp"), 'r').readlines()
        wr_dat = open(osp.join(os.getcwd(), "pwscf.mfp_data"), 'w')
        [nk, nbnd, nT] = [int(item.split()[0]) for item in data[3].split(':')[1:]]
        wr_dat.write(data[3])
        for i in range(len(data)):
            if not "#" in data[i]:
                wr_dat.write(data[i])
        wr_dat.close()
        return nk, nbnd, nT
    def mfp(self):
        data = np.loadtxt(osp.join(os.getcwd(), 'pwscf.mfp_data'), skiprows=1, dtype=float)
        ik = set(data[:, 0])
        ik_len = len(ik)
        iband = set(data[:, 2])
        iband_len = len(iband)
        modes_ik_E = data[:, 3].reshape(ik_len, -1)
        modes_ik_rt = data[:, 4].reshape(ik_len, -1)
        modes_ik_mfp = data[:, 5].reshape(ik_len, -1)
        if not osp.exists(osp.join(os.getcwd(), 'mfp_data')):
            os.makedirs(osp.join(os.getcwd(), 'mfp_data'))
        self.write_ITemp_mfp("RT", modes_ik_E[0], modes_ik_rt)
        self.write_ITemp_mfp("MFP", modes_ik_E[0], modes_ik_mfp)
        self.write_modes_resolved_data_MFP(iband_len, 'RT', modes_ik_E, modes_ik_rt)
        self.write_modes_resolved_data_MFP(iband_len, 'MFP', modes_ik_E, modes_ik_mfp)
    def plot_mode_MFP(self, fermi, eng, label_RTorMFP='MFP', ylabel="MFP", figname="MFP"):
        data = np.loadtxt(osp.join(os.getcwd(), 'pwscf.mfp_data'), skiprows=1, dtype=float)
        iband = set(data[:, 2])
        iband_len = len(iband)
        plt.figure(figsize=self.fig_sz)
        # SR
        ax1 = plt.subplot(111)
        ax1.spines['top'].set_linewidth(self.bwith)
        ax1.spines['right'].set_linewidth(self.bwith)
        ax1.spines['left'].set_linewidth(self.bwith)
        ax1.spines['bottom'].set_linewidth(self.bwith)
        min_x1 = []
        max_x1 = []
        min_y1 = []
        max_y1 = []
        for ii in range(1, iband_len+1):
            tmp_data1 = np.loadtxt(osp.join(os.getcwd(), 'mfp_data', label_RTorMFP + '_modes_'+str(ii)+'.dat'), skiprows=1, dtype=float)
            min_x1.append(np.min(tmp_data1[:, 0]))
            max_x1.append(np.max(tmp_data1[:, 0]))
            min_y1.append(np.min(tmp_data1[:, 1:]))
            max_y1.append(np.max(tmp_data1[:, 1:]))
            for jj in range(1, len(tmp_data1[0, :])):
                ax1.scatter(tmp_data1[:, 0]-fermi, tmp_data1[:, jj], s=5, color=self.pointcolor[ii-1])
        if self.log:
            ax1.set_yscale('log')
        else:
            ax1.set_ylim(min(min_y1), max(max_y1))
        ax1.set_xlim(eng)
        ax1.set_ylabel(ylabel, fontdict=font)
        ax1.set_xlabel("Energy (eV)", fontdict=font)
        labels = ax1.get_xticklabels() + ax1.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-self.font_diff) for label in labels]
        for kk in range(1, iband_len+1):
            ax1.scatter([], [], color=self.pointcolor[kk-1], s=10, label="modes: "+str(kk))
        leg1 = ax1.legend(frameon=True, loc='lower right')
        labels = leg1.get_texts()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-6) for label in labels]   
        plt.savefig(osp.join(os.getcwd(), 'mfp_data', figname + "_modes.png"), dpi=300, bbox_inches='tight')
        plt.close()


def Pertubo_welcome():
    print("===================================================================")
    print("================= Perturbo-Post Processing ========================")
    print("       P1) bands                      P2) bands (,PBE bands)       ")
    print("       P3) phonon                     P4) phonon (,PBE phonon)     ")
    print("       P5) imsigma                    P6) Relaxation Time          ")
    print("       P7) RT & MFP                                                ")
    print("------------------->>")
    select = input("Input : ")
    return select

def manipulate():
    sel_perturbo = Pertubo_welcome()
    if sel_perturbo == 'P1':
        eng = [float(item) for item in input("Please input energy range: ").split()]
        perturbop = perturboPost(eng)
        perturbop.plot_bands('bands')

    elif sel_perturbo == 'P2':
        eng = [float(item) for item in input("Please input energy range: ").split()]
        perturbop = perturboPost(eng)
        perturbop.plot_bands_compare('band.dat', 'MLWFs')
        sys.dont_write_bytecode = True 

    elif sel_perturbo == 'P3':
        eng = [float(item) for item in input("Please input energy range: ").split()]
        perturbop = perturboPost(eng)
        perturbop.plot_phonon_single('phonon')
        sys.dont_write_bytecode = True  

    elif sel_perturbo == 'P4':
        eng = [float(item) for item in input("Please input energy range: ").split()]
        perturbop = perturboPost(eng)
        perturbop.plot_phonon()
        sys.dont_write_bytecode = True 
        
    elif sel_perturbo == 'P5':
        import perturbo_pp as ppp
        ppp.imsigma_pp(os.getcwd())

    elif sel_perturbo == 'P6':
        import configfile.config as conf
        path_perturbo = conf.perturbo_path
        eng = [float(item) for item in input("Please input energy range: ").split()]
        perturbop = perturboPost(eng)
        os.system(osp.join(path_perturbo, '../utils/relaxation_time.py'))
        perturbop.relaxation_time()
        perturbop.plot_relaxationTime_RT_SR()
        perturbop.plot_mode_RT_ST()
        perturbop.plot_iband_SR()
        sys.dont_write_bytecode = True 

    elif sel_perturbo == 'P7':
        fermi = float(input('please input the Fermi level: '))
        eng_range = input('please input the plot energy range: ').split()
        eng = [float(item) for item in eng_range]
        perturbop = perturboPost(eng)
        perturbop.get_data_mfp()
        perturbop.mfp()
        perturbop.plot_mode_MFP(fermi, eng, 'MFP', 'MFP (nm)', 'MFP')
        perturbop.plot_mode_MFP(fermi, eng, 'RT', 'Relaxation time(fs)', 'RT')
        sys.dont_write_bytecode = True 

    else:
        print("Please input right parameters !")

if __name__ == '__main__':
    manipulate()
