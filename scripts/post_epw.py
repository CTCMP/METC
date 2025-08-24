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
import module.data_output as IO
import AutoinputQE as inpp
curPath = os.path.abspath(os.path.dirname(__file__))
fontpath = osp.join(curPath, "../font/arial.ttf")
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

class epwPost:
    def __init__(self):
        self.file_path = os.getcwd()
        self.scatterpoints_size = 8
        self.pointcolor = ['#B8860B', '#00FFFF', '#808000', '#00008B', '#006400', '#FF00FF', '#D2691E', '#DC143C', \
                        '#B8860B', '#00FFFF', '#00FF00', "#8B4513", "#FFB6C1", "#6A5ACD", "#1E90FF", "#40E0D0", "#008080", \
                        "#3CB371", "#00FF7F", "#6B8E23", "DAA520", "#DC143C", "#FF4500"] 
        self.font_size=20
        self.mev2freq = 1/(4.1356676969)
        self.lw = 2.0
        self.pl2tau = 10/6.582119
        self.cmTomeV = 1/8.06554429
        self.efermi = bas.get_fermi_level(self.file_path, 'scf.log')
        self.fig_sz = (5, 4)
        self.fig_sz_large = (10, 4)
        self.mag = bas.mag_judge()
    
    def epw_write_linewidth(self, freq, phononLinewidth, itemFile):
        ff = open(osp.join(self.file_path, "BTE.sr_phel_metal"+itemFile.split('.')[2]), 'w')
        for i in range(len(freq[0, :])):
            for j in range(len(freq)):
                ff.write(" " + str(np.round(freq[j, i] * self.mev2freq, 10)).rjust(12, ' ') + "    " + 
                         str(np.round(phononLinewidth[j, i] * self.mev2freq * 2 * np.pi, 10)).rjust(12, ' ') + '\n')
                # ff.write(" " + str(np.round(freq[j, i], 10)).rjust(12, ' ') + "    " + 
                #          str(np.round(phononLinewidth[j, i], 10)).rjust(12, ' ') + '\n')
        ff.close()

    def epw_linewidth_BTE(self):            # Transform epw linewidth file to BTE file
        file_list = os.listdir(self.file_path)
        # branches = int(input("Please input phonon branches: "))
        for itemFile in file_list:
            if 'linewidth.phself.' in itemFile:
                data = np.loadtxt(osp.join(self.file_path, itemFile), skiprows=2, dtype=float)
                branches = int(np.max(data[:, 1]))
                print("Phonon branches: ", branches)
                freq = data[:, 2].reshape(-1, branches)
                phononLinewidth = data[:, 3].reshape(-1, branches)
                self.epw_write_linewidth(freq, phononLinewidth, itemFile)

    def epw_scattering_plot(self):
        fig = plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        try:
            files = os.listdir()
            file_use = ''
            for itt in range(len(files)):
                if 'linewidth.phself' in files[itt]:
                    file_use = files[itt]
            data = np.loadtxt(osp.join(os.getcwd(), file_use), skiprows=2, dtype=float)
            ax.scatter(data[:, 2]*self.mev2freq, data[:, 3] * self.mev2freq * 2 * np.pi, s=self.scatterpoints_size, marker='o', c=self.pointcolor[0])
        except:
            print("Something Wrong !")
        ax.set_xlabel("Frequency (THz)", fontdict=font)
        ax.set_ylabel("Scattering Rate (1/ps)", fontdict=font)
        ax.set_ylim(10**-5, 10**1)
        # ax.set_ylim(0, 130)
        # ax.set_xlim(0, 10)

        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-4) for label in labels]
        ax.set_yscale("log")
        plt.savefig(osp.join(os.getcwd(), "EPW_Scattering.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
        plt.close()

    def epw_scattering_plot_elec(self):
        fig = plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        try:
            files = os.listdir()
            file_use = ''
            for itt in range(len(files)):
                if 'linewidth.elself' in files[itt]:
                    file_use = files[itt]
            data = np.loadtxt(osp.join(os.getcwd(), file_use), skiprows=2, dtype=float)
            ax.scatter(data[:, 2], data[:, 3] * self.mev2freq, s=self.scatterpoints_size, marker='o', c=self.pointcolor[0])
        except:
            print("Something Wrong !")
        ax.set_xlabel("Energy (eV)", fontdict=font)
        ax.set_ylabel("Scattering Rate (1/ps)", fontdict=font)
        # ax.set_ylim(10**-5, 10**2)
        # ax.set_ylim(0, 130)
        # ax.set_xlim(0, 10)

        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-4) for label in labels]
        # ax.set_yscale("log")
        plt.savefig(osp.join(os.getcwd(), "EPW_Scattering_ele.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
        plt.close()

    def epw_compare_kgrid(self):
        filelist = os.listdir()
        ori_path = os.getcwd()
        fig = plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        for it in range(len(filelist)):
            try:
                if 'epw_' in filelist[it]: 
                    markerlabel = filelist[it].split('_')[1]
                    os.chdir(osp.join(os.getcwd(), filelist[it], "epw_nscf"))
                    files = os.listdir()
                    file_use = ''
                    for itt in range(len(files)):
                        if 'linewidth.phself' in files[itt]:
                            file_use = files[itt]
                    data = np.loadtxt(osp.join(os.getcwd(), file_use), skiprows=2, dtype=float)
                    ax.scatter(data[:, 2]/(2*np.pi), data[:, 3] * self.mev2freq, s=self.scatterpoints_size, marker='o', c=self.pointcolor[it], label="kgrid: " + markerlabel)
                    print("End processing kgrid: ", markerlabel)
                    os.chdir(ori_path)
            except:
                os.chdir(ori_path)
                continue
            os.chdir(ori_path)
        ax.set_xlabel("Frequency (THz)", fontdict=font)
        ax.set_ylabel("Scattering Rate (1/ps)", fontdict=font)
        # ax.set_ylim(10**-3, 10)
        leg = ax.legend(frameon=True, loc='upper right').get_texts()
        [label.set_fontproperties(font_properties) for label in leg]
        [label.set_fontsize(self.font_size-8) for label in leg]

        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-4) for label in labels]
        ax.set_yscale("log")
        plt.savefig(osp.join(ori_path, "EPW_TestKgrid.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
        plt.close()

    def epw_scattering_compare(self, filepathlist, filenamelist, labellist):
        filelist = os.listdir()
        ori_path = os.getcwd()
        fig = plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        for it in range(len(filelist)):
            try:
                if 'epw_' in filelist[it]: 
                    markerlabel = filelist[it].split('_')[1]
                    os.chdir(osp.join(os.getcwd(), filelist[it], "epw_nscf"))
                    files = os.listdir()
                    file_use = ''
                    for itt in range(len(files)):
                        if 'linewidth.phself' in files[itt]:
                            file_use = files[itt]
                    data = np.loadtxt(osp.join(os.getcwd(), file_use), skiprows=2, dtype=float)
                    ax.scatter(data[:, 2]/(2*np.pi), data[:, 3] * self.mev2freq, s=self.scatterpoints_size, marker='o', c=self.pointcolor[it], label=labellist[it])
                    print("End processing kgrid: ", markerlabel)
                    os.chdir(ori_path)
            except:
                continue
            os.chdir(ori_path)
        ax.set_xlabel("Frequency (THz)", fontdict=font)
        ax.set_ylabel("Scattering Rate (1/ps)", fontdict=font)
        ax.set_ylim(10**-3, 10)
        leg = ax.legend(frameon=True, loc='upper right').get_texts()
        [label.set_fontproperties(font_properties) for label in leg]
        [label.set_fontsize(self.font_size-8) for label in leg]

        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-4) for label in labels]
        ax.set_yscale("log")
        plt.savefig(osp.join(ori_path, "linewidth_epw_kgrid.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
        plt.close()

    def band_plot(self, energyrange):
        fig = plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        ax.plot()
        linewidth = 2
        banddata = np.loadtxt(osp.join(os.getcwd(), 'band.eig.gnu'), skiprows=0, dtype=float)
        import scripts.band as bd
        if self.mag == 2:
            nbands, nks, hsp, band_k, kk, energy  = bd.band_data_pro(os.getcwd(), "band_up.dat")
            for js in range(nbands-2):
                ax.plot(band_k, energy[js*nks:(js+1)*nks, 0], color="blue", lw=linewidth)
            ax.plot(band_k, energy[(nbands-1)*nks:((nbands-1)+1)*nks, 0], color="blue", lw=linewidth, label="Spin-up: PBE")
            nbands1, nks1, hsp1, band_k1, kk1, energy1  = bd.band_data_pro(os.getcwd(), "band_dw.dat")
            for jd in range(nbands1-2):
                ax.plot(band_k1, energy1[jd*nks1:(jd+1)*nks1, 0], color="green", lw=linewidth)
            ax.plot(band_k1, energy1[(nbands1-1)*nks1:((nbands1-1)+1)*nks1, 0], color="green", lw=linewidth, label="Spin-down: PBE")
        else:
            nbands, nks, hsp, band_k, kk, energy  = bd.band_data_pro(os.getcwd(), "band.dat")
            for j in range(nbands-1):
                ax.plot(band_k, energy[j*nks:(j+1)*nks, 0], color="black", lw=linewidth)
            ax.plot(band_k, energy[(nbands-1)*nks:((nbands-1)+1)*nks, 0], color="black", lw=linewidth, label='PBE: band')
        for i in range(0,len(hsp['klabels'])):
            ax.vlines(hsp['hspdata'][i], ymin=energyrange[0], ymax=energyrange[-1],linewidth=1,linestyle="--", color='gray')
        ####### EPW plot
        for bp in range(1, len(banddata[0, :])-1):
            ax.plot(banddata[:, 0], banddata[:, bp], color='red', ls='--')
        ax.plot(banddata[:, 0], banddata[:, -1], color='red', ls='--', label='EPW')
        ####### EPW plot
        ax.set_title("BAND", fontdict=font)
        ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1], linewidth=1.0, linestyle="--", color='gray')
        ax.set_xticks(hsp['hspdata'])
        ax.set_xticklabels(hsp['klabels'])
        ax.set_ylabel("Energy (eV)", fontdict=font)
        ax.set_ylim(energyrange)
        ax.set_xlim(0, np.max(band_k))
        leg = ax.legend(frameon=True, loc='upper right')
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size) for label in labels]
        labels = leg.get_texts()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size-6) for label in labels]
        plt.savefig("band_epw.png", dpi=400, bbox_inches = 'tight', pad_inches=0.2)
        plt.close()

    def phon_plot(self, filepath, filename):
        import scripts.phonon_band as pb
        nbands, nks, hsp, band_k, kk, energy  = pb.band_data_pro(filepath, filename)
        fig = plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        phonondata = np.loadtxt(osp.join(os.getcwd(), 'phband.freq.gnu'), skiprows=0, dtype=float)
        for kk in range(0, len(energy[0, :])-1):
            ax.plot(band_k, energy[:, kk], color='black', lw=2.0)
        ax.plot(band_k, energy[:, -1], color='black', lw=2.0, label="PBE")
        max_value = np.max(energy)
        ####### EPW plot
        for bp in range(1, len(phonondata[0, :])-1):
            ax.plot(phonondata[:, 0], phonondata[:, bp], color='red', ls='--')
        ax.plot(phonondata[:, 0], phonondata[:, -1], color='red', ls='--', label='EPW')
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
        plt.savefig("phonon_epw.png", dpi=400, bbox_inches = 'tight', pad_inches=0.1)
        plt.close()

    def bandphon_compare(self):
        import module.funcs as modFunc
        modFunc.epwband_phonon(os.getcwd(), "band.eig", self.efermi)
        modFunc.epwband_phonon(os.getcwd(), "phband.freq", 0)
        energy_range = [float(item) for item in input("Please input energy for band plot: ").split()]
        self.band_plot(energy_range)
        self.phon_plot(os.getcwd(), "../phonon/pwscf.freq")

    def get_data_ph_el_lw(self, filename):
        linew = np.loadtxt(osp.join(self.file_path, filename), skiprows=2, dtype=float)
        branche_min = int(np.min(linew[:, 1]))
        branche_max = int(np.max(linew[:, 1]))
        len_branch = branche_max - branche_min + 1
        len_x = int(len(linew[:, 0])/len_branch)
        x = linew[:, 0].reshape(len_x, -1).T
        y = linew[:, -1].reshape(len_x, -1).T
        branches = [it for it in range(branche_min, branche_max+1)]
        return branches, x[0], y

    def get_hsp(self, x):
        hsp_labels = [str(items.split("!")[1]).strip() for items in inpp.kpoints_band.split("\n")[1:-1]]
        len_x = int(x/(len(hsp_labels)-1))
        hsp_data = [item*len_x for item in range(len(hsp_labels))]
        hsp = {
            'klabels' : hsp_labels,
            'hspdata' : hsp_data,
        }
        return hsp

    def epw_lambda_plot(self):
        fig = plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        files = os.listdir()
        file_use = ''
        for itt in range(len(files)):
            if 'lambda.phself' in files[itt]:
                file_use = files[itt]
        print("Lambda File Name: ", file_use)
        data_lines = open(osp.join(os.getcwd(), file_use), 'r').readlines()
        lines = 0
        for i in range(3, 7):
            if data_lines[i].split()[0].strip() == '1':
                lines = i 
                break
        print("Skip lines: ", lines)
        data = np.loadtxt(osp.join(os.getcwd(), file_use), skiprows=lines, dtype=float)
        nbands, nks, hsp, band_k, kk, energy  = bas.band_data_pro_phon(osp.join(self.file_path, '../phonon'), 'pwscf.freq')

        for k in range(len(data[0, :])-1):
            ax.plot(band_k, data[:, k+1], lw=self.lw, label="branch: " + str(k+1))
        max_value = np.max(data[:, 1:])
        try:
            for i in range(0,len(hsp['klabels'])):
                ax.vlines(hsp['hspdata'][i], ymin=0, ymax=max_value*1.5,linewidth=1.0,linestyle="--", color='gray')
            ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
            ax.set_xticks(hsp['hspdata'])
            ax.set_xticklabels(hsp['klabels'])
        except:
            print("Warning: Not specify hsp labels !")
        ax.set_ylabel("$\lambda_{qv}$", fontdict=font)
        ax.set_ylim(0, max_value*1.03)
        ax.set_xlim(0, np.max(band_k))

        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size) for label in labels]
        leg = ax.legend(loc='upper right')
        labelleg = leg.get_texts()
        [label.set_fontproperties(font_properties) for label in labelleg]
        [label.set_fontsize(self.font_size-10) for label in labelleg]
        plt.savefig(osp.join(os.getcwd(), "Lambda.png"), dpi=300, bbox_inches='tight', pad_inches=0.1)
        plt.close()

    def plot_phon_linewidth(self):
        files = os.listdir()
        filename = ''
        for i in range(len(files)):
            if 'linewidth.phself' in files[i]:
                filename = files[i]
                break
        print("Linewidth File Name: ", filename)
        try:
            nbands, nks, hsp, band_k, kk, energy  = bas.band_data_pro_phon(osp.join(self.file_path, '../phonon'), 'pwscf.freq')
            phon_linew = np.loadtxt(osp.join(self.file_path, filename), skiprows=2, dtype=float)
            phon_linew_data = phon_linew[:, -1].reshape(len(band_k), -1).T
            branche_min = int(np.min(phon_linew[:, 1]))
            branche_max = int(np.max(phon_linew[:, 1]))
            label_phonon = [it for it in range(branche_min, branche_max+1)]
        except:
            print("Use linewidth data !")
            label_phonon, band_k, phon_linew_data = self.get_data_ph_el_lw(filename)
            try:
                hsp = self.get_hsp(band_k)
            except:
                print("No high symmetry data found !")
        plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        for k in range(len(phon_linew_data)):
            ax.plot(band_k, phon_linew_data[k], lw=self.lw, label="branch: " + str(int(label_phonon[k])))
        IO.output_linewidth(self.file_path, 'PHON_LW.dat', band_k, phon_linew_data)
        print("Save phonon linewidth data to PHON_LW.dat")
        max_value = np.max(phon_linew_data)
        try:
            for i in range(0,len(hsp['klabels'])):
                ax.vlines(hsp['hspdata'][i], ymin=0, ymax=max_value*1.5,linewidth=1.0,linestyle="--", color='gray')
            ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
            ax.set_xticks(hsp['hspdata'])
            ax.set_xticklabels(hsp['klabels'])
        except:
            print("Warning: Not specify hsp labels !")
        ax.set_ylabel("$\gamma_{qv}$ (meV)", fontdict=font)
        ax.set_title('Phonon linewidth')
        if max_value >= 50:
            ax.set_ylim(0, 50)
        else:
            ax.set_ylim(0, max_value*1.03)
        ax.set_xlim(0, np.max(band_k))
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size) for label in labels]
        leg = ax.legend(loc='upper right')
        labelleg = leg.get_texts()
        [label.set_fontproperties(font_properties) for label in labelleg]
        [label.set_fontsize(self.font_size-10) for label in labelleg]
        plt.savefig("PhonLineWidth.png", dpi=300, bbox_inches = 'tight', pad_inches=0.1)
        plt.close()

    def plot_elec_linewidth(self):
        files = os.listdir()
        filename = ''
        for i in range(len(files)):
            if 'linewidth.elself' in files[i]:
                filename = files[i]
                break
        print("Linewidth File Name: ", filename)
        if osp.exists(osp.join(self.file_path, 'band.dat')):
            nbands, nks, hsp, band_k, kk, energy = bas.band_data_pro_band(self.file_path, 'band.dat')
            elec_linew = np.loadtxt(osp.join(self.file_path, filename), skiprows=2, dtype=float)
            elec_linew_data = elec_linew[:, -1].reshape(len(band_k), -1).T
            branche_min = int(np.min(elec_linew[:, 1]))
            branche_max = int(np.max(elec_linew[:, 1]))
            label_band = [it for it in range(branche_min, branche_max+1)]
        else:
            label_band, band_k, elec_linew_data = self.get_data_ph_el_lw(filename)
            try:
                hsp = self.get_hsp(band_k)
            except:
                print("No high symmetry data found !")
        plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        for k in range(len(elec_linew_data)):
            ax.plot(band_k, elec_linew_data[k], lw=self.lw, label='band index: ' + str(int(label_band[k])))
        IO.output_linewidth(self.file_path, 'ELEC_LW.dat', band_k, elec_linew_data)
        print("Save electron linewidth data to ELEC_LW.dat")
        max_value = np.max(elec_linew_data)
        try:
            for i in range(0,len(hsp['klabels'])):
                ax.vlines(hsp['hspdata'][i], ymin=0, ymax=max_value*1.5,linewidth=1.0,linestyle="--", color='gray')
            ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
            ax.set_xticks(hsp['hspdata'])
            ax.set_xticklabels(hsp['klabels'])
        except:
            print("Warning: Not specify hsp labels !")
        ax.set_ylabel("Im$\Sigma$ (meV)", fontdict=font)
        ax.set_title('Electron linewidth')
        ax.set_ylim(0, max_value*1.03)
        ax.set_xlim(0, np.max(band_k))
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size) for label in labels]
        leg = ax.legend(loc='upper right')
        labelleg = leg.get_texts()
        [label.set_fontproperties(font_properties) for label in labelleg]
        [label.set_fontsize(self.font_size-10) for label in labelleg]
        plt.savefig("ElecLineWidth.png", dpi=300, bbox_inches = 'tight', pad_inches=0.1)
        plt.close()

    def phon_linewidth_weight(self):
        files = os.listdir()
        filename = ''
        for i in range(len(files)):
            if 'linewidth.phself' in files[i]:
                filename = files[i]
                break
        print("Filename: ", filename)
        nbands, nks, hsp, band_k, kk, energy  = bas.band_data_pro_phon(osp.join(self.file_path, '../phonon'), 'pwscf.freq')
        phon_linew = np.loadtxt(osp.join(self.file_path, filename), skiprows=2, dtype=float)
        phon_linew_data = phon_linew[:, -1].reshape(len(band_k), -1)
        label_phonon = phon_linew[:list(phon_linew_data.shape)[1], 1]
        plt.figure(figsize=self.fig_sz)
        ax = plt.subplot(111)
        for k in range(len(phon_linew_data[0])):
            ax.plot(band_k, phon_linew_data[:, k], lw=self.lw, label="ph branch: " + str(int(label_phonon[k])))
        max_value = np.max(phon_linew_data)
        for i in range(0,len(hsp['klabels'])):
            ax.vlines(hsp['hspdata'][i], ymin=0, ymax=max_value*1.5,linewidth=1.0,linestyle="--", color='gray')
        ax.axhline(0, xmin=0, xmax=hsp['hspdata'][-1],linewidth=1.0,linestyle="--", color='gray')
        ax.set_xticks(hsp['hspdata'])
        ax.set_xticklabels(hsp['klabels'])
        ax.set_ylabel("$\gamma_{qv}$ (meV)", fontdict=font)
        ax.set_title('Phonon linewidth')
        ax.set_ylim(0, max_value*1.03)
        ax.set_xlim(0, np.max(band_k))
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size) for label in labels]
        leg = ax.legend(loc='upper right')
        labelleg = leg.get_texts()
        [label.set_fontproperties(font_properties) for label in labelleg]
        [label.set_fontsize(self.font_size-10) for label in labelleg]
        plt.savefig("PhonLineWidth.png", dpi=300, bbox_inches = 'tight', pad_inches=0.1)
        plt.close()

    def plot_super_SC_gap(self, temp):
        plt.figure(figsize=self.fig_sz_large)
        ######### imag
        ax = plt.subplot(121)
        ######### real imag
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size) for label in labels]
        leg = ax.legend(loc='upper right')
        labelleg = leg.get_texts()
        [label.set_fontproperties(font_properties) for label in labelleg]
        [label.set_fontsize(self.font_size-10) for label in labelleg]
        plt.savefig("SC_gap_"+ str(temp) + ".png", dpi=300, bbox_inches = 'tight', pad_inches=0.1)
        plt.close

    def super_conductivity_SC_gap(self):
        files = os.listdir()
        filename_imag = []
        filename_pade = []
        filename_acon = []
        for i in range(len(files)):
            if 'imag_iso' in files[i]:
                filename_imag.append(files[i])
            elif 'pade_iso' in files[i]:
                filename_pade.append(files[i])
            elif 'acon_iso' in files[i]:
                filename_acon.append(files[i])
            else:
                continue

def EPW_welcome():
    print("===================================================================")
    print("=================== EPW-Post Processing ===========================")
    print("======================== checking =================================")
    print("       C1) Check band phonon          C2) plot decay matrix        ")
    print("===================== Phonon related ==============================")
    print("       P1) phonon linewidth           P2) Phonon scattering        ")
    print("       P3) phon lw weighted phonon    P4) phonon spectra fucn      ")
    print("       P5) phonon specfun             P6) renormal phonon          ")
    print("       P7) Lambda                                                  ")
    print("==================== Electron related =============================")
    print("       E1) electron linewidth         E2) electron scattering      ")
    print("       E3) electron specfun           E4) a2f                      ")
    print("================ Superconductivity related ========================")
    print("       S1) iso                        S2)  aniso                   ")
    print("====================== Other revelant =============================")
    print("       T1) EPW2ShengBTE               T2) comp DFPT and EPW (|g|)  ")
    print("       T3) EPW renormal phonon                                     ")
    print("------------------->>")
    select = input("Input : ")
    return select

def manipulate():
    sel_EPW = EPW_welcome()
    if sel_EPW == 'C1':
        epwp = epwPost()
        epwp.bandphon_compare()

    elif sel_EPW == 'C2':
        print(os.getcwd())
        bas.epmat_decay(os.getcwd(), 'decay.H', ['R_e', 'H_{nm}(R_e)'])
        bas.epmat_decay(os.getcwd(), 'decay.dynmat', ['R_p', 'D_{nm}(R_p)'])
        bas.epmat_decay(os.getcwd(), 'decay.epmate', ['R_e', 'g_{nm}(R_e)'])
        bas.epmat_decay(os.getcwd(), 'decay.epmatp', ['R_p', 'g_{nm}(R_p)']) 

    elif sel_EPW == 'P1':
        epwp = epwPost()
        epwp.plot_phon_linewidth() 

    elif sel_EPW == 'P2':
        epwp = epwPost()
        epwp.epw_scattering_plot()

    elif sel_EPW == 'P3':
        import scripts.epw_phonon_linewidth_weight as epwImag
        epwImag.manipulate()

    elif sel_EPW == 'P4':
        os.system("grep 'Omega' epw_phon_spectra.out | awk -F ' ' '{print $7, $10}' > omega.dat")

    elif sel_EPW == 'P5':
        import scripts.epwPhon_renormal as eprenormal
        eprenormal.specfun_phon()
        
    elif sel_EPW == 'P6':
        import scripts.epwPhon_renormal as eprenormal
        eprenormal.specfun_phon_renormal()

    elif sel_EPW == 'P7':
        epwp = epwPost()
        epwp.epw_lambda_plot()

    elif sel_EPW == 'E1':
        epwp = epwPost()
        epwp.plot_elec_linewidth()

    elif sel_EPW == 'E2':
        epwp = epwPost()
        epwp.epw_scattering_plot_elec()

    # elif sel_EPW == 'E3':
    #     specfun_elec()

    elif sel_EPW == 'E4':
        import scripts.epw_superconductivity as epws
        epws.plot_a2f(os.getcwd(), 'pwscf.a2f')

    elif sel_EPW == 'PK':
        epwp = epwPost()
        epwp.epw_compare_kgrid()

    elif sel_EPW == 'PE':
        filepathlist = input("Input filepath list: ").split()
        filenamelist = input("Input filename list: ").split()
        labelist = input("Input labelist list: ").split()
        epwp = epwPost()
        epwp.epw_scattering_compare(filepathlist, filenamelist, labelist)

    elif sel_EPW == 'S1':
        import scripts.epw_superconductivity as epws
        epws.epw_superconductivity_iso_manipulate()

    elif sel_EPW == 'S2':
        import scripts.epw_superconductivity as epws
        epws.epw_superconductivity_aniso_manipulate()

    elif sel_EPW == 'T1':
        epwp = epwPost()
        epwp.epw_linewidth_BTE()

    elif sel_EPW == 'T2':
        import scripts.epwPhon_renormal as eprenormal
        eprenormal.comp_dfpt_epw_g(os.getcwd(), 'dfpt.dat', 'epw.dat')

    elif sel_EPW == 'T3':
        import scripts.epw_real as epwr
        epwr.manipulate()

    else:
        print("Please input right parameters !")

if __name__ == '__main__':
    manipulate()
