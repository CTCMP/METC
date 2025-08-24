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
import matplotlib as mpl
from matplotlib import font_manager
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
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
		'size': 20, 
        }

def get_qpoints(filepath, filename="BTE.qpoints"):
    data = np.loadtxt(osp.join(filepath, filename), skiprows=0, dtype=float)
    return len(data[:, 0])

color_phonon = ['black', '#B8860B', '#00FFFF', '#00008B', '#006400', '#FF00FF', '#D2691E', '#DC143C', \
                        '#B8860B', '#00FFFF', '#00FF00', "#8B4513", "#FFB6C1", "#6A5ACD", "#1E90FF", "#40E0D0", "#008080", \
                        "#3CB371", "#00FF7F", "#6B8E23", "#DAA520", "#DC143C", "#FF4500", '#808000'] 
cor_ph = []
for pp in range(60):
    if pp < len(color_phonon):
        cor_ph.append(color_phonon[pp])
    else:
        cor_ph.append('red')
ori_path = os.getcwd()
qpts = get_qpoints(ori_path, "BTE.qpoints") 
print(qpts)

font_sz = 20
fig_sz =  (5, 4)

def plot_thermal_conductivity(filename, Figname):
    data = np.loadtxt(osp.join(os.getcwd(), filename), skiprows=0, dtype=float)
    plt.figure(figsize=fig_sz)
    ax = plt.subplot(111)
    ax.plot(data[:, 0], data[:, 1], '-o', label='a-axis', color='red', lw=2)
    ax.plot(data[:, 0], data[:, 5], '-o', label='b-axis', color='green', lw=2)
    ax.plot(data[:, 0], data[:, 9], '-o', label='c-axis', color='blue', lw=2)
    max_y = []
    min_y = []
    for ii in [1, 5, 9]:
        max_y.append(np.max(data[:, ii]))
        min_y.append(np.min(data[:, ii]))
    ax.set_xlabel("Temperature (K)", fontdict=font)
    ax.set_ylabel("$\kappa$ (W/mK)", fontdict=font)
    ax.set_xlim(np.min(data[:, 0]), np.max(data[:, 0]))
    ax.set_ylim(min(min_y)*0.95, max(max_y) * 1.1)
    leg = ax.legend(frameon=False)
    labels = ax.get_yticklabels() + ax.get_xticklabels() + leg.get_texts()
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(font_sz - 4) for label in labels]
    plt.savefig(osp.join(os.getcwd(), Figname + ".png"), dpi=300, bbox_inches='tight')
    plt.close()

class shengbte:
    def __init__(self):
        self.filepath = os.getcwd()
        self.alldirs = [item for item in os.listdir() if os.path.isdir(item) and "T" in item]
        self.fig_size = fig_sz
        self.scatterpoints_size = 5
        self.font_size = font_sz
        self.pointcolor = '#008B8B'
        self.pcolor = cor_ph
        font['size'] = self.font_size
        self.num_qpoints = qpts
        self.bwith = 1.5

    def write_data(self, x, y, fpath, fname, marklabel):
        file = osp.join(fpath, fname)
        conf = open(file, 'w')
        conf.write(marklabel + '\n')
        for i in range(len(x)):
            conf.write(str(np.round(x[i], 6)).rjust(10, " ") + "  " + str(np.round(y[i], 6)).rjust(10, " ") + '\n')
        conf.close()
    def write_data_branches(self, x, y, fpath, fname, marklabel):
        file = osp.join(fpath, fname)
        conf = open(file, 'w')
        conf.write(marklabel + '\n')
        for mk in range(len(x)):
            conf.write(str("qpoints: "+str(mk+1)).rjust(12, " ") + "  " + str("branches: " + str(mk+1)).rjust(12, " "))
        conf.write("\n")
        for i in range(len(x[0])):
            for j in range(len(x)):
                conf.write(str(np.round(x[j, i], 6)).rjust(12, " ") + "  " + str(np.round(y[j, i], 6)).rjust(12, " "))
            conf.write('\n')
        conf.close()

    def plot_anharmonic_branches(self):
        data = np.loadtxt(osp.join(os.getcwd(), "Anharmonic_branches.dat"), skiprows=2, dtype=float)
        plt.figure(figsize=self.fig_size)
        ax = plt.subplot(111)
        ax.spines['top'].set_linewidth(self.bwith)
        ax.spines['bottom'].set_linewidth(self.bwith)
        ax.spines['left'].set_linewidth(self.bwith)
        ax.spines['right'].set_linewidth(self.bwith)
        x_max = []
        for i in range(int(len(data[0, :])/2)):
            ax.scatter(data[:, 2*i], data[:, 2*i+1], s=self.scatterpoints_size, marker='o', color='none', \
                       edgecolors=cor_ph[i], label="mode " + str(i+1), alpha=0.6)
            x_max.append(np.max(data[:, 2*i]))
        ax.set_xlabel("Frequency (THz)", fontdict=font)
        ax.set_ylabel("P3 scattering rate (1/ps)", fontdict=font)
        ax.set_xlim(0, max(x_max))
        leg = ax.legend(frameon=True, loc='upper right').get_texts()
        [label.set_fontproperties(font_properties) for label in leg]
        [label.set_fontsize(self.font_size-10) for label in leg]
        labels = ax.get_yticklabels() + ax.get_xticklabels()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size - 6) for label in labels]
        plt.savefig(osp.join(os.getcwd(), "anharmonic_branches.png"), dpi=300, bbox_inches='tight')
        plt.close()

    def plot_anharmonic_branches_4ph(self):
        data = np.loadtxt(osp.join(os.getcwd(), "Anharmonic_branches_4ph.dat"), skiprows=2, dtype=float)
        plt.figure(figsize=self.fig_size)
        ax = plt.subplot(111)
        ax.spines['top'].set_linewidth(self.bwith)
        ax.spines['bottom'].set_linewidth(self.bwith)
        ax.spines['left'].set_linewidth(self.bwith)
        ax.spines['right'].set_linewidth(self.bwith)
        x_max = []
        for i in range(int(len(data[0, :])/2)):
            ax.scatter(data[:, 2*i], data[:, 2*i+1], s=self.scatterpoints_size, marker='o', color='none', \
                       edgecolors=cor_ph[i], label="mode " + str(i+1), alpha=0.6)
            x_max.append(np.max(data[:, 2*i]))
        ax.set_xlabel("Frequency (THz)", fontdict=font)
        ax.set_ylabel("P4 scattering rate (1/ps)", fontdict=font)
        ax.set_xlim(0, max(x_max))
        leg = ax.legend(frameon=True, loc='upper right').get_texts()
        [label.set_fontproperties(font_properties) for label in leg]
        [label.set_fontsize(self.font_size-10) for label in leg]
        labels = ax.get_yticklabels() + ax.get_xticklabels()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size - 6) for label in labels]
        plt.savefig(osp.join(os.getcwd(), "anharmonic_branches_4ph.png"), dpi=300, bbox_inches='tight')
        plt.close()

    def plotT_BTE_w_anharmonic(self):
        for items in self.alldirs:
            os.chdir(osp.join(self.filepath, items))
            print("Filepath: " + os.getcwd())
            try:
                print("Begin process " + str(items) + " catalogue !")
                ff = 'BTE.w_anharmonic'
                if osp.exists(osp.join(os.getcwd(), 'BTE.w_3ph')):
                    ff = 'BTE.w_3ph'
                data = np.loadtxt(osp.join(os.getcwd(), ff), skiprows=0, dtype=float)
                plt.figure(figsize=self.fig_size)
                ax = plt.subplot(111)
                ax.spines['top'].set_linewidth(self.bwith)
                ax.spines['bottom'].set_linewidth(self.bwith)
                ax.spines['left'].set_linewidth(self.bwith)
                ax.spines['right'].set_linewidth(self.bwith)
                tmp_x = data[:,0]/(2 * np.pi)
                tmp_y = data[:, 1]
                ax.scatter(tmp_x, tmp_y, s=self.scatterpoints_size, marker='o', c='blue')
                self.write_data(tmp_x, tmp_y, os.getcwd(), 'Anharmonic.dat', 'BTE.w_anharmonic')
                anhar_branches_q = tmp_x.reshape(int(len(data[:, 0])/self.num_qpoints), self.num_qpoints)
                anhar_branches_scatter = tmp_y.reshape(int(len(data[:, 0])/self.num_qpoints), self.num_qpoints)
                self.write_data_branches(anhar_branches_q, anhar_branches_scatter, os.getcwd(), 
                                            "Anharmonic_branches.dat", "BTE.w_anharmonic phonon branches")
                self.plot_anharmonic_branches()
                # try:
                #     self.plot_anharmonic_branches()
                # except:
                #     print("Not successfully plot anharmonic_branches, you can plot it by use origin !")
                plt.yscale("log")
                # ax.set_xlim(0, np.max(data[:, 0]/(2 * np.pi))*1.05)
                # ax.set_ylim(0, np.max(data[:, 1])*2)
                ax.set_xlabel("Frequency (THz)", fontdict=font)
                ax.set_ylabel("P3 scattering rate (1/ps)", fontdict=font)
                labels = ax.get_yticklabels() + ax.get_xticklabels()
                [label.set_fontproperties(font_properties) for label in labels]
                [label.set_fontsize(self.font_size - 6) for label in labels]
                plt.savefig(osp.join(os.getcwd(), "anharmonic.png"), dpi=300, bbox_inches='tight')
                plt.close()
                os.chdir(ori_path)
            except:
                os.chdir(ori_path)
                print("No BTE.w_anharmonic data found !")
            os.chdir(ori_path)

    def plotT_BTE_w_anharmonic_4th(self):
        for items in self.alldirs:
            os.chdir(osp.join(self.filepath, items))
            try:
                # print("Begin process " + str(items) + " catalogue !")
                ff = 'BTE.w_4ph'
                data = np.loadtxt(osp.join(os.getcwd(), ff), skiprows=0, dtype=float)
                plt.figure(figsize=self.fig_size)
                ax = plt.subplot(111)
                ax.spines['top'].set_linewidth(self.bwith)
                ax.spines['bottom'].set_linewidth(self.bwith)
                ax.spines['left'].set_linewidth(self.bwith)
                ax.spines['right'].set_linewidth(self.bwith)
                tmp_x = data[:,0]/(2 * np.pi)
                tmp_y = data[:, 1]
                ax.scatter(tmp_x, tmp_y, s=self.scatterpoints_size, marker='o', c='blue')
                self.write_data(tmp_x, tmp_y, os.getcwd(), 'Anharmonic_4ph.dat', 'BTE.w_4ph')
                anhar_branches_q = tmp_x.reshape(int(len(data[:, 0])/self.num_qpoints), self.num_qpoints)
                anhar_branches_scatter = tmp_y.reshape(int(len(data[:, 0])/self.num_qpoints), self.num_qpoints)
                self.write_data_branches(anhar_branches_q, anhar_branches_scatter, os.getcwd(), 
                                            "Anharmonic_branches_4ph.dat", "BTE.w_anharmonic phonon branches")
                self.plot_anharmonic_branches_4ph()
                plt.yscale("log")
                ax.set_xlabel("Frequency (THz)", fontdict=font)
                ax.set_ylabel("P4 scattering rate (1/ps)", fontdict=font)
                labels = ax.get_yticklabels() + ax.get_xticklabels()
                [label.set_fontproperties(font_properties) for label in labels]
                [label.set_fontsize(self.font_size - 6) for label in labels]
                plt.savefig(osp.join(os.getcwd(), "anharmonic_4ph.png"), dpi=300, bbox_inches='tight')
                plt.close()
                os.chdir(ori_path)
            except:
                os.chdir(ori_path)
                print("No BTE.w_4ph data found !")
            os.chdir(ori_path)

    def plotT_BTE_w_anharmonic_nolog_4ph(self):
        for items in self.alldirs:
            os.chdir(osp.join(self.filepath, items))
            try:
                # print("Begin process " + str(items) + " catalogue !")
                ff = 'BTE.w_4ph'
                data = np.loadtxt(osp.join(os.getcwd(), ff), skiprows=0, dtype=float)
                plt.figure(figsize=self.fig_size)
                ax = plt.subplot(111)
                ax.spines['top'].set_linewidth(self.bwith)
                ax.spines['bottom'].set_linewidth(self.bwith)
                ax.spines['left'].set_linewidth(self.bwith)
                ax.spines['right'].set_linewidth(self.bwith)
                ax.scatter(data[:, 0]/(2 * np.pi), data[:, 1], s=self.scatterpoints_size, marker='o', c=self.pointcolor)
                self.write_data(data[:, 0]/(2 * np.pi), data[:, 1], os.getcwd(), 'Anharmonic_4ph.dat', 'w_anharmonic')
                # plt.yscale("log")
                # ax.set_xlim(0, np.max(data[:, 0]/(2 * np.pi))*1.05)
                # ax.set_ylim(0, np.max(data[:, 1])*2)
                ax.set_xlabel("Frequency (THz)", fontdict=font)
                ax.set_ylabel("P4 scattering rate (1/ps)", fontdict=font)
                labels = ax.get_yticklabels() + ax.get_xticklabels()
                [label.set_fontproperties(font_properties) for label in labels]
                [label.set_fontsize(self.font_size - 6) for label in labels]
                plt.savefig(osp.join(os.getcwd(), "anharmonic_notlogy_4ph.png"), dpi=300, bbox_inches='tight')
                plt.close()
                os.chdir(self.filepath)
            except:
                os.chdir(self.filepath)
                print("No BTE.w_4ph data found !")
            os.chdir(self.filepath)

    def plotT_BTE_w_anharmonic_nolog(self):
        for items in self.alldirs:
            os.chdir(osp.join(self.filepath, items))
            try:
                # print("Begin process " + str(items) + " catalogue !")
                ff = 'BTE.w_anharmonic'
                if osp.exists(osp.join(os.getcwd(), 'BTE.w_3ph')):
                    ff = 'BTE.w_3ph'
                data = np.loadtxt(osp.join(os.getcwd(), ff), skiprows=0, dtype=float)
                plt.figure(figsize=self.fig_size)
                ax = plt.subplot(111)
                ax.spines['top'].set_linewidth(self.bwith)
                ax.spines['bottom'].set_linewidth(self.bwith)
                ax.spines['left'].set_linewidth(self.bwith)
                ax.spines['right'].set_linewidth(self.bwith)
                ax.scatter(data[:, 0]/(2 * np.pi), data[:, 1], s=self.scatterpoints_size, marker='o', c=self.pointcolor)
                self.write_data(data[:, 0]/(2 * np.pi), data[:, 1], os.getcwd(), 'Anharmonic.dat', 'w_anharmonic')
                # plt.yscale("log")
                # ax.set_xlim(0, np.max(data[:, 0]/(2 * np.pi))*1.05)
                # ax.set_ylim(0, np.max(data[:, 1])*2)
                ax.set_xlabel("Frequency (THz)", fontdict=font)
                ax.set_ylabel("P3 scattering rate (1/ps)", fontdict=font)
                labels = ax.get_yticklabels() + ax.get_xticklabels()
                [label.set_fontproperties(font_properties) for label in labels]
                [label.set_fontsize(self.font_size - 6) for label in labels]
                plt.savefig(osp.join(os.getcwd(), "anharmonic_notlogy.png"), dpi=300, bbox_inches='tight')
                plt.close()
                os.chdir(self.filepath)
            except:
                os.chdir(self.filepath)
                print("No BTE.anharmonic data found !")
            os.chdir(self.filepath)

    def plotT_BTEWP4(self):
        for items in self.alldirs:
            os.chdir(osp.join(self.filepath, items))
            try:
                # print("Begin process " + str(items) + " catalogue !")
                data = np.loadtxt(osp.join(os.getcwd(), 'BTE.WP4'), skiprows=0, dtype=float)
                plt.figure(figsize=self.fig_size)
                ax = plt.subplot(111)
                ax.spines['top'].set_linewidth(self.bwith)
                ax.spines['bottom'].set_linewidth(self.bwith)
                ax.spines['left'].set_linewidth(self.bwith)
                ax.spines['right'].set_linewidth(self.bwith)
                ax.scatter(data[:, 0]/(2 * np.pi), data[:, 1], s=self.scatterpoints_size, marker='o', c=self.pointcolor)
                self.write_data(data[:, 0]/(2 * np.pi), data[:, 1], os.getcwd(), 'WP4.dat', 'BTE.WP4')
                # ax.scatter(data[:, 0], data[:, 1], s=self.scatterpoints_size, marker='o', c=self.pointcolor)
                # self.write_data(data[:, 0], data[:, 1], os.getcwd(), 'WP3.dat', 'BTE.WP3')
                plt.yscale("log")
                # ax.set_xlim(0, np.max(data[:, 0]/(2 * np.pi))*1.05)
                # ax.set_ylim(np.min(data[:, 1]), np.max(data[:, 1])*2)
                ax.set_xlabel("Frequency (THz)", fontdict=font)
                ax.set_ylabel("P4 Phase space (ps$^{4}$/rad$^{4}$)", fontdict=font)
                labels = ax.get_yticklabels() + ax.get_xticklabels()
                [label.set_fontproperties(font_properties) for label in labels]
                [label.set_fontsize(self.font_size - 6) for label in labels]
                plt.savefig(osp.join(os.getcwd(), "w_4ph.png"), dpi=300, bbox_inches='tight')
                plt.close()
                os.chdir(self.filepath)
            except:
                print("No BTE.WP4 data found !")
                os.chdir(self.filepath)
            os.chdir(self.filepath)

    def plotT_BTEWP3(self):
        for items in self.alldirs:
            os.chdir(osp.join(self.filepath, items))
            try:
                # print("Begin process " + str(items) + " catalogue !")
                data = np.loadtxt(osp.join(os.getcwd(), 'BTE.WP3'), skiprows=0, dtype=float)
                plt.figure(figsize=self.fig_size)
                ax = plt.subplot(111)
                ax.spines['top'].set_linewidth(self.bwith)
                ax.spines['bottom'].set_linewidth(self.bwith)
                ax.spines['left'].set_linewidth(self.bwith)
                ax.spines['right'].set_linewidth(self.bwith)
                ax.scatter(data[:, 0]/(2 * np.pi), data[:, 1], s=self.scatterpoints_size, marker='o', c=self.pointcolor)
                self.write_data(data[:, 0]/(2 * np.pi), data[:, 1], os.getcwd(), 'WP3.dat', 'BTE.WP3')
                # ax.scatter(data[:, 0], data[:, 1], s=self.scatterpoints_size, marker='o', c=self.pointcolor)
                # self.write_data(data[:, 0], data[:, 1], os.getcwd(), 'WP3.dat', 'BTE.WP3')
                plt.yscale("log")
                # ax.set_xlim(0, np.max(data[:, 0]/(2 * np.pi))*1.05)
                # ax.set_ylim(np.min(data[:, 1]), np.max(data[:, 1])*2)
                ax.set_xlabel("Frequency (THz)", fontdict=font)
                ax.set_ylabel("Phase space (ps$^{4}$/rad$^{4}$)", fontdict=font)
                labels = ax.get_yticklabels() + ax.get_xticklabels()
                [label.set_fontproperties(font_properties) for label in labels]
                [label.set_fontsize(self.font_size - 6) for label in labels]
                plt.savefig(osp.join(os.getcwd(), "w_3ph.png"), dpi=300, bbox_inches='tight')
                plt.close()
                os.chdir(self.filepath)
            except:
                print("No BTE.w_3ph data found !")
                os.chdir(self.filepath)
            os.chdir(self.filepath)

    def plotT_BTE_NU_4ph(self):
        for items in self.alldirs:
            os.chdir(osp.join(self.filepath, items))
            try:
                # print("Begin process " + str(items) + " catalogue !")
                data_N = np.loadtxt(osp.join(os.getcwd(), 'BTE.w_4ph_normal'), skiprows=0, dtype=float)
                data_U = np.loadtxt(osp.join(os.getcwd(), 'BTE.w_4ph_Umklapp'), skiprows=0, dtype=float)
                plt.figure(figsize=self.fig_size)
                ax = plt.subplot(111)
                ax.spines['top'].set_linewidth(self.bwith)
                ax.spines['bottom'].set_linewidth(self.bwith)
                ax.spines['left'].set_linewidth(self.bwith)
                ax.spines['right'].set_linewidth(self.bwith)
                ax.scatter(data_N[:, 0]/(2 * np.pi), data_N[:, 3], s=self.scatterpoints_size, marker='s', c='red', label='Normal')
                ax.scatter(data_U[:, 0]/(2 * np.pi), data_U[:, 3], s=self.scatterpoints_size, marker='o', c='black', label="Umklapp")
                plt.yscale("log")
                ax.set_xlabel("Frequency (THz)", fontdict=font)
                ax.set_ylabel("P4 Scattering rate (ps$^{-1}$)", fontdict=font)

                leg = ax.legend(frameon=False, loc='lower right').get_texts()
                [label.set_fontproperties(font_properties) for label in leg]
                [label.set_fontsize(self.font_size - 6) for label in leg]

                labels = ax.get_yticklabels() + ax.get_xticklabels()
                [label.set_fontproperties(font_properties) for label in labels]
                [label.set_fontsize(self.font_size - 6) for label in labels]
                plt.savefig(osp.join(os.getcwd(), "Normal_Umklapp_4ph.png"), dpi=300, bbox_inches='tight')
                plt.close()
                os.chdir(self.filepath)
            except:
                print("No normal or umklapp data found for 4ph, go on !")
                os.chdir(self.filepath)
            os.chdir(self.filepath)

    def plotT_BTE_NU(self):
        for items in self.alldirs:
            os.chdir(osp.join(self.filepath, items))
            try:
                # print("Begin process " + str(items) + " catalogue !")
                data_N = np.loadtxt(osp.join(os.getcwd(), 'BTE.w_anharmonic_normal'), skiprows=0, dtype=float)
                data_U = np.loadtxt(osp.join(os.getcwd(), 'BTE.w_anharmonic_Umklapp'), skiprows=0, dtype=float)
                plt.figure(figsize=self.fig_size)
                ax = plt.subplot(111)
                ax.spines['top'].set_linewidth(self.bwith)
                ax.spines['bottom'].set_linewidth(self.bwith)
                ax.spines['left'].set_linewidth(self.bwith)
                ax.spines['right'].set_linewidth(self.bwith)
                ax.scatter(data_N[:, 0]/(2 * np.pi), data_N[:, 3], s=self.scatterpoints_size, marker='s', c='red', label='Normal')
                ax.scatter(data_U[:, 0]/(2 * np.pi), data_U[:, 3], s=self.scatterpoints_size, marker='o', c='black', label="Umklapp")
                plt.yscale("log")
                ax.set_xlabel("Frequency (THz)", fontdict=font)
                ax.set_ylabel("Scattering rate (ps$^{-1}$)", fontdict=font)

                leg = ax.legend(frameon=False, loc='lower right').get_texts()
                [label.set_fontproperties(font_properties) for label in leg]
                [label.set_fontsize(self.font_size - 6) for label in leg]

                labels = ax.get_yticklabels() + ax.get_xticklabels()
                [label.set_fontproperties(font_properties) for label in labels]
                [label.set_fontsize(self.font_size - 6) for label in labels]
                plt.savefig(osp.join(os.getcwd(), "Normal_Umklapp.png"), dpi=300, bbox_inches='tight')
                plt.close()
                os.chdir(self.filepath)
            except:
                print("No normal or umklapp data found for 3ph, go on !")
                os.chdir(self.filepath)
            os.chdir(self.filepath)

    def plotT_BTEv(self):
        for items in self.alldirs:
            os.chdir(osp.join(self.filepath, items))
            try:
                # print("Begin process " + str(items) + " catalogue !")
                ff = 'BTE.w_anharmonic'
                if osp.exists(osp.join(os.getcwd(), 'BTE.w_3ph')):
                    ff = 'BTE.w_3ph'
                data_x = np.loadtxt(osp.join(os.getcwd(), ff), skiprows=0, dtype=float)
                data = np.loadtxt(osp.join(os.getcwd(), '../BTE.v'), skiprows=0, dtype=float)
                plt.figure(figsize=self.fig_size)
                ax = plt.subplot(111)
                ax.spines['top'].set_linewidth(self.bwith)
                ax.spines['bottom'].set_linewidth(self.bwith)
                ax.spines['left'].set_linewidth(self.bwith)
                ax.spines['right'].set_linewidth(self.bwith)
                ax.scatter(data_x[:, 0]/(2 * np.pi), np.abs(data[:, 0]), s=self.scatterpoints_size, marker='o', c=self.pointcolor)
                self.write_data(data_x[:, 0]/(2 * np.pi), np.abs(data[:, 0]), os.getcwd(), 'v.dat', 'group velocity')
                # ax.scatter(data_x[:, 0], np.abs(data[:, 0]), s=self.scatterpoints_size, marker='o', c=self.pointcolor)
                # self.write_data(data_x[:, 0], np.abs(data[:, 0]), os.getcwd(), 'v.dat', 'group velocity')
                if np.max(np.abs(data[:, 0])) > 10:
                    plt.yscale("log")
                # ax.set_xlim(0, np.max(data_x[:, 0]/(2 * np.pi))*1.05)
                # ax.set_ylim(0, np.max(data[:, 0])*2)
                ax.set_xlabel("Frequency (THz)", fontdict=font)
                ax.set_ylabel("Group Velocity (km/s)", fontdict=font)
                labels = ax.get_yticklabels() + ax.get_xticklabels()
                [label.set_fontproperties(font_properties) for label in labels]
                [label.set_fontsize(self.font_size - 6) for label in labels]
                plt.savefig(osp.join(os.getcwd(), "v.png"), dpi=300, bbox_inches='tight')
                plt.close()
                os.chdir(self.filepath)
            except:
                os.chdir(self.filepath)
                print("No BTE.v data found !")
            os.chdir(self.filepath)

    def plotT_BTE_gruneisen(self):
            os.chdir(self.filepath)
            try:
                # print("Begin process gruneisen catalogue !")
                data_x = np.loadtxt(osp.join(os.getcwd(), 'BTE.omega'), skiprows=0, dtype=float)
                data = np.loadtxt(osp.join(os.getcwd(), 'BTE.gruneisen'), skiprows=0, dtype=float)
                plt.figure(figsize=self.fig_size)
                ax = plt.subplot(111)
                ax.spines['top'].set_linewidth(self.bwith)
                ax.spines['bottom'].set_linewidth(self.bwith)
                ax.spines['left'].set_linewidth(self.bwith)
                ax.spines['right'].set_linewidth(self.bwith)
                omega = []
                greuneisen = []
                for i in range(len(data_x[0, :])):
                    for j in range(len(data_x[:, 0])):
                        ax.scatter(data_x[j, i]/(2 * np.pi), data[j, i], s=self.scatterpoints_size, marker='o', c=self.pointcolor)
                        omega.append(data_x[j, i]/(2 * np.pi))
                        # ax.scatter(data_x[j, i], data[j, i], s=self.scatterpoints_size, marker='o', c=self.pointcolor)
                        # omega.append(data_x[j, i])
                        greuneisen.append(data[j, i])
                self.write_data(np.array(omega), np.array(greuneisen), os.getcwd(), 'gruneisen.dat', 'gruneisen')
                # ax.set_xlim(np.min(data_x/(2 * np.pi)), np.max(data_x/(2 * np.pi))*1.05)
                # ax.set_xlim(0, np.max(data_x/(2 * np.pi))*1.05)
                # ax.set_ylim(0, np.max(data)*1.1)
                ax.set_xlabel("Frequency (THz)", fontdict=font)
                ax.set_ylabel("Gruneisen parameter", fontdict=font)
                # plt.yscale("log")
                labels = ax.get_yticklabels() + ax.get_xticklabels()
                [label.set_fontproperties(font_properties) for label in labels]
                [label.set_fontsize(self.font_size - 6) for label in labels]
                plt.savefig(osp.join(os.getcwd(), "gruneisen.png"), dpi=300, bbox_inches='tight')
                plt.close()
                os.chdir(self.filepath)
            except:
                print("No BTE.gruneisen data found !")
                os.chdir(self.filepath)

    def plotT_BTE_cumulative(self):
        for items in self.alldirs:
            os.chdir(osp.join(self.filepath, items))
            try:
                # print("Begin process " + str(items) + " catalogue !")
                ff = 'BTE.cumulative_kappaVsOmega_tensor'
                data = np.loadtxt(osp.join(os.getcwd(), ff), skiprows=0, dtype=float)
                plt.figure(figsize=self.fig_size)
                ax = plt.subplot(111)
                ax.spines['top'].set_linewidth(self.bwith)
                ax.spines['bottom'].set_linewidth(self.bwith)
                ax.spines['left'].set_linewidth(self.bwith)
                ax.spines['right'].set_linewidth(self.bwith)
                axis =['a-axis', 'b-axis', 'c-axis']
                index = [1, 5, 9]
                max_ylim = []
                for ii in range(3):
                    ax.plot(data[:, 0]/(2 * np.pi), np.abs(data[:, index[ii]]), c=self.pcolor[ii], label=axis[ii])
                    max_ylim.append(np.max(data[:, index[ii]]))
                    # ax.plot(data[:, 0], np.abs(data[:, index[ii]]), c=self.pcolor[ii], label=axis[ii])
                    # max_ylim.append(np.max(data[:, index[ii]]))
                max_y = max(max_ylim)
                ax.set_ylim(0, max_y*1.05)
                ax.set_xlabel("Frequency (THz)", fontdict=font)
                ax.set_ylabel("Cumulative $\kappa$ (W/mK)", fontdict=font)

                leg = ax.legend(frameon=True, loc='lower right').get_texts()
                [label.set_fontproperties(font_properties) for label in leg]
                [label.set_fontsize(self.font_size-6) for label in leg]

                labels = ax.get_yticklabels() + ax.get_xticklabels()
                [label.set_fontproperties(font_properties) for label in labels]
                [label.set_fontsize(self.font_size - 6) for label in labels]
                plt.savefig(osp.join(os.getcwd(), "cumulative.png"), dpi=300, bbox_inches='tight', pad_inches=0.05)
                plt.close()
                os.chdir(self.filepath)
            except:
                os.chdir(self.filepath)
                print("No BTE.cumulative_kappaVsOmega_tensor data found !")
            os.chdir(self.filepath)

    def compoare_eph_3ph(self, temps='300K', mark=0):
        plt.figure(figsize=self.fig_size)
        ax = plt.subplot(111)
        ax.spines['top'].set_linewidth(self.bwith)
        ax.spines['bottom'].set_linewidth(self.bwith)
        ax.spines['left'].set_linewidth(self.bwith)
        ax.spines['right'].set_linewidth(self.bwith)
        data_eph = np.loadtxt(osp.join(self.filepath, "BTE.sr_phel_metal"), skiprows=0, dtype=float)
        ax.scatter(data_eph[:, 0], data_eph[:, 1], s=self.scatterpoints_size, marker='o', c=self.pcolor[0], label="e-ph")
        count = 0
        if mark:
            ff = 'BTE.w_anharmonic'
            if osp.exists(osp.join(osp.join(self.filepath, 'T' + temps), 'BTE.w_3ph')):
                ff = 'BTE.w_3ph'
            data = np.loadtxt(osp.join(osp.join(self.filepath, 'T' + temps), ff), skiprows=0, dtype=float)
            ax.scatter(data[:, 0]/(2 * np.pi), data[:, 1], s=self.scatterpoints_size, marker='o', c=self.pcolor[1], label="3ph: " + temps)
        else:
            for items in self.alldirs:
                try:
                    ff = 'BTE.w_anharmonic'
                    if osp.exists(osp.join(osp.join(self.filepath, items), 'BTE.w_3ph')):
                        ff = 'BTE.w_3ph'
                    data = np.loadtxt(osp.join(osp.join(self.filepath, items), ff), skiprows=0, dtype=float)
                    count = count + 1
                    ax.scatter(data[:, 0]/(2 * np.pi), data[:, 1], s=self.scatterpoints_size, marker='o', c=self.pcolor[count], label="3ph: " + items)
                except:
                    print("No BTE.anharmonic data found !")
                    continue
        plt.yscale("log")
        # ax.set_xlim(0, np.max(data[:, 0]/(2 * np.pi))*1.05)
        # ax.set_ylim(0, np.max(data[:, 1])*2)
        ax.set_xlabel("Frequency (THz)", fontdict=font)
        ax.set_ylabel("Scattering Rate (1/ps)", fontdict=font)
        ax.set_ylim(10**(-3), 10)
        labels = ax.get_yticklabels() + ax.get_xticklabels()
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(self.font_size - 6) for label in labels]

        leg = ax.legend(frameon=True, loc='upper right').get_texts()
        [label.set_fontproperties(font_properties) for label in leg]
        [label.set_fontsize(self.font_size-6) for label in leg]
        plt.savefig(osp.join(self.filepath, "Eph-3ph.png"), dpi=300, bbox_inches='tight')
        plt.close()

def shengbteFig():
    print("===================================================================")
    print("================== PostProcess of ShengBTE ========================")
    print("             (S1) generate all general figs                        ")
    print("             (S2) Scattering Rate for all temperatures             ")
    print("             (S3) Thermal conductivity kappa                       ")
    print("             (S4) gruneisen constant                               ")
    print("             (S5) Compare e-ph and 3ph scattering rate (T)         ")
    print("             (S6) Compare e-ph and 3ph scattering rate             ")
    print("------------------->>")
    select = input("Input : ")
    return select

def manipulate():
    sel = shengbteFig()
    print("Begin processing !")
    if sel == 'S1':
        try:
            plot_thermal_conductivity('BTE.KappaTensorVsT_CONV', 'kappa_CONV')
            plot_thermal_conductivity('BTE.KappaTensorVsT_RTA', 'kappa_RTA')
        except:
            print("Not plot thermal conductivity, go on !")
            # plot_thermal_conductivity('BTE.KappaTensorVsT_RTA', 'kappa_RTA')
        btefig = shengbte()
        btefig.plotT_BTE_w_anharmonic()
        btefig.plotT_BTE_w_anharmonic_nolog()
        btefig.plotT_BTEWP3()
        btefig.plotT_BTEv()
        btefig.plotT_BTE_gruneisen()
        btefig.plotT_BTE_cumulative()
        btefig.plotT_BTE_NU()
        # 4ph
        btefig.plotT_BTE_w_anharmonic_4th()
        btefig.plotT_BTE_w_anharmonic_nolog_4ph()
        btefig.plotT_BTE_NU_4ph()
        btefig.plotT_BTEWP4()
        # try:
        #     btefig.plotT_BTE_w_anharmonic_4th()
        #     btefig.plotT_BTE_w_anharmonic_nolog_4ph()
        #     btefig.plot_anharmonic_branches_4ph()
        #     btefig.plotT_BTE_NU_4ph()
        #     btefig.plotT_BTEWP4()
        # except:
        #     print("No 4ph data found, go on !")

    elif sel == 'S2':
        btefig = shengbte()
        btefig.plotT_BTE_w_anharmonic()
        btefig.plotT_BTE_w_anharmonic_nolog()
    elif sel == 'S3':
        plot_thermal_conductivity()
    elif sel == 'S4':
        btefig = shengbte()
        btefig.plotT_BTE_gruneisen()
    elif sel == 'S5':
        btefig = shengbte()
        temps = '300K'
        temps = input("Input temperature (default 300K): ")
        btefig.compoare_eph_3ph(temps, 1)
    elif sel == 'S6':
        btefig = shengbte()
        btefig.compoare_eph_3ph('300K', 0)
    else:
        print("\n")
        print("Please input write parametes !")

    print("End process !")



