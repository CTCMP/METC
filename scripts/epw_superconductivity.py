
import numpy as np 
from matplotlib import pyplot as plt 
import os, sys
import os.path as osp
import AutoinputQE as inpp
from matplotlib import font_manager
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
		'size': 22, 
        }
sys.dont_write_bytecode = True       ## not generate __pycache__ files

me2VToTHz =  0.241799050402417

font_size = font['size']
font_size_axis = font_size - 8
lwidth = 1.5

def plot_a2f(filepath, filename):
    x = []
    y1 = []
    y2 = []
    af = open(osp.join(filepath, filename), 'r').readlines()
    for i in range(1, len(af)):
        if 'Integrated' in af[i]:
            break
        else:
            x.append(float(af[i].split()[0]))
            y1.append(float(af[i].split()[10]))
            y2.append(float(af[i].split()[-1]))
    xx = np.array(x)
    yy1 = np.array(y1)
    yy2 = np.array(y2)
    plt.figure(figsize=(5, 4))
    ax = plt.subplot(111)
    ax.plot(xx, yy1, color='blue', lw=1.5)
    ax.plot(xx, yy2, color='red', ls='--', lw=1.5)
    ax.set_xlim(0, np.max(xx))
    ax.set_xlabel(r'$\omega$ (meV)', fontdict=font)
    ax.set_ylabel(r'$\alpha_{2}$F ($\omega$)', fontdict=font)
    labels = ax.get_xticklabels() + ax.get_yticklabels() 
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(font['size']-4) for label in labels]
    plt.savefig(osp.join(os.getcwd(), "a2f.png"), dpi=150, bbox_inches = 'tight', pad_inches=0.1)
    plt.close()
##################  superconductivity anisotropic ME eqs ###############
def lambda_pairs(filepath, filename):
    lpirs = np.loadtxt(osp.join(filepath, filename[0]), skiprows=1, dtype=float)
    plt.figure(figsize=(5, 4))
    ax = plt.subplot(211)
    ax.plot(lpirs[:, 0], lpirs[:, 1], lw=lwidth)
    ax.set_xlim(np.min(lpirs[:, 0]), np.max(lpirs[:, 0]))
    ax.set_xlabel('$\lambda_{nk}$, mk+q', fontdict=font)
    labels = ax.get_xticklabels() + ax.get_yticklabels() 
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(font_size_axis) for label in labels]
    ax2 = plt.subplot(212)
    klpirs = np.loadtxt(osp.join(filepath, filename[1]), skiprows=1, dtype=float)
    ax2.plot(klpirs[:, 0], klpirs[:, 1], lw=lwidth)
    ax2.set_xlim(np.min(klpirs[:, 0]), np.max(klpirs[:, 0]))
    ax2.set_xlabel('$\lambda_{nk}$', fontdict=font)
    labels = ax2.get_xticklabels() + ax2.get_yticklabels() 
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(font_size_axis) for label in labels]
    plt.savefig(osp.join(os.getcwd(), "Lambda_kpairs.png"), dpi=200, bbox_inches = 'tight', pad_inches=0.2)
    plt.close()
def Imag_aniso(filepath, flist, figdir, index_column, label_x):
    if not osp.exists(osp.join(os.getcwd(), figdir)):
        os.makedirs(osp.join(os.getcwd(), figdir))
    for i in range(len(flist)):
        data = np.loadtxt(osp.join(filepath, flist[i]), skiprows=1, dtype=float)
        tmp = flist[i].split('.')[1].split('_')[-1]
        plt.clf()
        plt.figure(figsize=(5, 4))
        ax = plt.subplot(111)
        ax.scatter(data[:, 0]*1000, data[:, index_column]*1000, s=8, color='black')
        ax.set_xlabel(label_x, fontdict=font)
        ax.set_ylabel(r'$\Delta_{nk}$ (meV)', fontdict=font)
        ax.set_xlim(np.min(data[:, 0])*1000, np.max(data[:, 0])*1000)
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(font_size_axis) for label in labels]
        plt.savefig(osp.join(os.getcwd(), figdir, "Aniso_imag_"+ tmp +".png"), dpi=200, bbox_inches = 'tight', pad_inches=0.2)
        plt.close()
def manipulate_temps(filepath, label_files, figdir, index_column, label_x):
    listdirs =  os.listdir(filepath)
    dirfiles = []
    for i in range(len(listdirs)):
        if label_files in listdirs[i]:
            dirfiles.append(listdirs[i])
    Imag_aniso(filepath, dirfiles, figdir, index_column, label_x)

###################  superconductivity isotropic ME eqs ################
def imag_iso_single(filepath, flist, figdir, index_column, label_x):
    if not osp.exists(osp.join(os.getcwd(), figdir)):
        os.makedirs(osp.join(os.getcwd(), figdir))
    for i in range(len(flist)):
        data = np.loadtxt(osp.join(filepath, flist[i]), skiprows=1, dtype=float)
        tmp = flist[i].split('.')[1].split('_')[-1]
        plt.clf()
        plt.figure(figsize=(5, 4))
        ax = plt.subplot(111)
        ax.plot(data[:, 0]*1000, data[:, index_column]*1000, color='black')
        ax.set_xlabel(label_x, fontdict=font)
        ax.set_ylabel(r'$\Delta_{nk}$ (meV)', fontdict=font)
        ax.set_xlim(np.min(data[:, 0])*1000, np.max(data[:, 0])*1000)
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(font_size_axis) for label in labels]
        plt.savefig(osp.join(os.getcwd(), figdir, "Iso_imag_"+ tmp +".png"), dpi=200, bbox_inches = 'tight', pad_inches=0.2)
        plt.close()
def imag_iso_mul(filepath, dirfiles_pade, dirfiles_acon, figdir, label_x):
    if not osp.exists(osp.join(os.getcwd(), figdir)):
        os.makedirs(osp.join(os.getcwd(), figdir))
    for i in range(len(dirfiles_pade)):
        data_pade = np.loadtxt(osp.join(filepath, dirfiles_pade[i]), skiprows=1, dtype=float)
        data_acon = np.loadtxt(osp.join(filepath, dirfiles_acon[i]), skiprows=1, dtype=float)
        tmp = dirfiles_pade[i].split('.')[1].split('_')[-1]
        plt.figure(figsize=(5, 4))
        ax = plt.subplot(111)
        ax.plot(data_pade[:, 0]*1000, data_pade[:, 3]*1000,  color='red', label=r'Re $\Delta$-Pade approx')
        ax.plot(data_pade[:, 0]*1000, data_pade[:, 4]*1000,  color='green', ls='--', label=r'Im $\Delta$-Pade approx')
        ax.plot(data_acon[:, 0]*1000, data_acon[:, 3]*1000,  color='blue', label=r'Re $\Delta$-analytic approx')
        ax.plot(data_acon[:, 0]*1000, data_acon[:, 4]*1000,  color='black', ls='--', label=r'Im $\Delta$-analytic approx')
        ax.set_xlabel(label_x, fontdict=font)
        ax.set_ylabel(r'$\Delta_{nk}$ (meV)', fontdict=font)
        ax.set_xlim(np.min(data_pade[:, 0])*1000, np.max(data_pade[:, 0])*1000)
        leg = ax.legend(frameon=True, loc='upper right').get_texts()
        [label.set_fontproperties(font_properties) for label in leg]
        [label.set_fontsize(font_size_axis) for label in leg]
        labels = ax.get_xticklabels() + ax.get_yticklabels() 
        [label.set_fontproperties(font_properties) for label in labels]
        [label.set_fontsize(font_size_axis) for label in labels]
        plt.savefig(osp.join(os.getcwd(), figdir, "PadAcon_imag_"+ tmp +".png"), dpi=200, bbox_inches = 'tight', pad_inches=0.2)
        plt.close()
def manipulate_temps_single(filepath, label_files, figdir, index_column, label_x):
    listdirs =  os.listdir(filepath)
    dirfiles = []
    for i in range(len(listdirs)):
        if label_files in listdirs[i]:
            dirfiles.append(listdirs[i])
    imag_iso_single(filepath, dirfiles, figdir, index_column, label_x)

def manipulate_temps_mul(filepath, label_files, figdir, label_x):
    listdirs =  os.listdir(filepath)
    dirfiles_pade = []
    dirfiles_acon = []
    for i in range(len(listdirs)):
        if label_files in listdirs[i]:
            headfiles = listdirs[i].split('.')[1].split('_')[-1] + '.' + listdirs[i].split('.')[-1]
            dirfiles_pade.append('pwscf.pade_iso_' + headfiles)
            dirfiles_acon.append('pwscf.acon_iso_' + headfiles)
    imag_iso_mul(filepath, dirfiles_pade, dirfiles_acon, figdir, label_x)

def epw_superconductivity_aniso_manipulate():
    print('--------- plotting lambda pairs ----------')
    lambda_pairs(os.getcwd(), ['pwscf.lambda_pairs', 'pwscf.lambda_k_pairs'])
    print('******************************************')
    print('-------------- plotting a2f --------------')
    plot_a2f(os.getcwd(), 'pwscf.a2f')
    print('******************************************')
    print('---------- plotting imag_aniso -----------')
    manipulate_temps(os.getcwd(), 'pwscf.imag_aniso_0', 'Imag_aniso', 3, r'i$\omega$ (meV)')
    print('******************************************')
    print('------- plotting pade imag_aniso ---------')
    manipulate_temps(os.getcwd(), 'pwscf.pade_aniso_0', 'Imag_pade_aniso', 4, r'$\omega$ (meV)')

def epw_superconductivity_iso_manipulate():
    print('----------- plotting imag_iso ------------')
    manipulate_temps_single(os.getcwd(), 'pwscf.imag_iso_0', 'Imag_iso', 2, r'i$\omega$ (meV)')
    print('******************************************')
    print('-------- plotting pade imag_iso ----------')
    manipulate_temps_mul(os.getcwd(), 'pwscf.pade_iso_0', 'Imag_pade_iso', r'$\omega$ (meV)')
    print('******************************************')
    print('-------------- plotting a2f --------------')
    plot_a2f(os.getcwd(), 'pwscf.a2f')
    print('******************************************')
