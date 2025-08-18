
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


def get_data_specfun(filepath, filename1, filename2, skiprow1, skiprow2):
    try:
        specfun = np.loadtxt(osp.join(os.getcwd(), filename1), skiprows=skiprow1, dtype=float)
        specfun_sup = np.loadtxt(osp.join(os.getcwd(), filename2), skiprows=skiprow2, dtype=float)
    except:
        print("Exception: too small data exist, substituting it to 0.0 !")
        specfun_tmp = open(osp.join(filepath, filename1), 'r').readlines()
        specfun_sup_tmp = open(osp.join(filepath, filename2), 'r').readlines()
        row_specfun = len(specfun_tmp[skiprow1].split())
        row_specfun_tmp = len(specfun_sup_tmp[skiprow2].split())
        specfun = np.zeros(shape=(len(specfun_tmp)-skiprow1, row_specfun))
        specfun_sup = np.zeros(shape=(len(specfun_sup_tmp)-skiprow2, row_specfun_tmp))
        for i in range(0, len(specfun_tmp)-skiprow1):
            tmp = specfun_tmp[i+skiprow1].split()
            for j in range(row_specfun):
                if '.' in tmp[j] and 'E' not in tmp[j].split('.')[1] and '-' in tmp[j].split('.')[1]:
                        specfun[i, j] = float(tmp[j].split('.')[0]+tmp[j].split('.')[1].replace('-', 'E-'))
                else:
                    specfun[i, j] = float(tmp[j])
        for ii in range(0, len(specfun_sup_tmp)-skiprow2):
            tmp1 = specfun_sup_tmp[i+skiprow2].split()
            for jj in range(row_specfun_tmp):
                if '.' in tmp1[jj] and 'E' not in tmp1[jj].split('.')[1] and '-' in tmp1[jj].split('.')[1]:
                        specfun_sup[i, j] = float(tmp1[jj].split('.')[0]+tmp1[jj].split('.')[1].replace('-', 'E-'))
                else:
                    specfun_sup[ii, jj] = float(tmp1[jj])
    return specfun, specfun_sup


def specfun_elec():
    files = os.listdir()
    filename = ''
    for i in range(len(files)):
        if 'specfun.elself' in files[i]:
            filename = files[i].replace('specfun', '')
            break
    filename1 = "specfun" + filename
    filename2 = "specfun_sup" + filename
    print("Filename: ", filename1)
    print("Filename: ", filename2)
    specfun, specfun_sup = get_data_specfun(os.getcwd(), filename1, filename2, 6, 6)
    # cm = plt.cm.jet
    # cm = plt.cm.gist_earth
    # cm = plt.cm.cubehelix
    # cm = plt.cm.seismic
    # cm = plt.cm.seismic
    # cm = plt.cm.viridis
    # cm = plt.cm.twilight_shifted
    cm = plt.cm.hot
    branches = len(set(specfun_sup[:, 1]))
    plt.figure(figsize=(6, 4))
    ax = plt.subplot(111)
    x = specfun[:, 0]
    y = specfun[:, 1]
    z = specfun[:, 2]
    z1 = specfun_sup[:, 2].reshape(branches, -1)
    cm = plt.cm.rainbow
    p = ax.scatter(x, y, s=8, c=z, cmap=cm, marker="o", vmin=np.min(z), vmax=np.max(z), alpha=0.8)
    cb = plt.colorbar(p)
    ax.set_ylabel("$\omega$ (meV)")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y), np.max(y))
    sp = 0.5
    ax.plot(x, z1[0, :], color='red', linewidth=sp)
    ax.plot(x, z1[1, :], color='red', linewidth=sp)
    ax.plot(x, z1[2, :], color='red', linewidth=sp)
    ax.plot([np.min(x), np.max(x)], [0, 0], color='black', ls='--', linewidth=0.8)
    plt.savefig(osp.join(os.getcwd(), "spec_elec.png"), dpi=300, bbox_inches = 'tight', pad_inches=0.1)

def specfun_phon():
    files = os.listdir()
    filename = ''
    for i in range(len(files)):
        if 'specfun.phon' in files[i]:
            filename = files[i].replace('specfun', '')
            break
    filename1 = "specfun" + filename
    filename2 = "specfun_sup.phon"
    print("Filename: ", filename1)
    print("Filename: ", filename2)
    specfun, specfun_sup = get_data_specfun(os.getcwd(), filename1, filename2, 4, 2)
    cm = plt.cm.hot
    branches = len(set(specfun_sup[:, 1]))
    plt.figure(figsize=(6, 4))
    ax = plt.subplot(111)
    x = specfun[:, 0]
    y = specfun[:, 1]
    z = specfun[:, 2]
    # z1 = specfun_sup[:, 2].reshape(branches, -1)
    cm = plt.cm.rainbow
    p = ax.scatter(x, y, s=8, c=z, cmap=cm, marker="o", vmin=np.min(z), vmax=np.max(z), alpha=0.8)
    cb = plt.colorbar(p)
    ax.set_ylabel("$\omega$ (meV)")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y), np.max(y))
    # sp = 0.5
    # ax.plot(x, z1[0, :], color='red', linewidth=sp)
    # ax.plot(x, z1[1, :], color='red', linewidth=sp)
    # ax.plot(x, z1[2, :], color='red', linewidth=sp)
    # ax.plot([np.min(x), np.max(x)], [0, 0], color='black', ls='--', linewidth=0.8)
    plt.savefig(osp.join(os.getcwd(), "spec_phon.png"), dpi=300, bbox_inches = 'tight', pad_inches=0.1)

