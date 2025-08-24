
import numpy as np 
import matplotlib as mpl
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
            tmp1 = specfun_sup_tmp[ii+skiprow2].split()
            for jj in range(row_specfun_tmp):
                if '.' in tmp1[jj] and 'E' not in tmp1[jj].split('.')[1] and '-' in tmp1[jj].split('.')[1]:
                        specfun_sup[ii, jj] = float(tmp1[jj].split('.')[0]+tmp1[jj].split('.')[1].replace('-', 'E-'))
                else:
                    specfun_sup[ii, jj] = float(tmp1[jj])
    return specfun, specfun_sup

def extract_realOmega(omega, w, real, Qpoints, branches):  # extract real phonon self-energy
    omega_Qpoints_branches = omega.reshape(Qpoints, -1)
    w_Qpoints_branches = w.reshape(Qpoints, -1)
    real_Qpoints_branches = real.reshape(Qpoints, -1)
    real_new = []
    phonon_data = []
    for i in range(Qpoints):
        qp_wqv = omega_Qpoints_branches[i, :branches]
        phonon_data.append(qp_wqv)
        QP_phon_branches = np.zeros(branches)
        for k in range(len(qp_wqv)):
            judge_threshold = []
            for j in range(len(w_Qpoints_branches[i])):
                judge_threshold.append(np.abs(omega_Qpoints_branches[i, j] - qp_wqv[k]))
            min_threshold = min(judge_threshold)
            for kk in range(len(w_Qpoints_branches[i])):
                if np.abs(omega_Qpoints_branches[i, kk]-qp_wqv[k]) == min_threshold:
                    QP_phon_branches[k] = real_Qpoints_branches[i, kk]
                    break 
        real_new.append(QP_phon_branches)
    return np.array(phonon_data), np.array(real_new)

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
    plt.close()

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
    z = specfun[:, 3]
    # z1 = specfun_sup[:, 2].reshape(branches, -1)
    cm = plt.cm.rainbow
    p = ax.scatter(x, y, s=8, c=z, cmap=cm, marker="o", norm=mpl.colors.LogNorm(), vmin=10**-10, vmax=10**2, alpha=0.8)
    # p = ax.scatter(x, y, s=8, c=z, cmap=cm, marker="o", norm=mpl.colors.LogNorm(), alpha=0.8)
    cb = plt.colorbar(p)
    ax.set_ylabel("$\omega$ (meV)")
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y), np.max(y))
    # sp = 0.5
    # ax.plot(x, z1[0, :], color='red', linewidth=sp)
    # ax.plot(x, z1[1, :], color='red', linewidth=sp)
    # ax.plot(x, z1[2, :], color='red', linewidth=sp)
    # ax.plot([np.min(x), np.max(x)], [0, 0], color='black', ls='--', linewidth=0.8)
    plt.savefig(osp.join(os.getcwd(), "spec_phon2.png"), dpi=300, bbox_inches = 'tight', pad_inches=0.1)
    plt.close()

def specfun_phon_renormal():
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
    omega = specfun_sup[:, 3]
    w = specfun_sup[:, 4]
    real_ori = specfun_sup[:, 5]
    Qpoints = len(set(specfun_sup[:, 0]))
    branches = len(set(specfun_sup[:, 1]))
    x = np.arange(Qpoints)
    print("---------------------- begin processing specfun_sup real data -----------------------")
    print("Qpoints, branches:", Qpoints, branches)
    phon, real = extract_realOmega(omega, w, real_ori, Qpoints, branches)
    print("----------------------- end processing specfun_sup real data ------------------------")
    renormal_ph = np.zeros(shape=(Qpoints, branches))
    for i in range(Qpoints):
        for j in range(branches):
            tmp = np.abs(phon[i, j]**2 + 2*phon[i, j]*real[i, j]*0.001)
            if tmp >=0:
                renormal_ph[i, j] = np.sqrt(np.abs(phon[i, j]**2 + 2*phon[i, j]*real[i, j]*0.001))
            else:
                renormal_ph[i, j] = -np.sqrt(np.abs(phon[i, j]**2 + 2*phon[i, j]*real[i, j]*0.001))
    # return x, phon.T, np.array(renormal_ph).T
    x, y0, y1 = x, phon.T, np.array(renormal_ph).T
    ori_phon = np.vstack((x, phon.T))
    re_phon = np.vstack((x, np.array(renormal_ph).T))
    print("Data:", ori_phon.shape, re_phon.shape)
    try:
        np.savetxt(osp.join(os.getcwd(), "Ori_phon.dat"), ori_phon.T)
        np.savetxt(osp.join(os.getcwd(), "Renormal_phonon.dat"), re_phon.T)
    except:
        print("Not save phonon and renormalized data !")
    plt.figure(figsize=(5, 4))
    ax = plt.subplot(111)
    for i in range(len(y0)):
        ax.scatter(x, y0[i], color="black", s=3)
        ax.scatter(x, y1[i], color='red', s=3)
    ax.scatter([], [], color="black", label="phonon")
    ax.scatter([], [], color="red", label="renormal")
    ax.set_ylabel("Frequency (eV)")
    ax.set_xlim(0, max(x))
    ax.set_ylim(0, max([np.max(y0), np.max(y1)])*1.02)
    ax.legend(loc='upper left')
    plt.savefig(osp.join(os.getcwd(), "renormal.png"), dpi=300, bbox_inches = 'tight', pad_inches=0.2)
    plt.close()

def comp_dfpt_epw_g(filepath, filedfpt='dfpt.dat', fileepw='epw.dat'):
    font_size = font['size']
    kp = inpp.phon_points
    hsp_labels = [str(items.split("!")[1]).strip() for items in inpp.kpoints_band.split("\n")[1:-1]]
    hsp_data = [kp*it for it in range(len(hsp_labels))]
    print("HSP labels: ", hsp_labels)
    print("HSP   data: ", hsp_data)
    dfpt = np.loadtxt(osp.join(filepath, filedfpt), skiprows=0, dtype=float)
    epw = np.loadtxt(osp.join(filepath, fileepw), skiprows=0, dtype=float)
    plt.figure(figsize=(5, 4))
    ax = plt.subplot(111)
    len_x = len(dfpt[:, 0])
    x = [it for it in range(0, len_x)]
    ax.scatter(x, dfpt[:, -1], s=5, label="ph.x", color='none', marker="o", edgecolors='red')
    ax.scatter(x, epw[:, -1], s=4, label="EPW", color='gray')
    ax.set_xticks(hsp_data)
    max_dfpt = np.max(dfpt[:, -1])
    max_epw = np.max(dfpt[:, -1])
    max_y = max([max_dfpt, max_epw])
    for i in range(len(hsp_data)):
        ax.vlines(hsp_data, ymin=-10, ymax=2*max_y,linewidth=1,linestyle="--", color='gray')
    ax.set_xticklabels(hsp_labels)
    ax.set_ylim(0, 1.1*max_y)
    ax.set_xlim(0, len_x)
    ax.set_ylabel("|g|$_{avg}$ (meV)", fontdict=font)
    labels = ax.get_xticklabels() + ax.get_yticklabels() 
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(font_size-4) for label in labels]
    leg = ax.legend(frameon=True, loc='upper right').get_texts()
    [label.set_fontproperties(font_properties) for label in leg]
    [label.set_fontsize(font_size-6) for label in leg]
    plt.savefig("g_DFPT_EPW.png", dpi=400, bbox_inches = 'tight', pad_inches=0.2)
    plt.close()