# -*- coding: utf-8 -*-
"""
Created on 16:16 25-08-2021

@author: XY Ding
mail to: dxy_vasp@163.com
python3: gen_kp.py
"""
import numpy as np
import AutoinputQE as inpp
import os, sys
import os.path as osp 

sys.path.append(os.getcwd())
sys.dont_write_bytecode = True  

############################ get high symmetry points from KPATH.in file ##################################

def get_kpt():
    kpoints_bands = inpp.kpoints_band.split("\n")[1:-1]
    kpt = []
    hsp = []
    for kk in range(len(kpoints_bands)):
        tmp = []
        for ii in range(3):
            tmp.append(kpoints_bands[kk].split("!")[0].split()[ii])
        tmpp = [float(item) for item in tmp]
        kpt.append(tmpp)
        hsp.append(kpoints_bands[kk].split("!")[1])
    return kpt, hsp

################################ generate kpoints for all the kpt ###################################
def generate_kpoint(k_num):
    kpath, high_kpt_point = get_kpt()
    len_k = len(kpath)
    kpoints = np.zeros(shape=(k_num * (len_k-1), 3))
    for i in range(1, len_k):
        for j in range(0, k_num):
            kpoints[(i - 1) * k_num + j][0] = (kpath[i][0] - kpath[i-1][0]) * (j + 1) / k_num + kpath[i-1][0]
            kpoints[(i - 1) * k_num + j][1] = (kpath[i][1] - kpath[i - 1][1]) * (j + 1) / k_num + kpath[i - 1][1]
            kpoints[(i - 1) * k_num + j][2] = (kpath[i][2] - kpath[i - 1][2]) * (j + 1) / k_num + kpath[i - 1][2]
    return kpoints


################################### output all the kpoints #############################################
def output_kpoints(k_num=20):
    a, high_kpt_point = get_kpt()
    kpoints = generate_kpoint(k_num)
    if not osp.exists(osp.join(os.getcwd(), "meshes")):
        os.makedirs(osp.join(os.getcwd(), "meshes"))
    else:
        # os.system('rm -rf meshes')
        os.makedirs(osp.join(os.getcwd(), "meshes"), exist_ok=True)
    file = osp.join(os.getcwd(), "meshes/kpoints")
    with open(file, "w", encoding='utf-8') as kp:
        tot_k = int(k_num * (len(a) - 1) + 1)
        kp.write(str(int(k_num * (len(a) - 1) + 1)) + '  crystal \n')
        kp.write(str(a[0][0])[:6].ljust(8, " ") + str(a[0][1])[:6].ljust(8, " ") + str(a[0][2])[:6].ljust(8, " ") + \
                "  " + str(1.0).ljust(8, " ") + "\n")
        for flag_knum in range(0, k_num * (len(a)-1)):
            kp.write(str(kpoints[flag_knum][0])[:6].ljust(8, " ") + str(kpoints[flag_knum][1])[:6].ljust(8, " ") + \
                    str(kpoints[flag_knum][2])[:6].ljust(8, " ") + "  " + str(1.0).ljust(8, " ") + "\n")
    kp.close()

################################### output all the qpoints #############################################
def output_qpoints(k_num=20):
    a, high_kpt_point = get_kpt()
    kpoints = generate_kpoint(k_num)
    if not osp.exists(osp.join(os.getcwd(), "meshes")):
        os.makedirs(osp.join(os.getcwd(), "meshes"))
    else:
        os.system('rm -rf meshes')
        os.makedirs(osp.join(os.getcwd(), "meshes"))
    file = osp.join(os.getcwd(), "meshes/qpoints")
    with open(file, "w", encoding='utf-8') as kp:
        tot_k = int(k_num * (len(a) - 1) + 1)
        kp.write(str(int(k_num * (len(a) - 1) + 1)) + '  crystal \n')
        kp.write(str(a[0][0])[:6].ljust(8, " ") + str(a[0][1])[:6].ljust(8, " ") + str(a[0][2])[:6].ljust(8, " ") + \
                "  " + str(1.0).ljust(8, " ") + "\n")
        for flag_knum in range(0, k_num * (len(a)-1)):
            kp.write(str(kpoints[flag_knum][0])[:6].ljust(8, " ") + str(kpoints[flag_knum][1])[:6].ljust(8, " ") + \
                    str(kpoints[flag_knum][2])[:6].ljust(8, " ") + "  " + str(1.0).ljust(8, " ") + "\n")
    kp.close()

def manipulate():
    output_kpoints(20)
    
if __name__ == '__main__':
    manipulate()






