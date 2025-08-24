# -*- coding: utf-8 -*-
"""
Created on 10:14 13-11-2022

@author: Xian-Yong Ding
mail to: dxy_vasp@163.com
python3: funcs.py
"""
import numpy as np
import math
import os, sys 
import os.path as osp
import AutoinputQE as inpp
sys.dont_write_bytecode = True       ## not generate __pycache__ files

# get reciprocal unit cell
def ChaCheng(a1, a2):
    c = []
    c1 = float(a1[1] * a2[2] - a1[2] * a2[1])
    c.append(c1)
    c2 = float(a1[2] * a2[0] - a1[0] * a2[2])
    c.append(c2)
    c3 = float(a1[0] * a2[1] - a1[1] * a2[0])
    c.append(c3)
    return c
def DianCheng(b1, b2):
    d1 = float(b1[0] * b2[0])
    d2 = float(b1[1] * b2[1])
    d3 = float(b1[2] * b2[2])
    d = d1 + d2 + d3
    return d
def get_poscar(filepath="", filename='POSCAR'):
    filepath = osp.join(filepath, filename)
    poscar1 = open(filepath, 'r')
    poscar_lines = poscar1.readlines()
    scale = float(poscar_lines[1].strip())
    pos = []
    for flag_lines in range(2, 5):
        for i in range(0, 3):
            pos.append(float(poscar_lines[flag_lines].split()[i]))
    latt = np.array(pos).reshape(3, 3) * scale
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
    # print(elements)
    ############# get coordinate ###################
    # cor = np.loadtxt(filepath, skiprows=8, dtype=np.str_, encoding="UTF-8")
    cor_tmp = np.loadtxt(filepath, skiprows=8)
    if tot_ele == 1:
        tmp = [cor_tmp]
        cor = np.array(tmp)
    else:
        cor = cor_tmp[:tot_ele, 0:3]
    return np.round(latt, 10), element, num_element, tot_ele, np.round(cor, 10)

def get_nbnds(filepath, filename):
    file1 = open(osp.join(filepath, filename), 'r').readlines()
    nbands = int(file1[0].split(',')[0].split('=')[1])
    return nbands
    
def get_pos_all(filepath):
    path = osp.join(filepath, 'POSCAR')
    poscar = open(path, 'r')
    poscar_lines = poscar.readlines()
    pos = []
    for flag_lines in range(2, 5):
        for i in range(0, 3):
            pos.append(float(poscar_lines[flag_lines].split()[i]))   
    pos = get_poscar(filepath, 'POSCAR')
    a1 = [pos[0], pos[1], pos[2]]
    a2 = [pos[3], pos[4], pos[5]]
    a3 = [pos[6], pos[7], pos[8]]
    latt = [a1, a2, a3]
    return pos 
########## get elements informations
def get_element(filepath="", filename="POSCAR"):
    file = osp.join(filepath, filename)
    poscar = open(file, 'r').read().strip('\n').splitlines()
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
    return element, num_element, tot_ele

def get_reciprocal():
    pos = get_poscar()
    rep = []
    a1 = [pos[0], pos[1], pos[2]]
    a2 = [pos[3], pos[4], pos[5]]
    a3 = [pos[6], pos[7], pos[8]]
    volume = DianCheng(a1, ChaCheng(a2, a3))
    scalar = 2 * math.pi / volume
    for i in range(0, 3):
        b1 = scalar*ChaCheng(a2, a3)[i]
        rep.append(b1)
    for j in range(0, 3):
        b2 = scalar*ChaCheng(a3, a1)[j]
        rep.append(b2)
    for k in range(0, 3):
        b3 = scalar * ChaCheng(a1, a2)[k]
        rep.append(b3)
    return rep

# transform the qe output file of relax to POSCAR
# get number of atoms
def get_natom(filename = "relax.out"):
    relax = open(filename, 'r')
    data = relax.readlines()
    natom = 0
    for i in range(0, len(data)):
        if "number of atoms/cell" in data[i]:
            natom = data[i].split('=')[1]
            break
    return int(natom)

def trans_stdPOSCAR(filepath, filename="PRIMCELL.vasp"):
    file = osp.join(filepath, filename)
    data = open(file, 'r').readlines()
    posQE = open(osp.join(filepath, "POSCAR"), 'w')
    for i in range(0, 8):
        posQE.write(data[i])
    for j in range(8, len(data)):
        for k in range(3):
            posQE.write(data[j].split()[k].rjust(20, " ") + "    ")
        posQE.write('\n')
    posQE.close()    

def transPOSCAR_QE(filepath, filename="POSCAR"):
    file = osp.join(filepath, filename)
    data = open(file, 'r').readlines()
    posQE = open(osp.join(filepath, "POSCAR-QE"), 'w')
    for i in range(2, len(data)):
        posQE.write(data[i])
    posQE.close()

def get_ibrav(filepath, filename="cellinfo"):
    file = osp.join(filepath, filename)
    data = open(file, 'r').readlines()
    for i in range(len(data)):
        if 'ibrav' in data[i]:
            ibrav = int(data[i].split()[-1])
            return ibrav
    
def get_celldm(filepath, filename="cellinfo"):
    file = osp.join(filepath, filename)
    data = open(file, 'r').readlines()
    # 读取POSCAR文件
    celldm = {}
    for j in range(len(data)):
        if 'celldm' in data[j]:
            celldm[data[j].split()[0]] = np.round(float(data[j].split()[-1]), 6)
            # celldm[data[j].split()[0]] = np.round(float(data[j].split()[-1]) * 1.8897161646320724, 6) 
    celldm['celldm(1)'] = celldm['celldm(1)'] * 1.8897161646320724
    return celldm

def get_ispin_SOC(filepath, filename):
    data = open(osp.join(filepath, filename), 'r').readlines()
    nspin = 1
    soc = 0
    for i in range(len(data)):
        if 'nspin' in data[i]:
            nspin = int(data[i].split('=')[1].split(',')[0])
        if 'noncolin' in data[i]:
            if str(data[i].split('=')[1].split(',')[0]).lower() == 'true'.lower():
                soc = 1
    return [nspin, soc]

def get_fermi_level(filepath, filename="scf.log"):
    try:
        data = open(osp.join(filepath, filename), 'r').readlines()
        fermi_all = []
        for i in range(len(data)):
            if "the Fermi energy" in data[i]:
                fermi_all.append(float(data[i].split()[4]))
        # print(fermi)
        return fermi_all[-1]
    except:
        return 0.0

def get_nks_hsplabels(filepath, filename):
    import AutoinputQE as inpp
    file1 = open(osp.join(filepath, filename), 'r').readlines()
    nks = int(file1[0].split(',')[1].split("=")[1].split("/")[0])
    hsp_labels = [str(items.split("!")[1]).strip() for items in inpp.kpoints_band.split("\n")[1:-1]]
    hsp_points = [[float(str(items.split("!")[0]).split()[i]) for i in range(3) ] for items in inpp.kpoints_band.split("\n")[1:-1]]
    print(hsp_points)
    for i in range(len(hsp_labels)):
        # if str(hsp_labels[i]).lower() == "\Gamma".lower():
        #     hsp_labels[i] = u"Γ"
        if '\\' in hsp_labels[i].lower():
            hsp_labels[i] = "$" + hsp_labels[i] + "$"
    return nks, hsp_labels
def cal(b1, b2):
    l = np.sqrt((b2[0] - b1[0])**2 + (b2[1] - b1[1])**2 + (b2[2] - b1[2])**2)
    return l
def get_xk_phon(filepath):
    data = np.loadtxt(osp.join(filepath, '../phonon/pwscf.freq.gp'), skiprows=0, dtype=float)
    return data[:, 0]

def get_xk_band(filepath,):
    data = np.loadtxt(osp.join(filepath, "./band.dat.gnu"), skiprows=0, dtype=float)
    data_k = []
    data_k.append(data[0, 0])
    for i in range(1, len(data[:, 0])):
        if float(data[i, 0]) == 0:
            break 
        else:
            data_k.append(data[i, 0])
    return np.array(data_k)

def epwband_phonon(filepath, filename, fermi=0):                    # process the kpoints file for EPW band and phonon 
    data = open(osp.join(filepath, filename), 'r').readlines()
    nks = int(data[0].split(',')[1].split()[-2])
    kpt = []
    data2 = []
    for i in range(nks):
        kpt.append(np.array(data[2*i + 1].split(), dtype=float))
        data2.append(np.array(data[2*i + 2].split(), dtype=float))
    data3 = np.array(data2) - fermi
    if fermi:
        xk = get_xk_band(filepath)
    else:
        xk = get_xk_phon(filepath)
    band = open(osp.join(filepath, filename+'.gnu'), 'w')
    for k in range(len(xk)):
        band.write(str(np.round(xk[k], 8)).rjust(12, ' ') + "    " )
        for ee in range(len(data3[0])):
            band.write(str(np.round(data3[k, ee], 8)).rjust(12, ' ') + "    " )
        band.write('\n')
    band.close()

def epmat_decay(filepath, filename, label):
    if osp.exists(osp.join(filepath, filename)):
        import json
        decay = open(osp.join(filepath, filename + ".gnu"), 'w')
        decay.write('set encoding iso_8859_1 \n')
        decay.write('set terminal png truecolor enhanced font "Arial, 40" size 1200, 1000 \n')
        decay.write('set output \'' + filename.replace('.', '_') + '.png' + '\' \n')
        decay.write('set grid \n')
        decay.write('set style data linespoints \n')
        decay.write('set xlabel ')
        decay.write('\"|' + label[0] + '|  {\\305} \"')
        # decay.write('\"|' + label[0] + '| {/Symbol \\305} \"')
        decay.write('\n')
        decay.write('set ylabel ')
        json.dump('max ' + label[1] + ' [Ry] ', decay)
        decay.write('\n')
        decay.write('set logscale y \n')
        decay.write('set format y \"10^{%.0T}\" \n')
        decay.write('set pointsize 3 \n')
        decay.write('plot \"'+ filename +'\" u 1:2 w p pt 7 \n')
        decay.close()
        os.system('gnuplot ' + filename + '.gnu')
    else:
        print("No file: ", filename)

#------------------- Data manipulating ----------------------
#-------------------------  QE  -----------------------------
def band_data_pro_phon(filepath, filename):
    file1 = open(osp.join(filepath, filename), 'r').readlines()
    file2 = open(osp.join(filepath, filename + ".gp"), 'r').readlines()
    nbands = int(file1[0].split(',')[0].split('=')[1])
    nks = int(file1[0].split(',')[1].split("=")[1].split("/")[0])
    band_data = np.loadtxt(file2, skiprows=0, dtype=float)
    band_k = band_data[:nks, 0]
    kk = band_data[:, 0]
    energy = band_data[:, 1:len(band_data[0, :])]/33.35641
    hsp_labels = [str(items.split("!")[1]).strip() for items in inpp.kpoints_band.split("\n")[1:-1]]
    for i in range(len(hsp_labels)):
        if str(hsp_labels[i]).lower() == "\Gamma".lower():
            hsp_labels[i] = u"Γ"
    hsp_nums = len(hsp_labels)-1
    hsp_data = band_k[0:nks:(int(nks/(hsp_nums)))]
    hsp = {
        'klabels' : hsp_labels,
        'hspdata' : hsp_data,
    }
    # print(hsp['klabels'])
    return nbands, nks, hsp, band_k, kk, energy 

def band_data_pro_band(filepath, filename):
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
        if str(hsp_labels[i]).lower() == "\Gamma".lower():
            hsp_labels[i] = u"Γ"
    hsp_nums = len(hsp_labels)-1
    hsp_data = band_k[0:nks:(int(nks/(hsp_nums)))]
    if hsp_data[-1] != band_k[-1]:
        hsp_data = np.append(hsp_data, band_k[-1])
    hsp = {
        'klabels' : hsp_labels,
        'hspdata' : hsp_data,
    }
    return nbands, nks, hsp, band_k, kk, energy 

def write_band_k_perturbo(filepath, filename, nks):
    conf = open(osp.join(filepath, filename), 'w')
    # conf.write("\n")
    kpoints_bands = inpp.kpoints_band.split("\n")[1:-1]
    conf.write(str(len(kpoints_bands)) + '\n')
    for pp in range(len(kpoints_bands)):
        if pp != len(kpoints_bands) -1:
            conf.write(kpoints_bands[pp].split("!")[0].ljust(len(kpoints_bands[pp].split("!")[0])) + "  " + str(nks) + '\n')
        else:
            conf.write(kpoints_bands[pp].split("!")[0].ljust(len(kpoints_bands[pp].split("!")[0])) + "  " + str(1) + '\n')
    conf.close()

def write_perturbo_temper(filepath, filename, temp, fermi, ne):
    conf = open(osp.join(filepath, filename), 'w')
    logical_ne = False 
    logical_temp = False
    logical_fermi = False
    len_dat = 1
    if '[' in str(ne):
        logical_ne = True 
        len_dat = len(ne)
    if '[' in str(temp):
        logical_temp = True
        len_dat = len(temp)
    if '[' in str(fermi):
        logical_fermi = True
        len_dat = len(fermi)
    # conf.write(str(len_dat) + "    F \n")
    conf.write(str(len_dat) + " \n")
    if logical_fermi:
        if logical_temp:
            for i in range(len(fermi)):
                conf.write(str(temp[i]) + '    ' + str(fermi[i]) + '    ' + f"{ne: .2e}" + '\n')
        else:
            for i in range(len(fermi)):
                conf.write(str(temp) + '    ' + str(fermi[i]) + '    ' + f"{ne: .2e}" + '\n')
    else:
        if logical_temp and logical_ne:
            for i in range(len(temp)):
                conf.write(str(temp[i]) + '    ' + str(fermi) + '    ' + f"{ne[i]: .2e}" + '\n')
        elif logical_temp and not logical_ne:
            for i in range(len(temp)):
                conf.write(str(temp[i]) + '    ' + str(fermi) + '    ' + f"{ne: .2e}" + '\n')
        elif not logical_temp and logical_ne:
            for i in range(len(ne)):
                conf.write(str(temp) + '    ' + str(fermi) + '    ' + f"{ne[i]: .2e}" + '\n')
        else:
            conf.write(str(temp) + '    ' + str(fermi) + '    ' + str(ne) + '\n')
    conf.close()

def write_perturbo_kpt(filepath, filename, kdim):
    conf = open(osp.join(filepath, filename), 'w')
    conf.write('\n')
    kpoints = inpp.kpt_nscf
    str_kpt = str(kpoints[0]) + " " + str(kpoints[1]) + " " + str(kpoints[2])
    if osp.exists(osp.join(filepath, "kpoints_kmesh")):
        os.system('rm kpoints_kmesh')
    os.system('kmesh.pl ' + str_kpt +' >> kpoints_kmesh')
    data = open(osp.join(filepath, "kpoints_kmesh"), 'r').readlines()
    conf.write(str(len(data)) + '  ' + 'crystal')
    count = 0
    for keys, vals in kdim.items():
        conf.write(' ' + str(kdim['boltz_kdim('+str(count+1) + ')']))
        count = count + 1
    conf.write('  ')
    conf.write('\n')
    for kkk in range(2, len(data)):
        conf.write(data[kkk])
    conf.close()

def write_perturbo_kpt_hsp(filepath, filename):
    conf = open(osp.join(filepath, filename), 'w')
    kpoints = inpp.kpoints_band.split("\n")[1:-1]
    conf.write(str(len(kpoints)) + '\n')
    kpoints = inpp.kpoints_band.split("\n")[1:-1]
    for i in range(len(kpoints)):
        # conf.write(kpoints[i].split("!")[0].ljust(len(kpoints[i].split("!")[0])) + ' ' + str(1/len(kpoints)) + '\n')
        conf.write(kpoints[i].split("!")[0].ljust(len(kpoints[i].split("!")[0])) + ' ' + str(1.0) + '\n')
    conf.close()

def write_kmesh_pl_homogenous(filepath, filename):
    conf = open(osp.join(filepath, filename), 'w')
    conf.write('\n')
    kpoints = inpp.kpt_nscf
    str_kpt = str(kpoints[0]) + " " + str(kpoints[1]) + " " + str(kpoints[2])
    if osp.exists(osp.join(filepath, "kpoints_kmesh")):
        os.system('rm kpoints_kmesh')
    os.system('kmesh.pl ' + str_kpt +' >> kpoints_kmesh')
    data = open(osp.join(filepath, "kpoints_kmesh"), 'r').readlines()
    for kkk in range(len(data)):
        conf.write(data[kkk])
    return conf

def mag_judge():
    '''
    mag:  
    1: ispin 2 and noncolin
    2: ispin 2
    3: ispin 1 and SOC 
    4: ispin 1
    '''
    mag = 4
    ############################# magmom ##################
    if inpp.SOC or inpp.noncolin:
        logical = True 
    else:
        logical = False
    if inpp.ISPIN == 2 and inpp.noncolin:
        mag = 1
    elif inpp.ISPIN == 2 and not logical:
        mag = 2
    elif inpp.ISPIN == 1 and inpp.SOC:
        mag = 3
    else:
        mag = 4
    ############################# magmom ##################
    return mag

def mag_set(system):
    mag = mag_judge()
    ############################# magmom ##################
    if mag == 1:
        del system['nspin']
        system['noncolin'] = True
        for key_mag, val_mag in inpp.Magmom.items():
            system['starting_magnetization'+"("+ key_mag + ")"] = val_mag[0]
            system['angle1(' + key_mag + ')'] = val_mag[1]
            system['angle2(' + key_mag + ')'] = val_mag[2]
    elif mag == 2:
        system['nspin'] = 2
        for key_mag, val_mag in inpp.Magmom.items():
            system['starting_magnetization'+"("+ key_mag + ")"] = val_mag[0]
    elif mag == 3:
        system['noncolin'] = True
        system['lspinorb'] = True
    else:
        system['nspin'] = 1
    ############################# magmom ##################
    return system


def superconductivity_TC(filename):    # Tc calculation.
    import json
    with open(osp.join(os.getcwd(), filename), 'w') as tc:
        tc.write('-- \n')
        input_epw_sup = {}
        input_epw_sup['epwread'] = True
        input_epw_sup['fila2f'] = 'pwscf.a2f_iso'
        input_epw_sup['liso'] = True 
        input_epw_sup['limag'] = True 
        input_epw_sup['lpade'] = False 
        input_epw_sup['lacon'] = False 
        input_epw_sup['tc_linear'] = True 
        input_epw_sup['tc_linear_solver'] = 'power' 
        # input_epw_sup['mp_mesh_k'] = True
        for key, val in input_epw_sup.items():
            tc.write(key.ljust(20, ' ') + '  =  ')
            json.dump(val, tc)
            tc.write('\n')
    tc.close()

def superconductivity_TC_aniso(filename):    # Tc calculation.
    import json
    with open(osp.join(os.getcwd(), filename), 'w') as tc:
        tc.write('-- \n')
        input_epw_sup = {}
        input_epw_sup['ep_coupling'] = True
        input_epw_sup['elph'] = True
        input_epw_sup['epwwrite'] = False 
        input_epw_sup['epwread'] = True 
        input_epw_sup['fermi_plot'] = False 
        input_epw_sup['wannierize'] = False 
        input_epw_sup['ephwrite'] = False 
        input_epw_sup['imag_read'] = True
        temp = inpp.epw_input['temps']
        # input_epw_sup['mp_mesh_k'] = True
        for key, val in input_epw_sup.items():
            tc.write(key.ljust(20, ' ') + '  =  ')
            json.dump(val, tc)
            tc.write('\n')
        insert = int((temp[1]-temp[0])/10)
        tc.write('temps   =  ')
        for i in range(10):
            tc.write(str(temp[0] + i * insert) + '  ')
        tc.write('\n')
    tc.close()

def KPATH_Phonopy(filepath, filename):
    data = inpp.kpoints_band.split("\n")[1:-1]
    kpath = []
    klabel = []
    nums = inpp.phon_points
    for i in range(len(data)):
        kpath.append([item for item in data[i].split('!')[0].split()])
        klabel.append(data[i].split('!')[1].strip())
    ff = open(osp.join(filepath, filename), 'w')
    ff.write("SYSTEM=Phonopy \n")
    ff.write("DIM = " + str(inpp.qpoints[0]) + ' ' + str(inpp.qpoints[1]) + 
             " " + str(inpp.qpoints[2]) +  " \n")
    ff.write("BAND = ")
    for i in range(len(kpath)):
        for j in range(3):
            ff.write(kpath[i][j] + ' ')
    ff.write('\n')
    ff.write('BAND_LABELS = ')
    for i in range(len(klabel)):
        if '\\' in klabel[i]: 
            ff.write('$'+klabel[i]+'$' + ' ')
        else:
            ff.write(klabel[i] + ' ')
    ff.write('\n')
    ff.write('FORCE_CONSTANTS = READ \n')
    ff.close()

def get_data_relaxationT(filepath, filename):
    temp = open(osp.join(filepath, 'pwscf.temper'), 'r').readlines()
    nums = int(temp[0].split()[0])
    print("Temperature: ", nums)
    rt = open(osp.join(filepath, filename), 'r').readlines()
    firstLine = rt[0]
    fname = ["relaxation_time_"+str(i+1) for i in range(nums)]
    num_lines = int((len(rt))/nums) + 1
    num_data = num_lines - 5 
    print("Line nums: ", num_lines, num_data)
    count = 1
    for i in range(nums):
        print("File: ", i+1)
        write_data = open(osp.join(filepath, fname[i]), 'w')
        write_data.write(firstLine)
        for j in range(count, count+num_data): 
            write_data.write(rt[j])
        count = count + num_lines
        write_data.close()

def get_data_SR(filepath):
    data = open(osp.join(filepath, "pwscf.imsigma_mode"), 'r').readlines()
    wr_dat = open(osp.join(filepath, "pwscf.imsigma_mode_data"), 'w')
    [nk, nbnd, nT, nmodes] = [int(item.split()[0]) for item in data[3].split(':')[1:]]
    wr_dat.write(data[3])
    for i in range(len(data)):
        if not "#" in data[i]:
            wr_dat.write(data[i])
    wr_dat.close()
    return nk, nbnd, nT, nmodes

def get_ScatteringRate(filepath, filename):
    get_data_SR(filepath)
    temp = open(osp.join(filepath, 'pwscf.temper'), 'r').readlines()
    nums = int(temp[0].split()[0])
    print("Temperature: ", nums)
    sr = open(osp.join(filepath, filename), 'r').readlines()
    firstLine = sr[0]
    fname = ["SR_"+str(i+1) for i in range(nums)]
    num_lines = int((len(sr))/nums)
    print("Line nums: ", num_lines)
    count = 1
    for i in range(nums):
        print("File: ", i+1)
        write_data = open(osp.join(filepath, fname[i]), 'w')
        write_data.write(firstLine)
        for j in range(count, count+num_lines): 
            write_data.write(sr[j])
        count = count + num_lines
        write_data.close()