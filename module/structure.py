# -*- coding: utf-8 -*-
"""
Created on 2:06 AM 28-11-2022 

@author: X-Y Ding
mail to: dxy_vasp@163.com
python3: structure.py
"""
import os, sys
import os.path as osp
import numpy as np
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
import funcs as bas
import AutoinputQE as inpp
sys.dont_write_bytecode = True 

def vc_relax(filepath=os.getcwd(), filename='vc-relax.log'):
    file = osp.join(filepath, filename)
    data = open(file, 'r').readlines()
    latt = []
    cor = []
    tmp = 0
    for i in range(len(data)):
        if 'Begin final coordinates' in data[i]:
            tmp = i+10
            tmp1 = np.array([float(item) for item in data[i+5].split()])
            tmp2 = np.array([float(item) for item in data[i+6].split()])
            tmp3 = np.array([float(item) for item in data[i+7].split()])
            latt.append([tmp1, tmp2, tmp3])
        elif 'End final coordinates' in data[i]:
            for j in range(tmp, i):
                cor.append([float(item) for item in data[j].split()[1:]])
            return np.around(np.array(latt).reshape(3, 3), 8), np.around(np.array(cor), 8)
        else:
            continue

def vc_relax_structure(filepath=os.getcwd(), filename='vc-relax.log'):
    file = osp.join(filepath, filename)
    data = open(file, 'r').readlines()
    latt = []
    cor = []
    len_cell_para = 0
    for j in range(len(data)):
        if 'CELL_PARAMETERS' in data[j]:
            len_cell_para = len_cell_para + 1
        else:
            continue
    count_i = 0
    tmp = 0
    print("Last step: ", len_cell_para)
    for i in range(len(data)):
        if 'CELL_PARAMETERS' in data[i]:
            if count_i == len_cell_para-1:
                tmp1 = np.array([float(item) for item in data[i+1].split()])
                tmp2 = np.array([float(item) for item in data[i+2].split()])
                tmp3 = np.array([float(item) for item in data[i+3].split()])
                latt.append(np.array([tmp1, tmp2, tmp3]))
                tmp = i + 6
                while(data[tmp].strip() !=''):
                    # print(data[tmp])
                    cor.append(np.array([float(item) for item in data[tmp].split()[1:]]))
                    tmp = tmp + 1
                return np.around(np.array(latt).reshape(3, 3), 8), np.around(np.array(cor), 8)
            else:
                count_i = count_i + 1
                continue
        else:
            continue
    
def relax(filepath=os.getcwd(), filename='relax.log'):
    file = osp.join(filepath, filename)
    data = open(file, 'r').readlines()
    cor = []
    tmp = 0
    for i in range(len(data)):
        if 'Begin final coordinates' in data[i]:
            tmp = i+3
        elif 'End final coordinates' in data[i]:
            for j in range(tmp, i):
                cor.append([float(item) for item in data[j].split()[1:]])
            return np.around(np.array(cor), 8)
        else:
            continue
def get_saclar_factor(filepath, filename):
    file = osp.join(filepath, filename)
    data = open(file, 'r').readlines()
    for i in range(len(data)):
        if "Begin final coordinates" in data[i]:
            try:
                if 'vc-relax.log' in filename and 'CELL_PARAMETERS' in data[i+4] and 'alat' in data[i+4]:
                    scalar_factor = float(data[i+4].split()[-1].split(')')[0])
                    return scalar_factor
                else:
                    return 1.0
            except:
                return 1.0
        elif 'vc-relax.log' in filename and 'CELL_PARAMETERS' in data[i] and 'alat' in data[i]:
                scalar_factor = float(data[i].split()[-1].split(')')[0])
                return scalar_factor

def write_contcar(filepath=os.getcwd(), filename='vc-relax.log'):
    latt, element, num_element, tot_ele, cor = bas.get_poscar(filepath, 'POSCAR')
    f = open(osp.join(filepath, "CONTCAR"), 'w')
    if filename == 'vc-relax.log':
        # latt_opted, cor_opted = vc_relax(filepath, filename)
        try:
            latt_opted, cor_opted = vc_relax(filepath, filename)
        except:
            latt_opted, cor_opted = vc_relax_structure(filepath, filename)
        if inpp.ibrav != 0:
            scalar = get_saclar_factor(filepath, filename)
            latt_opted = np.round(latt_opted * scalar * 0.52918, 5)
        f.write("system-opted \n")
        f.write("1.0 \n")
        for i in range(len(latt_opted)):
            f.write(str(latt_opted[i, 0]).rjust(15) + "  " + str(latt_opted[i, 1]).rjust(15) + "  " + str(latt_opted[i, 2]).rjust(15) + "\n")
        for k in range(len(element)):
            f.write(str(element[k]).ljust(3) + " ")
        f.write('\n')
        for k in range(len(element)):
            f.write(str(num_element[k]).ljust(3) + " ")
        f.write('\n')
        f.write("Direct \n")
        for j in range(len(cor_opted)):
            f.write(str(cor_opted[j, 0]).rjust(15) + "  " + str(cor_opted[j, 1]).rjust(15) + "  " + str(cor_opted[j, 2]).rjust(15) + "\n")  
    elif filename == 'relax.log':
        # cor_opted = relax(filepath, filename)
        try:
            latt_opted, cor_opted = vc_relax(filepath, filename)
        except:
            latt_opted, cor_opted = vc_relax_structure(filepath, filename)
        f.write("system-opted \n")
        f.write("1.0 \n")
        for i in range(len(latt)):
            f.write(str(latt[i, 0]).rjust(15) + "  " + str(latt[i, 1]).rjust(15) + "  " + str(latt[i, 2]).rjust(15) + "\n")
        for k in range(len(element)):
            f.write(str(element[k]).ljust(3) + " ")
        f.write('\n')
        for k in range(len(element)):
            f.write(str(num_element[k]).ljust(3) + " ")
        f.write('\n')
        f.write("Direct \n")
        for j in range(len(cor_opted)):
            f.write(str(cor_opted[j, 0]).rjust(15) + "  " + str(cor_opted[j, 1]).rjust(15) + "  " + str(cor_opted[j, 2]).rjust(15) + "\n")  
    else:
        print("No result files: vc-relax.log or relax.log")

def get_stru_from_QE_input(filepath, filename):
    from ase.io import read 
    structure = read(osp.join(filepath, filename), format='espresso-in')
    # 初始化一个字典来跟踪每个元素的数量
    element_counts = {}
    # 遍历每个原子，统计元素的数量
    for atom in structure:
        element = atom.symbol
        if element in element_counts:
            element_counts[element] += 1
        else:
            element_counts[element] = 1
    print("###############################################")
    cell = structure.get_cell()
    for vector in cell:
        print(vector)
    # 输出每个元素的原子个数
    f = ''
    for element, count in element_counts.items():
        f = f + atom.symbol
        print(f"Element: {element}, Count: {count}")
    print("Total atoms: ", len(structure))
    print("###############################################")
    for atom in structure:
        print(f"Atom: {atom.symbol}, Coordinates: {atom.scaled_position}")
    print("###############################################")
    with open(osp.join(filepath, 'str.poscar'), 'w') as stru:
        stru.write(f + '\n')
        stru.write(str(1.0) + '\n')
        for vector in cell:
            for i in range(len(vector)):
                stru.write("    " + str(np.round(vector[i], 12)).rjust(15, ' '))
            stru.write('\n')
        for element, count in element_counts.items():
            stru.write("    " + element.rjust(4, ' '))
        stru.write('\n')
        for element, count in element_counts.items():
            stru.write("    " + str(count).rjust(4, ' '))
        stru.write('\n')
        stru.write('Direct \n')
        for atom in structure:
            for j in range(len(atom.scaled_position)):
                stru.write("    " + str(np.round(atom.scaled_position[j], 12)).rjust(15, ' '))
            stru.write('\n')
    stru.close()

if __name__ == '__main__':
    import os.path as osp 
    if osp.exists(osp.join(os.getcwd(), 'vc-relax.log')):
        write_contcar(filepath=os.getcwd(), filename='vc-relax.log')
    else:
        write_contcar(filepath=os.getcwd(), filename='relax.log')
