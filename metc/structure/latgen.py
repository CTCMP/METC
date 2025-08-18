# -*- coding: utf-8 -*-
import os, sys
import os.path as osp
import numpy as np 
import utilities as util
from utilities import atom_mass

#######################################
prec = 5
Ang2bohr= 1.8897161646320724
#######################################


def reciprocal_cell(cell):
    a1 = np.array(cell[0])
    a2 = np.array(cell[1])
    a3 = np.array(cell[2])
    volume = np.dot(a1, np.cross(a2, a3))
    scalar = 2*np.pi/volume
    b1 = scalar * np.cross(a2, a3)
    b2 = scalar * np.cross(a3, a1)
    b3 = scalar * np.cross(a1, a2)
    return np.array([b1, b2, b3])

def get_poscar(filepath="", filename='POSCAR'):
    filepath = osp.join(filepath, filename)
    poscar1 = open(filepath, 'r')
    poscar_lines = poscar1.readlines()
    scale = float(poscar_lines[1].strip())
    pos = []
    for flag_lines in range(2, 5):
        for i in range(0, 3):
            pos.append(float(poscar_lines[flag_lines].split()[i]))
    latt = np.array(pos).reshape(3, 3) * scale * Ang2bohr
    ############# get element ###################
    poscar = open(filepath, 'r').read().strip('\n').splitlines()
    elements = poscar[5].lstrip().split()
    element = []
    masses = {}
    masses_ifc = []
    num_element = []
    numbers = poscar[6].lstrip().split()
    for flag_element in range(0, len(elements)):
        ele = elements[flag_element]
        element.append(ele)
        masses[ele] = atom_mass[ele]
        masses_ifc.append(atom_mass[ele] * util.QE_atom_mass_factor)
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
    
    mass = {}
    count = 0
    for ii in range(len(num_element)):
        for jj in range(num_element[ii]):
            count = count + 1
            mass[str(count)] = masses[element[ii]]
    return np.round(latt, 6), element, num_element, tot_ele, np.round(cor, 6), scale, mass, masses_ifc

def volume(a1, a2, a3):
    return np.dot(a1, np.cross(a2, a3))

def latgen_lib(ibrav, celldm, cell):
    # -----------------------------------------------------------------------
    #      sets up the crystallographic vectors a1, a2, and a3.
    # 
    #      ibrav is the structure index:
    #        1  cubic P (sc)                8  orthorhombic P
    #        2  cubic F (fcc)               9  1-face (C) centered orthorhombic
    #        3  cubic I (bcc)              10  all face centered orthorhombic
    #        4  hexagonal and trigonal P   11  body centered orthorhombic
    #        5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
    #        6  tetragonal P (st)          13  one face (base) centered monoclinic
    #        7  tetragonal I (bct)         14  triclinic P
    #      Also accepted:
    #        0  "free" structure          -12  monoclinic P (unique axis: b)
    #       -3  cubic bcc with a more symmetric choice of axis
    #       -5  trigonal R, threefold axis along (111)
    #       -9  alternate description for base centered orthorhombic
    #      -13  one face (base) centered monoclinic (unique axis: b)
    #       91  1-face (A) centered orthorombic
    # 
    #      celldm are parameters which fix the shape of the unit cell
    #      omega is the unit-cell volume
    # 
    #      NOTA BENE: all axis sets are right-handed
    #      Boxes for US PPs do not work properly with left-handed axis
    a1 = np.array(cell[0])
    a2 = np.array(cell[1])
    a3 = np.array(cell[2])
    sr2 = np.sqrt(2)
    sr3 = np.sqrt(3)
    Omega = 0

    if ibrav == 0:
        if np.sqrt( a1[0]**2 + a1[1]**2 + a1[2]**2 ) == 0:
            print('wrong at for ibrav=0')
        if np.sqrt( a2[0]**2 + a2[1]**2 + a2[2]**2 ) == 0:
            print('wrong at for ibrav=0')
        if np.sqrt( a3[0]**2 + a3[1]**2 + a3[2]**2 ) == 0:
            print('wrong at for ibrav=0')

        if celldm[0] != 1 :
            # ... input at are in units of alat => convert them to a.u.
            a1 = a1 * celldm[0]
            a2 = a2 * celldm[0]
            a3 = a3 * celldm[0]
        else:
            # ... input at are in atomic units: define celldm(1) from a1
            celldm[0] = np.sqrt(a1[0]**2 + a1[1]**2 + a1[2]**2)
    else:
        a1[:] = 0
        a2[:] = 0
        a3[:] = 0
    if celldm[0] <= 0:
        print('wrong celldm(1) !')

    if ibrav == 1:
        # simple cubic lattice
        a1[0] = celldm[0]
        a2[1] = celldm[0]
        a3[2] = celldm[0]
    elif ibrav == 2:
        # fcc lattice
        term = celldm[0] / 2
        a1[0] = -term 
        a1[2] = term 
        a2[1] = term 
        a2[2] = term
        a3[0] = -term 
        a3[1] = term 
    elif np.abs(ibrav) == 3:
        # bcc lattice
        term = celldm[0] / 2.0 
        for i in range(1, 4):
            a1[i] = term 
            a2[i] = term 
            a3[i] = term 
        if ibrav < 0: 
            a1[0] = -a1[0]
            a2[1] = -a2[1]
            a3[2] = -a3[2]
        else:
            a2[0] = -a2[0]
            a3[0] = -a3[0]
            a3[1] = -a3[1]
    elif ibrav == 4:
        # hexagonal lattice
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        cbya = celldm[2]
        a1[0] = celldm[0]
        a2[0] = -celldm[0] / 2.0
        a2[1] = celldm[0] * sr3 / 2.0
        a3[2] = celldm[0] * cbya
    elif np.abs(ibrav) == 5:
        # trigonal lattice
        if celldm[3] <= -0.5 or celldm[3] >= 1.0:
            print('wrong celldm(4)')
        term1 = np.sqrt(1.0 + 2.0 * celldm[3])
        term2 = np.sqrt(1.0 - celldm[3])
        if ibrav == 5:
            # threefold axis along c (001)
            a2[1]=sr2*celldm[0]*term2/sr3
            a2[2]=celldm[0]*term1/sr3
            a1[0]=celldm[0]*term2/sr2
            a1[1]=-a1[0]/sr3
            a1[2]= a2[2]
            a3[0]=-a1[0]
            a3[1]= a1[1]
            a3[2]= a2[2]
        elif ibrav == -5:
            '''
                threefold axis along (111)
            Notice that in the cubic limit (alpha=90, celldm(4)=0, term1=term2=1)
            does not yield the x,y,z axis, but an equivalent rotated triplet:
                a/3 (-1,2,2), a/3 (2,-1,2), a/3 (2,2,-1)
            If you prefer the x,y,z axis as cubic limit, you should modify the
            definitions of a1(1) and a1(2) as follows:'
                a1(1) = celldm(1)*(term1+2.0_dp*term2)/3.0_dp
                a1(2) = celldm(1)*(term1-term2)/3.0_dp
            (info by G. Pizzi and A. Cepellotti)
            '''
            a1[0] = celldm[0]*(term1-2.0*term2)/3.0
            a1[1] = celldm[0]*(term1+term2)/3.0
            a1[2] = a1[1]
            a2[0] = a1[2]
            a2[1] = a1[0]
            a2[2] = a1[1]
            a3[0] = a1[1]
            a3[1] = a1[2]
            a3[2] = a1[0]
    elif ibrav == 6:
        # tetragonal lattice
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        cbya=celldm[2]
        a1[0]=celldm[0]
        a2[1]=celldm[0]
        a3[2]=celldm[0]*cbya  
    elif ibrav == 7:
        # body centered tetragonal lattice
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        cbya=celldm[2]
        a2[0]=celldm[0]/2.0
        a2[1]=a2[0]
        a2[2]=cbya*celldm[0]/2.0
        a1[0]= a2[0]
        a1[1]=-a2[0]
        a1[2]= a2[2]
        a3[0]=-a2[0]
        a3[1]=-a2[0]
        a3[2]= a2[2] 
    elif ibrav == 8:
        # Simple orthorhombic lattice
        if celldm[1] <= 0:
            print('wrong celldm(2)')
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        a1[0]=celldm[0]
        a2[1]=celldm[0]*celldm[1]
        a3[2]=celldm[0]*celldm[2]
    elif np.abs(ibrav) == 9:
        # One face (base) centered orthorhombic lattice  (C type)
        if celldm[1] <= 0:
            print('wrong celldm(2)')
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        if ibrav == 9:
            a1[0] = 0.5 * celldm[0]
            a1[1] = a1[0] * celldm[1]
            a2[0] = -a1[0]
            a2[1] = a1[1]
        else:
            # alternate description
            a1[0] = 0.5 * celldm[0]
            a1[1] = -a1[0] * celldm[1]
            a2[0] = a1[0]
            a2[1] = -a1[1] 
        a3[2] = celldm[0] * celldm[2]
    elif ibrav == 91:
        # One face (base) centered orthorhombic lattice  (A type)
        if celldm[1] <= 0:
            print('wrong celldm(2)')
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        a1[0] = celldm[0]
        a2[1] = celldm[0] * celldm[1] * 0.5
        a2[2] = - celldm[0] * celldm[2] * 0.5
        a3[1] = a2[1]
        a3[2] = - a2[2]
    elif ibrav == 10:
        # All face centered orthorhombic lattice
        if celldm[1] <= 0:
            print('wrong celldm(2)')
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        a2[0] = 0.5 * celldm[0]
        a2[1] = a2[0] * celldm[1]
        a1[0] = a2[0]
        a1[2] = a2[0] * celldm[2]
        a3[1] = a2[0] * celldm[1]
        a3[2] = a1[2]
    elif ibrav == 11:
        #  Body centered orthorhombic lattice
        if celldm[1] <= 0:
            print('wrong celldm(2)')
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        a1[0] = 0.5 * celldm[0]
        a1[1] = a1[0] * celldm[1]
        a1[2] = a1[0] * celldm[2]
        a2[0] = - a1[0]
        a2[1] = a1[1]
        a2[2] = a1[2]
        a3[0] = - a1[0]
        a3[1] = - a1[1]
        a3[2] = a1[2]
    elif ibrav == 12:
        # Simple monoclinic lattice, unique (i.e. orthogonal to a) axis: c
        if celldm[1] <= 0:
            print('wrong celldm(2)')
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        if np.abs(celldm[3]) >= 1.0:
            print('wrong celldm(4)')
        sen=np.sqrt(1.0-celldm[3]**2)
        a1[0]=celldm[0]
        a2[0]=celldm[0]*celldm[1]*celldm[3]
        a2[1]=celldm[0]*celldm[1]*sen
        a3[2]=celldm[0]*celldm[2]
    elif ibrav == -12:
        # Simple monoclinic lattice, unique axis: b (more common)
        if celldm[1] <= 0:
            print('wrong celldm(2)')
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        if np.abs(celldm[4]) >= 1.0:
            print('wrong celldm(5)')
        sen=np.sqrt(1.0-celldm[4]**2)
        a1[0]=celldm[0]
        a2[1]=celldm[0]*celldm[1]
        a3[0]=celldm[0]*celldm[2]*celldm[4]
        a3[2]=celldm[0]*celldm[2]*sen
    elif ibrav == 13:
        # One face centered monoclinic lattice unique axis c
        if celldm[1] <= 0:
            print('wrong celldm(2)')
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        if np.abs(celldm[3]) >= 1.0:
            print('wrong celldm(4)')
        sen = sqrt( 1.0 - celldm[3] ** 2 )
        a1[0] = 0.5 * celldm[0]
        a1[2] =-a1[0] * celldm[2]
        a2[0] = celldm[0] * celldm[1] * celldm[3]
        a2[1] = celldm[0] * celldm[1] * sen
        a3[0] = a1[0]
        a3[2] =-a1[2]
    elif ibrav == -13:
        # One face centered monoclinic lattice unique axis b
        if celldm[1] <= 0:
            print('wrong celldm(2)')
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        if np.abs(celldm[4]) >= 1.0:
            print('wrong celldm(5)')
        sen = np.sqrt( 1.0 - celldm[4] ** 2 )
        a1[0] = 0.5 * celldm[0]
        a1[1] = a1[0] * celldm[1]
        a2[0] =-a1[0]
        a2[1] = a1[1]
        a3[0] = celldm[0] * celldm[2] * celldm[4]
        a3[2] = celldm[0] * celldm[2] * sen
    elif ibrav == 14:
        # Triclinic lattice
        if celldm[1] <= 0:
            print('wrong celldm(2)')
        if celldm[2] <= 0:
            print('wrong celldm(3)')
        if np.abs(celldm[3]) >= 1.0:
            print('wrong celldm(4)')
        if np.abs(celldm[4]) >= 1.0:
            print('wrong celldm(5)')
        if np.abs(celldm[5]) >= 1.0:
            print('wrong celldm(6)')
        singam=np.sqrt(1.0-celldm[5]**2)
        term= (1.0 + 2.0 * celldm[3]*celldm[4]*celldm[5] - celldm[3]**2-celldm[4]**2-celldm[5]**2)
        if (term < 0.0):
            print('celldm do not make sense, check your data')
        term= np.sqrt(term/(1.0-celldm[5]**2))
        a1[0]=celldm[0]
        a2[0]=celldm[0]*celldm[1]*celldm[5]
        a2[1]=celldm[0]*celldm[1]*singam
        a3[0]=celldm[0]*celldm[2]*celldm[4]
        a3[1]=celldm[0]*celldm[2]*(celldm[3]-celldm[4]*celldm[5])/singam
        a3[2]=celldm[0]*celldm[2]*term
    elif ibrav != 0:
        print('nonexistent bravais lattice')
    omega = volume(a1, a2, a3)
    return ibrav, celldm, a1,a2,a3, omega
    
def eqq(a, b):
    return np.round(a, prec) == np.round(b, prec)

def at2celldm(ibrav, alat, cell):
    a1 = np.array(cell[0])
    a2 = np.array(cell[1])
    a3 = np.array(cell[2])
    celldm = np.zeros(6)
    if ibrav == 0:
        celldm[0] = 1.0
        # celldm[0] = np.sqrt(np.dot(a1, a1))
    elif ibrav == 1:
        celldm[0] = np.sqrt(np.dot(a1, a1))
    elif ibrav == 2:
        celldm[0] = np.sqrt(np.dot(a1, a1) * 2.0)
    elif ibrav == 3 or ibrav == -3:
        celldm[0] = np.sqrt(np.dot(a1, a1)/3) * 2
    elif ibrav == 4:
        celldm[0] = np.sqrt(np.dot(a1, a1))
        celldm[2] = np.sqrt(np.dot(a3, a3)) / celldm[0]
    elif ibrav == -5 or ibrav == 5:
        celldm[0] = np.sqrt(np.dot(a1, a1))
        celldm[2] = np.dot(a1, a2) / celldm[0] / np.sqrt(np.dot(a2, a2))
    elif ibrav == 6:
        celldm[0] = np.sqrt(np.dot(a1, a1))
        celldm[2] = np.sqrt(np.dot(a3, a3)) / celldm[0]
    elif ibrav == 7:
        celldm[0] = np.abs(a1[0]) * 2.0
        celldm[2] = np.abs(a1[2]/a1[0])
    elif ibrav == 8:
        celldm[0] = np.sqrt(np.dot(a1, a1))
        celldm[1] = np.sqrt(np.dot(a2, a2))/celldm[0]
        celldm[2] = np.sqrt(np.dot(a3, a3))/celldm[0]
    elif ibrav == 9 or ibrav == -9:
        celldm[0] = np.abs(a1[0]) * 2.0
        celldm[1] = np.abs(a2[1]) * 2.0 / celldm[0]
        celldm[2] = np.abs(a3[2]) / celldm[0]
    elif ibrav == 91:
        celldm[0] = np.sqrt(np.dot(a1,a1))
        celldm[1] = abs(a2[1])*2.0 / celldm[0]
        celldm[2] = abs(a3[2])*2.0 / celldm[0]
    elif ibrav == 10:
        celldm[0] = np.abs(a1[0])*2
        celldm[1] = np.abs(a2[1])*2.0/celldm[0]
        celldm[2] = np.abs(a3[2])*2.0/celldm[0]
    elif ibrav == 11:
        celldm[0] = abs(a1[0])*2.0
        celldm[1] = abs(a1[1])*2.0/celldm[0]
        celldm[2] = abs(a1[2])*2.0/celldm[0]
    elif ibrav == 12 or ibrav == -12:
        celldm[0] = np.sqrt( np.dot(a1,a1) )
        celldm[1] = np.sqrt( np.dot(a2,a2) ) / celldm[0]
        celldm[2] = np.sqrt( np.dot(a3,a3) ) / celldm[0]
        if ibrav == 12:
            celldm[3] = np.dot(a1,a2) / celldm[0] / np.sqrt(np.dot(a2,a2))
        else:
            celldm[4] = np.dot(a1,a3) / celldm[0] / np.sqrt(np.dot(a3,a3)) 
    elif ibrav == 13:
        celldm[0] = np.abs(a1[0])*2.0
        celldm[1] = np.sqrt(np.dot(a2, a2)) / celldm[0]
        celldm[2] = np.abs(a1[2]/a1[0])
        celldm[3] = a2[0]/a1[0]/celldm[1]/2.0
    elif ibrav == -13:
        celldm[0] = np.abs(a1[0])*2.0
        celldm[1] = np.abs(a2[1]/a2[0])
        celldm[2] = np.sqrt(np.dot(a3,a3)) / celldm[0]
        celldm[4] = a3[0]/a1[0]/celldm[2]/2.0
    elif ibrav == 14:
        celldm[0] = np.sqrt(np.dot(a1,a1))
        celldm[1] = np.sqrt( np.dot(a2,a2)) / celldm[0]
        celldm[2] = np.sqrt( np.dot(a3,a3)) / celldm[0]
        celldm[3] = np.dot(a3,a2)/np.sqrt(np.dot(a2,a2) * np.dot(a3,a3))
        celldm[4] = np.dot(a3,a1) / celldm[0] / np.sqrt( np.dot(a3,a3))
        celldm[5] = np.dot(a1,a2) / celldm[0] / np.sqrt( np.dot(a2,a2))
    else:
        print('wrong ibrav ?')
    celldm[0] = celldm[0] * alat
    return celldm 

def at2ibrav(cell):
    '''
    Returns ibrav from lattice vectors if recognized, 0 otherwise
    '''
    a1 = np.array(cell[0])
    a2 = np.array(cell[1])
    a3 = np.array(cell[2])
    v1 = np.sqrt(np.dot(a1, a1))
    v2 = np.sqrt(np.dot(a2, a2))
    v3 = np.sqrt(np.dot(a3, a3))
    cosbc = np.dot(a2, a3)/v2/v3 
    cosac = np.dot(a1, a3)/v1/v3 
    cosab = np.dot(a1, a2)/v1/v2 
    ibrav = 14
    if eqq(v1, v2) and eqq(v1, v3):
        # Case: a=b=c
        if eqq(cosab, cosac) and eqq(cosab, cosbc):
            # Case: alpha = beta = gamma
            if eqq(cosab, 0):
                # Cubic P - ibrav=1
                ibrav = 1
            elif eqq(cosab, 0.5):
                # Cubic F - ibrav=2
                ibrav = 2
            elif eqq(cosab, -1/3):
                # Cubic I - ibrav=-3
                ibrav = -3
            else:
                if eqq(np.abs(a1[2]), np.abs(a2[2])) and eqq(np.abs(a2[2]), np.abs(a3[2])):
                    # Trigonal 001 axis
                    ibrav = 5
                else:
                    # Trigonal, 111 axis
                    ibrav = -5 
        elif eqq(cosab, cosac) and not eqq(cosab, cosbc):
            if eqq(np.abs(a1[0]), np.abs(a1[1])) and eqq(np.abs(a2[0]), np.abs(a2[1])):
                # Tetragonal I
                ibrav = 7
            else:
                # Cubic I - ibrav=3
                ibrav = 3
        elif eqq(cosab, -cosac) and eqq(cosab, cosbc) and eqq(cosab, 1/3):
            # Cubic I - ibrav=3
            ibrav = 3
        elif eqq(np.abs(a1[0]), np.abs(a2[0])) and eqq(np.abs(a1[1]), np.abs(a2[1])):
            ibrav = 11

    elif eqq(v1, v2) and not eqq(v1, v3):
        # Case: a=b/=c
        if eqq(cosab, 0) and eqq(cosac, 0) and eqq(cosbc, 0):
            # Case: alpha = beta = gamma = 90
            # Simple tetragonal
            ibrav = 6
        elif eqq(cosab, -0.5) and eqq(cosac, 0) and eqq(cosbc, 0):
            # Case: alpha = 120, beta = gamma = 90 => simple hexagonal
            # Simple hexagonal
            ibrav = 4
        elif eqq(cosac, 0) and eqq(cosbc, 0):
            # Orthorhombic bco
            if eqq(a1[0], a2[0]) and eqq(a1[1], -a2[1]):
                ibrav = -9
            elif eqq(a1[0], -a2[0]) and eqq(a1[1], a2[1]):
                ibrav = 9
        elif eqq(cosac, -cosbc):
            # bco (unique axis b)
            ibrav = -13

    elif eqq(v1, v3) and not eqq(v1, v2):
        # Case: a=c/=b
        # Monoclinic bco (unique axis c)
        ibrav = 13

    elif eqq(v2, v3) and not eqq(v1, v2):
        # Case: a/=b=c
        # Orthorhombic 1-face bco
        ibrav = 91
    
    elif not eqq(v1, v2) and not eqq(v1, v3) and not eqq(v2, v3):
        # Case: a/=b/=c
        if eqq(cosab, 0) and eqq(cosac, 0) and eqq(cosbc, 0):
            # Case: alpha = beta = gamma = 90
            # Orthorhombic P
            ibrav = 8
        elif not eqq(cosab, 0) and eqq(cosac, 0) and eqq(cosbc, 0):
            # Case: alpha /= 90,  beta = gamma = 90
            # Monoclinic P, unique axis c
            ibrav = 12
        elif eqq(cosab, 0) and not eqq(cosac, 0) and eqq(cosbc, 0):
            # Case: beta /= 90, alpha = gamma = 90
            # Monoclinic P, unique axis b
            ibrav = -12
        elif not eqq(cosab, 0) and not eqq(cosac, 0) and eqq(cosbc, 0):
            # Case: beta /= 90, alpha = gamma = 90
            if eqq(np.abs(a1[0]), np.abs(a2[0])) and eqq(np.abs(a1[2]), np.abs(a3[2])) and eqq(np.abs(a2[1]), np.abs(a3[1])):
                # Orthorhombic F
                ibrav = 10
            else:
                ibrav = 14
    return ibrav 

def abc2celldm(ibrav, a, b, c, cosab, cosac, cosbc):
    celldm = np.zeros(6)
    if (a <= 0.0):
        print('abc2celldm',' incorrect lattice parameter (a)')
    if (b <  0.0):
        print('abc2celldm',' incorrect lattice parameter (b)')
    if (c <  0.0):
        print('abc2celldm',' incorrect lattice parameter (c)')
    if ( abs (cosab) > 1.0):
        print('abc2celldm', ' incorrect lattice parameter (cosab)')
    if ( abs (cosac) > 1.0):
        print('abc2celldm', '  incorrect lattice parameter (cosac)')
    if ( abs (cosbc) > 1.0):
        print('abc2celldm', ' incorrect lattice parameter (cosbc)')
    
    celldm[0] = a / BOHR_RADIUS_ANGS
    celldm[1] = b/a 
    celldm[2] = c/a 

    if ibarv == 14 or ibrav == 0:
        celldm[3] = cosbc 
        celldm[4] = cosac 
        celldm[5] = cosab 
    elif ibrav == -12 or ibrav == -13:
        celldm[3] = 0 
        celldm[4] = cosac 
        celldm[5] = 0
    elif ibarv == -5 or ibrav == 5 or ibrav == 12 or ibrav == 13:
        celldm[3] = cosab
        celldm[4] = 0
        celldm[5] = 0
    else:
        celldm[3] = 0
        celldm[4] = 0
        celldm[5] = 0 
    return celldm

def celldm2abc(ibrav, celldm, a, b, c):
    a = celldm[0] * BOHR_RADIUS_ANGS
    b = celldm[0] * celldm[1] * BOHR_RADIUS_ANGS
    c = celldm[0] * celldm[2] * BOHR_RADIUS_ANGS
    cosab = 0 
    cosac = 0 
    cosbc = 0
    if ibrav == 14 or ibrav == 0:
        cosbc = celldm[3]
        cosac = celldm[4]
        cosab = celldm[5]
    elif ibrav == -12 or ibrav == -13:
        cosab = 0 
        cosac = celldm[0]
        cosbc = 0
    elif ibrav == -5 or ibrav == 5 or ibrav == 12 or ibrav == 13:
        cosab = celldm[3]
        cosac = 0
        cosbc = 0
    else:
        cosab = 0
        cosac = 0
        cosbc = 0
    return a, b, c, cosab, cosac, cosbc

def remake_cell(ibrav, alat, cell):
    a1 = np.array(cell[0])
    a2 = np.array(cell[1])
    a3 = np.array(cell[2])
    if ibrav == 0:
        print('WARNING! With ibrav=0, cell_dofree= ibrav has no effect. ')
        # alat = np.sqrt(np.dot(a1, a1))
    celldm_internal = at2celldm(ibrav, alat, cell)
    e1 = np.array(a1) 
    e2 = np.array(a2) 
    e3 = np.array(a3) 
    ibrav, celldm_internal, a1, a2, a3, omega = latgen_lib(ibrav, celldm_internal, [a1, a2, a3])
    print("============================= structure info ===========================")
    print(f"ibrav = {ibrav}")
    print(f"celldm: {celldm_internal[0]}")
    if celldm_internal[1] != 0 :
        print(f'celldm(2) = {celldm_internal[1]}')
    if celldm_internal[2] != 0 :
        print(f'celldm(3) = {celldm_internal[2]}')
    if celldm_internal[3] != 0 :
        print(f'celldm(4) = {celldm_internal[3]}')
    if celldm_internal[4] != 0 :
        print(f'celldm(5) = {celldm_internal[4]}')
    if celldm_internal[5] != 0 :
        print(f'celldm(6) = {celldm_internal[5]}')
    print('-------------------------------------------------------------')
    print('Input lattice vectors:')
    print(e1)
    print(e2)
    print(e3)
    print('-------------------------------------------------------------')
    print('New lattice vectors in INITIAL alat:')
    print(a1/alat)
    print(a2/alat)
    print(a3/alat)
    print('omega: ', omega)
    print('-------------------------------------------------------------')
    print('lattice vectors of reciprocal axes:')
    rep = reciprocal_cell(np.array([e1, e2, e3]))
    print(rep[0])
    print(rep[1])
    print(rep[2])
    print('-------------------------------------------------------------')
    print('========================================================================')
    print('New lattice vectors in NEW alat (for information only):')
    print(a1/celldm_internal[0])
    print(a2/celldm_internal[0])
    print(a3/celldm_internal[0])
    a1 = a1/celldm_internal[0]
    a2 = a2/celldm_internal[0]
    a3 = a3/celldm_internal[0]
    print('-------------------------------------------------------------')
    print(f'(Discrepancy in bohr = {np.sqrt(np.sum(a1-e1)**2)}  {np.sqrt(np.sum(a2-e2)**2)}  {np.sqrt(np.sum(a3-e3)**2)})', )
    print('-------------------------------------------------------------')
    print('New lattice vectors of reciprocal axes (in units 2 pi/alat):')
    rep = reciprocal_cell(np.array([a1, a2, a3])) / 2 / np.pi
    print(rep[0])
    print(rep[1])
    print(rep[2])
    print('-------------------------------------------------------------')
    print('========================================================================')
    new_alat = celldm_internal[0]

def re_cell(ibrav, alat, cell):
    a1 = np.array(cell[0])
    a2 = np.array(cell[1])
    a3 = np.array(cell[2])
    celldm_internal = at2celldm(ibrav, alat, cell)
    e1 = np.array(a1) 
    e2 = np.array(a2) 
    e3 = np.array(a3) 
    ibrav, celldm_internal, a1, a2, a3, omega = latgen_lib(ibrav, celldm_internal, [a1, a2, a3])
    a11 = a1/celldm_internal[0]
    a22 = a2/celldm_internal[0]
    a33 = a3/celldm_internal[0]
    rep = reciprocal_cell(np.array([a11, a22, a33])) / 2 / np.pi
    new_alat = celldm_internal[0]
    return np.array([a11, a22, a33]), new_alat 

class QElattice:
    """ Class to get strucutre info """
    def __init__(self, symbols=None, ntype=None, numbers=None, natoms=None, alat=None, masses=None,
                mass_ifc=None, scale=None, ibrav=None, celldm=None, cell=None, cell_QE=None, cor=None, cor_cart=None, bg_QE=None, bg=None):
        self._symbols = symbols 
        self._numbers = numbers 
        self._masses = masses
        self._ibrav = ibrav
        self._celldm = celldm
        self._scale = scale
        self._cell = cell 
        self._cell_QE = cell_QE
        self._coordinate = cor
        self._cor_cart = cor_cart 
        self._natoms = natoms
        self._ntype = ntype
        self._mass_ifc = mass_ifc
        self._bg = bg
        self._bg_QE = bg_QE 
        self._alat = alat
        self._set_parameter()

    def _set_parameter(self):
        latt, element, num_element, tot_ele, cor, scale, mass, masses_ifc = \
            get_poscar(os.getcwd(), 'POSCAR')
        if self.ibrav != 0:
            self.ibrav = at2ibrav(latt)
            celldm = at2celldm(self.ibrav, 1.0, latt)
            ########
            alat_tmp = np.sqrt(np.dot(latt[0], latt[0]))
            cell_new, new_alat = re_cell(self.ibrav, 1.0, latt)
            self.alat = new_alat
            ########
            self.celldm = celldm
            self.scale = scale 
            self.coordinate = cor 
            self.cor_cart = np.dot(cor, latt)
            self.symbols = element
            self.numbers = num_element 
            self.masses = mass
            self.natoms = len(cor)
            self.ntype = len(num_element)
            self.mass_ifc = masses_ifc
            self.cell = latt 
            self.cell_QE = cell_new
            self.bg = reciprocal_cell(latt)
            self.bg_QE = reciprocal_cell(self.cell_QE) / 2 / np.pi
        else:
            self.ibrav = 0
            a1 = np.array(latt[0])
            alat_tmp = np.sqrt(np.dot(a1, a1))
            celldm = at2celldm(self.ibrav, alat_tmp, latt)
            ########
            cell_new, new_alat = re_cell(self.ibrav, 1.0, latt)
            print("new alat: ", new_alat)
            self.alat = new_alat
            ########
            self.celldm = celldm
            self.scale = scale 
            self.coordinate = cor 
            self.cor_cart = np.dot(cor, latt)
            self.symbols = element
            self.numbers = num_element 
            self.masses = mass
            self.natoms = len(cor)
            self.ntype = len(num_element)
            self.mass_ifc = masses_ifc
            self.cell = latt 
            self.cell_QE = cell_new
            self.bg = reciprocal_cell(latt)
            self.bg_QE = reciprocal_cell(self.cell_QE) / 2 / np.pi

    def print_structure_info(self):
        print("============================= structure info ===========================")
        print(f"ibrav = {self.ibrav}")
        print(f"celldm: {np.round(self.celldm[0], 8):< 12.8f}  {np.round(self.celldm[1], 8):< 12.8f}  {np.round(self.celldm[2], 8):< 12.8f}  {np.round(self.celldm[3], 8):< 12.8f}  {np.round(self.celldm[4], 8):< 12.8f}  {np.round(self.celldm[5], 8):< 12.8f}")
        print('alat: ', self.alat)
        print('-------------------------------------------------------------')
        print('cell:')
        print(f' {np.round(self.cell[0, 0], 8):< 12.8f}   {np.round(self.cell[0, 1], 8):< 12.8f}   {np.round(self.cell[0, 2], 8):< 12.8f}')
        print(f' {np.round(self.cell[1, 0], 8):< 12.8f}   {np.round(self.cell[1, 1], 8):< 12.8f}   {np.round(self.cell[1, 2], 8):< 12.8f}')
        print(f' {np.round(self.cell[2, 0], 8):< 12.8f}   {np.round(self.cell[2, 1], 8):< 12.8f}   {np.round(self.cell[2, 2], 8):< 12.8f}')
        print('-------------------------------------------------------------')
        print('reciprocal cell:')
        print(f' {np.round(self.bg[0, 0], 8):< 12.8f}   {np.round(self.bg[0, 1], 8):< 12.8f}   {np.round(self.bg[0, 2], 8):< 12.8f}')
        print(f' {np.round(self.bg[1, 0], 8):< 12.8f}   {np.round(self.bg[1, 1], 8):< 12.8f}   {np.round(self.bg[1, 2], 8):< 12.8f}')
        print(f' {np.round(self.bg[2, 0], 8):< 12.8f}   {np.round(self.bg[2, 1], 8):< 12.8f}   {np.round(self.bg[2, 2], 8):< 12.8f}')
        print('========================================================================')
        print('cart. coord. in units of alat:')
        print(f' {np.round(self.cell_QE[0, 0], 8):< 12.8f}   {np.round(self.cell_QE[0, 1], 8):< 12.8f}   {np.round(self.cell_QE[0, 2], 8):< 12.8f}')
        print(f' {np.round(self.cell_QE[1, 0], 8):< 12.8f}   {np.round(self.cell_QE[1, 1], 8):< 12.8f}   {np.round(self.cell_QE[1, 2], 8):< 12.8f}')
        print(f' {np.round(self.cell_QE[2, 0], 8):< 12.8f}   {np.round(self.cell_QE[2, 1], 8):< 12.8f}   {np.round(self.cell_QE[2, 2], 8):< 12.8f}')
        print('-------------------------------------------------------------')
        print('cart. coord. in units 2 pi/alat:')
        print(f' {np.round(self.bg_QE[0, 0], 8):< 12.8f}   {np.round(self.bg_QE[0, 1], 8):< 12.8f}   {np.round(self.bg_QE[0, 2], 8):< 12.8f}')
        print(f' {np.round(self.bg_QE[1, 0], 8):< 12.8f}   {np.round(self.bg_QE[1, 1], 8):< 12.8f}   {np.round(self.bg_QE[1, 2], 8):< 12.8f}')
        print(f' {np.round(self.bg_QE[2, 0], 8):< 12.8f}   {np.round(self.bg_QE[2, 1], 8):< 12.8f}   {np.round(self.bg_QE[2, 2], 8):< 12.8f}')
        print('========================================================================')
    
    @property
    def bg(self):
        return self._bg
    @bg.setter
    def bg(self, bg):
        self._bg = bg

    @property
    def alat(self):
        return self._alat
    @alat.setter
    def alat(self, alat):
        self._alat = alat

    @property
    def bg_QE(self):
        return self._bg_QE
    @bg_QE.setter
    def bg_QE(self, bg_QE):
        self._bg_QE = bg_QE

    @property
    def cor_cart(self):
        return self._cor_cart
    @cor_cart.setter
    def cor_cart(self, cor_cart):
        self._cor_cart = cor_cart

    @property
    def cell_QE(self):
        return self._cell_QE
    @cell_QE.setter
    def cell_QE(self, cell_QE):
        self._cell_QE = cell_QE

    @property
    def scale(self):
        return self._scale
    @scale.setter
    def scale(self, scale):
        self._scale = scale

    @property
    def mass_ifc(self):
        return self._mass_ifc
    @mass_ifc.setter
    def mass_ifc(self, mass_ifc):
        self._mass_ifc = mass_ifc

    @property
    def ntype(self):
        return self._ntype
    @ntype.setter
    def ntype(self, ntype):
        self._ntype = ntype

    @property
    def natoms(self):
        return self._natoms
    @natoms.setter
    def natoms(self, natoms):
        self._natoms = natoms

    @property
    def coordinate(self):
        return self._coordinate
    @coordinate.setter
    def coordinate(self, coordinate):
        self._coordinate = coordinate

    @property
    def celldm(self):
        return self._celldm
    @celldm.setter
    def celldm(self, celldm):
        self._celldm = celldm

    @property
    def ibrav(self):
        return self._ibrav
    @ibrav.setter
    def ibrav(self, ibrav):
        self._ibrav = ibrav

    @property
    def masses(self):
        return self._masses
    @masses.setter
    def masses(self, masses):
        self._masses = masses

    @property
    def numbers(self):
        return self._numbers 
    @numbers.setter
    def numbers(self, numbers):
        self._numbers = numbers 
    
    @property
    def symbols(self):
        return self._symbols 
    @symbols.setter
    def symbols(self, symbols):
        self._symbols = symbols 

    @property
    def cell(self):
        """ cell parameters """
        return self._cell 
    @cell.setter
    def cell(self, cell):
        """ set cell parameters """
        self._cell = cell 

    def __len__(self):
        """Return number of atoms."""
        return len(self._coordinate)



