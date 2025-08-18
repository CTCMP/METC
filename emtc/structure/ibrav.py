# -*- coding: utf-8 -*-

from __future__ import annotations

import numpy as np 
import math
import qskit.utilities as util
from qskit.utilities import atom_mass

# def QEAtoms(*args, **kwargs):
#     """Atoms class for Quantum espresso.
#     """
#     warnings.warn(
#         "phonopy.atoms.Atoms is deprecated. Please use PhonopyAtoms instead of Atoms.",
#         DeprecationWarning,
#         stacklevel=2,
#     )
#     return PhonopyAtoms(*args, **kwargs)

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
    return np.round(latt, 10), element, num_element, tot_ele, np.round(cor, 10), scale, mass, masses_ifc

def get_reciprocal(latt):
    rep = []
    a1 = latt[0]
    a2 = latt[1]
    a3 = latt[2]
    volume = np.dot(a1, np.cross(a2, a3))
    scalar = 2*np.pi/volume
    for i in range(0, 3):
        b1 = scalar*np.cross(a2, a3)[i]
        rep.append(b1)
    for j in range(0, 3):
        b2 = scalar*np.cross(a3, a1)[j]
        rep.append(b2)
    for k in range(0, 3):
        b3 = scalar * np.cross(a1, a2)[k]
        rep.append(b3)
    rep_latt = np.reshape(np.array(rep), (3, 3))
    return rep_latt 

def read_poscar(filepath, filename):
    with open(osp.join(filepath, filename), 'r') as f:
        lines = f.readlines()
    
    # 读取缩放因子
    scale_factor = float(lines[1].strip())
    
    # 读取晶格矢量
    lattice_vectors = []
    for i in range(2, 5):
        vector = list(map(float, lines[i].strip().split()))
        lattice_vectors.append(vector)
    
    lattice_vectors = np.array(lattice_vectors) * scale_factor
    return np.array(lattice_vectors)

def calculate_lattice_parameters(lattice_vectors):
    a_vec = lattice_vectors[0]
    b_vec = lattice_vectors[1]
    c_vec = lattice_vectors[2]
    
    a = np.linalg.norm(a_vec)
    b = np.linalg.norm(b_vec)
    c = np.linalg.norm(c_vec)
    
    alpha = np.arccos(np.dot(b_vec, c_vec) / (b * c)) * 180 / np.pi
    beta = np.arccos(np.dot(a_vec, c_vec) / (a * c)) * 180 / np.pi
    gamma = np.arccos(np.dot(a_vec, b_vec) / (a * b)) * 180 / np.pi
    
    return a, b, c, alpha, beta, gamma

def identify_ibrav(a, b, c, alpha, beta, gamma, tolerance=1e-2):
    # 定义一个函数来比较浮点数是否相等，考虑到数值误差
    def is_equal(x, y):
        return abs(x - y) < tolerance

    # 简单立方 (SC)
    if is_equal(a, b) and is_equal(b, c) and \
       is_equal(alpha, 90) and is_equal(beta, 90) and is_equal(gamma, 90):
        return 1
    # 面心立方 (FCC)
    elif is_equal(a, b) and is_equal(b, c) and \
         is_equal(alpha, 60) and is_equal(beta, 60) and is_equal(gamma, 60):
        return 2
    # 体心立方 (BCC)
    elif is_equal(a, b) and is_equal(b, c) and \
         is_equal(alpha, 109.4712206) and is_equal(beta, 109.4712206) and is_equal(gamma, 109.4712206):
        return 3
    # 简单六方 (Hexagonal)
    elif is_equal(a, b) and not is_equal(c, a) and \
         is_equal(alpha, 90) and is_equal(beta, 90) and is_equal(gamma, 120):
        return 4
    # 菱方 (Trigonal R)
    elif is_equal(a, b) and is_equal(b, c) and \
         is_equal(alpha, beta) and is_equal(beta, gamma) and not is_equal(alpha, 90) and not is_equal(alpha, 60):
        return 5
    # 简单四方
    elif is_equal(a, b) and not is_equal(c, a) and \
         is_equal(alpha, 90) and is_equal(beta, 90) and is_equal(gamma, 90):
        return 6
    # 体心四方
    elif is_equal(a, b) and not is_equal(c, a) and \
         is_equal(alpha, 90) and is_equal(beta, 90) and is_equal(gamma, 90):
        return 7
    # 简单正交
    elif not is_equal(a, b) and not is_equal(b, c) and \
         is_equal(alpha, 90) and is_equal(beta, 90) and is_equal(gamma, 90):
        return 8
    # 基面中心正交
    elif not is_equal(a, b) and not is_equal(b, c) and \
         is_equal(alpha, 90) and is_equal(beta, 90) and is_equal(gamma, 90):
        return 9
    # 面心正交
    elif not is_equal(a, b) and not is_equal(b, c) and \
         is_equal(alpha, 90) and is_equal(beta, 90) and is_equal(gamma, 90):
        return 10
    # 体心正交
    elif not is_equal(a, b) and not is_equal(b, c) and \
         is_equal(alpha, 90) and is_equal(beta, 90) and is_equal(gamma, 90):
        return 11
    # 单斜 P
    elif not is_equal(a, b) and not is_equal(b, c) and \
         is_equal(alpha, 90) and not is_equal(beta, 90) and is_equal(gamma, 90):
        return 12
    # 单斜 B
    elif not is_equal(a, b) and not is_equal(b, c) and \
         is_equal(alpha, 90) and not is_equal(beta, 90) and is_equal(gamma, 90):
        return 13
    # 三斜
    else:
        return 14

def calculate_celldm(a, b, c, alpha, beta, gamma, ibrav):
    celldm = [0] * 6
    if ibrav in [1, 2, 3]:
        celldm[0] = a * util.AngToBohr
    elif ibrav in [4]:
        celldm[0] = a * util.AngToBohr
        celldm[2] = c / a
    elif ibrav in [6, 7]:
        celldm[0] = a * util.AngToBohr
        celldm[2] = c / a
    elif ibrav in [8, 9, 10, 11]:
        celldm[0] = a * util.AngToBohr
        celldm[1] = b / a
        celldm[2] = c / a
    elif ibrav in [12, 13]:
        celldm[0] = a * util.AngToBohr
        celldm[1] = b / a
        celldm[2] = c / a
        celldm[3] = math.cos(beta * math.pi / 180)
    elif ibrav == 14:
        celldm[0] = a * util.AngToBohr
        celldm[1] = b / a
        celldm[2] = c / a
        celldm[3] = math.cos(alpha * math.pi / 180)
        celldm[4] = math.cos(beta * math.pi / 180)
        celldm[5] = math.cos(gamma * math.pi / 180)
    return celldm

def calculate_lattice_vectors_from_celldm(celldm, ibrav):
    if ibrav == 1:
        # Simple cubic
        a = celldm[0]
        lattice_vectors = np.array([
            [a, 0, 0],
            [0, a, 0],
            [0, 0, a]
        ])
    
    elif ibrav == 2:
        # Face-centered cubic (FCC)
        a = celldm[0]
        lattice_vectors = np.array([
            [0, 0.5 * a, 0.5 * a],
            [0.5 * a, 0, 0.5 * a],
            [0.5 * a, 0.5 * a, 0]
        ])
    
    elif ibrav == 3:
        # Body-centered cubic (BCC)
        a = celldm[0]
        lattice_vectors = np.array([
            [-0.5 * a, 0.5 * a, 0.5 * a],
            [0.5 * a, -0.5 * a, 0.5 * a],
            [0.5 * a, 0.5 * a, -0.5 * a]
        ])
    
    elif ibrav == 4:
        # Hexagonal
        a = celldm[0]
        c = celldm[2] * a
        lattice_vectors = np.array([
            [a, 0, 0],
            [-0.5 * a, np.sqrt(3) / 2 * a, 0],
            [0, 0, c]
        ])
    
    elif ibrav == 6:
        # Simple tetragonal
        a = celldm[0]
        c = celldm[2] * a
        lattice_vectors = np.array([
            [a, 0, 0],
            [0, a, 0],
            [0, 0, c]
        ])
    
    elif ibrav == 7:
        # Body-centered tetragonal
        a = celldm[0]
        c = celldm[2] * a
        lattice_vectors = np.array([
            [0.5 * a, 0.5 * a, 0],
            [0.5 * a, -0.5 * a, 0],
            [0, 0, c]
        ])
    
    elif ibrav == 8:
        # Simple orthorhombic
        a = celldm[0]
        b = celldm[1] * a
        c = celldm[2] * a
        lattice_vectors = np.array([
            [a, 0, 0],
            [0, b, 0],
            [0, 0, c]
        ])
    
    elif ibrav == 9:
        # Body-centered orthorhombic
        a = celldm[0]
        b = celldm[1] * a
        c = celldm[2] * a
        lattice_vectors = np.array([
            [0.5 * a, 0.5 * b, 0],
            [-0.5 * a, 0.5 * b, 0],
            [0, 0, c]
        ])
    
    elif ibrav == 10:
        # Face-centered orthorhombic
        a = celldm[0]
        b = celldm[1] * a
        c = celldm[2] * a
        lattice_vectors = np.array([
            [0.5 * a, 0.5 * b, 0],
            [-0.5 * a, 0.5 * b, 0],
            [0, 0, c]
        ])
    
    elif ibrav == 11:
        # Base-centered orthorhombic
        a = celldm[0]
        b = celldm[1] * a
        c = celldm[2] * a
        lattice_vectors = np.array([
            [0.5 * a, 0.5 * b, 0],
            [-0.5 * a, 0.5 * b, 0],
            [0, 0, c]
        ])
    
    elif ibrav == 12:
        # Simple triclinic
        a = celldm[0]
        b = celldm[1] * a
        c = celldm[2] * a
        alpha = np.arccos(celldm[3])  # cos(alpha)
        beta = np.arccos(celldm[4])   # cos(beta)
        gamma = np.arccos(celldm[5])  # cos(gamma)
        
        lattice_vectors = np.array([
            [a, 0, 0],
            [b * np.cos(gamma), b * np.sin(gamma), 0],
            [c * np.cos(beta), 
             c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma), 
             c * np.sqrt(1 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma) - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2) / np.sin(gamma)]
        ])
    
    elif ibrav == 13:
        # Body-centered triclinic (Example: For reference, not standard)
        a = celldm[0]
        b = celldm[1] * a
        c = celldm[2] * a
        alpha = np.arccos(celldm[3])  # cos(alpha)
        beta = np.arccos(celldm[4])   # cos(beta)
        gamma = np.arccos(celldm[5])  # cos(gamma)
        
        lattice_vectors = np.array([
            [0.5 * a, 0.5 * b, 0],
            [-0.5 * a, 0.5 * b, 0],
            [0, 0, c]
        ])
    
    elif ibrav == 14:
        # Triclinic
        a = celldm[0]
        b = celldm[1] * a
        c = celldm[2] * a
        alpha = np.arccos(celldm[3])  # cos(alpha)
        beta = np.arccos(celldm[4])   # cos(beta)
        gamma = np.arccos(celldm[5])  # cos(gamma)
        
        lattice_vectors = np.array([
            [a, 0, 0],
            [b * np.cos(gamma), b * np.sin(gamma), 0],
            [c * np.cos(beta), 
             c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma), 
             c * np.sqrt(1 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma) - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2) / np.sin(gamma)]
        ])
    
    else:
        raise ValueError(f"ibrav={ibrav} is not supported in this example.")
    
    return lattice_vectors

def ibrav_celldm(filepath, filename):
    ibrav_user = int(input("Please input ibrav (0 or 1):"))
    lattice_vectors, element, num_element, tot_ele, cor, scale, mass, mass_ifc = get_poscar(os.getcwd(), filename)
    a, b, c, alpha, beta, gamma = calculate_lattice_parameters(lattice_vectors)
    ibrav = identify_ibrav(a, b, c, alpha, beta, gamma)
    celldm = np.array(calculate_celldm(a, b, c, alpha, beta, gamma, ibrav))
    if ibrav_user == 0:
        print("IBRAV: ", 0)
        ibrav = 0
        lattice_vec = lattice_vectors / celldm[0] * util.AngToBohr
    else:
        lattice_vec = calculate_lattice_vectors_from_celldm(celldm, ibrav)
    return ibrav, celldm, lattice_vec

class QEAtoms:
    """ Class to get Quantum espresso strucutre info.
    ----------
    cell : np.ndarray, lattice vectors
    symbols : list(str), List of chemical symbols of atoms.
    natoms : np.ndarray, numbers of each element
    masses : np.ndarray, atomic masses of each element
    scale : scale of lattice vector
    ibrav : Briav lattice type for Quantum espresso 
    celldm : celldm for Quantum espresso
    positions : np.ndarray, Positions of atoms in Cartesian coordinates
    frac_positons : np.ndarray, Positions of atoms in fraction coordinates
    bg : reciprocal lattice vector.
    volume : cell volume

    """
    def __init__(
            self, 
            cell=None,
            symbols=None, 
            natoms=None, 
            masses=None,
            scale=None, 
            ibrav=None, 
            celldm=None, 
            positions=None, 
            frac_positons = None,
            bg=None,
            volume=None,
        ):
        pass
        # self._symbols = symbols 
        # self._numbers = numbers 
        # self._masses = masses
        # self._ibrav = ibrav
        # self._celldm = celldm
        # self._scale = scale
        # self._cell = cell 
        # self._coordinate = cor
        # self._natoms = natoms
        # self._ntype = ntype
        # self._mass_ifc = mass_ifc
        # self._bg = bg

    @property
    def bg(self):
        return self._bg
    @bg.setter
    def bg(self, bg):
        self._bg = bg

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
    
def StructQEClass(filename, ibrav_user):
    bohrToAng = 1.8897161646320724
    # ibrav_user = int(input("Please input ibrav (0 or 1):"))
    lattice_vectors, element, num_element, tot_ele, cor, scale, masses, mass_ifc = get_poscar(os.getcwd(), filename)
    bg = get_reciprocal(lattice_vectors)
    a, b, c, alpha, beta, gamma = calculate_lattice_parameters(lattice_vectors)
    ibrav = identify_ibrav(a, b, c, alpha, beta, gamma)
    celldm = np.array(calculate_celldm(a, b, c, alpha, beta, gamma, ibrav))
    cor_qe = np.dot(cor, lattice_vectors)
    structc = StructQE()
    if ibrav_user == 0:
        ibrav = 0
        lattice_vec = lattice_vectors / celldm[0] * util.AngToBohr
        structc.ibrav = ibrav
    else:
        lattice_vec = calculate_lattice_vectors_from_celldm(celldm, ibrav)
        structc.ibrav = ibrav
    # structc.cell = lattice_vec  
    structc.cell = lattice_vectors * util.AngToBohr
    structc.symbols = element 
    structc.numbers = num_element 
    structc.celldm = celldm
    structc.coordinate = cor_qe / celldm[0] * bohrToAng
    structc.masses = masses
    structc.natoms = len(cor_qe)
    structc.ntype = len(num_element)
    structc.mass_ifc = mass_ifc
    structc.bg = bg
    return structc 