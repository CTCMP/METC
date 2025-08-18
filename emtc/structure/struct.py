# -*- coding: utf-8 -*-
import os, sys 
import os.path as osp
import numpy as np 
from src.structure.utils import atom_mass
import spglib

def length_l(a):
    return np.sqrt(a[0]**2 + a[1]**2 + a[2]**2)

def angle(a, b):
    dot_product = np.dot(a, b)
    norma = np.linalg.norm(a)
    normb = np.linalg.norm(b)
    cos_theta = dot_product / (norma * normb)
    angle = np.arccos(cos_theta)
    return np.round(np.degrees(angle))

def get_xyz_label(rot, tran):
    symbols = ['x', 'y', 'z']
    operation = []
    
    # 处理旋转矩阵部分
    for i in range(3):
        terms = []
        for j in range(3):
            coeff = rot[i, j]
            if coeff == 1:
                terms.append(f"{symbols[j]}")
            elif coeff == -1:
                terms.append(f"-{symbols[j]}")
        
        # 处理平移向量部分
        translation_term = tran[i]
        if translation_term != 0:
            fraction = f"{translation_term:.3f}".rstrip('0').rstrip('.')
            terms.append(f"{fraction}")
        
        # 连接所有项
        if not terms:
            terms.append("0")
        operation.append(" + ".join(terms))
    
    return ", ".join(operation)

def get_xyz_label_all(rots, trans):
    operations = []
    for i in range(len(rots)):
        operation = get_xyz_label(rots[i], trans[i])
        operations.append(operation)
    return operations

def get_crystal_system(spg):
    space_group_number = int(spg)
    if 1 <= space_group_number <= 2:
        return "Triclinic"
    elif 3 <= space_group_number <= 15:
        return "Monoclinic"
    elif 16 <= space_group_number <= 74:
        return "Orthorhombic"
    elif 75 <= space_group_number <= 142:
        return "Tetragonal"
    elif 143 <= space_group_number <= 167:
        return "Trigonal"
    elif 168 <= space_group_number <= 194:
        return "Hexagonal"
    elif 195 <= space_group_number <= 230:
        return "Cubic"
    else:
        return "P1"

def save_cif_vesta(filepath, filename):
    struct = StructAtoms('POSCAR')
    ff = open(osp.join(filepath, filename), 'w')
    ff.write('\n')
    ff.write('#====================================================================== \n')
    ff.write("# CRYSTAL DATA \n")
    ff.write("#---------------------------------------------------------------------- \n")
    ff.write('data_VESTA_phase_1 \n')
    ff.write('\n')
    
    ff.write(f"{'_chemical_name_common':<28}  '{'Xianyong Ding':<20}'\n")
    ff.write(f"{'_cell_length_a':28}  {struct.la} \n")
    ff.write(f"{'_cell_length_b':28}  {struct.lb} \n")
    ff.write(f"{'_cell_length_c':28}  {struct.lc} \n")

    ff.write(f"{'_cell_angle_alpha':28}  {struct.alpha} \n")
    ff.write(f"{'_cell_angle_beta':28}  {struct.beta} \n")
    ff.write(f"{'_cell_angle_gamma':28}  {struct.gamma} \n")
    ff.write(f"{'_cell_volume':<28}  {struct.omega} \n")
    # ff.write(f"{'_space_group_name_H-M_alt':<28}  '{struct.spg['international_short']:<20}'\n")
    # ff.write(f"{'_space_group_IT_number':<28}  {struct.spg['number']:<20} \n")
    ff.write(f"{'_space_group_name_H-M_alt':<28}  'P 1'\n")
    ff.write(f"{'_space_group_IT_number':<28}  1 \n")
    
    ff.write('\n')

    ff.write("loop_\n")
    ff.write('_space_group_symop_operation_xyz \n')
    ff.write('x, y, z \n')
    # for i in range(len(struct.sym_operation_xyz)):
    #     ff.write(f"  {struct.sym_operation_xyz[i]} \n")

    ff.write('\n')
    ff.write('loop_ \n')
    ff.write(f"  {'_atom_site_label':<30} \n")
    ff.write(f"  {'_atom_site_occupancy':<30} \n")
    ff.write(f"  {'_atom_site_fract_x':<30} \n")
    ff.write(f"  {'_atom_site_fract_y':<30} \n")
    ff.write(f"  {'_atom_site_fract_z':<30} \n")
    # ff.write(f"  {'_atom_site_adp_type':<30} \n")
    # ff.write(f"  {'_atom_site_U_iso_or_equiv':<30} \n")
    # ff.write(f"  {'_atom_site_type_symbol':<30} \n")
    count = 0
    for i in range(len(struct.numbers)):
        for j in range(struct.numbers[i]):
            # ff.write(f"  {struct.symbols[i]+str(j+1):<10}  {1.0:6}  {struct.coordinate[count][0]:< 12.10f}  {struct.coordinate[count][1]:< 12.10f}  {struct.coordinate[count][2]:< 12.10f}  Uiso  ?  {struct.symbols[i]} \n")
            ff.write(f"  {struct.symbols[i]+str(j+1):<10}  {1.0:6}  {struct.coordinate[count][0]:< 12.10f}  {struct.coordinate[count][1]:< 12.10f}  {struct.coordinate[count][2]:< 12.10f}  \n")
            count = count + 1
    ff.write('\n')
    ff.close()

class StructAtoms:
    """ Class to get strucutre info """
    def __init__(self, filename=None, symbols = None,numbers = None, masses = None, scale = None,cell = None, cor=None, atomic_numbers=None, atoms_labels=None,
                 spg=None, symmetry_operations=None, sym_operation_xyz=None, crystal_type=None, bg=None, alpha=None, beta=None, gamma=None, la=None, lb=None, lc=None, omega=None):
        # self._symbols = symbols 
        # self._numbers = numbers 
        # self._masses = masses
        # self._scale = scale
        # self._cell = cell 
        # self._coordinate = cor
        # self._atomic_numbers = atomic_numbers
        # self._spg = spg
        # self._symmetry_operations = symmetry_operations
        # self._bg = bg
        # self._alpha = alpha 
        # self._beta = beta 
        # self._gamma = gamma 
        # self._la = la 
        # self._lb = lb 
        # self._lc = lc
        self._set_parameter(filename)

    def _set_parameter(self, filename="POSCAR"):
        latt, element, num_element, tot_ele, cor, scale, masses, atomic_numbers, atoms_labels = get_poscar(os.getcwd(), filename)
        cell = (latt, cor, atomic_numbers)
        bg = get_reciprocal(os.getcwd(), filename)
        symmetry_operations = spglib.get_symmetry(cell)
        spg = spglib.get_spacegroup_type_from_symmetry(symmetry_operations['rotations'], symmetry_operations['translations'])
        self.cell = latt  
        self.bg = bg
        self.symbols = element 
        self.numbers = num_element 
        self.scale = scale
        self.coordinate = cor 
        self.masses = masses
        self.atomic_numbers = atomic_numbers
        self.atoms_labels = atoms_labels
        self.spg = spg
        self.symmetry_operations = symmetry_operations
        self.la = length_l(self.cell[0])
        self.lb = length_l(self.cell[1])
        self.lc = length_l(self.cell[2])
        self.alpha = angle(self.cell[1], self.cell[2]) 
        self.beta = angle(self.cell[0], self.cell[2]) 
        self.gamma = angle(self.cell[0], self.cell[1]) 
        self.omega = np.dot(self.cell[0], np.cross(self.cell[1], self.cell[2]))
        self.sym_operation_xyz = get_xyz_label_all(self.symmetry_operations['rotations'], self.symmetry_operations['translations'])
        self.crystal_type = get_crystal_system(self.spg['number'])

    @property
    def atoms_labels(self):
        return self._atoms_labels
    @atoms_labels.setter
    def atoms_labels(self, atoms_labels):
        self._atoms_labels = atoms_labels

    @property
    def crystal_type(self):
        return self._crystal_type
    @crystal_type.setter
    def crystal_type(self, crystal_type):
        self._crystal_type = crystal_type

    @property
    def sym_operation_xyz(self):
        return self._sym_operation_xyz
    @sym_operation_xyz.setter
    def sym_operation_xyz(self, sym_operation_xyz):
        self._sym_operation_xyz = sym_operation_xyz

    @property
    def omega(self):
        return self._omega
    @omega.setter
    def omega(self, omega):
        self._omega = omega

    @property
    def alpha(self):
        return self._alpha
    @alpha.setter
    def alpha(self, alpha):
        self._alpha = alpha

    @property
    def la(self):
        return self._la
    @la.setter
    def la(self, la):
        self._la = la

    @property
    def lb(self):
        return self._lb
    @lb.setter
    def lb(self, lb):
        self._lb = lb

    @property
    def lc(self):
        return self._lc
    @lc.setter
    def lc(self, lc):
        self._lc = lc

    @property
    def beta(self):
        return self._beta
    @beta.setter
    def beta(self, beta):
        self._beta = beta

    @property
    def gamma(self):
        return self._gamma
    @gamma.setter
    def gamma(self, gamma):
        self._gamma = gamma

    @property
    def coordinate(self):
        return self._coordinate
    @coordinate.setter
    def coordinate(self, coordinate):
        self._coordinate = coordinate

    @property
    def spg(self):
        return self._spg
    @spg.setter
    def spg(self, spg):
        self._spg = spg

    @property
    def symmetry_operations(self):
        return self._symmetry_operations
    @symmetry_operations.setter
    def symmetry_operations(self, symmetry_operations):
        self._symmetry_operations = symmetry_operations

    @property
    def scale(self):
        return self._scale
    @scale.setter
    def scale(self, scale):
        self._scale = scale

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
    def atomic_numbers(self):
        return self._atomic_numbers
    @atomic_numbers.setter
    def atomic_numbers(self, atomic_numbers):
        self._atomic_numbers = atomic_numbers 

    @property
    def cell(self):
        """ cell parameters """
        return self._cell 
    @cell.setter
    def cell(self, cell):
        """ set cell parameters """
        self._cell = cell 

    @property
    def bg(self):
        """ reciprocal cell parameters """
        return self._bg 
    @bg.setter
    def bg(self, bg):
        """ set reciprocal cell parameters """
        self._bg = bg 

    def __len__(self):
        """Return number of atoms."""
        return len(self._coordinate)

# def StructClass(filename):
#     latt, element, num_element, tot_ele, cor, scale, masses, atomic_numbers = get_poscar(os.getcwd(), filename)
#     cell = (latt, cor, atomic_numbers)
#     bg = get_reciprocal(os.getcwd(), filename)
#     symmetry_operations = spglib.get_symmetry(cell)
#     spg = spglib.get_spacegroup_type_from_symmetry(symmetry_operations['rotations'], symmetry_operations['translations'])
#     structc = StructAtoms()
#     structc.cell = latt  
#     structc.bg = bg
#     structc.symbols = element 
#     structc.numbers = num_element 
#     structc.scale = scale
#     structc.coordinate = cor 
#     structc.masses = masses
#     structc.atomic_numbers = atomic_numbers
#     structc.spg = spg
#     structc.symmetry_operations = symmetry_operations
#     return structc 

