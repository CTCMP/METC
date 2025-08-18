import os, sys
import os.path as osp
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.getcwd())
import AutoinputQE as inpp
import configfile.config as configg
sys.dont_write_bytecode = True
####################################################
submit_pbs = configg.submit_order

# path_QE = '/work/wangr/data/dxy/soft/perturbo/QE/qe-7.0/bin'
# path_QE ='/work/wangr/data/dxy/soft/perturbo/QE/qe-7.0/bin'
# path_QE = '/public1/soft/qe/7.0.0-oneAPI.2022.1/bin'
path_QE = configg.QE_path
python_path = configg.python_path
# def get_pseudo():
#   if inpp.user_provide:
#     return osp.join(os.getcwd(), inpp.path_pp)
#   else:
#     path_presudo = osp.join(configg.path_presudo, 'pbe')
#     if inpp.SOC:
#       path_presudo = osp.join(configg.path_presudo, 'rel-'+inpp.functional)
#     else:
#       path_presudo = osp.join(configg.path_presudo, inpp.functional)
#     return path_presudo

def get_path_pseudo():
    if inpp.user_provide:
      return osp.join(os.getcwd(), inpp.path_pp)
    else:
      path_p = osp.join(configg.path_presudo, 'pbe')
      if inpp.SOC:
          path_p = osp.join(configg.path_presudo, 'rel-'+inpp.functional)
          if inpp.pseudo_potential == 'sssp':
              print("sssp pseduo is not support for SOC !")
          elif inpp.pseudo_potential == 'oncv':
              print("sssp pseduo is not support for oncv !")
      else:
          path_p = osp.join(configg.path_presudo, inpp.functional)
          if inpp.pseudo_potential == 'sssp':
              path_p = osp.join(configg.path_presudo, inpp.pseudo_potential)
          elif inpp.pseudo_potential == 'oncv':
              path_p = osp.join(configg.path_presudo, inpp.pseudo_potential, inpp.functional)
      print("Pseudo path: ", path_p)
      return path_p

path_presudo = get_path_pseudo()


try:
    perturbo_path = configg.perturbo_path
except:
    print("Not specify perturbo software path, Continue !")

try:
    thirdorder_path = configg.third_order_path
    # shengbte_path = '/public1/soft/ShengBTE/1.1.1/ShengBTE'
    shengbte_path = configg.shengbte_path
except:
    print("Not specify shengbte software path, Continue !")

try:
    fourthorder_path = configg.fourthorder_path
except:
    print("Not specify Fourthorder_espresso.py path, Continue !")

try:
  python2_path = configg.python2_path
except:
  python2 = "python2"

SOC = {
    'noncolin' : True,
    'lspinorb' : True,
}

CONTROL = {
    'calculation'   : 'scf',
    'prefix'        : 'pwscf',
    #'nstep'         : 600,
    'restart_mode'  : 'from_scratch',
    'pseudo_dir'    : path_presudo,
	  'outdir'        : './tmp',
	  #'verbosity'     : 'low',
    #'tprnfor'       : True,
	  #'tstress'       : True,
	  #'etot_conv_thr' : '1.0d-5',        
	  #'forc_conv_thr' : '1.0d-5', 
    # 'iprint'        :   2,          
}
SYSTEM = {
    'ibrav'       :  0,
    'ntyp'        :  1,
    'nat'         :  1,
    #'ecutwfc'     :  50,
    #'ecutrho'     :  600,
    #'occupations' : 'smearing',          
    #'smearing'    : 'gaussian',
    #'degauss'     :  0.05,
    'nspin'       : 1,   
}
ELECTRONS = {
    'conv_thr'         :  '1.0d-10',       
    'electron_maxstep' :  300,
    #'mixing_beta'      :  '0.7d0',
    #'mixing_mode'      : 'plain',
    #'diagonalization'  : 'david',
}
IONS = {
    # 'ion_dynamics' : 'bfgs'
}

CELL = {
    #'cell_dynamics'  : 'bfgs',                     
    # 'press_conv_thr' : 0.1,          
}

PHONON = {
  # 'verbosity': 'default',
  'tr2_ph'   : '1.0d-12',
  'prefix'   : 'pwscf',
  'outdir'   : './tmp',
  'fildyn'   : 'pwscf.dyn',
  'fildvscf' : 'dvscf',
#   'epsil'    :  True,
#   'lqdir'    :  True,
  'nq1'      : 6, 
  'nq2'      : 6, 
  'nq3'      : 6,
}
PHONON_q2r = {
    'fildyn'    : 'pwscf.dyn',
    'zasr'      : 'simple',
    'flfrc'     : 'pwscf.fc',
}
PHONON_mat = {
  'asr'            : 'simple',
  'flfrc'          : 'pwscf.fc',
  'flfrq'          : 'pwscf.freq',
  'q_in_band_form' : True,
  'q_in_cryst_coord':True,
}

PHONON_dos = {
 'asr'             : 'simple', 
 'flfrc'           : 'pwscf.fc', 
 'dos'             : True,
 'fldos'           : "phonon.dos",
  # 'nk1'            : 6, 
  # 'nk2'            : 6, 
  # 'nk3'            : 6,
}

band_plot = {
'prefix'     : 'pwscf',
'outdir'     : './tmp',
'filband'    :'band.dat',
'lp'         : True,
}

pro_band = {
'prefix'     : 'pwscf',
'outdir'     : './tmp',
'filproj'    :'proband',
'lsym'       : False,
'ngauss'     : 0,
'DeltaE'     : 0.05,
'degauss'    : 0.005,
}

dos_input = {
  'prefix'   :  'pwscf',
  'outdir'   : './tmp',
#   'bz_sum'    : 'tetrahedra_opt',
  'ngauss'   : 1,
  'degauss'  : 0.005,
  'deltaE'   : 0.01,
  'lsym'     : True,
  'filpdos'  : 'pwscf',
  'filproj'  : 'pwscf',
}

w90 = {
'outdir'     : './tmp',
'prefix'     : 'pwscf',
'seedname'   :'pwscf',
# 'spin_component' : 'none',
'write_mmn'  : True,
'write_amn'  : True,
'write_unk'  : False,
'write_dmn'  : False,
# 'write_sHu'  : True,
# 'write_sIu'  : True,
'write_spn'  : True,
'wan_mode'   : 'standalone',
}

w90_SCDM = {
'outdir'     : './tmp',
'prefix'     : 'pwscf',
'seedname'   :'pwscf',
'spin_component' : 'none',
'write_mmn'  : True,
'write_amn'  : True,
'write_unk'  : False,
'write_dmn'  : False,
# 'write_sHu'  : True,
# 'write_sIu'  : True,
'write_spn'  : True,
################ SCDM algorithm ##############################
'scdm_proj'  : True,
'scdm_entanglement' : 'erfc',   # erfc_smooth, amn_smooth
# 'scdm_mu'    :,
'scdm_sigma' : 0.01,
################ SCDM algorithm ##############################
'wan_mode'   : 'standalone',
}
try:
    nbnd = inpp.nbnd       # default: get nbnd from input_system tags.
except: 
    nbnd = inpp.input_system['nbnd']
wannier90_para = {
'num_wann'            : 18,                         #  set to NBANDS by VASP(numbers of wannier orbital)
'num_bands'           : nbnd,          # numbers of total bands
'####### energy windows'             : '###########',
'dis_froz_max'        : 3.0,
'dis_froz_min'        : -3.0,
'#************'             : '************',
# 'dis_win_max'         : 4.0,
# 'dis_win_min'         : -10,
'######################'             : '###########',
'spinors'             : False,
'bands_plot'          : True,
'bands_num_points'    : 101,
'num_iter'            : 10,
'conv_tol'            : 1.0e-10,
'dis_num_iter'        : 500,
'dis_conv_tol'        : 1.0e-10,
'write_hr'            : True,
'write_u_matrices'    : True,
'write_xyz'           : True,
# 'kmesh_tol'           : 0.0001,
# 'guiding_centres'     : True,
}

epw = {
  'prefix'            : 'pwscf',
  'outdir'            : './tmp',
  'dvscf_dir'         : '../phonon/save',
  'iverbosity'        : 0,
  #################################
  # 'elph'              : True,         # calculation of electron-phonon matrix
  # 'epbwrite'          : True,
  # 'epbread'           : False,
  # 'epwwrite'          : True,
  # 'epwread'           : False,
  # 'etf_mem'           : 1,
  ################################
  # 'lifc'              : True,
  # 'asr_typ'           : 'crystal',
  ################################
  # 'wannierize'        : True,
  # 'num_iter'          : 300,
  # 'iprint'            : 2,
  # 'elecselfen'        : False,
  # 'phonselfen'        : True,
  # 'a2f'               : True,
  # 'efermi_read'       : True,
  # 'fermi_energy'      : 0.0,
}
shengbte_allocations = {
  #   'ngrid(:)'          : inpp.kpt,
  # 'norientations'     : 0,          # number of orientations along which to study nanowires
}

shengbte_crystals = {
  'lfactor'           : 0.1,          # unit nm 
  'scell(:)'          : '3 3 3',
  'born'              : 0,            # add born effective matrix, 0: not, 1: yes and give born parameters.
  #'orientations'      : ,            # terns of integer indices defining the crystallographic directions along which to study nanowires
}

shengbte_parameters = {
  # 'T'                 : 300,        # temperature
  'T_min'             : 100,
  'T_max'             : 800,
  'T_step'            : 50,
  'scalebroad'        : 1.0,          # default 1.0 
#  'maxiter'           : 1000,         # default 1000, maximum number of iterations allowed in the BTE convergence process
#  'nticks'            : 100,          # number of different values of the mean free path at which to compute the cumulative thermal conductivity
}

shengbte_flags = {
  'nonanalytic '      : True,         # default true
#  'isotopes'          : True,         # default true
#  'nanowires'         : False,        # default false
  'espresso'          : True,         # default false
  # 'nthreads'          : 1,            # default 1, number of OpenMP threads each MPI process will use.
}

pert_qe2pert = {
  'prefix'           : 'pwscf',
  'outdir'           : './tmp',
  'phdir'            : '../phonon/pert_save',
  'lwannier'         : True,
}
pert_pertubo = {
'prefix'             : 'pwscf',
# 'tmp_dir'            : './tmp',
}

atom_mass = {
  "H"   : 1.00794,
  "He"  : 4.0026,
  "Li"  : 6.941,
  "Be"  : 9.01218,
  "B"   : 10.811,
  "C"   : 12.01070,
  "N"   : 14.0067,
  "O"   : 15.9994,
  "F"   : 18.9984,
  "Ne"  : 20.1797, 
  "Na"  : 22.9898, 
  "Mg"  : 24.305,
  "Al"  : 26.9815, 
  "Si"  : 28.0850, 
  "P"   : 30.9738,
  "S"   : 32.066, 
  "Cl"  : 35.4527, 
  "Ar"  : 39.948, 
  "K"   : 39.0983,
  "Ca"  : 40.078,
  "Sc"  : 44.9559,
  "Ti"  : 47.88,
  "V"   : 50.9415,
  "Cr"  : 51.996,
  "Mn"  : 54.938,
  "Fe"  : 55.847,
  "Co"  : 58.9332,
  "Ni"  : 58.6934, 
  "Cu"  : 63.546,
  "Zn"  : 65.39, 
  "Ga"  : 69.723, 
  "Ge"  : 72.61,
  "As"  : 74.9216,
  "Se"  : 78.96, 
  "Br"  : 79.904,
  "Kr"  : 83.8, 
  "Rb"  : 85.4678,
  "Sr"  : 87.62, 
  "Y"   : 88.9059,
  "Zr"  : 91.224, 
  "Nb"  : 92.9064,
  "Mo"  : 95.94,
  "Tc"  : 98, 
  "Ru"  : 101.07,
  "Rh"  : 102.906, 
  "Pd"  : 106.42,
  "Ag"  : 107.868, 
  "Cd"  : 112.41,
  "In"  : 114.82,
  "Sn"  : 118.71, 
  "Sb"  : 121.757, 
  "Te"  : 127.6, 
  "I"   : 126.904, 
  "Xe"  : 131.29, 
  "Cs"  : 132.905, 
  "Ba"  : 137.33, 
  "La"  : 138.905, 
  "Ce"  : 140.12, 
  "Pr"  : 140.908, 
  "Nd"  : 144.24, 
  "Pm"  : 145, 
  "Sm"  : 150.36, 
  "Eu"  : 151.965, 
  "Gd"  : 157.25, 
  "Tb"  : 158.925, 
  "Dy"  : 162.5, 
  "Ho"  : 164.93, 
  "Er"  : 167.26, 
  "Tm"  : 168.934, 
  "Yb"  : 173.04, 
  "Lu"  : 174.967, 
  "Hf"  : 178.49, 
  "Ta"  : 180.948, 
  "W"   : 183.85, 
  "Re"  : 186.207, 
  "Os"  : 190.2, 
  "Ir"  : 192.22, 
  "Pt"  : 195.08, 
  "Au"  : 196.966, 
  "Hg"  : 200.59, 
  "Tl"  : 204.383, 
  "Pb"  : 207.2, 
  "Bi"  : 208.98, 
  "Po"  : 209, 
  "At"  : 210, 
  "Rn"  : 222, 
  "Fr"  : 223, 
  "Ra"  : 226.025, 
  "Ac"  : 227, 
  "Th"  : 232.038, 
  "Pa"  : 231.036, 
  "U"   : 238.029, 
  "Np"  : 237.048, 
  "Pu"  : 244, 
  "Am"  : 243, 
  "Cm"  : 247, 
  "Bk"  : 247, 
  "Cf"  : 251, 
  "Es"  : 252, 
  "Fm"  : 257, 
  "Md"  : 258, 
  "No"  : 259, 
  "Lr"  : 262, 
  "Rf"  : 261, 
  "Db"  : 262, 
  "Sg"  : 263,
  "Bh"  : 262, 
  "Hs"  : 265, 
  "Mt"  : 266, 
  "Uun" : 269, 
  "Uuu" : 272,
  "Uub" : 277,
}