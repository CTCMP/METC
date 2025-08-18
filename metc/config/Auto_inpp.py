#****************** serve info *************************
node = 'node05'          # sr850
num_nodes = 1
num_process = 40
mpi_par = '-npool 4'     # '-nk 4', '-npool 4', or '0': do not use parallel par.
#****************** pseudopential ***********************
user_provide = False  # user define pseudo potential
path_pp = "../pp"     # you need to rename you pseudo file to element.UPF (eg. C.UPF, Si.UPF)

# select options for pseudo potential 
pseudo_potential = 'rrkj'       # rrkj, rrkjus, kjpaw, bpaw, nc, ae, mt, bhs, vbc, van             
functional  = 'pbe'             # pbe, bp, pz, vwn, blyp, pw91, tpss, coulomb 
nl_corr = ''                    # spdfn(nonlinear core-correction)
#*********************************************************
# LDA + U calculations, For QE version > 7.1
SOC = False                                      # SOC or not
noncolin = False                                 # noncolin calculations.
ISPIN = 1                                        # 2 : spin polarized
LDAU = False 
Hubbard_projectors = 'atomic'              # 'atomic', 'ortho-atomic', 'norm-atomic', 'wf', 'pseudo'
LDAU_par = {'1':['3d', 2.0], '3': ['3d', 2.0], } # 1 denote the first tyep atoms, '3d' means 3d orbital + U.                               
Magmom = {'1':[0.8], '2':[1], '3':[1], }           # setting initial magmom, more big is better, the examples represents 3 types of atoms.
# if noncolin of magmom : Magmom = {'1': [0.8, '90.0', '0.0']}, 0.8: initial magmom, 90.0: angle1(1) = 90.0, angle2(1) = 0.0.
#************** input parameters for QE ******************
ibrav = 0                                       # if 0, all 0 for all calculation, if not, get info from POSCAR
input_control = {
'calculation'      : 'vc-relax',                 # calculate type: 'vc-relax', 'relax', 'nscf', 'band', 'dos', 'phonon', 'wannier90'
#------------ opt parameters -----------------
'etot_conv_thr'    : '1.0d-4',        
'forc_conv_thr'    : '1.0d-3', 
'verbosity'        : 'low',                      # high, medium, or nowf: no wave function and charge density.
'disk_io'          : 'low',
'tprnfor'          :  True,
'wf_collect'       : False,
}
#---------- main input parameters -------------    
input_system ={      
'ecutwfc'          :  40,
'ecutrho'          :  400,                   ## INCUT for QE calculation
'nbnd'             :  40,
'occupations'      : 'smearing',             ## smearing: metals, fixed: for gapped materials, tetrahedra_opt.
'smearing'         : 'gaussian',   # 'marzari-vanderbilt', 'gaussian', mp: conductor
'degauss'          :  0.05,                  ## gap: ->0,  metal: 0.01-0.1
#'nosym'            : True,
# 'vdw_corr'         : 'grimme-d3',  # DFT-D3, DFT-D (D2)
# 'dftd3_version'    : 3,
}
# setting nbnd for the following calculations.
try:
    nbnd = input_system['nbnd']       # default: get nbnd from input_system tags.
except: 
    nbnd = 40
input_electrons = {
##--------- electron congervence ---------------
'mixing_beta'      :  '0.7d0',
'mixing_ndim'      : 10,                     # default: 8
#'mixing_ndim'      :   4,
'mixing_mode'      :'plain',
'electron_maxstep' : 300,
'conv_thr'         :  '1.0d-12',             ## default: 
'diagonalization'  :  "david",               # "david",  "cg"
'diago_full_acc'   :  True,
#'tot_magnetization':  '2',
}
input_ions = {
  'ion_dynamics' : 'bfgs',
  # 'trust_radius_ini': 0.5,   # initial ion displacement
  # 'trust_radius_min': 0.001, # minimum ion displacement
  # 'trust_radius_max' : 0.8,  # maximum ion displacement
}
input_cell = {
  'cell_dynamics'  : 'bfgs', 
  # 'press'          : 0.0,
  # 'cell_factor'    : 3.0,                    
  # 'press_conv_thr' : 0.1,    
  # 'cell_dofree'    : '2Dxy', # 'all', '2Dxy', 'ibrav', 'fixa/b/c', 'volume', 'epitaxial_ab/bc/ac/'  # For 2D systems
}
#********************************************************
#---------------- kpoints -----------------------
kpt = [8, 8, 8]
kpt_nscf = [8, 8, 8]
kpt_shift = [0, 0, 0]
#--------- bands and phonon path ----------------
band_points = 51                                
kpoints_band = '''
   0.0000   0.0000   0.0000   ! \Gamma
   0.5000   0.0000   0.5000   ! X
   0.5000   0.2500   0.7500   ! W
   0.3750   0.3750   0.7500   ! K
   0.0000   0.0000   0.0000   ! \Gamma
   0.5000   0.5000   0.5000   ! L
   0.6250   0.2500   0.6250   ! U
   0.5000   0.2500   0.7500   ! W
   0.5000   0.5000   0.5000   ! L
   0.3750   0.3750   0.7500   ! K
'''
#------------------ dos --------------------------
dos = {
#   'bz_sum'  : 'tetrahedra_opt',
  'ngauss'   : 1,
  'degauss'  : 0.005,
  'deltaE'   : 0.01,
  'lsym'     : True,
}
#*****************************************************************************
#**************************** phonon *****************************************
phon_points = 51  
qpoints = [6, 6, 6]
eph_k = [0, 0, 0]         # electron phonon coupling for defined k grid, default is gamma.
phonon_input = {
  'verbosity': 'default',   # debug, high, medium, low, default, minimal
  'tr2_ph'   : '1.0d-12',
  # 'reduce_io': False,
  'ldisp'    : True,
  # 'epsil'    : True,            # calculate the macroscopic dielectric constant.
  'recover'  : True,            # restart from an interrupted run.
# 'alpha_mix(1)' : 0.5,
# 'nmix_ph'   : 10,             # speed up convergence, but use more momery.
  'nq1'      : qpoints[0], 
  'nq2'      : qpoints[1],
  'nq3'      : qpoints[2],         
  'diagonalization' : 'david',  # 'cg'
  # 'start_q'  : 1,
  # 'last_q'   : 1,
}

phon_q2r = {
    'zasr'      : 'simple',
}
phonon_mat_input = {
  'asr'            : 'simple',
  'q_in_band_form' : True,
  # 'q_in_cryst_coord': True,
}

phon_dos = {
 'asr'             : 'simple', 
 'flfrc'           : 'pwscf.fc', 
 'dos'             : True,
 'fldos'           : "phonon.dos",
 'nk1'             : kpt[0], 
 'nk2'             : kpt[1], 
 'nk3'             : kpt[2],
}
#*****************************************************************************
#**************************** Wannier90 *****************************************
porjection_w90 = {
  'Si' : 'sp3',
}
wannier90 = {
'num_wann'            : 8,                         #  set to NBANDS by VASP(numbers of wannier orbital)
'num_bands'           : nbnd,          # numbers of total bands
#'exclude_bands'       : '21-40',
#----------- max energy --------------------------
'dis_froz_max'        : 3.0,
'dis_froz_min'        : -3.0,
#----------- min energy --------------------------
# 'dis_win_max'         : 4.0,
# 'dis_win_min'         : -10,
#-------------------------------------------------
'spinors'             : False,
'num_iter'            : 1000,
'dis_num_iter'        : 1000,
'bands_num_points'    : 101,
##------------------------------------------------
'write_hr'            : True,
'write_u_matrices'    : True,
'write_xyz'           : True,
# 'use_ws_distance'     : True,
# 'guiding_centres'     : True,
#'use_bloch_phases' : True,
}
#*****************************************************************************
#****************************** EPW ******************************************
nktonkf = 2
nqtonqf = 4
nk = kpt_nscf 
nq = qpoints
use_projection_w90 = True      # use w90 projection or not 
if use_projection_w90:
  porjection_epw = porjection_w90
else:
  porjection_epw = {
    'Si' : 'sp3',
}
epw_input = {
  'prefix'            : 'pwscf',
  'outdir'            : './tmp/',
  'dvscf_dir'         : '../phonon/save',
  ################### memory and format of output #####
  'iverbosity'        : 0,
  'max_memlt'         : 5.0,
  # 'etf_mem'           : 1,
  'lifc'              : True,
  'asr_typ'           : 'crystal',
  #--------------------------------
  # 'lpolar'            : True,       # polar material Lo mode is important
  # 'system_2d'           : False,    # for 2d system
  #--------------------------------
  #################################
  'elph'              : True,         # calculation of electron-phonon matrix
  'epbwrite'          : True,
  'epbread'           : False,
  'epwwrite'          : True,
  'epwread'           : False,
  ########### confirm the band and phonon dispersion
  'band_plot'         : False,
  #'filkf'             : './meshes/kpoints',
  #'filqf'             : './meshes/kpoints',
  ################################
  'efermi_read'       : True,
  'fermi_energy'      : 0,
  ###################### calculation #########################
  'scattering'        : True, 
  'phonselfen'        : True,
  'elecselfen'        : True, 
  # 'eliashberg'        : True,
  # 'ephwrite'          : True,
  # 'specfun'           : True,         # specfun_el, specfun_ph, sepcfun_pl
  # 'cumulant'          : True,         # electron spectral function : cal band in different temps
  #'prtgkk'            : False,          # allow to print the el-ph vertex |g| (in meV), would slow down the calculation
  #############################################################
  ############# calculation : superconductivity ###############
  # 'ep_coupling'       : True,
  # 'liso'              : True,
  # 'limag'             : True, 
  # 'lpade'             : True, 
  # 'lacon'             : True, 
  # 'laniso             : True,            # solve anisotropic ME eps.
  # 'wscut'             : 0.1,             # default: 1.0d0
  # 'muc'               : 0.1,             # default: 0.0d0
  # 'degaussq'          : 0.05,
  # 'nsiter'            : 500,
  # 'conv_thr_iaxis'    : '1.0d-3', # ! convergence threshold for solving ME eqs. on imaginary axis
  # 'conv_thr_racon'    : '1.0d-3', # convergence threshold for solving ME eqs. on real axis
  # 'fermi_plot'        : True, # write files to plot fermi surface.
  #############################################################
  'fsthick'           : 3.0,
  # 'nstemp'            : 1,              # default: 1
  'temps'             : 300,          # list for temps: [100, 200, 300]
  'degaussw'          : 0.005,
  # 'nsiter'            : 500,            # default: 40
  # 'npade'             : 40,             # default: 90
  # 'conv_thr_iaxis'    : 0.001,          # default: 1.0d0-5
  # 'conv_thr_racon'    : 0.001,          # default: 1.0d0-5
  ################################
  # 'specfun_ph'        : True,           # Calculate the phonon spectral function from the e-ph interaction. 
  # 'wmin_specfun'      : 0.0,            # The lower boundary for the frequency in the electron spectral function in [eV].
  # 'wmax_specfun'      : 0.05,           # The upper boundary for the frequency in the electron spectral function in [eV]
  # 'nw_specfun'        : 100,            # Number of bins for frequency in electron spectral function.
  ################################
  'wannierize'        : True,
  'nbndsub'           : wannier90['num_wann'],
  # 'bands_skipped'     : 'exclude_bands = 1-10',
  'num_iter'          : wannier90['num_iter'],
  '!####### energy windows'             : '###########',
  'dis_froz_max'      : wannier90['dis_froz_max'],
  'dis_froz_min'      : wannier90['dis_froz_min'],
  '!######################'             : '###########',
  # 'wannier_plot'      : True,
  'wdata(1)'          : 'bands_plot = .true.',
  'wdata(2)'          : 'bands_num_points = 51',
  'wdata(3)'          : 'bands_plot_format = gnuplot',
  'wdata(4)'          : 'conv_tol      = 1.0e-10 ',
  'wdata(5)'          : 'dis_conv_tol  = 1.0e-10 ',
  'wdata(6)'          : 'dis_num_iter=' + str(wannier90['dis_num_iter']),
  'wdata(7)'          : 'kmesh_tol = 0.0001',
  # 'wdata(8)'          : 'guiding_centres = .true.',
}
nkf = [nk[0]*nktonkf, nk[1]*nktonkf, nk[2]*nktonkf]
nqf = [nq[0]*nqtonqf, nq[1]*nqtonqf, nq[2]*nqtonqf]
epw_nkq = {
  'nk1'               : nk[0],
  'nk2'               : nk[1],
  'nk3'               : nk[2],
  'nq1'               : nq[0],
  'nq2'               : nq[1],
  'nq3'               : nq[2], 

  'nkf1'              : nkf[0],
  'nkf2'              : nkf[1],
  'nkf3'              : nkf[2],
  'nqf1'              : nqf[0],
  'nqf2'              : nqf[1],
  'nqf3'              : nqf[2], 
}
#*****************************************************************************
#**************************** ShengBTE ***************************************
kpt_supercell = [3, 3, 3]              # kpoints for supercell calculation of 3rd
third_pp = '3 3 3 -5'
third_order_par = 'thirdorder_espresso.py scf.in sow ' + third_pp + ' scf_sc.in'
third_order_collect = 'find . -name \'DISP.scf_sc.in*log\' | sort -n | thirdorder_espresso.py scf.in reap ' + third_pp
fourth_pp = "3 3 3 -2"
fourth_order_par = 'Fourthorder_espresso.py scf.in sow ' + fourth_pp + ' scf_sc.in'
fourth_order_collect = 'find . -name \'DISP.scf_sc.in*log\' | sort -n | Fourthorder_espresso.py scf.in reap ' + fourth_pp
vasp = False                    # Use vasp or QE 
shengbte_epw =  False                 # 
shengbte_allocation = {
  'ngrid(:)'          : kpt,
  # 'norientations'     : 0,          # number of orientations along which to study nanowires
}
shengbte_crystal = {
  'lfactor'           : 0.1,          # unit nm 
  'scell(:)'          : str(phonon_input['nq1']) + " " + str(phonon_input['nq2']) + " " + str(phonon_input['nq3']),
  'born'              : 0,            # add born effective matrix, 0: not, 1: yes and give born parameters.
  #'orientations'      : ,            # terns of integer indices defining the crystallographic directions along which to study nanowires
}
shengbte_parameter = {
  # 'T'                 : 300,        # temperature
  'T_min'             : 200,
  'T_max'             : 500,
  'T_step'            : 50,
  'scalebroad'        : 1.0,          # default 1.0 
#  'maxiter'           : 1000,         # default 1000, maximum number of iterations allowed in the BTE convergence process
#  'nticks'            : 100,          # number of different values of the mean free path at which to compute the cumulative thermal conductivity
}
shengbte_flag = {
  'nonanalytic '      : True,         # default true
#  'isotopes'          : True,         # default true
#  'nanowires'         : False,        # default false
  'espresso'          : True,         # default false
  # 'onlyharmonic'      : False,
  # 'nthreads'          : 1,            # default 1, number of OpenMP threads each MPI process will use.
}
born_data = '''
  epsilon(:,1)=16.1676793 0.00000000 0.00000000,
  epsilon(:,2)=0.00000000 16.1676793 0.00000000,
  epsilon(:,3)=0.00000000 0.00000000 16.1676793,
  born(:,1,1)=2.1587779 0.0000000 0.0000000,
  born(:,2,1)=0.0000000 2.1587779 0.0000000,
  born(:,3,1)=0.0000000 0.0000000 2.1587779,
  born(:,1,2)=-2.1587779  0.0000000  0.0000000,
  born(:,2,2)= 0.0000000 -2.1587779  0.0000000,
  born(:,3,2)= 0.0000000  0.0000000 -2.1587779,
'''
#*****************************************************************************
#**************************** Perturbo ***************************************
qe2pert_in = {
  'prefix'           : 'pwscf',
  'outdir'           : './tmp',
# ‘phdir’ 存放处理后的声子的文件。计算完声子后请使用QEkit-> M4->P1处理。
  'phdir'            : '../phonon/pert_save',     
  'lwannier'         : True,     # 默认 True (通过wannier90计算e-ph matrix), 若改成False，则计算量很大。
# nk1, nk2, nk3, 设置与 nscf k点一致。
  'nk1'              : kpt_nscf[0],
  'nk2'              : kpt_nscf[1],
  'nk3'              : kpt_nscf[2],
# 'dft_band_min' and 'dft_band_max' determine the range of bands we are interested in, 
# and should be the same as the values used in the Wannierization process. 
# For example, if we used 40 bands in the nSCF calculation and we 
# excluded bands 1-4 and 31-40 in the Wannierization, then dft_band_min=5 and dft_band_max=30.
# 结论：可以理解为num_bands包含的能带数目，如果没有设置excludes_bands,则从1开始到nbnds结束。
  'dft_band_min'     : 1,                         
  'dft_band_max'     : nbnd,      
  'num_wann'         : wannier90['num_wann'],     # 同 wannier90里的num_wann.
  'load_ephmat'      : False,                     
  'system_2d'        : False,                     # 是否二维材料。
}
n_carriers = 0.9945847 * 10**18  # 0.9945847E+18
temperature = [200, 250, 300, 350, 400, 450, 500]    # temperature list
energy_fermi = 12.0626
perturbo_in = {
'tmp_dir'            : './tmp',
'hole'               : False,

'boltz_kdim(1)'      : 60, 
'boltz_kdim(2)'      : 60, 
'boltz_kdim(3)'      : 60,

# 'boltz_qdim(1)'      : 80,  # if not set, equal to kdim.
# 'boltz_qdim(2)'      : 80, 
# 'boltz_qdim(3)'      : 80,

'boltz_emin'         : 6.4,  
'boltz_emax'         : 6.9,
'band_min'           : 5,
'band_max'           : 6,
'boltz_nstep'        : 50,   # 0: RTA; > 0: ita
# 'phfreq_cutoff'      : 1,    # meV
# 'delta_smear'        : 10,   # meV
# 'trans_thr'          : 0.002,
'sampling'           : 'uniform',  # 'uniform', 'cauchy'
#'load_scatter_eph'   : True,
#'cauchy_scale'       :   ,  # if sampling = cauchy, set it.
# 'nsamples'           : 1000000,
}
