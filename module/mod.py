# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: input.py
"""
import os, sys
import numpy as np 
import configfile.inputpara as parameters
import module.funcs as bas
from module.writter import WriCls
from configfile.jobber import job_submit
import copy
sys.path.append(os.getcwd())
import AutoinputQE as inpp
sys.dont_write_bytecode = True       ## not generate __pycache__ files

class moduleslib():
    def __init__(self, filepath=os.getcwd()):
        self.filepath = filepath
        self.latt, self.element, self.num_element, self.tot_ele, self.cor = bas.get_poscar(self.filepath, "POSCAR")
        self.writter = WriCls(self.filepath)
        self.mag = bas.mag_judge()

    def ibrav_confirm_vaspkit(self):
        if inpp.ibrav != 0:
            try:
                os.system('cp POSCAR POSCAR_ori')
                os.system('echo -e "602\n " | vaspkit > kit')
                bas.trans_stdPOSCAR(self.filepath, "PRIMCELL.vasp")
            except:
                print("Not get standard primitive cell !")
            bas.transPOSCAR_QE(self.filepath, "POSCAR")
            os.system(parameters.path_QE + "/cell2ibrav.x < POSCAR-QE > cellinfo")
            os.system('rm PRIMCELL.vasp kit POSCAR-QE ') 
            self.ibrav = bas.get_ibrav(self.filepath, "cellinfo")
            self.celldm = bas.get_celldm(self.filepath, "cellinfo")
    def ibrav_confirm(self):
        if inpp.ibrav != 0:
            # not use vaspkit
            bas.transPOSCAR_QE(self.filepath, "POSCAR")
            os.system(parameters.path_QE + "/cell2ibrav.x < POSCAR-QE > cellinfo")      
            # os.system('rm PRIMCELL.vasp kit POSCAR-QE ') 
            self.ibrav = bas.get_ibrav(self.filepath, "cellinfo")
            self.celldm = bas.get_celldm(self.filepath, "cellinfo")
    def single(self, cal_type=inpp.input_control['calculation']):
        paras = {}
        control = copy.deepcopy(parameters.CONTROL)                ########### control parameters
        for keys_control, val_control in inpp.input_control.items():
            control[keys_control] = val_control
        control['calculation'] = cal_type
        paras['CONTROL'] = control 

        system = copy.deepcopy(parameters.SYSTEM)                   ########### system parameters
        for keys_system, val_system in inpp.input_system.items():
            system[keys_system] = val_system
        system['nat'] = self.tot_ele
        system['ntyp'] = len(self.element)
        try:
            if 'relax' in cal_type:
                del system['nbnd']
        except:
            print("Calculation: " + str(cal_type))
        # if inpp.SOC:
        #     system['noncolin'] = True 
        #     system['lspinorb'] = True
        ######################### ibrav confirm ###############
        # if 'relax' in inpp.input_control['calculation'] and inpp.ibrav==1:
        if inpp.ibrav==1:
            self.ibrav_confirm()
            system["ibrav"] = self.ibrav
            for keys, vals in self.celldm.items():
                system[keys] = vals
        else:
            system["ibrav"] = inpp.ibrav
        ######################### ibrav confirm ###############
        ############################# magmom ##################
        system = bas.mag_set(system)
        ############################# magmom ##################
        paras['SYSTEM'] = system

        electrons = copy.deepcopy(parameters.ELECTRONS)              ########### electron parameters
        for keys_electrons, val_electrons in inpp.input_electrons.items():
            electrons[keys_electrons] = val_electrons
        paras['ELECTRONS'] = electrons

        ions = copy.deepcopy(parameters.IONS)                        ########### ions parameters
        for keys_ions, val_ions in inpp.input_ions.items():
            ions[keys_ions] = val_ions 
        paras['IONS'] = ions

        cell = copy.deepcopy(parameters.CELL)                        ########### cell parameters
        for keys_cell, val_cell in inpp.input_cell.items():
            cell[keys_cell] = val_cell
        paras['CELL'] = cell
        return paras
    def single_modify(self, ibrav, cal_type):
        paras = {}
        control = copy.deepcopy(parameters.CONTROL)                ########### control parameters
        for keys_control, val_control in inpp.input_control.items():
            control[keys_control] = val_control
        control['calculation'] = cal_type
        paras['CONTROL'] = control 

        system = copy.deepcopy(parameters.SYSTEM)                   ########### system parameters
        for keys_system, val_system in inpp.input_system.items():
            system[keys_system] = val_system
        system['nat'] = self.tot_ele
        system['ntyp'] = len(self.element)
        ######################### ibrav confirm ###############
        # if 'relax' in inpp.input_control['calculation'] and inpp.ibrav==1:
        if ibrav==1:
            self.ibrav_confirm()
            system["ibrav"] = self.ibrav
            for keys, vals in self.celldm.items():
                system[keys] = vals
        else:
            system["ibrav"] = ibrav
        ######################### ibrav confirm ###############
        ############################# magmom ##################
        system = bas.mag_set(system)
        ############################# magmom ##################
        paras['SYSTEM'] = system

        electrons = parameters.ELECTRONS              ########### electron parameters
        for keys_electrons, val_electrons in inpp.input_electrons.items():
            electrons[keys_electrons] = val_electrons
        paras['ELECTRONS'] = electrons

        ions = parameters.IONS                        ########### ions parameters
        for keys_ions, val_ions in inpp.input_ions.items():
            ions[keys_ions] = val_ions 
        paras['IONS'] = ions

        cell = parameters.CELL                        ########### cell parameters
        for keys_cell, val_cell in inpp.input_cell.items():
            cell[keys_cell] = val_cell
        paras['CELL'] = cell
        return paras
    def config_single_modify(self, cal_type, ibrav, filepath=os.getcwd()):
        paras = copy.deepcopy(self.single_modify(ibrav, cal_type))
        self.writter.QE_input_RelaxScfNscf_ibrav0(cal_type, paras)
        job_submit(filepath).job_notAuto(cal_type) 
    def config_single(self, cal_type=inpp.input_control['calculation'], filepath=os.getcwd()):
        paras = copy.deepcopy(self.single(cal_type))
        self.writter.QE_input_RelaxScfNscf(cal_type, paras)
        job_submit(filepath).job_notAuto(cal_type)
    def config_AutoBands(self):
        self.config_single('scf')
        self.config_single('bands')
        self.writter.QE_input_postband(parameters)
        job_submit(self.filepath).job_AutoBand()
        
    def config_AutoNscf(self):
        self.config_single('scf')
        self.config_single('nscf')
        job_submit(self.filepath).job_AutoNscf()
    def config_AutoNscf_EPW(self, filepath):
        self.config_single('scf')
        paras = copy.deepcopy(self.single('nscf'))
        self.writter.QE_input_RelaxScfNscf2(paras, filepath)
        job_submit(filepath).job_AutoNscf()
    def config_AutoBandNscf_EPW(self, filepath):
        self.config_single('scf')
        self.config_single('bands')
        self.writter.QE_input_postband(parameters)
        paras = copy.deepcopy(self.single('nscf'))
        self.writter.QE_input_RelaxScfNscf2(paras, filepath)
        job_submit(filepath).job_AutoBandsNscf()
    def single_epw_redefined(self, cal_type=inpp.input_control['calculation']):
        paras = {}
        control = copy.deepcopy(parameters.CONTROL)                ########### control parameters
        for keys_control, val_control in inpp.input_control.items():
            control[keys_control] = val_control
        control['calculation'] = cal_type
        paras['CONTROL'] = control 

        system = copy.deepcopy(parameters.SYSTEM)                   ########### system parameters
        for keys_system, val_system in inpp.input_system.items():
            system[keys_system] = val_system
        system['nat'] = self.tot_ele
        system['ntyp'] = len(self.element)
        try:
            if 'relax' in cal_type:
                del system['nbnd']
        except:
            print("Calculation: " + str(cal_type))
        # if inpp.SOC:
        #     system['noncolin'] = True 
        #     system['lspinorb'] = True
        ######################### ibrav confirm ###############
        # if 'relax' in inpp.input_control['calculation'] and inpp.ibrav==1:
        if inpp.ibrav==1:
            self.ibrav_confirm()
            system["ibrav"] = self.ibrav
            for keys, vals in self.celldm.items():
                system[keys] = vals
        else:
            system["ibrav"] = inpp.ibrav
        ######################### ibrav confirm ###############
        ############################# magmom ##################
        system = bas.mag_set(system)
        ############################# magmom ##################

        paras['SYSTEM'] = system

        electrons = copy.deepcopy(parameters.ELECTRONS)              ########### electron parameters
        for keys_electrons, val_electrons in inpp.input_electrons.items():
            electrons[keys_electrons] = val_electrons
        paras['ELECTRONS'] = electrons

        ions = copy.deepcopy(parameters.IONS)                        ########### ions parameters
        for keys_ions, val_ions in inpp.input_ions.items():
            ions[keys_ions] = val_ions 
        paras['IONS'] = ions

        cell = copy.deepcopy(parameters.CELL)                        ########### cell parameters
        for keys_cell, val_cell in inpp.input_cell.items():
            cell[keys_cell] = val_cell
        paras['CELL'] = cell
        return paras
    def config_AutoBandNscf_EPW_SC(self, filepath):
        self.config_single('scf')
        self.config_single('bands')
        self.writter.QE_input_postband(parameters)
        paras = copy.deepcopy(self.single('nscf'))
        self.writter.QE_input_RelaxScfNscf2(paras, filepath)
        job_submit(filepath).job_AutoBandsNscf()
    def config_AutoDos(self):
        self.config_single('scf')
        self.config_single('nscf')
        paras = copy.deepcopy(parameters.dos_input)
        for keys_control, val_control in inpp.dos.items():
            paras[keys_control] = val_control
        self.writter.QE_dos(paras, self.filepath)
        job_submit(self.filepath).job_AutoDos()
    def config_AutoW90(self):
        self.config_single('scf')
        paras = copy.deepcopy(self.single('nscf'))
        self.writter.QE_input_RelaxScfNscf2(paras, os.getcwd())
        self.config_single('bands')
        self.writter.QE_input_postband(parameters)
        self.writter.w90(parameters.w90, self.filepath)
        paras = copy.deepcopy(parameters.wannier90_para)
        for keys_control, val_control in inpp.wannier90.items():
            paras[keys_control] = val_control
        if inpp.SOC or inpp.noncolin:
            paras['spinors'] = True
        if self.mag == 2:
            para = {}
            para['up'] = paras
            para['dw'] = paras
            self.writter.w90_win(para, self.filepath)
        else:
            self.writter.w90_win(paras, self.filepath)
        job_submit(self.filepath).job_AutoW90()
    def config_w90(self):
        self.writter.w90(parameters.w90, self.filepath)
        paras = copy.deepcopy(parameters.wannier90_para)
        for keys_control, val_control in inpp.wannier90.items():
            paras[keys_control] = val_control
        if inpp.SOC or inpp.noncolin:
            paras['spinors'] = True
        if self.mag == 2:
            para = {}
            para['up'] = paras
            para['dw'] = paras
            self.writter.w90_win(para, self.filepath)
        else:
            self.writter.w90_win(paras, self.filepath)
        job_submit(self.filepath).job_AutoW90()
    def config_distangle(self):
        os.system('cp ../POSCAR .')
        if self.mag == 2:
        # os.system('ln -s ../AutoinputQE.py .')
            os.system('cp ../pwscf_up.win .')
            os.system('cp ../pwscf_dw.win .')
        else:
            os.system('cp ../pwscf.win .')
        job_submit(os.getcwd()).job_distangle()
    def config_AutoPhon(self, filepath=os.getcwd()):
        self.config_single('scf')
        # self.config_single('nscf')
        ########## phonon 
        paras_phon = copy.deepcopy(parameters.PHONON)      
        for keys, vals in inpp.phonon_input.items():
            paras_phon[keys] = vals
        for i in range(len(self.element)):
            paras_phon['amass('+ str(i+1) +')'] = parameters.atom_mass[self.element[i]]
        self.writter.QE_phonon(paras_phon, filepath)
        ######### q2r
        paras_q2r = copy.deepcopy(parameters.PHONON_q2r)
        for keys_q2r, vals_q2r in inpp.phon_q2r.items():
            paras_phon[keys_q2r] = vals_q2r
        if inpp.SOC or inpp.noncolin:
            paras_phon['flfrc'] = 'pwscf.fc.xml'
        self.writter.QE_q2r(paras_q2r, filepath)
        ######### matdyn
        paras_matdyn = copy.deepcopy(parameters.PHONON_mat)
        for keys_mat, val_mat in inpp.phonon_mat_input.items():
            paras_matdyn[keys_mat] = val_mat
        if inpp.SOC or inpp.noncolin:
            paras_matdyn['flfrc'] = 'pwscf.fc.xml'
        self.writter.QE_matdyn(paras_matdyn, filepath)
        ######### 
        paras_dos = copy.deepcopy(parameters.PHONON_dos)
        for key_dos, val_dos in inpp.phon_dos.items():
            paras_dos[key_dos] = val_dos 
        self.writter.QE_phdos(paras_dos, filepath)
        job_submit(self.filepath).job_AutoPhon()
    def config_AutoPhon_TEST(self, filepath=os.getcwd()):
        paras = copy.deepcopy(self.single('scf'))
        self.writter.QE_input_RelaxScfNscf_kpt('scf', paras)
        # self.config_single('nscf')
        ########## phonon 
        paras_phon = copy.deepcopy(parameters.PHONON)      
        for keys, vals in inpp.phonon_input.items():
            paras_phon[keys] = vals
        for i in range(len(self.element)):
            paras_phon['amass('+ str(i+1) +')'] = parameters.atom_mass[self.element[i]]
        self.writter.QE_phonon(paras_phon, filepath)
        ########## phonon2 
        paras_phon2 = copy.deepcopy(parameters.PHONON)      
        for keys, vals in inpp.phonon_input.items():
            paras_phon2[keys] = vals
        # for i in range(len(self.element)):
        #     paras_phon2['amass('+ str(i+1) +')'] = parameters.atom_mass[self.element[i]]
        try:
            del paras_phon2['verbosity']
            del paras_phon2['recover']
            del paras_phon2['diagonalization']
        except:
            print('Go ON !')
        paras_phon2['electron_phonon'] = 'epw'
        paras_phon2['q_in_band_form'] = True
        paras_phon2['qplot'] = True
        paras_phon2['fildyn'] = 'eph.dyn'
        try:
            paras_phon2['kx'] = inpp.eph_k[0]
            paras_phon2['ky'] = inpp.eph_k[1]
            paras_phon2['kz'] = inpp.eph_k[2]
        except:
            print("Not specify k grid for electron phonon calculation !")
        self.writter.QE_phononTEST(paras_phon2, filepath)
        ######### q2r
        paras_q2r = copy.deepcopy(parameters.PHONON_q2r)
        for keys_q2r, vals_q2r in inpp.phon_q2r.items():
            paras_phon[keys_q2r] = vals_q2r
        if inpp.SOC or inpp.noncolin:
            paras_phon['flfrc'] = 'pwscf.fc.xml'
        self.writter.QE_q2r(paras_q2r, filepath)
        ######### matdyn
        paras_matdyn = copy.deepcopy(parameters.PHONON_mat)
        for keys_mat, val_mat in inpp.phonon_mat_input.items():
            paras_matdyn[keys_mat] = val_mat
        if inpp.SOC or inpp.noncolin:
            paras_matdyn['flfrc'] = 'pwscf.fc.xml'
        self.writter.QE_matdyn(paras_matdyn, filepath)
        job_submit(self.filepath).job_AutoPhonTest()
    def config_AutoPhon_Pertubo(self, filepath=os.getcwd()):
        self.config_single('scf')
        # self.config_single('nscf')
        ########## phonon 
        paras_phon = copy.deepcopy(parameters.PHONON)      
        for keys, vals in inpp.phonon_input.items():
            paras_phon[keys] = vals
        for i in range(len(self.element)):
            paras_phon['amass('+ str(i+1) +')'] = parameters.atom_mass[self.element[i]]
        # paras_phon['verbosity'] = 'debug'
        paras_phon['fildyn'] = 'pwscf.dyn.xml'
        self.writter.QE_phonon(paras_phon, filepath)
        ######### q2r
        paras_q2r = copy.deepcopy(parameters.PHONON_q2r)
        for keys_q2r, vals_q2r in inpp.phon_q2r.items():
            paras_phon[keys_q2r] = vals_q2r
        paras_q2r['fildyn'] = "pwscf.dyn.xml"
        paras_q2r['flfrc'] = 'pwscf.fc'
        self.writter.QE_q2r(paras_q2r, filepath)
        ######### matdyn
        paras_matdyn = copy.deepcopy(parameters.PHONON_mat)
        for keys_mat, val_mat in inpp.phonon_mat_input.items():
            paras_matdyn[keys_mat] = val_mat
        paras_matdyn['flfrc'] = 'pwscf.fc.xml'
        self.writter.QE_matdyn(paras_matdyn, filepath)
        ######### 
        paras_dos = copy.deepcopy(parameters.PHONON_dos)
        for key_dos, val_dos in inpp.phon_dos.items():
            paras_dos[key_dos] = val_dos 
        paras_dos['flfrc'] = 'pwscf.fc.xml'
        self.writter.QE_phdos(paras_dos, filepath)
        job_submit(self.filepath).job_AutoPhon_perturbo()

    def config_Epw_prepare(self, paras):
        epw = copy.deepcopy(parameters.epw)
        ''' configure the atoms 
            mass for all atoms'''
        for atom_flag in range(len(self.element)):
            epw['amass'+'(' + str(atom_flag+1) + ')'] = parameters.atom_mass[self.element[atom_flag]]
        i = 1
        for ele in range(len(self.element)):
            keys_ele = 'amass(' + str(ele+1) + ')'
            vals_ele = parameters.atom_mass[self.element[ele]]
            epw[keys_ele] = vals_ele
        ################ SOC
        if inpp.SOC or inpp.noncolin:
            epw['wdata(' + str(i+1) + ')'] = 'spinors = true'
            i = i + 1
        ################
        for keys, vals in paras['epw_input'].items():
            if 'wdata' in keys:
                keys = 'wdata(' + str(i) + ')'
                i = i + 1
            epw[keys] = vals
        i = i - 1
        if inpp.use_projection_w90:
            ''' configure high symmetry points '''
            kpoints_band = [items.split("!")[0] for items in inpp.kpoints_band.split("\n")[1:-1]]
            hsp_labels = [str(items.split("!")[1]).strip() for items in inpp.kpoints_band.split("\n")[1:-1]]
            kpoints_band_tmp = [[] for i in range(len(kpoints_band))]
            for kpp in range(len(kpoints_band)):
                for pp in range(3):
                    kpoints_band_tmp[kpp].append(np.round(float(kpoints_band[kpp].split()[pp]), 5))
            kpoints_bands = []
            for kkp in range(len(kpoints_band)):
                ktmp = " "
                for kp in range(3):
                    ktmp = ktmp + "  " + str(kpoints_band_tmp[kkp][kp])
                kpoints_bands.append(ktmp)
            for ii in range(len(hsp_labels)):
                if str(hsp_labels[ii]).lower() == "\Gamma".lower():
                    hsp_labels[ii] = "G"
            epw['wdata(' + str(i+1)+')'] = 'begin kpoint_path'      
            for j in range(0, len(kpoints_bands)-1):
                hsp_keys = 'wdata(' + str(i+j+2) + ')'
                hsp_vals = hsp_labels[j] + "  " + kpoints_bands[j] + "  " + hsp_labels[j+1] + "  " + kpoints_bands[j+1]
                epw[hsp_keys] = hsp_vals
            epw['wdata(' + str(i+len(kpoints_bands)+2)+')'] = 'end kpoint_path'
        else:
            print("User specifie wannier90 parameters !")
        ###### projections for wannier90
        proj = inpp.porjection_epw
        pp = 0
        for keyss, valss in proj.items():
            key_proj = 'proj(' + str(pp + 1) + ')'
            val_proj = keyss + ' : ' + valss
            epw[key_proj] = val_proj
            pp = pp + 1
        for keys_nkq, val_nkq in paras['epw_nkq'].items():
            epw[keys_nkq] = val_nkq
        try:
            if paras['epw_input']['efermi_read']:
                fermi_eng = bas.get_fermi_level(os.getcwd(), "scf.log")
                print("Fermi Energy: " + str(fermi_eng))
                epw['fermi_energy'] = fermi_eng
            else:
                del epw['efermi_read']
                if 'fermi_energy' in epw:
                    del epw['fermi_energy']
        except:
            print("Not obtain Fermi energy, need set fermi_energy by yourself !")
        return epw
    def config_Epw(self, paras, filepath, filename):
        epwpar = copy.deepcopy(self.config_Epw_prepare(paras))
        self.writter.epw(epwpar, filepath, filename)
        job_submit(self.filepath).job_epw(filepath, filename)

    def config_Epw_single(self, paras, filepath, filename):
        epwpar = copy.deepcopy(self.config_Epw_prepare(paras))
        self.writter.epw(epwpar, filepath, filename)

    def config_qe2pert(self, paras, filepath, filename):
        self.writter.qe2pert(paras, filepath, filename)
        job_submit(filepath).job_qe2pert(filepath, filename)

    def config_pert(self, paras, filepath, filename):
        self.writter.perturbo(paras, filepath, filename)
        job_submit(filepath).job_pert(filepath, filename)

    def config_pert_ephmat(self, paras, filepath, filename):
        self.writter.perturbo(paras, filepath, filename)
        job_submit(filepath).job_pert_ephmat(filepath, filename)

    def config_pert_trans(self, paras, filepath, filename):
        self.writter.perturbo(paras, filepath, filename)
        job_submit(filepath).job_pert_trans(filepath, filename)

    def config_AutoPhon_phonopy(self, filepath):
        self.config_single_modify('scf', 0)
        paras = {}
        control = parameters.CONTROL                ########### control parameters
        for keys_control, val_control in inpp.input_control.items():
            control[keys_control] = val_control
        control['calculation'] = 'scf'
        control['disk_io'] = 'nowf'
        paras['CONTROL'] = control 
        system = parameters.SYSTEM                   ########### system parameters
        for keys_system, val_system in inpp.input_system.items():
            system[keys_system] = val_system
        system['nat'] = self.tot_ele * inpp.phonon_input['nq1']*inpp.phonon_input['nq2']*inpp.phonon_input['nq3']
        system['ntyp'] = len(self.element)
        
        system = bas.mag_set(system)
        ############################# magmom ##################
        try:
            del system['nbnd']
        except:
            print("No nbnd parameter in QE, go on !")
        system['ibrav'] = 0
        paras['SYSTEM'] = system

        electrons = parameters.ELECTRONS              ########### electron parameters
        for keys_electrons, val_electrons in inpp.input_electrons.items():
            electrons[keys_electrons] = val_electrons
        paras['ELECTRONS'] = electrons

        ions = parameters.IONS                        ########### ions parameters
        for keys_ions, val_ions in inpp.input_ions.items():
            ions[keys_ions] = val_ions 
        paras['IONS'] = ions

        cell = parameters.CELL                        ########### cell parameters
        for keys_cell, val_cell in inpp.input_cell.items():
            cell[keys_cell] = val_cell
        paras['CELL'] = cell
        self.writter.phonopy_2nd(paras, filepath)
        bas.KPATH_Phonopy(os.getcwd(), "KPATH.phonopy")
        os.system("phonopy --qe -d --dim=\"" + str(inpp.phonon_input['nq1']) + ' ' + 
                    str(inpp.phonon_input['nq2']) + ' ' + str(inpp.phonon_input['nq3']) + "\" -c scf.in")
    def shengbte_3rd(self, filepath):
        self.config_single_modify('scf', 0)
        paras = {}
        control = parameters.CONTROL                ########### control parameters
        for keys_control, val_control in inpp.input_control.items():
            control[keys_control] = val_control
        control['calculation'] = 'scf'
        control['disk_io'] = 'nowf'
        paras['CONTROL'] = control 

        system = parameters.SYSTEM                   ########### system parameters
        for keys_system, val_system in inpp.input_system.items():
            system[keys_system] = val_system
        system['nat'] = '##NATOMS##'
        system['ntyp'] = len(self.element)
        
        ############################# magmom ##################
        # if inpp.ISPIN == 2 and inpp.SOC:
        #     system['nspin'] = 2
        #     system['noncolin'] = True
        #     system['lspinorb'] = True
        #     for key_mag, val_mag in inpp.Magmom.items():
        #         system['starting_magnetization'+"("+ key_mag + ")"] = val_mag[0]
        #         system['angle1(' + key_mag + ')'] = val_mag[1]
        #         system['angle2(' + key_mag + ')'] = val_mag[2]
        ############################# magmom ##################
        system = bas.mag_set(system)
        ############################# magmom ##################
        try:
            del system['nbnd']
        except:
            print("No nbnd parameter in QE, go on !")
        system['ibrav'] = 0
        paras['SYSTEM'] = system

        electrons = parameters.ELECTRONS              ########### electron parameters
        for keys_electrons, val_electrons in inpp.input_electrons.items():
            electrons[keys_electrons] = val_electrons
        paras['ELECTRONS'] = electrons

        ions = parameters.IONS                        ########### ions parameters
        for keys_ions, val_ions in inpp.input_ions.items():
            ions[keys_ions] = val_ions 
        paras['IONS'] = ions

        cell = parameters.CELL                        ########### cell parameters
        for keys_cell, val_cell in inpp.input_cell.items():
            cell[keys_cell] = val_cell
        paras['CELL'] = cell
        self.writter.ShengBTE_3rd(paras, filepath)

    def shengbte_4th(self, filepath):
        self.config_single_modify('scf', 0)
        paras = {}
        control = parameters.CONTROL                ########### control parameters
        for keys_control, val_control in inpp.input_control.items():
            control[keys_control] = val_control
        control['calculation'] = 'scf'
        control['disk_io'] = 'nowf'
        paras['CONTROL'] = control 

        system = parameters.SYSTEM                   ########### system parameters
        for keys_system, val_system in inpp.input_system.items():
            system[keys_system] = val_system
        system['nat'] = '##NATOMS##'
        system['ntyp'] = len(self.element)
        
        ############################# magmom ##################
        # if inpp.ISPIN == 2 and inpp.SOC:
        #     system['nspin'] = 2
        #     system['noncolin'] = True
        #     system['lspinorb'] = True
        #     for key_mag, val_mag in inpp.Magmom.items():
        #         system['starting_magnetization'+"("+ key_mag + ")"] = val_mag[0]
        #         system['angle1(' + key_mag + ')'] = val_mag[1]
        #         system['angle2(' + key_mag + ')'] = val_mag[2]
        ############################# magmom ##################
        system = bas.mag_set(system)
        ############################# magmom ##################
        try:
            del system['nbnd']
        except:
            print("No nbnd parameter in QE, go on !")
        system['ibrav'] = 0
        paras['SYSTEM'] = system

        electrons = parameters.ELECTRONS              ########### electron parameters
        for keys_electrons, val_electrons in inpp.input_electrons.items():
            electrons[keys_electrons] = val_electrons
        paras['ELECTRONS'] = electrons

        ions = parameters.IONS                        ########### ions parameters
        for keys_ions, val_ions in inpp.input_ions.items():
            ions[keys_ions] = val_ions 
        paras['IONS'] = ions

        cell = parameters.CELL                        ########### cell parameters
        for keys_cell, val_cell in inpp.input_cell.items():
            cell[keys_cell] = val_cell
        paras['CELL'] = cell
        self.writter.ShengBTE_3rd(paras, filepath)

    def shengbte(self, filepath):
        paras = {}
        allocation = parameters.shengbte_allocations
        allocation['nelements'] = len(self.element)
        allocation['natoms'] = self.tot_ele
        for key_allo, val_allo in inpp.shengbte_allocation.items():
            allocation[key_allo] = val_allo
        ngrid = str(inpp.shengbte_allocation['ngrid(:)'][0])  + " " + str(inpp.shengbte_allocation['ngrid(:)'][1]) + " " + str(inpp.shengbte_allocation['ngrid(:)'][2])
        allocation['ngrid(:)'] = ngrid 
        paras['allocations'] = allocation

        crystal = parameters.shengbte_crystals
        for key_cry, val_cry in inpp.shengbte_crystal.items():
            crystal[key_cry] = val_cry
        del crystal['scell(:)']
        for i in range(3):
            crystal['lattvec(:,' + str(i+1) + ')'] = str(self.latt[i, 0]) + "  " + str(self.latt[i, 1]) + "  " + str(self.latt[i, 2])
        tmp = []
        typ = ''
        for j in range(len(self.element)):
            # tmp = tmp + self.element[j] + ' '
            tmp.append(self.element[j])
            for k in range(self.num_element[j]):
                typ = typ + str(j+1) + " "
        crystal['elements'] = tmp
        crystal['types'] = typ
            
        for ii in range(self.tot_ele):
            crystal['positions(:,'+str(ii+1)+')'] = str(self.cor[ii, 0]) + '  ' + str(self.cor[ii, 1]) + ' ' + str(self.cor[ii, 2])
        if inpp.shengbte_crystal['born'] != 0:
            crystal['borndata'] = inpp.born_data
        del crystal['born']
        crystal['scell(:)'] = inpp.shengbte_crystal['scell(:)']
        paras['crystals'] = crystal

        par = parameters.shengbte_parameters
        for key_par, val_par in inpp.shengbte_parameter.items():
            par[key_par] = val_par
        paras['parameter'] = par

        flags = parameters.shengbte_flags
        for key_flg, val_flg in inpp.shengbte_flag.items():
            flags[key_flg] = val_flg
        if inpp.shengbte_epw:
            flags['electrons'] = True
        try:
            if inpp.vasp:
                flags['espresso'] = False 
        except:
            print("CONTROL file for QE")
        paras['flags'] = flags
        self.writter.ShengBTE(paras, filepath)
        job_submit(filepath).job_ShengBTE()

