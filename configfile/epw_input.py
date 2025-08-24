# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: input.py
"""
import os.path as osp 
import os, sys
import configfile.inputpara as parameters
import module.mod as modd
import copy
sys.path.append(os.getcwd())
import AutoinputQE as inpp
sys.dont_write_bytecode = True        ## not generate __pycache__ files

class epwinput():
    def __init__(self, filepath=os.getcwd()):
        self.filepath = filepath
        self.cal = set(['elph', 'epbwrite', 'epbread', 'epwwrite', 'epwread', 'scattering',
                        'ep_coupling', 'phonselfen', 'elecselfen', 'a2f', 
                        'cumulant', 'wannierize', 'ephwrite', 'eliashberg', 'specfun_el', 'specfun_ph',
                        'specfun_pl'])
        self.nks = inpp.phon_points
    def cal_control(self, cal_use, input_epw):
        cal = copy.deepcopy(self.cal)
        calculation_exclude = list(cal - cal_use)
        epwio = copy.deepcopy(inpp.epw_input)
        for i in range(len(calculation_exclude)):
            for keys, vals in epwio.items():
                if calculation_exclude[i] in keys:
                    del input_epw[calculation_exclude[i]]
        return input_epw
    
    def check_bands(self):    # the first step for epw calculation.
        paras_bd = {}
        cal_use_bands = set(['elph', 'epwwrite', 'wannierize'])
        calculation_bd = list(cal_use_bands)
        input_epw_bands = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            # if 'wdata' in key_input:
            #     continue
            # else:
            #     input_epw[key_input] = val_input
            input_epw_bands[key_input] = val_input
        for i in range(len(calculation_bd)):
            input_epw_bands[calculation_bd[i]] = True
        input_epw_bands['band_plot'] = True
        # input_epw_bands['wannier_plot'] = True
        input_epw_bands['filkf'] = './meshes/kpoints'
        input_epw_bands['filqf'] = './meshes/kpoints'
        input_epw_bands.pop('fsthick', 1)
        input_epw_bands.pop('temps', 1)
        input_epw_bands.pop('degaussw', 1)
        input_epw_bands.pop('efermi_read', 1)
        input_epw_bands.pop('fermi_energy', 1)
        try:
            input_epw_bands = copy.deepcopy(self.cal_control(cal_use_bands, input_epw_bands))
        except:
            print("Go on !")
        # try:
        #     # del input_epw_bands['fsthick']
        #     # del input_epw_bands['temps']
        #     del input_epw_bands['degaussw']
        # except:
        #     print("Go on !")
        paras_bd['epw_input'] = input_epw_bands
        epwnkfq = copy.deepcopy(inpp.epw_nkq) 
        epwnkfq.pop('nkf1', 1)
        epwnkfq.pop('nkf2', 1)
        epwnkfq.pop('nkf3', 1)
        epwnkfq.pop('nqf1', 1)
        epwnkfq.pop('nqf2', 1)
        epwnkfq.pop('nqf3', 1)
        paras_bd['epw_nkq'] = epwnkfq
        return paras_bd

    def Test_gavg(self):    # Test |g| for \Gamma point
        paras_lp = {}
        cal_use_ph = set(['elph', 'epwwrite', 'wannierize'])
        calculation_lp = list(cal_use_ph)
        input_epw_lp = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            input_epw_lp[key_input] = val_input
        for i in range(len(calculation_lp)):
            input_epw_lp[calculation_lp[i]] = True
        try:
            input_epw_lp = copy.deepcopy(self.cal_control(cal_use_ph, input_epw_lp))
        except:
            print("Go on !")
        # input_epw_lp['band_plot'] = False
        input_epw_lp['filqf'] = './meshes/kpoints'
        input_epw_lp['prtgkk'] = True
        # input_epw['delta_approx'] = True 
        try:
            del input_epw_lp['filkf']
        except:
            print("Go on !")
        paras_lp['epw_input'] = input_epw_lp
        epwnkfq = copy.deepcopy(inpp.epw_nkq) 
        epwnkfq.pop('nqf1', 1)
        epwnkfq.pop('nqf2', 1)
        epwnkfq.pop('nqf3', 1)
        epwnkfq['nkf1'] = 1
        epwnkfq['nkf2'] = 1
        epwnkfq['nkf3'] = 1
        paras_lp['epw_nkq'] = epwnkfq
        return paras_lp

    def epw_linewidth_phonon(self):    # phonon linewidth calculation
        paras_lp = {}
        cal_use_ph = set(['elph', 'epwread', 'phonselfen'])
        calculation_lp = list(cal_use_ph)
        input_epw_lp = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            input_epw_lp[key_input] = val_input
        for i in range(len(calculation_lp)):
            input_epw_lp[calculation_lp[i]] = True
        try:
            input_epw_lp = copy.deepcopy(self.cal_control(cal_use_ph, input_epw_lp))
        except:
            print("Go on !")
        input_epw_lp['band_plot'] = False
        input_epw_lp['epwwrite'] = False 
        input_epw_lp['filqf'] = './meshes/kpoints'
        # input_epw['delta_approx'] = True 
        try:
            del input_epw_lp['filkf']
        except:
            print("Go on !")
        paras_lp['epw_input'] = input_epw_lp
        epwnkfq = copy.deepcopy(inpp.epw_nkq) 
        epwnkfq.pop('nqf1', 1)
        epwnkfq.pop('nqf2', 1)
        epwnkfq.pop('nqf3', 1)
        paras_lp['epw_nkq'] = epwnkfq
        return paras_lp

    def epw_phonon_spectraFunc(self):    # phonon spectra function calculation.
        paras_lp = {}
        cal_use_ph = set(['elph', 'epwread', 'phonselfen', 'specfun_ph'])
        calculation_lp = list(cal_use_ph)
        input_epw_lp = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            input_epw_lp[key_input] = val_input
        for i in range(len(calculation_lp)):
            input_epw_lp[calculation_lp[i]] = True
        try:
            input_epw_lp = copy.deepcopy(self.cal_control(cal_use_ph, input_epw_lp))
        except:
            print("Go on !")
        input_epw_lp['band_plot'] = False
        input_epw_lp['epwwrite'] = False 
        input_epw_lp['filqf'] = './meshes/kpoints'
        # input_epw['delta_approx'] = True 
        try:
            del input_epw_lp['filkf']
        except:
            print("Go on !")
        paras_lp['epw_input'] = input_epw_lp
        epwnkfq = copy.deepcopy(inpp.epw_nkq) 
        epwnkfq.pop('nqf1', 1)
        epwnkfq.pop('nqf2', 1)
        epwnkfq.pop('nqf3', 1)
        paras_lp['epw_nkq'] = epwnkfq
        return paras_lp

    def epw_linewidth_electron(self):    # electron linewidth calculation.
        paras_le = {}
        cal_use_el = set(['elph', 'epwread', 'elecselfen'])
        calculation_le = list(cal_use_el)
        input_epw_le = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            input_epw_le[key_input] = val_input
        for i in range(len(calculation_le)):
            input_epw_le[calculation_le[i]] = True
        try:
            input_epw_le = copy.deepcopy(self.cal_control(cal_use_el, input_epw_le))
        except:
            print("Go on !")
        input_epw_le['band_plot'] = False
        input_epw_le['epwwrite'] = False 
        input_epw_le['filkf'] = './meshes/kpoints'
        # input_epw['delta_approx'] = True 
        try:
            del input_epw_le['filqf']
        except:
            print("Go on !")
        paras_le['epw_input'] = input_epw_le
        epwnkfq = copy.deepcopy(inpp.epw_nkq) 
        epwnkfq.pop('nkf1', 1)
        epwnkfq.pop('nkf2', 1)
        epwnkfq.pop('nkf3', 1)
        paras_le['epw_nkq'] = epwnkfq
        return paras_le

    def epw_electron_spectraFun(self):    # electron spectral function calculation.
        paras_le = {}
        cal_use_el = set(['elph', 'epwread', 'elecselfen', 'specfun_el'])
        calculation_le = list(cal_use_el)
        input_epw_le = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            input_epw_le[key_input] = val_input
        for i in range(len(calculation_le)):
            input_epw_le[calculation_le[i]] = True
        try:
            input_epw_le = copy.deepcopy(self.cal_control(cal_use_el, input_epw_le))
        except:
            print("Go on !")
        input_epw_le['band_plot'] = False
        input_epw_le['epwwrite'] = False 
        input_epw_le['filkf'] = './meshes/kpoints'
        # input_epw['delta_approx'] = True 
        try:
            del input_epw_le['filqf']
        except:
            print("Go on !")
        paras_le['epw_input'] = input_epw_le
        epwnkfq = copy.deepcopy(inpp.epw_nkq) 
        epwnkfq.pop('nkf1', 1)
        epwnkfq.pop('nkf2', 1)
        epwnkfq.pop('nkf3', 1)
        paras_le['epw_nkq'] = epwnkfq
        return paras_le

    def epw_Auto_cal(self):
        paras_checkbands = self.check_bands()
        paras_lp = self.epw_linewidth_phonon()
        paras_le = self.epw_linewidth_electron()
        paras_all = {
            'epw_band.in'           : paras_checkbands,
            'epw_phonon.in'         : paras_lp,
            'epw_electron.in'       : paras_le,
        }
        return paras_all

    def epw_cal_user_defined(self):    # user defined epw calculation.
        paras = {}
        input_epw = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            input_epw[key_input] = val_input
        # input_epw['delta_approx'] = True 
        paras['epw_input'] = input_epw
        epwnkfq = copy.deepcopy(inpp.epw_nkq)            #### only nqf
        paras['epw_nkq'] = epwnkfq
        return paras

    ##################### resistivity ##################################
    def resistivity_ziman(self):    # electron resistivity calculation.
        paras_resistivity = {}
        cal_use_resistivity = set(['elph', 'epwwrite', 'epbwrite', 'wannierize', 'phonselfen', 'a2f'])
        calculation_resistivity = list(cal_use_resistivity)
        input_epw_res = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            input_epw_res[key_input] = val_input
        for i in range(len(calculation_resistivity)):
            input_epw_res[calculation_resistivity[i]] = True
        try:
            input_epw_res = copy.deepcopy(self.cal_control(cal_use_resistivity, input_epw_res))
        except:
            print("Go on !")
        input_epw_res['vme'] = 'wannier'
        input_epw_res['band_plot'] = False
        # input_epw_res['epwwrite'] = False 
        input_epw_res['nc'] = 4.0
        input_epw_res['assume_metal'] = True
        input_epw_res['degaussq'] = 0.05 
        input_epw_res['ngaussw'] = -99
        input_epw_res['vme'] = 'wannier' 
        input_epw_res['delta_approx'] = True 
        paras_resistivity['epw_input'] = input_epw_res
        epwnkfq = copy.deepcopy(inpp.epw_nkq) 
        paras_resistivity['epw_nkq'] = epwnkfq
        return paras_resistivity
    def resistivity_BTE(self):    # electron resistivity calculation.
        paras_resistivity = {}
        cal_use_resistivity = set(['elph', 'epwwrite', 'epbwrite', 'scattering', 'wannierize'])
        calculation_resistivity = list(cal_use_resistivity)
        input_epw_res = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            input_epw_res[key_input] = val_input
        for i in range(len(calculation_resistivity)):
            input_epw_res[calculation_resistivity[i]] = True
        try:
            input_epw_res = copy.deepcopy(self.cal_control(cal_use_resistivity, input_epw_res))
        except:
            print("Go on !")
        input_epw_res['band_plot'] = False
        input_epw_res['vme'] = 'wannier' 
        # input_epw_res['epwwrite'] = False 
        input_epw_res['nstemp'] = 9
        input_epw_res['temps'] = [100, 500] 
        input_epw_res['assume_metal'] = False
        input_epw_res['degaussw'] = 0.0
        input_epw_res['ngaussw'] = -99
        input_epw_res['int_mob'] = True
        input_epw_res['iterative_bte'] = True 
        input_epw_res['carrier'] = False 
        input_epw_res['mp_mesh_k'] = True  
        input_epw_res['epmatkqread'] = False 
        input_epw_res['mob_maxiter'] = 300 
        input_epw_res['broyden_beta'] = 0.7
        # input_epw_res['restart'] = True    
        # input_epw_res['restrat_step'] = 50 
        input_epw_res['selecqread'] = False 
        # input_epw['delta_approx'] = True 
        try:
            del input_epw_res['nc']
        except:
            print("No nc parameters !")
        paras_resistivity['epw_input'] = input_epw_res
        epwnkfq = copy.deepcopy(inpp.epw_nkq) 
        paras_resistivity['epw_nkq'] = epwnkfq
        return paras_resistivity
    ##################### superconductivity ##################################
    def superconductivity_iso(self):    # superconductivity calculation.
        paras_superconductivity = {}
        cal_use_superconductivity = set(['elph', 'epwwrite', 'ephwrite', 'wannierize', 'ep_coupling', 'ephwrite', 'eliashberg', 'wannierize'])
        calculation_superconductivity = list(cal_use_superconductivity)
        input_epw_sup = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            input_epw_sup[key_input] = val_input
        for i in range(len(calculation_superconductivity)):
            input_epw_sup[calculation_superconductivity[i]] = True
        try:
            input_epw_sup = copy.deepcopy(self.cal_control(cal_use_superconductivity, input_epw_sup))
        except:
            print("Go on !")
        input_epw_sup['iverbosity'] = 2
        input_epw_sup['epwread'] = False 
        input_epw_sup['mp_mesh_k'] = True
        input_epw_sup['etf_mem'] = 1
        # SC
        input_epw_sup['eliashberg'] = True
        input_epw_sup['ephwrite'] = True
        input_epw_sup['fermi_plot'] = True

        input_epw_sup['liso'] = True 
        # input_epw_sup['liso'] = True 
        input_epw_sup['limag'] = True 
        input_epw_sup['lpade'] = True 
        input_epw_sup['lacon'] = True 
        input_epw_sup['nsiter'] = 500 
        input_epw_sup['npade'] = 40 
        input_epw_sup['conv_thr_iaxis'] = '1.0d-3'
        input_epw_sup['conv_thr_racon'] = '1.0d-3'
        input_epw_sup['wscut'] = 0.1
        input_epw_sup['muc'] = 0.1
        if 'phonselfen' in input_epw_sup:
            del input_epw_sup['phonselfen']
        if 'elecselfen' in input_epw_sup:
            del input_epw_sup['elecselfen']
        if 'etf_mem' in input_epw_sup:
            del input_epw_sup['etf_mem']
        # input_epw_sup['lacon'] = True 
        # input_epw['delta_approx'] = True 
        paras_superconductivity['epw_input'] = input_epw_sup
        epwnkfq = copy.deepcopy(inpp.epw_nkq) 
        paras_superconductivity['epw_nkq'] = epwnkfq
        return paras_superconductivity
    def superconductivity_TC_iso(self):    # Tc calculation.
        paras_superconductivity = {}
        cal_use_superconductivity = set(['elph', 'epwwrite', 'ephwrite', 'wannierize', 'ep_coupling', 'ephwrite', 'eliashberg', 'wannierize'])
        calculation_superconductivity = list(cal_use_superconductivity)
        input_epw_sup = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            input_epw_sup[key_input] = val_input
        for i in range(len(calculation_superconductivity)):
            input_epw_sup[calculation_superconductivity[i]] = True
        try:
            input_epw_sup = copy.deepcopy(self.cal_control(cal_use_superconductivity, input_epw_sup))
        except:
            print("Go on !")
        input_epw_sup['iverbosity'] = 2
        input_epw_sup['epwread'] = False 
        input_epw_sup['mp_mesh_k'] = True
        input_epw_sup['etf_mem'] = 1
        # SC
        input_epw_sup['eliashberg'] = True
        input_epw_sup['ephwrite'] = True
        input_epw_sup['fermi_plot'] = True

        input_epw_sup['liso'] = True 
        # input_epw_sup['liso'] = True 
        input_epw_sup['limag'] = True 
        input_epw_sup['lpade'] = True 
        input_epw_sup['lacon'] = True 
        input_epw_sup['nsiter'] = 500 
        input_epw_sup['npade'] = 40 
        input_epw_sup['conv_thr_iaxis'] = '1.0d-3'
        input_epw_sup['conv_thr_racon'] = '1.0d-3'
        input_epw_sup['wscut'] = 0.1
        input_epw_sup['muc'] = 0.1
        if 'phonselfen' in input_epw_sup:
            del input_epw_sup['phonselfen']
        if 'elecselfen' in input_epw_sup:
            del input_epw_sup['elecselfen']
        if 'etf_mem' in input_epw_sup:
            del input_epw_sup['etf_mem']
        input_epw_sup['epwread'] = True
        input_epw_sup['fila2f'] = 'pwscf.a2f_iso'
        input_epw_sup['liso'] = True 
        input_epw_sup['limag'] = True 
        input_epw_sup['lpade'] = False 
        input_epw_sup['lacon'] = False 
        input_epw_sup['tc_linear'] = True 
        input_epw_sup['tc_linear_solver'] = 'power' 
        paras_superconductivity['epw_input'] = input_epw_sup
        epwnkfq = copy.deepcopy(inpp.epw_nkq) 
        paras_superconductivity['epw_nkq'] = epwnkfq
        return paras_superconductivity

    def superconductivity_aniso(self):    # superconductivity calculation.
        paras_superconductivity = {}
        cal_use_superconductivity = set(['elph', 'epwwrite', 'ephwrite', 'wannierize', 'ep_coupling', 'ephwrite', 'eliashberg', 'wannierize'])
        calculation_superconductivity = list(cal_use_superconductivity)
        input_epw_sup = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            input_epw_sup[key_input] = val_input
        for i in range(len(calculation_superconductivity)):
            input_epw_sup[calculation_superconductivity[i]] = True
        try:
            input_epw_sup = copy.deepcopy(self.cal_control(cal_use_superconductivity, input_epw_sup))
        except:
            print("Go on !")
        input_epw_sup['iverbosity'] = 2
        input_epw_sup['epwread'] = False 
        input_epw_sup['mp_mesh_k'] = True
        input_epw_sup['etf_mem'] = 1
        # input_epw_sup['ivervosity'] = 2
        # SC
        input_epw_sup['eliashberg'] = True
        input_epw_sup['ephwrite'] = True
        input_epw_sup['fermi_plot'] = True

        input_epw_sup['laniso'] = True 
        # input_epw_sup['liso'] = True 
        input_epw_sup['limag'] = True 
        input_epw_sup['lpade'] = True 
        # input_epw_sup['lacon'] = True 
        input_epw_sup['nsiter'] = 500 
        # input_epw_sup['npade'] = 40 
        input_epw_sup['conv_thr_iaxis'] = '1.0d-3'
        # input_epw_sup['conv_thr_racon'] = '1.0d-3'
        input_epw_sup['wscut'] = 0.1
        input_epw_sup['muc'] = 0.1
        if 'phonselfen' in input_epw_sup:
            del input_epw_sup['phonselfen']
        elif 'elecselfen' in input_epw_sup:
            del input_epw_sup['elecselfen']
        # input_epw_sup['lacon'] = True 
        # input_epw['delta_approx'] = True 
        paras_superconductivity['epw_input'] = input_epw_sup
        epwnkfq = copy.deepcopy(inpp.epw_nkq) 
        paras_superconductivity['epw_nkq'] = epwnkfq
        return paras_superconductivity

    def superconductivity_TC_aniso(self):    # Tc calculation.
        paras_superconductivity = {}
        cal_use_superconductivity = set(['elph', 'epwwrite', 'ephwrite', 'wannierize', 'ep_coupling', 'ephwrite', 'eliashberg', 'wannierize'])
        calculation_superconductivity = list(cal_use_superconductivity)
        input_epw_sup = {}
        epwio = copy.deepcopy(inpp.epw_input)
        for key_input, val_input in epwio.items():
            input_epw_sup[key_input] = val_input
        for i in range(len(calculation_superconductivity)):
            input_epw_sup[calculation_superconductivity[i]] = True
        try:
            input_epw_sup = copy.deepcopy(self.cal_control(cal_use_superconductivity, input_epw_sup))
        except:
            print("Go on !")
        input_epw_sup['epwread'] = False 
        input_epw_sup['mp_mesh_k'] = True
        input_epw_sup['etf_mem'] = 1
        # SC
        input_epw_sup['eliashberg'] = True
        input_epw_sup['ephwrite'] = True
        input_epw_sup['fermi_plot'] = True

        input_epw_sup['laniso'] = True 
        # input_epw_sup['liso'] = True 
        input_epw_sup['limag'] = True 
        input_epw_sup['lpade'] = True 
        # input_epw_sup['lacon'] = True 
        input_epw_sup['nsiter'] = 500 
        # input_epw_sup['npade'] = 40 
        input_epw_sup['conv_thr_iaxis'] = '1.0d-3'
        # input_epw_sup['conv_thr_racon'] = '1.0d-3'
        input_epw_sup['wscut'] = 0.1
        input_epw_sup['muc'] = 0.1
        if 'phonselfen' in input_epw_sup:
            del input_epw_sup['phonselfen']
        elif 'elecselfen' in input_epw_sup:
            del input_epw_sup['elecselfen']
        input_epw_sup['ep_coupling'] = True
        input_epw_sup['elph'] = True
        input_epw_sup['epwwrite'] = False 
        input_epw_sup['epwread'] = True 
        input_epw_sup['fermi_plot'] = False 
        input_epw_sup['wannierize'] = False 
        input_epw_sup['ephwrite'] = False 
        input_epw_sup['imag_read'] = True
        paras_superconductivity['epw_input'] = input_epw_sup
        epwnkfq = copy.deepcopy(inpp.epw_nkq) 
        paras_superconductivity['epw_nkq'] = epwnkfq
        return paras_superconductivity

def EPW_welcome():
    print("====================================================================")
    print("============= EPW calculation file preparation =====================")
    print("      E1) PP-collect                 E2) check bands                ")
    print("      E3) Test Gamma |g|                                            ")
    print("========================== modules =================================")
    print("      P1) phonon linewidth           P2) phon spectral function     ")
    print("      P3) electron linewidth         P4) ele spectral function      ")
    print("      P5) resistivity (Ziman formula)                               ")
    print("      P6) resistivity (BTE)                                         ")
    print("      P7) superconductivity (iso)     P8) Tc iso                    ")
    print("      P9) superconductivity (aniso)  P10) TC aniso                  ")
    print("      U ) User defined                                              ")
    print("====================================================================")
    print("====================== Post processing =============================")
    print("      P ) EPW-Post Processing modules                               ")
    print("------------------->>")
    select = input("Input : ")
    return select

def manipulate():
    sel_EPW = EPW_welcome()
    back = True 
    nks = inpp.phon_points
    while back:
        if sel_EPW == 'E1':
            import scripts.phCollect as collect
            collect.manipulate()
            back = False
            sys.dont_write_bytecode = True
        elif sel_EPW == 'E2':
            epi = epwinput() 
            paras = epi.check_bands()
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), 'epw_band.in')
            import scripts.kgen_epw as kgenepw
            kgenepw.output_kpoints(nks)
            back = False
        elif sel_EPW == 'E3':
            epi = epwinput() 
            paras = epi.Test_gavg()
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), 'epw_test.in')
            import scripts.kgen_epw as kgenepw
            kgenepw.output_kpoints(nks)
            back = False
        elif sel_EPW == 'P1':
            epi = epwinput() 
            paras = epi.epw_linewidth_phonon()
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), 'epw_phonon.in')
            import scripts.kgen_epw as kgenepw
            kgenepw.output_kpoints(nks)
            back = False
        elif sel_EPW == 'P2':
            epi = epwinput() 
            paras = epi.epw_phonon_spectraFunc()
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), 'epw_phon_spectra.in')
            import scripts.kgen_epw as kgenepw
            kgenepw.output_kpoints(nks)
            back = False
        elif sel_EPW == 'P3':
            epi = epwinput() 
            paras = epi.epw_linewidth_electron()
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), 'epw_electron.in')
            import scripts.kgen_epw as kgenepw
            kgenepw.output_kpoints(nks)
            back = False
        elif sel_EPW == 'P4':
            epi = epwinput() 
            paras = epi.epw_electron_spectraFun()
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), 'epw_elec_spectra.in')
            import scripts.kgen_epw as kgenepw
            kgenepw.output_kpoints(nks)
            back = False
        elif sel_EPW == 'P5':
            epi = epwinput() 
            paras = epi.resistivity_ziman()
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), 'epw_res_Ziman.in')
            back = False
            sys.dont_write_bytecode = True
        elif sel_EPW == 'P6':
            epi = epwinput() 
            paras = epi.resistivity_BTE()
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), 'epw_res_BTE.in')
            back = False
            sys.dont_write_bytecode = True
        elif sel_EPW == 'A1':
            from configfile.jobber import job_submit
            epi = epwinput() 
            paras = epi.epw_Auto_cal()
            filename = []
            for keys, vals in paras.items():
                modd.moduleslib(os.getcwd()).config_Epw_single(vals, os.getcwd(), keys)
                filename.append(keys)
            job_submit(os.getcwd()).job_AutoEpw(os.getcwd(), filename)
            import scripts.kgen_epw as kgenepw
            kgenepw.output_kpoints(nks)
            back = False
            sys.dont_write_bytecode = True
        elif sel_EPW == 'U':
            epi = epwinput() 
            paras = epi.epw_cal_user_defined()
            userfname = input("User defined the filename: ")
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), userfname)
            import scripts.kgen_epw as kgenepw
            kgenepw.output_kpoints(nks)
            back = False
            sys.dont_write_bytecode = True
        elif sel_EPW == 'P7':
            epi = epwinput() 
            paras = epi.superconductivity_iso()
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), 'epw_sc_iso.in')
            back = False
            sys.dont_write_bytecode = True
        elif sel_EPW == 'P8':
            import module.funcs as bas 
            epi = epwinput() 
            paras = epi.superconductivity_TC_iso()
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), 'epw_TC_iso.in')
            back = False
            sys.dont_write_bytecode = True

        elif sel_EPW == 'P9':
            epi = epwinput() 
            paras = epi.superconductivity_aniso()
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), 'epw_sc_aniso.in')
            back = False
            sys.dont_write_bytecode = True

        elif sel_EPW == 'P10':
            import module.funcs as bas 
            epi = epwinput() 
            paras = epi.superconductivity_TC_aniso()
            modd.moduleslib(os.getcwd()).config_Epw(paras, os.getcwd(), 'epw_TC_aniso.in')
            back = False
            sys.dont_write_bytecode = True

        elif sel_EPW == 'P':
            import scripts.post_epw as pe
            pe.manipulate()
            back = False
            sys.dont_write_bytecode = True
        else:
            print("Please input right parameters !")
def manipulate_AutoNscfPhon():
    ori_path = os.getcwd()
    if not osp.exists(osp.join(os.getcwd(), 'epw_nscf')):
        os.makedirs(osp.join(os.getcwd(), 'epw_nscf'))
        os.chdir(osp.join(os.getcwd(), 'epw_nscf'))
        import module.mod as modd
        os.system('ln -s ../AutoinputQE.py .')
        os.system('ln -s ../POSCAR .')
        modd.moduleslib(os.getcwd()).config_AutoNscf_EPW(os.getcwd())
        os.system(parameters.submit_pbs + "  submit_QE.sh")
        os.chdir(ori_path)
    os.chdir(ori_path)

    if not osp.exists(osp.join(os.getcwd(), 'phonon')):
        os.makedirs(osp.join(os.getcwd(), 'phonon'))
        os.chdir(osp.join(os.getcwd(), 'phonon'))
        import module.mod as modd1
        os.system('ln -s ../AutoinputQE.py .')
        os.system('ln -s ../POSCAR .')
        modd1.moduleslib(os.getcwd()).config_AutoPhon(os.getcwd())
        os.system(parameters.submit_pbs + "  submit_QE.sh")
        os.chdir(ori_path) 
    os.chdir(ori_path)
def manipulate_AutoBandNscfPhon():
    ori_path = os.getcwd()
    if not osp.exists(osp.join(os.getcwd(), 'epw_nscf')):
        os.makedirs(osp.join(os.getcwd(), 'epw_nscf'))
        os.chdir(osp.join(os.getcwd(), 'epw_nscf'))
        import module.mod as modd
        os.system('ln -s ../AutoinputQE.py .')
        os.system('ln -s ../POSCAR .')
        modd.moduleslib(os.getcwd()).config_AutoBandNscf_EPW(os.getcwd())
        os.system(parameters.submit_pbs + "  submit_QE.sh")
        os.chdir(ori_path)
    os.chdir(ori_path)

    if not osp.exists(osp.join(os.getcwd(), 'phonon')):
        os.makedirs(osp.join(os.getcwd(), 'phonon'))
        os.chdir(osp.join(os.getcwd(), 'phonon'))
        import module.mod as modd1
        os.system('ln -s ../AutoinputQE.py .')
        os.system('ln -s ../POSCAR .')
        modd1.moduleslib(os.getcwd()).config_AutoPhon(os.getcwd())
        os.system(parameters.submit_pbs + "  submit_QE.sh")
        os.chdir(ori_path) 
    os.chdir(ori_path)

def manipulate_Auto():
    epwinfo = epwinput()
    epwinfo.check_bands()