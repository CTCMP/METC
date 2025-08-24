# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: input.py
"""
import os.path as osp 
import os, sys
import module.funcs as basFucn
import configfile.inputpara as parameters
import AutoinputQE as inpp
import module.mod as modd
sys.path.append(os.getcwd())
import copy
sys.dont_write_bytecode = True        ## not generate __pycache__ files

ori_path = os.getcwd()

def check_file(fpath):
    if not osp.exists(osp.join(fpath, 'nscf')):
        print("PATH: %s does not exists, exit !" %(osp.join(fpath, 'nscf')))
        exit()
    if not osp.exists(osp.join(osp.join(fpath, 'w90'))):
        print("PATH: %s does not exists, exit !" %(osp.join(fpath, 'w90')))
        exit()
    if not osp.exists(osp.join(fpath, 'phonon')):
        print("PATH: %s does not exists, exit !" %(osp.join(fpath, 'phonon')))
        exit()

class perturboinput():
    def __init__(self, filepath=os.getcwd()):
        self.filepath = filepath 
        self.nks = 30
        self.nks2 = 20
        try:
            self.fermi = basFucn.get_fermi_level(osp.join(filepath, '../../nscf'), 'scf.log')
        except:
            self.fermi = basFucn.get_fermi_level(osp.join(filepath, '../nscf'), 'scf.log')
        self.kb = 8.6173324*10**(-5)  # ev/k
        self.boltzMin = self.fermi - 8*self.kb*300
        self.boltzMax = self.fermi + 8*self.kb*300
        print('Recomand Boltz_emin/ax:(%s, %s)' %(str(self.boltzMin), str(self.boltzMax)))

    def qe2pertinput(self):
        pars = {}
        qe2pertin = copy.deepcopy(parameters.pert_qe2pert)
        paras = copy.deepcopy(inpp.qe2pert_in)  
        for keys, vals in paras.items():
            qe2pertin[keys] = vals
        pars['qe2pert'] = qe2pertin
        return pars
    
    def get_pwscf_setup(self):
        if not osp.exists(osp.join(os.getcwd(), 'pwscf.temper')):
            if osp.exists(osp.join(os.getcwd(), '../setup/pwscf.temper')):
                os.system('ln -s ' + osp.join(os.getcwd(), '../setup/pwscf.temper') + ' .')
                os.system('ln -s ' + osp.join(os.getcwd(), '../setup/pwscf_tet.kpt') + ' .')
            else:
                try:
                    fermi_define = inpp.energy_fermi
                except:
                    fermi_define = float(input("Input energy for fermi: "))
                basFucn.write_perturbo_temper(self.filepath, 'pwscf.temper', \
                                        inpp.temperature, fermi_define, inpp.n_carriers)
        # try:
        #     fermi_define = inpp.energy_fermi
        # except:
        #     fermi_define = float(input("Input energy for fermi: "))
        # basFucn.write_perturbo_temper(self.filepath, 'pwscf.temper', \
        #                         inpp.temperature, fermi_define, inpp.n_carriers)

    # def get_pwscf_setup_2(self):
    #     if not osp.exists(osp.join(os.getcwd(), 'pwscf.temper')):
    #         if osp.exists(osp.join(os.getcwd(), '../setup/pwscf.temper')):
    #             os.system('ln -s ' + osp.join(os.getcwd(), '../setup/pwscf.temper') + ' .')
    #             basFucn.write_band_k_perturbo(self.filepath, 'bands.qpt', self.nks2)
    #             # os.system('ln -s ' + osp.join(os.getcwd(), '../setup/pwscf_tet.kpt') + ' .')
    #         else:
    #             basFucn.write_perturbo_temper(self.filepath, 'pwscf.temper', \
    #                                 inpp.temperature, self.fermi, inpp.n_carriers)
    #             basFucn.write_band_k_perturbo(self.filepath, 'bands.qpt', self.nks2)

    def perturbo_bands(self):
        pars = {}
        perturboin = copy.deepcopy(parameters.pert_pertubo)
        basFucn.write_band_k_perturbo(self.filepath, 'bands.kpt', self.nks)
        perturboin['fklist'] = 'bands.kpt'
        perturboin['calc_mode'] = 'bands'
        pars['perturbo'] = perturboin
        return pars
    
    def perturbo_phonon(self):
        pars = {}
        perturboin = copy.deepcopy(parameters.pert_pertubo)
        basFucn.write_band_k_perturbo(self.filepath, 'phdisp.qpt', self.nks)
        perturboin['fqlist'] = 'phdisp.qpt'
        perturboin['calc_mode'] = 'phdisp'
        pars['perturbo'] = perturboin
        return pars
    def perturbo_ephmat(self):
        pars = {}
        perturboin = copy.deepcopy(parameters.pert_pertubo)
        perturboin['calc_mode'] = 'ephmat'
        basFucn.write_band_k_perturbo(self.filepath, 'eph.qpt', self.nks)
        basFucn.write_perturbo_kpt_hsp(self.filepath, 'eph.kpt')
        perturboin['fklist'] = 'eph.kpt'
        perturboin['fqlist'] = 'eph.qpt'
        perturboin['band_min'] = inpp.perturbo_in['band_min']
        perturboin['band_max'] = inpp.perturbo_in['band_max']
        perturboin['phfreq_cutoff'] = inpp.perturbo_in['phfreq_cutoff']
        pars['perturbo'] = perturboin
        return pars
    def perturbo_setup(self):
        pars = {}
        perturboin = copy.deepcopy(parameters.pert_pertubo)
        perturboin['calc_mode'] = 'setup'
        try:
            if inpp.perturbo_in['hole']:
                perturboin['hole'] = True
            else:
                perturboin['hole'] = False 
        except:
            perturboin['hole'] = False 
        for keys, vals in inpp.perturbo_in.items():
            if 'boltz_' in keys:
                perturboin[keys] = vals 
        try:
            del perturboin['boltz_nstep']
        except:
            print('go on !')
        perturboin['band_min'] = inpp.perturbo_in['band_min']
        perturboin['band_max'] = inpp.perturbo_in['band_max']
        # basFucn.write_perturbo_temper(self.filepath, 'pwscf.temper', \
        #                     inpp.temperature, self.fermi, inpp.n_carriers)
        self.get_pwscf_setup()
        perturboin['ftemper'] = 'pwscf.temper'
        try:
            perturboin['find_efermi'] = inpp.perturbo_in['find_efermi']
        except:
            perturboin['find_efermi'] = False
        pars['perturbo'] = perturboin
        return pars
    
    def perturbo_imsigma(self):
        pars = {}
        perturboin = copy.deepcopy(parameters.pert_pertubo)
        perturboin['calc_mode'] = 'imsigma'
        perturboin['band_min'] = inpp.perturbo_in['band_min']
        perturboin['band_max'] = inpp.perturbo_in['band_max']
        perturboin['phfreq_cutoff'] = inpp.perturbo_in['phfreq_cutoff']
        perturboin['delta_smear'] = inpp.perturbo_in['delta_smear']
        try:
            perturboin['sampling'] = inpp.perturbo_in['sampling']
            perturboin['nsamples'] = inpp.perturbo_in['nsamples']
            if inpp.perturbo_in['sampling'] == 'cauchy':
                perturboin['cauchy_scale'] = inpp.perturbo_in['cauchy_scale']
        except:
            print("Something Wrong happend with sampling !")
        # basFucn.write_band_k_perturbo(self.filepath, 'pwscf_tet.kpt', self.nks)
        perturboin['fklist'] = 'pwscf_tet.kpt'
        # basFucn.write_perturbo_temper(self.filepath, 'pwscf.temper', \
        #                     inpp.temperature, self.fermi, inpp.n_carriers)
        self.get_pwscf_setup()
        perturboin['ftemper'] = 'pwscf.temper'
        pars['perturbo'] = perturboin
        return pars

    def perturbo_imsigma_hsp(self):
        pars = {}
        perturboin = copy.deepcopy(parameters.pert_pertubo)
        perturboin['calc_mode'] = 'imsigma'
        perturboin['band_min'] = inpp.perturbo_in['band_min']
        perturboin['band_max'] = inpp.perturbo_in['band_max']
        perturboin['phfreq_cutoff'] = inpp.perturbo_in['phfreq_cutoff']
        perturboin['delta_smear'] = inpp.perturbo_in['delta_smear']
        try:
            perturboin['sampling'] = inpp.perturbo_in['sampling']
            perturboin['nsamples'] = inpp.perturbo_in['nsamples']
            if inpp.perturbo_in['sampling'] == 'cauchy':
                perturboin['cauchy_scale'] = inpp.perturbo_in['cauchy_scale']
        except:
            basFucn.write_band_k_perturbo(self.filepath, 'bands.qpt', self.nks2)
            perturboin['fqlist'] = 'bands.qpt'
        # basFucn.write_band_k_perturbo(self.filepath, 'pwscf_tet.kpt', self.nks)
        perturboin['fklist'] = 'pwscf_tet.kpt'
        # basFucn.write_perturbo_temper(self.filepath, 'pwscf.temper', \
        #                     inpp.temperature, self.fermi, inpp.n_carriers)
        self.get_pwscf_setup()
        perturboin['ftemper'] = 'pwscf.temper'
        pars['perturbo'] = perturboin
        return pars

    def perturbo_meanfp(self):
        pars = {}
        perturboin = copy.deepcopy(parameters.pert_pertubo)
        perturboin['calc_mode'] = 'meanfp'
        perturboin['band_min'] = inpp.perturbo_in['band_min']
        perturboin['band_max'] = inpp.perturbo_in['band_max']
        perturboin['phfreq_cutoff'] = inpp.perturbo_in['phfreq_cutoff']
        perturboin['delta_smear'] = inpp.perturbo_in['delta_smear']
        perturboin['sampling'] = inpp.perturbo_in['sampling']
        if inpp.perturbo_in['sampling'] == 'cauchy':
            perturboin['cauchy_scale'] = inpp.perturbo_in['cauchy_scale']
        perturboin['nsamples'] = inpp.perturbo_in['nsamples']
        # basFucn.write_band_k_perturbo(self.filepath, 'eph.kpt', self.nks)
        perturboin['fklist'] = 'pwscf_tet.kpt'
        # basFucn.write_perturbo_temper(self.filepath, 'pwscf.temper', \
        #                     inpp.temperature, self.fermi, inpp.n_carriers)
        os.system('ln -s ' + osp.join(os.getcwd(), '../imsigma/pwscf.imsigma') + ' .')
        self.get_pwscf_setup()
        perturboin['ftemper'] = 'pwscf.temper'
        pars['perturbo'] = perturboin
        return pars
    
    def perturbo_trans_ITA(self):
        pars = {}
        perturboin = copy.deepcopy(parameters.pert_pertubo)
        perturboin['calc_mode'] = 'trans-ita'
        perturboin['band_min'] = inpp.perturbo_in['band_min']
        perturboin['band_max'] = inpp.perturbo_in['band_max']
        perturboin['delta_smear'] = inpp.perturbo_in['delta_smear']
        perturboin['phfreq_cutoff'] = inpp.perturbo_in['phfreq_cutoff']
        perturboin['tmp_dir'] = './tmp'
        try:
            perturboin['hole'] = inpp.perturbo_in['hole']
        except:
            perturboin['hole'] = False 
        try:
            perturboin['load_scatter_eph'] = inpp.perturbo_in['load_scatter_eph']
        except:
            print('Not specify load_scatter_eph parameter, go on !')
        for keys, vals in inpp.perturbo_in.items():
            if 'boltz_' in keys:
                perturboin[keys] = vals 
        # basFucn.write_perturbo_temper(self.filepath, 'pwscf.temper', \
        #                     inpp.temperature, self.fermi, inpp.n_carriers)
        os.system('ln -s ' + osp.join(os.getcwd(), '../setup/pwscf_tet.h5') + ' .')
        # os.system('ln -s ' + osp.join(os.getcwd(), '../setup/pwscf.temper') + ' .')
        # os.system('ln -s ' + osp.join(os.getcwd(), '../imsigma/pwscf.imsigma') + ' .')
        self.get_pwscf_setup()
        perturboin['ftemper'] = 'pwscf.temper'
        pars['perturbo'] = perturboin
        return pars
    
    def perturbo_trans_RTA(self):
        pars = {}
        perturboin = copy.deepcopy(parameters.pert_pertubo)
        perturboin['calc_mode'] = 'trans-rta'
        perturboin['band_min'] = inpp.perturbo_in['band_min']
        perturboin['band_max'] = inpp.perturbo_in['band_max']
        try:
            perturboin['hole'] = inpp.perturbo_in['hole']
        except:
            perturboin['hole'] = False 
        for keys, vals in inpp.perturbo_in.items():
            if 'boltz_' in keys:
                perturboin[keys] = vals 
        # basFucn.write_perturbo_temper(self.filepath, 'pwscf.temper', \
        #                     inpp.temperature, self.fermi, inpp.n_carriers)
        os.system('ln -s ' + osp.join(os.getcwd(), '../setup/pwscf_tet.h5') + ' .')
        os.system('ln -s ' + osp.join(os.getcwd(), '../imsigma/pwscf.imsigma') + ' .')
        self.get_pwscf_setup()
        perturboin['ftemper'] = 'pwscf.temper'
        pars['perturbo'] = perturboin
        return pars
    
class perturboinput_auto():
    def __init__(self, filepath=os.getcwd()):
        self.filepath = filepath 
        self.nks = 51

    def qe2pertinput(self):
        if not osp.exists(osp.join(self.filepath, 'qe2pert')):
            os.makedirs(osp.join(osp.join(self.filepath, 'qe2pert')))
            os.chdir(osp.join(self.filepath, 'qe2pert'))
        pars = {}
        qe2pertin = copy.deepcopy(parameters.pert_qe2pert)
        paras = copy.deepcopy(inpp.qe2pert_in)  
        for keys, vals in paras.items():
            qe2pertin[keys] = vals
        pars['qe2pert'] = qe2pertin
        return pars
    
    def pert_bands(self):
        pars = {}
        perturboin = copy.deepcopy(parameters.pert_qe2pert)
        paras = copy.deepcopy(inpp.perturbo_in)  
        for keys, vals in paras.items():
            perturboin[keys] = vals
        pars['perturbo'] = perturboin
        return pars

def Perturbo_welcome():
    print("===================================================================")
    print("========================= Perturbo ================================")
    print("       P1) pp-collect pert            P2) qe2pert prepare          ")
    print("       P3) bands                      P4) phdisp                   ")
    print("       P5) ephmat                     P6) setup                    ")
    print("       P7) imsigma                    P8) imsigma hsp              ")
    print("       P9) meanfp                                                  ")
    print("       T1) trans-ita                  T2) trans-rta                ")
    print("       T3) trans-mag-ita              T4) trans-mag-rta            ")
    print("===================================================================")
    print("       B1) pwscf_tet (HSP)                                         ")
    print("===================================================================")
    print("       PP) Perturbo-Post Processing Module                         ")
    print("------------------->>")
    select = input("Input : ")
    return select
        
def manipulate():
    sel_pert = Perturbo_welcome()
    # if sel_pert != 'P1':
    #     check_file(osp.join(os.getcwd(), '../'))
    if sel_pert == 'P1':
        import scripts.phCollectPerturbo as phcoll
        phcoll.manipulate()
        sys.dont_write_bytecode = True

    elif sel_pert == 'P2':
        check_file(osp.join(os.getcwd(), '../'))
        pert = perturboinput()
        paras = pert.qe2pertinput()
        modd.moduleslib(os.getcwd()).config_qe2pert(paras, osp.join(os.getcwd()), 'qe2pert.in')
        sys.dont_write_bytecode = True

    elif sel_pert == 'P3':
        pert = perturboinput()
        paras = pert.perturbo_bands()
        modd.moduleslib(os.getcwd()).config_pert(paras, osp.join(os.getcwd()), 'bands.in')
        sys.dont_write_bytecode = True
    
    elif sel_pert == 'P4':
        pert = perturboinput()
        paras = pert.perturbo_phonon()
        modd.moduleslib(os.getcwd()).config_pert(paras, osp.join(os.getcwd()), 'phdisp.in')
        sys.dont_write_bytecode = True

    elif sel_pert == 'P5':
        pert = perturboinput()
        paras = pert.perturbo_ephmat()
        modd.moduleslib(os.getcwd()).config_pert_ephmat(paras, osp.join(os.getcwd()), 'ephmat.in')
        sys.dont_write_bytecode = True

    elif sel_pert == 'P6':
        pert = perturboinput()
        paras = pert.perturbo_setup()
        modd.moduleslib(os.getcwd()).config_pert(paras, osp.join(os.getcwd()), 'setup.in')
        sys.dont_write_bytecode = True

    elif sel_pert == 'P7':
        pert = perturboinput()
        paras = pert.perturbo_imsigma()
        modd.moduleslib(os.getcwd()).config_pert_trans(paras, osp.join(os.getcwd()), 'imsigma.in')
        sys.dont_write_bytecode = True

    elif sel_pert == 'P8':
        pert = perturboinput()
        paras = pert.perturbo_imsigma_hsp()
        modd.moduleslib(os.getcwd()).config_pert_trans(paras, osp.join(os.getcwd()), 'imsigma.in')
        sys.dont_write_bytecode = True

    elif sel_pert == 'P9':
        pert = perturboinput()
        paras = pert.perturbo_meanfp()
        modd.moduleslib(os.getcwd()).config_pert(paras, osp.join(os.getcwd()), 'meanfp.in')
        sys.dont_write_bytecode = True

    elif sel_pert == 'T1':
        pert = perturboinput()
        paras = pert.perturbo_trans_ITA()
        modd.moduleslib(os.getcwd()).config_pert_trans(paras, osp.join(os.getcwd()), 'trans-ita.in')
        sys.dont_write_bytecode = True

    elif sel_pert == 'T2':
        pert = perturboinput()
        paras = pert.perturbo_trans_RTA()
        modd.moduleslib(os.getcwd()).config_pert_trans(paras, osp.join(os.getcwd()), 'trans-rta.in')
        sys.dont_write_bytecode = True

    elif sel_pert == 'B1':
        nks = int(input("Input nks for hsp cal: "))
        basFucn.write_band_k_perturbo(os.getcwd(), 'pwscf_tet.kpt', nks)

    elif sel_pert == 'PP':
        import scripts.post_perturbo as ppert
        ppert.manipulate()
        sys.dont_write_bytecode = True
    else:
        print("Please input right parameters !")
def manipulate_auto():
    sel_pert = Perturbo_welcome()
    check_file(osp.join(os.getcwd(), '../'))
    back = True 
    while back:
        if sel_pert == 'P1':
            import scripts.phCollectPerturbo as phcoll
            phcoll.manipulate()
            back = False
            sys.dont_write_bytecode = True

        elif sel_pert == 'P2':
            pert = perturboinput_auto()
            paras = pert.qe2pertinput()
            modd.moduleslib(os.getcwd()).config_qe2pert(paras, osp.join(os.getcwd(), 'qe2pert'), 'qe2pert.in')
            back = False
            sys.dont_write_bytecode = True
