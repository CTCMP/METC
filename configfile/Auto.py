# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: input.py
"""
import os, sys
sys.dont_write_bytecode = True        ## not generate __pycache__ files

def autocal():
    print("===================================================================")
    print("======================= AutoQE input ==============================")
    print("       A1)  scf-bands                 A2)  Nscf (scf-nscf)         ")
    print("       A3)  Nscf (homogeneous kpt)    A4)  dos (scf-nscf-dos)      ")
    print("======================= Auto-phonon ===============================")
    print("       AP1) Scf-phonon-band-dos       AP2) phonon (perturbo)       ")
    print("       AP3) scf-phonon (g-Gamma)      AP4) phonon (finite-diff)    ")
    print("====================== Auto Wannier90 =============================")
    print("       AW1) Auto W90                  AW2) Auto W90  (only w90)    ")
    print("       AW3) Distangle                                              ")
    print("========================= Auto EPW ================================")
    print("       AE1) Auto Bands Nscf           AE2) Auto band nscf: phonon  ")
    print("       AE3) Auto Nscf                 AE4) Auto nscf: phonon       ")
    # print("       AE5) Auto bands Nscf SC                                     ")
    print("                                                                   ")
    print("------------------->>")
    select = input("Input : ")
    return select

def manipulate():
    sel = autocal()
    back = True
    while back:
        if sel == 'A1':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoBands()
            back = False
        elif sel == 'A2':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoNscf()
            back = False
        elif sel == 'A3':
            import module.mod as modd
            modd.moduleslib(os.getcwd()).config_AutoNscf_EPW(os.getcwd())
            back = False  
        elif sel == 'A4':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoDos()
            back = False
        elif sel == 'AW1':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoW90()
            back = False
        elif sel == 'AW2':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_w90()
            back = False
        elif sel == 'AW3':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_distangle()
            back = False
        elif sel == 'AP1':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoPhon()
            back = False  
        elif sel == 'AP2':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoPhon_Pertubo()
            back = False   
        elif sel == 'AP3':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoPhon_TEST()
            back = False   
        elif sel == 'AP4':
            from module.mod import moduleslib
            import configfile.jobber as job
            moduleslib(os.getcwd()).config_AutoPhon_phonopy(os.getcwd())
            job.phonopy_2nd()
            back = False 
        elif sel == 'AE1':
            import module.mod as modd
            modd.moduleslib(os.getcwd()).config_AutoBandNscf_EPW(os.getcwd())
            back = False  
        elif sel == 'AE2':
            import configfile.epw_input as epwinp
            epwinp.manipulate_AutoBandNscfPhon()
            back = False  
        elif sel == 'AE3':
            import module.mod as modd
            modd.moduleslib(os.getcwd()).config_AutoNscf_EPW(os.getcwd())
            back = False  
        elif sel == 'AE4':
            import configfile.epw_input as epwinp
            epwinp.manipulate_AutoNscfPhon()
            back = False  
        # elif sel == 'AE5':
        #     import module.mod as modd
        #     modd.moduleslib(os.getcwd()).config_AutoBandNscf_EPW_SC(os.getcwd())
        #     back = False  
        # elif sel == 'AE6':
        #     from configfile.epw_input import manipulate_Auto
        #     manipulate_Auto()
        #     back = False  
        else:
            print("Please input right parameters !")

def manipulate_direct(sel):
    back = True
    while back:
        if sel == 'A1':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoBands()
            back = False
        elif sel == 'A2':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoNscf()
            back = False
        elif sel == 'A3':
            import module.mod as modd
            modd.moduleslib(os.getcwd()).config_AutoNscf_EPW(os.getcwd())
            back = False  
        elif sel == 'A4':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoDos()
            back = False
        elif sel == 'AW1':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoW90()
            back = False
        elif sel == 'AW2':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_w90()
            back = False
        elif sel == 'AW3':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_distangle()
            back = False
        elif sel == 'AP1':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoPhon()
            back = False  
        elif sel == 'AP2':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_AutoPhon_Pertubo()
            back = False  
        elif sel == 'AE1':
            import module.mod as modd
            modd.moduleslib(os.getcwd()).config_AutoBandNscf_EPW(os.getcwd())
            back = False  
        elif sel == 'AE2':
            import configfile.epw_input as epwinp
            epwinp.manipulate_AutoBandNscfPhon()
            back = False  
        elif sel == 'AE3':
            import module.mod as modd
            modd.moduleslib(os.getcwd()).config_AutoNscf_EPW(os.getcwd())
            back = False  
        elif sel == 'AE4':
            import configfile.epw_input as epwinp
            epwinp.manipulate_AutoNscfPhon()
            back = False  
        elif sel == 'AE5':
            from configfile.epw_input import manipulate_Auto
            manipulate_Auto()
            back = False  
        else:
            print("Please input right parameters !")
