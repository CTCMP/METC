# -*- coding: utf-8 -*-
"""
Created on 10:21 14-10-2022 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: qe.py
"""
import os, sys
import os.path as osp
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.getcwd())
sys.dont_write_bytecode = True       # not generate __pycache__ files

def welcome():
    print("===================================================================")
    print("============ Pre and post-preparation for QE calculation ==========")
    print("       1) Template files              2)  Get input files          ")
    print("       3) Get structure               4)  Get str from input       ")
    print("       5) Get structure (Not finished optimize)                    ")
    print("====================== Drawing Utilitites =========================")
    print("       D)  Plotting: post processing, band, dos, phonon ...        ")
    print("                                                                   ")
    print("======================= Other Modules =============================")
    print("*******************************************************************")
    print("       M1) Wannier90                  M2) EPW                      ")
    print("       M3) ShengBTE                   M4) Perturbo                 ")
    print("*******************************************************************")
    print("========================== Auto ===================================")
    print("       A)  Auto  bands, phonon, wannier90, EPW, Perturbo ...       ")
    print("------------------->>")
    select = input("Input : ")
    return select

def templatefiles():
    print("======================= Template files ============================")
    print("       T1) input template file                                     ")
    print("------------------->>")
    select = input("Input : ")
    return select

def checkfile_inputQE():
    from shutil import copyfile
    sys.dont_write_bytecode = True
    if not osp.exists(osp.join(os.getcwd(), "AutoinputQE.py")):
        copyfile(osp.join(osp.join(curPath, "configfile"), "Auto_inpp.py"), osp.join(os.getcwd(), "AutoinputQE.py"))
def checkfile_AutoinputQE(filename):
    from shutil import copyfile
    sys.dont_write_bytecode = True
    if not osp.exists(osp.join(os.getcwd(), "AutoinputQE.py")):
        copyfile(osp.join(osp.join(curPath, "configfile"), filename), osp.join(os.getcwd(), "AutoinputQE.py"))
def checkfile_relaxtype():
    import module.structure as struct
    sys.dont_write_bytecode = True
    if osp.exists(osp.join(os.getcwd(), "vc-relax.in")):
        struct.write_contcar(os.getcwd(), 'vc-relax.log')
    elif osp.exists(osp.join(os.getcwd(), "relax.in")):
        struct.write_contcar(os.getcwd(), 'relax.log')
    else:
        print("No output file: vc-relax.log or relax.log !")

def checkfile_relaxtype_NoOPT():
    import module.structure as struct
    sys.dont_write_bytecode = True
    if osp.exists(osp.join(os.getcwd(), "vc-relax.in")):
        struct.write_contcar(os.getcwd(), 'vc-relax.log')
    elif osp.exists(osp.join(os.getcwd(), "relax.in")):
        struct.write_contcar(os.getcwd(), 'relax.log')
    else:
        print("No output file: vc-relax.log or relax.log !")

def manipulate():
    select = welcome()
    if select == '1':
        checkfile_inputQE()
        os.system('chmod +x AutoinputQE.py')
    elif select == '2':
        checkfile_inputQE()
        from module.mod import moduleslib
        moduleslib(os.getcwd()).config_single()
        sys.dont_write_bytecode = True
    elif select == '3':
        print("Begin processing !")
        checkfile_relaxtype()
        print("End processing !")
    elif select == '4':
        from module.structure import get_stru_from_QE_input
        fname = input("Please input filename: ")
        get_stru_from_QE_input(os.getcwd(), fname)
    elif select == '5':
        print("Begin processing !")
        checkfile_relaxtype_NoOPT()
        print("End processing !")
    elif select == 'D':
        import scripts.post as postp
        postp.manipulate()
        sys.dont_write_bytecode = True
    elif select == 'M1':                           ################ Auto
        checkfile_inputQE()
        import configfile.w90_input as w90
        w90.manipulate()
    elif select == 'M2':                           ################ Auto
        checkfile_inputQE()
        import configfile.epw_input as epw
        epw.manipulate()
    elif select == 'M3':                           ################ Auto
        checkfile_inputQE()
        import configfile.shengbte_input as shengbte
        shengbte.manipulate()
    elif select == 'M4':                           ################ Auto
        checkfile_inputQE()
        import configfile.perturbo_input as perturbo
        perturbo.manipulate()
    elif select == 'A':                           ################ Auto
        checkfile_inputQE()
        import configfile.Auto as aut
        aut.manipulate()
    else:
        import configfile.Auto as ap
        # sel = ap.autocal()
        ap.manipulate_direct(select)
        # print("Please input right parameters !")
            
if __name__ == '__main__':
    # try:
    #     import configfile.config as configg
    #     configg.conf_intel()
    #     print("Intel environment config Success !")
    # except:
    #     print("Warning: Intel environment config failure !")
    manipulate()
