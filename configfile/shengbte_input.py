# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: input.py
"""
import os.path as osp 
import os, sys
import configfile.config as submit
sys.path.append(os.getcwd())
import configfile.inputpara as parameters
import AutoinputQE as inpp 

sys.dont_write_bytecode = True        ## not generate __pycache__ files

def three_rd_pre():
    import configfile.jobber as job
    if not osp.exists(osp.join(os.getcwd(), 'scf.in')) or not osp.exists(osp.join(os.getcwd(), 'scf_sc.in')):
        import module.mod as modd
        modd.moduleslib(os.getcwd()).shengbte_3rd(os.getcwd())
    os.system(parameters.thirdorder_path + '/' + inpp.third_order_par)
    
    job.shengbte_3rd()

def fourth_th_pre():
    import configfile.jobber as job
    if not osp.exists(osp.join(os.getcwd(), 'scf.in')) or not osp.exists(osp.join(os.getcwd(), 'scf_sc.in')):
        import module.mod as modd
        modd.moduleslib(os.getcwd()).shengbte_4th(os.getcwd())
    os.system(parameters.python2_path + " " + osp.join(parameters.fourthorder_path, inpp.fourth_order_par))
    job.shengbte_4th()

def three_rd_post():
    try:
        os.system(parameters.thirdorder_path + '/' + inpp.third_order_collect) 
        # os.system(inpp.third_order_collect)
    except:
        os.system(parameters.thirdorder_path + '/' + inpp.third_order_collect)    

def four_th_post():
    os.system(parameters.python2_path + " " + osp.join(parameters.fourthorder_path, inpp.fourth_order_collect)) 

def attention():
    import module.mod as modd 
    modd.moduleslib(os.getcwd()).shengbte(os.getcwd())

def ShengBTE_welcome():
    print("===================================================================")
    print("======================== ShengBTE =================================")
    print("   S1) 3rd_prepare                S2) 3rd_collect                  ") 
    print("   S3) 4th prepare                S4) 4th collect                  ")
    print("   S) configure input file                                         ")
    print("================= PostProcessing of ShengBTE ======================")
    print("   P) Output fig of shengBTE                                       ")
    print("------------------->>")
    select = input("Input : ")
    return select

def manipulate():
    sel_shengbte = ShengBTE_welcome()
    if sel_shengbte == 'S1':
        three_rd_pre()
        sys.dont_write_bytecode = True
    elif sel_shengbte == 'S2':
        three_rd_post()
        sys.dont_write_bytecode = True
    elif sel_shengbte == 'S3':
        fourth_th_pre()
        sys.dont_write_bytecode = True
    elif sel_shengbte == 'S4':
        four_th_post()
        sys.dont_write_bytecode = True
    elif sel_shengbte == 'S':
        attention()
        sys.dont_write_bytecode = True
    elif sel_shengbte == 'P':
        import scripts.post_shengBTE as psb
        sel2 = psb.manipulate()

