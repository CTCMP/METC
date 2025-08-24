# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: input.py
"""
import os, sys
sys.path.append(os.getcwd())
sys.dont_write_bytecode = True      

def wannier90_welcome():
    print("===================================================================")
    print("========================= Wannier90 ===============================")
    print("       W1)  Wannier90-input files                                  ")
    print("==================== W90-Post Processing ==========================")
    print("       P1) plot w90-pbe               P2) See occupied band        ")
    print("------------------->>")
    select = input("Input : ")
    return select

def manipulate():
    sel = wannier90_welcome()
    back = True
    while back:
        if sel == 'W1':
            from module.mod import moduleslib
            moduleslib(os.getcwd()).config_w90()
            back = False
        elif sel == 'W2':
            from scripts.nihe import manipulate_nihe
            manipulate_nihe()
            back = False
        elif sel == 'P1':
            # energy = [-20, 8]
            energy = [float(item) for item in input("Please input energy range: ").split()]
            from scripts.post_wannier90 import Pw90
            pw90 = Pw90()
            pw90.plot_nihe(os.getcwd(), "band.dat", energy, 'MLWFs')
            back = False
        elif sel == 'P2':
            # energy = [-20, 8]
            energy = [float(item) for item in input("Please input energy range: ").split()]
            from scripts.post_wannier90 import Pw90
            pw90 = Pw90()
            pw90.plot_occupiedBand(os.getcwd(), "band.dat", energy, 'MLWFs')
            back = False
        else:
            print("Please input right parameters !")
if __name__ == '__main__':
    manipulate()

