# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: band.py
"""
import os, sys
import os.path as osp
import AutoinputQE as inpp
import module.data_output as IO
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
sys.dont_write_bytecode = True       ## not generate __pycache__ files

def welcome():
    print("======================= Drawing Tools =============================")
    print("       D1) band                        D2) Projected band          ")
    print("       D3) dos                         D4) Projected dos           ")
    print("       D5) visual Brillouin Zone       D6) BZ and Lattice vectors  ")
    print("========================== Phonon =================================")
    print("       P1) Phonon band                 P2) Phonon dos              ")
    print("       P3) Phonon band dos             P4) Phonon band compare     ")
    print("------------------->>")
    select = input("Input : ")
    return select

def manipulate():
    select = welcome()
    if select == 'D1':
        import scripts.band as band 
        import module.funcs as bas 
        mag = bas.mag_judge()
        band_k, energy, hsp = band.plot_band(os.getcwd(), [-6, 6])
        if mag == 2:
            IO.output_bandf_data(os.getcwd(), "BAND_OUTPUT_UP.dat", band_k, energy, ['KPTs', 'BAND_index'])
            IO.output_bandf_data(os.getcwd(), "BAND_OUTPUT_DW.dat", band_k, energy, ['KPTs', 'BAND_index'])
        else:
            IO.output_bandf_data(os.getcwd(), "BAND_OUTPUT.dat", band_k, energy, ['KPTs', 'BAND_index'])
        IO.output_band_klines(os.getcwd(), "KLINES.dat", hsp['hspdata'], hsp['klabels'], ['hsp labels', 'hsp'], [-20, 20])
        from shutil import copyfile
        copyfile(osp.join(curPath, "../replot_tools/replot_band.py"), osp.join(os.getcwd(), "band.py"))
    elif select == 'D2':
        import scripts.proband as pb
        eng = [float(item) for item in input("Please input the energy range of y axis:").split()]
        pb.manipulate(eng)
        from shutil import copyfile
        # copyfile(osp.join(curPath, "dos.py"), osp.join(os.getcwd(), "replot_dos.py"))
    elif select == 'D4':
        import scripts.dos as dos
        dos.manipulate()
        from shutil import copyfile
        copyfile(osp.join(curPath, "dos.py"), osp.join(os.getcwd(), "replot_dos.py"))
    elif select == 'D5':
        import visualize_BZ
        visualize_BZ.manipulate()
    elif select == 'D6':
        import scripts.visualize_BZ as vbz
        vbz.manipulate_lattice()
    elif select == 'P1':
        import phonon_post as phpost
        kpts, phonon, hsp = phpost.manipulate_phonon_band()
        IO.output_bandf_data(os.getcwd(), "PHON_BAND.dat", kpts, phonon.T, ['KPTs', 'Branch'])
        import numpy as np
        IO.output_band_klines(os.getcwd(), "QLINES.dat", hsp['hspdata'], hsp['klabels'], ['hsp labels', 'hsp'], [0, np.max(phonon)*1.1])
        from shutil import copyfile
        copyfile(osp.join(curPath, "../replot_tools/replot_phonon.py"), osp.join(os.getcwd(), "phonon.py"))
    elif select == 'P2':
        import phonon_post as phpost
        phpost.manipulate_phonon_dos()
    elif select == 'P3':
        import phonon_post as phpost
        phpost.manipulate_phonon_band_dos()
    elif select == 'P4':
        import phonon_band_compare as pbc
        pbc.manipulate()
