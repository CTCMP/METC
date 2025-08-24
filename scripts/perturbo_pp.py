# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: band.py
"""
import os, sys
import os.path as osp
import numpy as np 
import pandas as pd 
from matplotlib import pyplot as plt 
from matplotlib import font_manager
import module.funcs as bas
import AutoinputQE as inpp
curPath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curPath)
fontpath = osp.join(curPath, "../font/arial.ttf")
# fontpath = "/work/wangr/data/dxy/scripts/program/QE/font/arial.ttf"
# ------------------ font setup ----------------------#
font_properties = font_manager.FontProperties(fname=fontpath)
styles = ['normal', 'italic', 'oblique']
weights = ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
font = {'fontproperties': font_properties,  # 'family': 'Liberation Sans'
		'style': styles[0],
		'color': 'black',
		'weight': weights[1],
		'size': 16, 
        }
sys.dont_write_bytecode = True       ## not generate __pycache__ files

def imsigma_data(filepath, iband, imode, df, select_plot, fermi):
    os.chdir(filepath)
    data = df[df['it'] == select_plot]
    hbar = 6.5821195 # meV
    if not osp.exists(osp.join(filepath, str(select_plot))):
        os.makedirs(osp.join(filepath, str(select_plot)))
    for i in range(iband):
        data_iband = data[data['ibnd']==i+1]
        for j in range(imode):
            data_imode = data_iband[data_iband['imode']==j+1]
            tmp_data = data_imode.copy()
            tmp_data.loc[:, 'Imsigma'] = tmp_data.loc[:, 'Imsigma'] * 20 / (hbar)
            tmp_data.loc[:, 'E'] = tmp_data.loc[:, 'E'] - fermi
            tmp_data.to_csv(osp.join(filepath, str(select_plot), 'it{}_ibnd{}_imode{}.dat'.format(select_plot, i+1, j+1)), 
                              sep='\t', index=False, float_format='%.8f')
            tmp_data.plot(kind='scatter', x='E', y='Imsigma', title='it{}_ibnd{}_imode{}.dat'.format(select_plot, i+1, j+1))
            plt.yscale('log')
            # plt.ylim(10**-6, 10**2)
            plt.xlabel("Energy (eV)", fontdict=font)
            plt.ylabel("Scattering rates (ps$^{-1}$)", fontdict=font)
            plt.savefig(osp.join(filepath, str(select_plot), 'it{}_ibnd{}_imode{}.png'.format(select_plot, i+1, j+1)), dpi=150)
            plt.close()
    os.chdir(filepath)

def imsigma_pp(filepath):
    imsigma = open(osp.join(filepath, 'pwscf.imsigma_mode'), 'r').readlines()
    if not osp.exists(osp.join(filepath, 'pwscf.imsigma_mode_tmp')):
        imsigma_tmp = open(osp.join(filepath, 'pwscf.imsigma_mode_tmp'), 'w')
        for i in range(len(imsigma)):
            if '#' not in imsigma[i]:
                imsigma_tmp.write(imsigma[i])
        imsigma_tmp.close()

    columns = ['it', 'ik', 'ibnd', 'E', 'imode', 'Imsigma']
    df = pd.read_csv(osp.join(filepath, 'pwscf.imsigma_mode_tmp'), sep='\s+', names=columns)
    #############
    it = df['it'].max()
    iband = df['ibnd'].max()
    imode = df['imode'].max()
    sel0 = input("which one temper you want to get (from 1-{}): ".format(it))
    fermi = float(input("Fermi Level: "))
    if '-' not in sel0:
        select_plot = int(sel0)
        imsigma_data(filepath, iband, imode, df, select_plot, fermi)
    else:
        sel = [int(item) for item in sel0.split('-')]
        select_plot = list(range(sel[0], sel[1]+1))
        print(select_plot)
        for sel_plot in select_plot:
            imsigma_data(filepath, iband, imode, df, sel_plot, fermi)

    

    


