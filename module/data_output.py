# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: input.py
"""
import os, sys
import os.path as osp
import numpy as np 
# import configfile.inputpara as parameters
# import module.funcs as bas
# import copy
# import AutoinputQE as inpp
sys.path.append(os.getcwd())
sys.dont_write_bytecode = True       ## not generate __pycache__ files

def output_band_klines(filepath, filename, kdata, klabel, labels, ylabel=[-50, 50]):
    output_klines = open(osp.join(filepath, filename), 'w')
    output_klines.write(labels[0]+': ')
    for k in range(len(klabel)-1):
        output_klines.write(klabel[k] + '-')
    output_klines.write(klabel[-1] + '\n')
    output_klines.write(labels[1] + ': ')
    for i in range(len(kdata)):
        output_klines.write(str(np.round(kdata[i], 5)).ljust(6, ' ')+ '  ')
    output_klines.write('\n')
    for l in range(len(kdata)-1):
        output_klines.write(str(np.round(kdata[l], 8)).rjust(10, ' ')+ '  ')
        output_klines.write(str(ylabel[1]) + '\n')
        output_klines.write(str(np.round(kdata[l], 8)).rjust(10, ' ')+ '  ')
        output_klines.write(str(ylabel[0]) + '\n')
        output_klines.write(str(np.round(kdata[l+1], 8)).rjust(10, ' ')+ '  ')
        output_klines.write(str(ylabel[0]) + '\n')
        output_klines.write(str(np.round(kdata[l+1], 8)).rjust(10, ' ')+ '  ')
        output_klines.write(str(ylabel[1]) + '\n')
    output_klines.close()

def output_bandf_data(filepath, filename, kpoints, data, labels):
    output_band = open(osp.join(filepath, filename), 'w')
    output_band.write(str(labels[0]).rjust(10, ' ') + "    ")
    for bd_index in range(list(data.shape)[0]):
        output_band.write(('Index: ' + str(bd_index+1)).rjust(10, ' ') + '    ')
    output_band.write('\n')
    for k_index in range(len(kpoints)):
        output_band.write(str(np.round(kpoints[k_index], 8)).rjust(10, ' ') + '    ')
        for bandindex in range(list(data.shape)[0]):
            output_band.write(str(np.round(data[bandindex, k_index], 8)).rjust(12, ' ') + '    ')
        output_band.write('\n')
    output_band.close()

def output_linewidth(filepath, filename, x, y):
    ff = open(osp.join(filepath, filename), 'w')
    ff.write("kpt".rjust(8, ' '))
    for j in range(len(y)):
        ff.write("   " + str("branches " + str(j+1)).rjust(10, ' '))
    ff.write('\n')
    for i in range(len(x)):
        ff.write(str(np.round(x[i], 6)).rjust(8, ' ') + '   ')
        for k in range(len(y)):
            ff.write(str(np.round(y[k][i], 8)).rjust(10, ' ') + '   ')
        ff.write('\n')
    ff.close()


def output_epw():
    pass