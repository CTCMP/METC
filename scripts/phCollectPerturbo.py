# -*- coding: utf-8 -*-
"""
Created on 19:54 10-11-2022 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: phCollectPerturbo.py
"""

import os
import os.path as osp 

class ppcolP:
    def __init__(self, prefix="pwscf", outdir='tmp'):
        self.filepath = os.getcwd()
        self.prefix = prefix
        self.outdir = outdir

    def manipulateFile(self):
        path_save = osp.join(self.filepath, 'pert_save')

        print("Creating a save dir...")
        os.system('mkdir -p ' + path_save + '/' + self.prefix + '.phsave')
        PHO_DIR = osp.join(self.filepath, 'tmp/_ph0')
        if not osp.exists(osp.join(os.getcwd(), 'pwscf.fc.xml')):
            os.system('cp pwscf.fc ' + path_save + '/ifc.q2r.xml')
        else:
            os.system('cp pwscf.fc.xml ' + path_save + '/ifc.q2r.xml') 

        print("Copying prefix.phsave...")
        os.system('cp ' + PHO_DIR + '/' + self.prefix + '.phsave/* ' + path_save + '/' + self.prefix + '.phsave')

        print("Copying dyn files...")
        os.system('cp ' + osp.join(os.getcwd(), self.prefix + '.dyn* ' + path_save))

        print("Copying the dvscf file for the first q-point...")
        os.system('cp ' + PHO_DIR + '/' + self.prefix + '.dvscf1 ' + path_save + '/' + self.prefix + '.dvscf_q1')

        print("Copy the dvscf for q-points > 1...")
        output = [dI for dI in os.listdir(PHO_DIR) if os.path.isdir(os.path.join(PHO_DIR,dI))]
        dirs = []
        for i in range(len(output)):
            if self.prefix + '.q' in output[i]:
                dirs.append(output[i])
        for q_folder in dirs:
            NQ = int(q_folder.split('_')[-1])
            os.system('cp ' + PHO_DIR + '/' + self.prefix + '.q_' + str(NQ) + '/' + self.prefix + '.dvscf1 '
                    + path_save + '/' + self.prefix + '.dvscf_q' + str(NQ))

    def manipulateFile_notxml(self):
        path_save = osp.join(self.filepath, 'pert_save')

        print("Creating a save dir...")
        os.system('mkdir -p ' + path_save + '/' + self.prefix + '.phsave')
        PHO_DIR = osp.join(self.filepath, 'tmp/_ph0')
        if not osp.exists(osp.join(os.getcwd(), 'pwscf.fc.xml')):
            os.system('cp pwscf.fc ' + path_save + '/ifc.q2r.xml')
        else:
            os.system('cp pwscf.fc.xml ' + path_save + '/ifc.q2r.xml') 

        print("Copying prefix.phsave...")
        os.system('cp ' + PHO_DIR + '/' + self.prefix + '.phsave/* ' + path_save + '/' + self.prefix + '.phsave')

        ff = os.listdir()
        fxml = False
        for x in range(len(ff)):
            if '.xml' in ff[x]:
                fxml = True
                break
            else:
                continue 
        print("fxml: ", fxml)
        if fxml:  
            print("Copying dyn files...")
            os.system('cp ' + osp.join(os.getcwd(), self.prefix + '.dyn* ' + path_save))
        else:
            print("Copying dyn files...")
            os.system('cp ' + osp.join(os.getcwd(), self.prefix + '.dyn* ' + path_save))
            os.system('for ffile in ' + self.prefix + '.dyn*' + '; do mv ' + osp.join(os.getcwd(), path_save, '$ffile ') 
                      + osp.join(os.getcwd(), path_save, '${ffile}.xml') + '; done')

        print("Copying the dvscf file for the first q-point...")
        os.system('cp ' + PHO_DIR + '/' + self.prefix + '.dvscf1 ' + path_save + '/' + self.prefix + '.dvscf_q1')

        print("Copy the dvscf for q-points > 1...")
        output = [dI for dI in os.listdir(PHO_DIR) if os.path.isdir(os.path.join(PHO_DIR,dI))]
        dirs = []
        for i in range(len(output)):
            if self.prefix + '.q' in output[i]:
                dirs.append(output[i])
        for q_folder in dirs:
            NQ = int(q_folder.split('_')[-1])
            os.system('cp ' + PHO_DIR + '/' + self.prefix + '.q_' + str(NQ) + '/' + self.prefix + '.dvscf1 '
                    + path_save + '/' + self.prefix + '.dvscf_q' + str(NQ))

def manipulate():
    prefix="pwscf"
    outdir='tmp'
    collect = ppcolP(prefix, outdir)
    collect.manipulateFile_notxml()

if __name__ == '__main__':
    manipulate()