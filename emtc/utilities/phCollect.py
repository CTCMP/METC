# -*- coding: utf-8 -*-
"""
Created on 19:54 10-11-2022 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: phCollect.py
"""

import os
import re
from xml.dom import minidom
import os.path as osp 

class ppcol:
    def __init__(self, prefix="pwscf", outdir='tmp'):
        self.filepath = os.getcwd()
        self.prefix = prefix
        self.outdir = outdir
    
    def hasXML(self):
        fname = osp.join(self.filepath, self.prefix + ".dyn1.xml")
        if os.path.isfile(fname):
            return True
        # check if the other without .xml extension exists
        # if not raise an error
        fname_no_xml = re.sub('\.xml$', '', fname)

        class FileNotFoundError(Exception):
            pass
        if not os.path.isfile(fname_no_xml):
            raise FileNotFoundError(
                    "No dyn0 file found cannot tell if xml format was used.")
        return False
    # Check if the calculation includes PAW
    def hasPAW(self):
        fname = osp.join(self.filepath, self.outdir + '/' + self.prefix + '.save/data-file-schema.xml')
        xmldoc = minidom.parse(fname)
        item = xmldoc.getElementsByTagName('paw')[0]
        lPAW = (item.childNodes[0].data == 'true')
        return lPAW

    # Check if the calculation used .fc or .fc.xml files
    def hasfc(self):
        fname = osp.join(self.filepath, self.prefix + '.fc.xml')
        if (os.path.isfile(fname)):
            lfc = True
        else:
            fname_no_xml = re.sub('\.xml$', '', fname)
            if (os.path.isfile(fname_no_xml)):
                lfc = True
            else:
                lfc = False
        return lfc
    # Check if the calculation include SOC
    def hasSOC(self):
        fname = osp.join(self.filepath, self.outdir + '/' + self.prefix+'.save/data-file-schema.xml')
        xmldoc = minidom.parse(fname)
        item = xmldoc.getElementsByTagName('spinorbit')[0]
        lSOC = item.childNodes[0].data
        return lSOC

    # Return the number of q-points in the IBZ
    def get_nqpt(self):
        fname = osp.join(self.filepath, self.outdir + '/' + '_ph0/' + self.prefix + '.phsave/control_ph.xml')
        fid = open(fname, 'r')
        lines = fid.readlines()
        # these files are relatively small so reading the whole thing shouldn't
        # be an issue
        fid.close()

        line_number_of_nqpt = 0
        while 'NUMBER_OF_Q_POINTS' not in lines[line_number_of_nqpt]:
            # increment to line of interest
            line_number_of_nqpt += 1
        line_number_of_nqpt += 1  # its on the next line after that text

        nqpt = int(lines[line_number_of_nqpt])

        return nqpt
    # Check if the calculation was done in sequential
    def isSEQ(self):
        fname = osp.join(self.filepath, self.outdir + '/' + '_ph0/'+str(self.prefix)+'.dvscf')
        if (os.path.isfile(fname)):
            lseq = True
        else:
            lseq = False

        return lseq

def manipulate():
    # user_input = input('Enter the prefix used for PH calculations (e.g. diam)\n')
    # prefix = str(user_input)
    prefix="pwscf"
    outdir='tmp'
    collect = ppcol(prefix, outdir)
    # # Test if SOC
    # SOC = hasSOC(prefix)
    # Test if '.xml' files are used
    XML = collect.hasXML()

    # Test if PAW
    PAW = collect.hasPAW()

    # Test if fc
    fc = collect.hasfc()

    # Test if seq. or parallel run
    SEQ = collect.isSEQ()

    if True:  # this gets the nqpt from the outputfiles
        nqpt = collect.get_nqpt()

    else:
        # Enter the number of irr. q-points
        user_input = input('Enter the number of irreducible q-points\n')
        nqpt = user_input
        try:
            nqpt = int(user_input)
        except ValueError:
            raise Exception('The value you enter is not an integer!')
    ori_file = os.getcwd() + '/save'
    os.system('mkdir ' + ori_file + ' 2>/dev/null')
    path_save = ori_file + '/'
    path_ori = os.getcwd() + '/'
    path_collect = osp.join(os.getcwd(), outdir) + '/_ph0/'
    for iqpt in range(1, nqpt+1):
        label = str(iqpt)

        # Case calculation in seq.
        if SEQ:
            # Case with XML files
            if XML:

                os.system('cp '+ path_ori + prefix +'.dyn0 '+ path_ori + prefix + '.dyn0.xml')
                os.system('cp '+ path_ori + prefix +'.dyn'+ str(iqpt) + '.xml ' + path_save + prefix
                        + '.dyn_q'+label+'.xml')
                if (iqpt == 1):
                    os.system('cp '+ path_collect + prefix+'.dvscf* ' + path_save + prefix + '.dvscf_q'
                            + label)
                    os.system('cp -r '+ path_collect + prefix+'.phsave ' + path_save)
                    if fc:
                        os.system('cp '+ path_ori + prefix +'.fc.xml ' + path_save + 'ifc.q2r.xml')
                    if PAW:
                        os.system('cp ' + path_collect + prefix+'.dvscf_paw* '+ path_save + prefix +
                                '.dvscf_paw_q'+label)
                else:
                    os.system('cp ' + path_collect + prefix+'.q_'+str(iqpt)+'/'+prefix +
                            '.dvscf* '+ path_save + prefix+'.dvscf_q'+label)
                    # os.system('rm ' + path_collect + prefix+'.q_'+str(iqpt)+'/*wfc*')
                    if PAW:
                        os.system('cp ' + path_collect + prefix+'.q_'+str(iqpt)+'/'+prefix +
                                '.dvscf_paw* ' + path_save + prefix+'.dvscf_paw_q'+label)
            # Case without XML files
            else:
                os.system('cp '+ path_ori + prefix+'.dyn'+str(iqpt)+' ' + path_save + prefix+'.dyn_q' +
                        label)
                if (iqpt == 1):
                    os.system('cp ' + path_collect + prefix+'.dvscf ' + path_save + prefix+'.dvscf_q' +
                            label)
                    os.system('cp -r ' + path_collect + prefix+'.phsave ' + path_save)
                    if fc:
                        os.system('cp ' + path_ori + prefix + '.fc ' + path_save + 'ifc.q2r')
                    if PAW:
                        os.system('cp ' + path_collect + prefix+'.dvscf_paw ' + path_save + prefix +
                                '.dvscf_paw_q'+label)
                else:
                    os.system('cp ' + path_collect + prefix+'.q_'+str(iqpt)+'/'+prefix +
                            '.dvscf ' + path_save + prefix+'.dvscf_q'+label)
                    # os.system('rm ' + path_collect + prefix+'.q_'+str(iqpt)+'/*wfc*')
                    if PAW:
                        os.system('cp ' + path_collect + prefix+'.q_'+str(iqpt)+'/'+prefix +
                                '.dvscf_paw ' + path_save + prefix+'.dvscf_paw_q'+label)
        else:
            # Case with XML format
            if XML:
                os.system('cp '+ path_ori + prefix+'.dyn0 '+ path_ori + prefix+'.dyn0.xml')
                os.system('cp '+path_ori + prefix+'.dyn'+str(iqpt)+'.xml ' + path_save + prefix +
                        '.dyn_q'+label+'.xml')
                if (iqpt == 1):
                    os.system('cp ' + path_collect + prefix+'.dvscf1 ' + path_save + prefix+'.dvscf_q' +
                            label)
                    os.system('cp -r ' + path_collect + prefix+'.phsave ' + path_save)
                    if fc:
                        os.system('cp '+ path_ori + prefix+'.fc.xml ' + path_save + 'ifc.q2r.xml')
                    if PAW:
                        os.system('cp ' + path_collect + prefix+'.dvscf_paw1 ' + path_save +prefix +
                                '.dvscf_paw_q'+label)
                else:
                    os.system('cp ' + path_collect + prefix+'.q_'+str(iqpt)+'/'+prefix +
                            '.dvscf1 ' + path_save + prefix+'.dvscf_q'+label)
                    # os.system('rm ' + path_collect + prefix+'.q_'+str(iqpt)+'/*wfc*')
                    if PAW:
                        os.system('cp ' + path_collect + prefix+'.q_'+str(iqpt)+'/'+prefix +
                                '.dvscf_paw1 ' + path_save + prefix+'.dvscf_paw_q'+label)
            # Case without XML format
            else:
                os.system('cp '+ path_ori + prefix+'.dyn'+str(iqpt)+' ' + path_save + prefix+'.dyn_q' +
                        label)
                if (iqpt == 1):
                    os.system('cp ' + path_collect + prefix + '.dvscf1 ' + path_save +prefix+'.dvscf_q' +
                            label)
                    os.system('cp -r ' + path_collect + prefix+'.phsave ' + path_save)
                    if fc:
                        os.system('cp '+ path_ori + prefix+'.fc ' + path_save + 'ifc.q2r')
                    if PAW:
                        os.system('cp ' + path_collect + prefix+'.dvscf_paw1 ' + path_save + prefix +
                                '.dvscf_paw_q'+label)
                else:
                    os.system('cp ' + path_collect + prefix+'.q_'+str(iqpt)+'/'+prefix +
                            '.dvscf1 ' + path_save + prefix+'.dvscf_q'+label)
                    # os.system('rm ' + path_collect + prefix+'.q_'+str(iqpt)+'/*wfc*')
                    if PAW:
                        os.system('cp ' + path_collect + prefix+'.q_'+str(iqpt)+'/'+prefix +
                                '.dvscf_paw1 ' + path_save +prefix+'.dvscf_paw_q'+label)

if __name__=="__main__":
    manipulate()
