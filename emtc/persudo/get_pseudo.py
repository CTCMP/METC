# -*- coding: utf-8 -*-
"""
Created on 19:54 04-11-2020 

@author: XY Ding
mail to: dxy_vasp@163.com
python3: pseudo.py
"""
import os.path as osp 
import os, sys
import configfile.inputpara as parameters
import module.funcs as bas
import AutoinputQE as inpp
import configfile.config as configg
sys.dont_write_bytecode = True        ## not generate __pycache__ files
############################ input parameters for QE #######################
class pseudo:
    def __init__(self):
        self.filepath = os.getcwd()
        self.latt, self.element, self.num_element, self.tot_ele, self.cor = bas.get_poscar(self.filepath, "POSCAR")

    def get_pseudo_sssp(self, label):
        filepath_pseudo = osp.join(configg.path_presudo_My_defined, label)
        name_pseudo = os.listdir(filepath_pseudo)
        pseudo_element = []
        for els in range(len(self.element)):
            for el in range(len(name_pseudo)):
                if str(self.element[els]+".").lower()  == str(name_pseudo[el][0:len(self.element[els])+1]).lower() or \
                    (str(self.element[els]+"_").lower() == str(name_pseudo[el][0:len(self.element[els])+1]).lower()):
                    pseudo_element.append(name_pseudo[el])
        return pseudo_element

    def get_pseudo_oncv(self):
        pseudo_element = []
        for els in range(len(self.element)):
            pseudo_element.append( self.element[els] + ".UPF")
        return pseudo_element

    def get_pseudo_oncv(self):
        pseudo_element = []
        for els in range(len(self.element)):
            pseudo_element.append( self.element[els] + ".UPF")
        return pseudo_element

    def get_pseudo_uspp(self, label):
        filepath_pseudo = osp.join(configg.path_presudo_My_defined, label)
        name_pseudo = os.listdir(filepath_pseudo)
        pseudo_element = []
        for els in range(len(self.element)):
            for el in range(len(name_pseudo)): 
                if (str(self.element[els]+".pbe").lower() ==str(name_pseudo[el][0:len(self.element[els])+4]).lower()) and ('us_psl.' in str(name_pseudo[el]).lower()):
                    pseudo_element.append(name_pseudo[el])
        if len(pseudo_element) == len(self.element):
            print("The pseudopential are rightly generated !")
        else:
            index_pseudo = {}
            for i in range(len(self.element)):
                tmp = []
                for j in range(len(pseudo_element)):
                    if (str(self.element[i]+".pbe").lower()==str(pseudo_element[j][0:len(self.element[i])+4]).lower()):
                        tmp.append(j)
                index_pseudo[self.element[i]] = tmp 
            ptmp = []
            for keys, vals in index_pseudo.items():
                ptmp.append(pseudo_element[vals[0]])
            pseudo_element = ptmp 
        return pseudo_element

def get_pseudo():
    pse = pseudo()
    try:
        if inpp.pseudo_potential == 'uspp':
            return pse.get_pseudo_uspp('uspp')
        elif inpp.pseudo_potential == 'oncv':
            return pse.get_pseudo_oncv('oncv')
        elif inpp.pseudo_potential == 'lda_oncv':
            return pse.get_pseudo_oncv('lda_oncv')
        elif inpp.pseudo_potential == 'pbesol_oncv':
            return pse.get_pseudo_oncv('pbesol_oncv')
        else:
            return pse.get_pseudo_sssp('sssp')
    except:
        return pse.get_pseudo_sssp('sssp')

pseudo_dict = {
    'sssp'          : 'sssp', 
    'uspp'          : 'rrkjus_psl', 
    'paw'           : 'kjpaw_psl', 
    # 'oncv'          : 'rrkj_psl', 
    'oncv'          : 'nc', 
    'van'           : 'van_psl',
}


class pseudoAll:
    def __init__(self):
        self.filepath = os.getcwd()
        self.latt, self.element, self.num_element, self.tot_ele, self.cor = bas.get_poscar(self.filepath, "POSCAR")

    def get_label(self):
        excorr = inpp.functional
        if inpp.SOC:
            excorr = 'rel-' + excorr 
        pseudoName = inpp.nl_corr + '-' + pseudo_dict[inpp.pseudo_potential]
        return pseudoName, excorr
    
    def path_pseudo(self):
        path_p = parameters.get_path_pseudo()
        return path_p
    
    def get_pseudo_oncv(self):
        # filepath_pseudo = osp.join(configg.path_presudo, inpp.pseudo_potential, inpp.functional)
        filepath_pseudo = self.path_pseudo()
        name_pseudo = os.listdir(filepath_pseudo)
        pseudo_element = []
        for els in range(len(self.element)):
            for el in range(len(name_pseudo)):
                if str(self.element[els]+".").lower()  == str(name_pseudo[el][0:len(self.element[els])+1]).lower():
                    pseudo_element.append(name_pseudo[el])
        print("pseudo oncv: ", pseudo_element)
        return pseudo_element

    def get_pseudo_pp(self, path):
        name_pseudo = os.listdir(osp.join(os.getcwd(), path))
        pseudo_element = []
        for els in range(len(self.element)):
            for el in range(len(name_pseudo)):
                if str(self.element[els]+".").lower()  == str(name_pseudo[el][0:len(self.element[els])+1]).lower():
                    pseudo_element.append(name_pseudo[el])
        print("pseudo sssp: ", pseudo_element)
        return pseudo_element

    def get_pseudo_sssp(self):
        filepath_pseudo = self.path_pseudo()
        # filepath_pseudo = osp.join(configg.path_presudo, inpp.pseudo_potential)
        name_pseudo = os.listdir(filepath_pseudo)
        pseudo_element = []
        for els in range(len(self.element)):
            for el in range(len(name_pseudo)):
                if str(self.element[els]+".").lower()  == str(name_pseudo[el][0:len(self.element[els])+1]).lower() or \
                    (str(self.element[els]+"_").lower() == str(name_pseudo[el][0:len(self.element[els])+1]).lower()):
                    pseudo_element.append(name_pseudo[el])
        print(pseudo_element)
        return pseudo_element
    
    def get_part_pseudo_ele(self, element, pseudoName, excorr): # get potential potential for single atoms
        filepath_pseudo = self.path_pseudo()
        name_pseudo = os.listdir(filepath_pseudo)
        pseudo_element = []                        
        for el in range(len(name_pseudo)):
            if str(element + ".").lower() == str(name_pseudo[el][0:len(element)+1]).lower() and \
                str(excorr).lower() in str(name_pseudo[el]).lower() and pseudoName.lower() in str(name_pseudo[el]).lower():
                pseudo_element.append(name_pseudo[el])
        print("#######################")
        print(pseudo_element)
        return pseudo_element
    
    def get_pseudo_others(self):                                # get all potential for all atoms 
        pseudoName, excorr = self.get_label()
        pseudo = []
        for ii in range(len(self.element)):
            tmp_ele = self.get_part_pseudo_ele(self.element[ii], pseudoName, excorr)
            if tmp_ele != []:
                pseudo.append(tmp_ele)
            else:
                pseudo.append(self.get_part_pseudo_ele(self.element[ii], '-' + pseudo_dict[inpp.pseudo_potential], excorr))
        pseudop = []
        for p in range(len(pseudo)):
            count = 0
            for pp in range(len(pseudo[p])):
                if len(pseudo[p]) == 1:
                    pseudop.append(pseudo[p][0])
                else:
                    if count > len(pseudo[p]):
                        pseudop.append(pseudo[p][0])
                    else:
                        # if '-n-' in pseudo[p][pp].lower() or '-spn-' in pseudo[p][pp].lower() or '-spdn-' in pseudo[p][pp].lower():
                        if 'n-' in pseudo[p][pp].lower() or '-spn-' in pseudo[p][pp].lower() or '-sl-' in pseudo[p][pp].lower() \
                            or '-spnl-' in pseudo[p][pp].lower():
                            pseudop.append(pseudo[p][pp])
                            break
                    count = count + 1
        print("************** potentials *****************") 
        if len(self.element) == len(pseudop):
            for k in range(len(pseudop)):
                print(f"%s : %s" % (self.element[k], pseudop[k]))
        else:
            print("Pseudo_potential are not rightly generated, please check the reasons. ")
        print("-------------------------------------------")
        return pseudop

def get_pseudo_all():
    if inpp.user_provide:
        pse = pseudoAll()
        return pse.get_pseudo_pp(inpp.path_pp)
    else:
        try:
            pse = pseudoAll()
            if inpp.pseudo_potential == 'sssp':
                return pse.get_pseudo_sssp()
            elif inpp.pseudo_potential == 'oncv':
                return pse.get_pseudo_oncv()
            else:
                return pse.get_pseudo_others()
        except:
            return get_pseudo()

class pseudoNew:
    def __init__(self):
        self.filepath = os.getcwd()
        self.latt, self.element, self.num_element, self.tot_ele, self.cor = bas.get_poscar(self.filepath, "POSCAR")

    def get_label(self):
        part1 = "."+inpp.functional + '-'
        part2 = "-" + inpp.pseudo_potential
        if inpp.SOC:
            part1 = 'rel-' + part1 
        return part1, part2
    
    def path_pseudo(self):
        path_p = parameters.get_path_pseudo()
        return path_p
    
    def get_pseudo_oncv(self):
        # filepath_pseudo = osp.join(configg.path_presudo, inpp.pseudo_potential, inpp.functional)
        filepath_pseudo = self.path_pseudo()
        name_pseudo = os.listdir(filepath_pseudo)
        pseudo_element = []
        for els in range(len(self.element)):
            for el in range(len(name_pseudo)):
                if str(self.element[els]+".").lower()  == str(name_pseudo[el][0:len(self.element[els])+1]).lower():
                    pseudo_element.append(name_pseudo[el])
        print("pseudo oncv: ", pseudo_element)
        return pseudo_element

    def get_pseudo_pp(self, path):
        name_pseudo = os.listdir(osp.join(os.getcwd(), path))
        pseudo_element = []
        for els in range(len(self.element)):
            for el in range(len(name_pseudo)):
                if str(self.element[els]+".").lower() == str(name_pseudo[el][0:len(self.element[els])+1]).lower():
                    pseudo_element.append(name_pseudo[el])
        print("pseudo sssp: ", pseudo_element)
        return pseudo_element

    def get_pseudo_sssp(self):
        filepath_pseudo = self.path_pseudo()
        # filepath_pseudo = osp.join(configg.path_presudo, inpp.pseudo_potential)
        name_pseudo = os.listdir(filepath_pseudo)
        pseudo_element = []
        for els in range(len(self.element)):
            for el in range(len(name_pseudo)):
                if str(self.element[els]+".").lower()  == str(name_pseudo[el][0:len(self.element[els])+1]).lower() or \
                    (str(self.element[els]+"_").lower() == str(name_pseudo[el][0:len(self.element[els])+1]).lower())or \
                    (str(self.element[els]+"-").lower() == str(name_pseudo[el][0:len(self.element[els])+1]).lower()):
                    pseudo_element.append(name_pseudo[el])
        print(pseudo_element)
        return pseudo_element
    
    def get_part_pseudo_ele(self, element, part1, part2): # get potential potential for single atoms
        filepath_pseudo = self.path_pseudo()
        name_pseudo = os.listdir(filepath_pseudo)
        pseudo_element = []                        
        for el in range(len(name_pseudo)):
            pesudo_ele = element + part1
            if pesudo_ele.lower() in name_pseudo[el].lower() and part2.lower() in name_pseudo[el].lower():
                pseudo_element.append(name_pseudo[el])
        # print(name_pseudo)
        return pseudo_element
    
    def get_pseudo_others(self):                                # get all potential for all atoms 
        part1, part2 = self.get_label()
        # print(part1, part2)
        pseudo = []
        for ii in range(len(self.element)):
            tmp_ele = self.get_part_pseudo_ele(self.element[ii], part1, part2)
            pseudo.append(tmp_ele)
        # print(pseudo)
        pseudop = []
        for p in range(len(pseudo)):
            if len(pseudo[p]) == 1:
                pseudop.append(pseudo[p][0])
            else:
                have_pse = True
                for pp in range(len(pseudo[p])):
                    if inpp.nl_corr  in pseudo[p][pp]:
                        pseudop.append(pseudo[p][pp])
                        have_pse = False 
                        break 
                if have_pse:
                    pseudop.append(pseudo[p][0])
        print("************** potentials *****************") 
        if len(self.element) == len(pseudop):
            for k in range(len(pseudop)):
                print(f"%s : %s" % (self.element[k], pseudop[k]))
        else:
            print("Pseudo_potential are not rightly generated, please check the reasons. ")
        print("-------------------------------------------")
        return pseudop

def get_pseudo_all_neW():
    if inpp.user_provide:
        pse = pseudoNew()
        return pse.get_pseudo_pp(inpp.path_pp)
    else:
        try:
            pse = pseudoNew()
            if inpp.pseudo_potential == 'sssp':
                return pse.get_pseudo_sssp()
            elif inpp.pseudo_potential == 'oncv':
                return pse.get_pseudo_oncv()
            else:
                return pse.get_pseudo_others()
        except:
            return get_pseudo()