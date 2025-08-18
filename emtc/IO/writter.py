import numpy as np
import os.path as osp 
import os, sys
import configfile.inputpara as parameters
import module.funcs as bas
import json
import copy
sys.path.append(os.getcwd())
import AutoinputQE as inpp
sys.dont_write_bytecode = True       ## not generate __pycache__ files

mag = bas.mag_judge()
def get_pseudopotential():
    import configfile.get_pseudo as pseu
    # return pseu.get_pseudo_all()
    return pseu.get_pseudo_all_neW()

def write_checkd0(conf, par):
    for keys, val in par.items():
        conf.write(str(keys).ljust(20, " ") + "  =  ")
        if '1.0d' in str(val).lower() or 'd0' in str(val).lower():
            conf.write(str(val).lower().ljust(len(str(val))))
        elif keys == 'exclude_bands':
            conf.write(str(val).lower().ljust(len(str(val))))
        else:
            json.dump(val, conf)
        conf.write(",")
        conf.write("\n")
    return conf
def write_LDAU(conf):
    conf.write('HUBBARD (' + inpp.Hubbard_projectors + ') \n')
    latt, element, num_element, tot_ele, cor = bas.get_poscar(os.getcwd(), "POSCAR")
    for keys, vals in inpp.LDAU_par.items():
        conf.write("U " + str(element[int(keys)-1]) + '-' + str(vals[0]) + ' ' + str(vals[1]) + '\n')
    conf.write('\n')
    return conf
def write_QEinput_CONTROL(conf, par):
    conf.write('&CONTROL \n')
    conf = write_checkd0(conf, par)
    conf.write('/ \n \n')
    return conf
def write_QEinput_SYSTEM(conf, par):
    conf.write('&SYSTEM \n')
    conf = write_checkd0(conf, par)
    conf.write('/ \n \n')
    return conf
def write_QEinput_ELECTRONS(conf, par):
    conf.write('&ELECTRONS \n')
    conf = write_checkd0(conf, par)
    conf.write('/ \n \n')
    return conf
def write_QEinput_IONS(conf, par):
    conf.write('&IONS \n')
    conf = write_checkd0(conf, par)
    conf.write('/ \n \n')
    return conf
def write_QEinput_CELL(conf, par):
    conf.write('&CELL \n')
    conf = write_checkd0(conf, par)
    conf.write('/ \n \n')
    return conf
def write_cellParameters(conf, latt):
    conf.write('CELL_PARAMETERS (angstrom)' + "\n")
    for j in range(0, 3):
        for k in range(0, 3):
            conf.write("  " + str(latt[j, k]).rjust(13, " ") + "  ")
        conf.write("\n")
    conf.write("\n")
    return conf
def write_atomP(conf):
    latt, element, num_element, tot_ele, cor = bas.get_poscar(os.getcwd(), "POSCAR")
    conf.write("ATOMIC_SPECIES \n")
    pseudo_allele = get_pseudopotential()   
    for ii in range(len(element)):
        conf.write(str(element[ii]).ljust(4, " ") + '  ')  
        try:
            conf.write(str(parameters.atom_mass[element[ii]]).rjust(10, " ") + "  ")  
        except:
            len_atom = len(element[ii])
            conf.write(str(parameters.atom_mass[element[ii][:len_atom-1]]).rjust(10, " ") + "  ")        
        conf.write(pseudo_allele[ii])
        conf.write("\n")
    conf.write("\n")                                  
    if inpp.ibrav == 0:
        write_cellParameters(conf, latt)
    conf.write('ATOMIC_POSITIONS (crystal)' + "\n")
    count = 0
    for jj in range(len(element)):
        for kk in range(num_element[jj]):
            conf.write(element[jj] + "    ")
            for h in range(3):
                conf.write(str(cor[count, h]).rjust(15, " ") + "   ")
            conf.write("\n")
            count = count  + 1
    conf.write("\n")
    if inpp.LDAU:
        write_LDAU(conf)
    return conf
def write_atomP_ibrav0(conf):
    latt, element, num_element, tot_ele, cor = bas.get_poscar(os.getcwd(), "POSCAR")
    conf.write("ATOMIC_SPECIES \n")
    pseudo_allele = get_pseudopotential()   
    for ii in range(len(element)):
        conf.write(str(element[ii]).ljust(4, " ") + '  ')  
        conf.write(str(parameters.atom_mass[element[ii]]).rjust(10, " ") + "  ")         
        conf.write(pseudo_allele[ii])
        conf.write("\n")
    conf.write("\n")                                  
    write_cellParameters(conf, latt)
    conf.write('ATOMIC_POSITIONS (crystal)' + "\n")
    count = 0
    for jj in range(len(element)):
        for kk in range(num_element[jj]):
            conf.write(element[jj] + "    ")
            for h in range(3):
                conf.write(str(cor[count, h]).rjust(15, " ") + "   ")
            conf.write("\n")
            count = count  + 1
    conf.write("\n")
    if inpp.LDAU:
        write_LDAU(conf)
    return conf

def write_kpoint_scf(conf):
    conf.write("\n")
    conf.write("K_POINTS {automatic} \n")
    kpoints1 = inpp.kpt_shift
    for i in range(len(inpp.kpt)):
        conf.write(str(inpp.kpt[i]).ljust(4, " "))
    for kpp in range(len(kpoints1)):
        conf.write(str(kpoints1[kpp]).ljust(4, " "))
    conf.write("\n \n")
    return conf
def write_kpoint_scf_kpt(conf, kpt):
    conf.write("\n")
    conf.write("K_POINTS {automatic} \n")
    kpoints1 = inpp.kpt_shift
    for i in range(len(kpt)):
        conf.write(str(kpt[i]).ljust(4, " "))
    for kpp in range(len(kpoints1)):
        conf.write(str(kpoints1[kpp]).ljust(4, " "))
    conf.write("\n \n")
    return conf
def write_kmesh_pl(conf, filepath):
    conf.write('\n')
    kpoints = inpp.kpt_nscf
    str_kpt = str(kpoints[0]) + " " + str(kpoints[1]) + " " + str(kpoints[2])
    if osp.exists(osp.join(filepath, "kpoints_kmesh")):
        os.system('rm kpoints_kmesh')
    os.system('kmesh.pl ' + str_kpt +' >> kpoints_kmesh')
    data = open(osp.join(filepath, "kpoints_kmesh"), 'r').readlines()
    for kkk in range(len(data)):
        conf.write(data[kkk])
    return conf
def write_kmesh_pl2(conf, filepath, kpt_nscf):
    conf.write('\n')
    kpoints = kpt_nscf
    str_kpt = str(kpoints[0]) + " " + str(kpoints[1]) + " " + str(kpoints[2])
    if osp.exists(osp.join(filepath, "kpoints_kmesh")):
        os.system('rm kpoints_kmesh')
    os.system('kmesh.pl ' + str_kpt +' >> kpoints_kmesh')
    data = open(osp.join(filepath, "kpoints_kmesh"), 'r').readlines()
    for kkk in range(len(data)):
        conf.write(data[kkk])
    return conf
def write_kpoints_band(conf):
    conf.write("\n")
    kpoints_bands = inpp.kpoints_band.split("\n")[1:-1]
    conf.write("K_POINTS {crystal_b}" + '\n')
    conf.write(str(len(kpoints_bands)) + '\n')
    for pp in range(len(kpoints_bands)):
        if pp != len(kpoints_bands) -1:
            conf.write(kpoints_bands[pp].split("!")[0].ljust(len(kpoints_bands[pp].split("!")[0])) + "  " + str(inpp.band_points) + '\n')
        else:
            conf.write(kpoints_bands[pp].split("!")[0].ljust(len(kpoints_bands[pp].split("!")[0])) + "  " + str(1) + '\n')
    return conf
def write_postband(filepath, fname, parameters_postband):
    with open(osp.join(filepath, fname), 'w') as postband:
        postband.write("&bands")
        postband.write('\n')
        for keys, vals in parameters_postband.items():
            postband.write(keys.ljust(10, " ") + " = ")
            json.dump(vals, postband)
            postband.write(',')
            postband.write('\n')
        postband.write('/')
        postband.write('\n')   
    postband.close() 
def write_pro_band(filepath, filename):
    with open(osp.join(filepath, filename), 'w') as pro:
        pro.write("&projwfc \n")
        for keys, vals in parameters.pro_band.items():
            pro.write(keys + ' = ' )
            json.dump(vals, pro)
            pro.write(',\n')
        pro.write('/\n')
    pro.close()
def write_band_post_input(parameters=parameters, filepath=os.getcwd()):
    if mag == 2:
        par_up = parameters.band_plot
        par_up['filband'] = "band_up.dat"
        par_up['spin_component'] = 1
        write_postband(filepath, 'postBand_up.in', par_up)
        par_dw = parameters.band_plot
        par_dw['filband'] = "band_dw.dat"
        par_dw['spin_component'] = 2
        write_postband(filepath, 'postBand_dw.in', par_dw)
        write_pro_band(filepath, 'proBand.in')
    else:
        write_postband(filepath, 'postBand.in', parameters.band_plot)
        write_pro_band(filepath, 'proBand.in')
def write_dos(par, filepath=os.getcwd()):
    conf = open(osp.join(filepath, "dos.in"), 'w')  
    conf.write(str('&PROJWFC') + "\n")    
    # conf.write("&inputph \n")               
    conf = write_checkd0(conf, par)
    conf.write("/" + "\n")
    conf.write("\n")
    conf.close() 
def write_phonon(par, filepath=os.getcwd()):
    conf = open(osp.join(filepath, "phonon.in"), 'w') 
    conf.write(str('phonons of systems') + "\n")    
    conf.write("&inputph \n")               
    conf = write_checkd0(conf, par)
    conf.write("/" + "\n")
    conf.write("\n")
    conf.close()  
def write_phononTEST(par, filepath=os.getcwd()):
    conf = open(osp.join(filepath, "phononEPW.in"), 'w') 
    for i in range(3):
        del par['nq'+str(i+1)]
    conf.write(str('phonons of systems') + "\n")    
    conf.write("&inputph \n")               
    conf = write_checkd0(conf, par)
    conf.write("/" + "\n")
    conf.write("\n")
    kpoints_bands = inpp.kpoints_band.split("\n")[1:-1]
    # conf.write("K_POINTS {crystal_b}" + '\n')
    conf.write(str(len(kpoints_bands)) + '\n')
    for pp in range(len(kpoints_bands)):
        if pp != len(kpoints_bands) -1:
            conf.write(kpoints_bands[pp].split("!")[0].ljust(len(kpoints_bands[pp].split("!")[0])) + "  " + str(inpp.phon_points) + '\n')
        else:
            conf.write(kpoints_bands[pp].split("!")[0].ljust(len(kpoints_bands[pp].split("!")[0])) + "  " + str(1) + '\n')
    conf.close() 
    conf.close() 
def write_q2r(par, filepath=os.getcwd()):
    conf = open(osp.join(filepath, "q2r.in"), 'w') 
    conf.write(str('&input') + "\n")                
    conf = write_checkd0(conf, par)
    conf.write("/" + "\n")
    conf.write("\n")
    conf.close() 
def write_matdyn(par, filepath=os.getcwd()):
    conf = open(osp.join(filepath, "matdyn.in"), 'w')
    conf.write(str('&input') + "\n")                
    conf = write_checkd0(conf, par)
    conf.write("/" + "\n")
    conf.write("\n")
    kpoints_bands = inpp.kpoints_band.split("\n")[1:-1]
    # conf.write("K_POINTS {crystal_b}" + '\n')
    conf.write(str(len(kpoints_bands)) + '\n')
    for pp in range(len(kpoints_bands)):
        if pp != len(kpoints_bands) -1:
            conf.write(kpoints_bands[pp].split("!")[0].ljust(len(kpoints_bands[pp].split("!")[0])) + "  " + str(inpp.phon_points) + '\n')
        else:
            conf.write(kpoints_bands[pp].split("!")[0].ljust(len(kpoints_bands[pp].split("!")[0])) + "  " + str(1) + '\n')
    conf.close() 
def write_phdos(par, filepath=os.getcwd()):
    conf = open(osp.join(filepath, "phdos.in"), 'w') 
    conf.write(str('&input') + "\n")                
    conf = write_checkd0(conf, par)
    conf.write("/" + "\n")
    conf.write("\n")
    conf.close() 
def write_w90(par, filepath):
    if mag == 2:
        with open(osp.join(filepath, "w90_up.in"), 'w') as w90upcal:
            w90upcal.write("&inputpp")
            w90upcal.write('\n')
            par['spin_component'] = 'up'
            par['seedname'] = 'pwscf_up'
            w90upcal = write_checkd0(w90upcal, par)
            w90upcal.write('/')
            w90upcal.write('\n')
        w90upcal.close()
        with open(osp.join(filepath, "w90_dw.in"), 'w') as w90dwcal:
            w90dwcal.write("&inputpp")
            w90dwcal.write('\n')
            par['spin_component'] = 'down'
            par['seedname'] = 'pwscf_dw'
            w90dwcal = write_checkd0(w90dwcal, par)
            w90dwcal.write('/')
            w90dwcal.write('\n')
        w90dwcal.close()
    else:
        with open(osp.join(filepath, "w90.in"), 'w') as w90cal:
            w90cal.write("&inputpp")
            w90cal.write('\n')
            w90cal = write_checkd0(w90cal, par)
            w90cal.write('/')
            w90cal.write('\n')
        w90cal.close()

def transform_cor(latt, cor):
    if np.max(cor) < 1:
        corr = np.zeros(shape=(len(cor), 3))
        for iii in range(len(cor)):
            corr[iii] = latt[0] * cor[iii, 0] + latt[1] * cor[iii, 1] + latt[2] * cor[iii, 2]
        return corr
    else:
        return cor
def write_w90_win_temp(para, filepath, filename):
    latt, element, num_element, tot_ele, cor = bas.get_poscar(os.getcwd(), "POSCAR")
    with open(osp.join(filepath, filename), 'w') as wannier:
        tmp = inpp.wannier90['num_wann']
        if mag == 2:
            para['num_wann'] = int(tmp/2)
        wannier = write_checkd0(wannier, para)
        ############### projection #########################
        wannier.write("\n")
        wannier.write("Begin Projections")
        wannier.write("\n")
        for keey, vaal in inpp.porjection_w90.items():
            wannier.write(keey)
            wannier.write(":")
            wannier.write(vaal)
            wannier.write("\n")
        wannier.write("End Projections")
        wannier.write("\n")
        wannier.write("\n")
        ############### KPOINTS #########################
        wannier.write("begin kpoint_path" + '\n')
        kpoints_bands = [items.split("!")[0] for items in inpp.kpoints_band.split("\n")[1:-1]]
        hsp_labels = [str(items.split("!")[1]).strip() for items in inpp.kpoints_band.split("\n")[1:-1]]
        for i in range(len(hsp_labels)):
            if str(hsp_labels[i]).lower() == "\Gamma".lower():
                hsp_labels[i] = "G"
        for pp in range(len(kpoints_bands)-1):
            wannier.write(hsp_labels[pp].ljust(2) + " ")
            wannier.write(kpoints_bands[pp] + " ")
            wannier.write(hsp_labels[pp+1].ljust(2) + " ")
            wannier.write(kpoints_bands[pp+1] + " ")
            wannier.write("\n")
        wannier.write("end kpoint_path" + '\n')
        wannier.write("\n")
        wannier.write("\n")
        wannier.write("begin unit_cell_cart" + '\n')
        # wannier.write("Ang" + '\n')
        for j in range(0, 3):
            for k in range(0, 3):
                wannier.write("  " + str(latt[j, k]).rjust(13, " ") + "  ")
            wannier.write("\n")
        wannier.write("end unit_cell_cart" + '\n')
        wannier.write("\n")
        wannier.write("\n")
        ############### structure position #########################
        wannier.write("begin atoms_cart" + '\n')
        # wannier.write("Ang" + '\n')
        corr = transform_cor(latt, cor)
        count = 0
        for jj in range(len(element)):
            for kk in range(num_element[jj]):
                wannier.write(element[jj] + "    ")
                for h in range(3):
                    wannier.write(str(corr[count, h]).rjust(20, " ") + "   ")
                wannier.write("\n")
                count = count  + 1
        wannier.write("end atoms_cart" + '\n')
    wannier.close()
def write_w90_win(par, filepath):
    if mag == 2:
        write_w90_win_temp(par['up'], filepath, 'pwscf_up.win')
        write_w90_win_temp(par['dw'], filepath, 'pwscf_dw.win')
    else:
        write_w90_win_temp(par, filepath, 'pwscf.win')
def epw_w(conf, par):
    for keys, val in par.items():
        conf.write(str(keys).ljust(20, " ") + "  =  ")
        if '1.0d' in str(val).lower() or 'd0' in str(val).lower():
            conf.write(str(val).lower().ljust(len(str(val))))
        elif isinstance(val,list):
            for pp in range(len(val)):
                conf.write(str(val[pp]).lower().ljust(len(str(val[pp])) + 1))
        else:
            json.dump(val, conf)
        conf.write(",")
        conf.write("\n")
    return conf
def write_epw(par, filepath, filename):
    conf = open(osp.join(filepath, filename), 'w')
    conf.write("epw calculation \n")
    conf.write("&inputepw \n")
    conf = epw_w(conf, par)
    conf.write("/" + "\n")
    conf.write("\n")
    conf.close()

def write_phonopy_2nd(para, filepath):
    latt, element, num_element, tot_ele, cor = bas.get_poscar(os.getcwd(), "POSCAR")
    with open(osp.join(filepath, "2nd.in"), 'w') as conf:
        for keys, values in para.items():                       ############### basic input parameters
            conf.write("&"+keys +"\n")
            for keyss, val in values.items():
                conf.write(str(keyss).ljust(20, " ") + "  =  ")
                if '1.0d' in str(val).lower() or 'd0' in str(val).lower():
                    conf.write(str(val).lower().ljust(len(str(val))))
                elif 'NATOMS' in str(val):
                    conf.write(str(val))
                else:
                    json.dump(val, conf)
                conf.write(", \n")
            conf.write("/" + "\n \n")
        conf.write("\n")
        ############### KPOINTS #########################
        conf.write("\n")
        conf.write("K_POINTS {automatic} \n")
        kpoints1 = inpp.kpt_shift
        for i in range(len(inpp.kpt_supercell)):
            conf.write(str(inpp.kpt_supercell[i]).ljust(4, " "))
        for kpp in range(len(kpoints1)):
            conf.write(str(kpoints1[kpp]).ljust(4, " "))
    conf.close()

def write_3rd(para, filepath):
    latt, element, num_element, tot_ele, cor = bas.get_poscar(os.getcwd(), "POSCAR")
    with open(osp.join(filepath, "scf_sc.in"), 'w') as conf:
        for keys, values in para.items():                       ############### basic input parameters
            conf.write("&"+keys +"\n")
            for keyss, val in values.items():
                conf.write(str(keyss).ljust(20, " ") + "  =  ")
                if '1.0d' in str(val).lower() or 'd0' in str(val).lower():
                    conf.write(str(val).lower().ljust(len(str(val))))
                elif 'NATOMS' in str(val):
                    conf.write(str(val))
                else:
                    json.dump(val, conf)
                conf.write(", \n")
            conf.write("/" + "\n \n")
        conf.write("\n")
        conf.write("ATOMIC_SPECIES \n")
        pseudo_allele = get_pseudopotential()                                       ############### configure pseudo file
        for ii in range(len(element)):
            conf.write(str(element[ii]).ljust(4, " ") + '  ')  
            conf.write(str(parameters.atom_mass[element[ii]]).rjust(10, " ") + "  ")         
            conf.write(pseudo_allele[ii])
            conf.write("\n")
        conf.write("\n")
        conf.write('##COORDINATES##')
        ############### KPOINTS #########################
        conf.write("\n")
        conf.write("K_POINTS {automatic} \n")
        kpoints1 = inpp.kpt_shift
        for i in range(len(inpp.kpt_supercell)):
            conf.write(str(inpp.kpt_supercell[i]).ljust(4, " "))
        for kpp in range(len(kpoints1)):
            conf.write(str(kpoints1[kpp]).ljust(4, " "))
        conf.write("\n \n")
        conf.write('##CELL##')
        conf.write("\n \n")
    conf.close()
def write_4th(para, filepath):
    latt, element, num_element, tot_ele, cor = bas.get_poscar(os.getcwd(), "POSCAR")
    with open(osp.join(filepath, "scf_sc.in"), 'w') as conf:
        for keys, values in para.items():                       ############### basic input parameters
            conf.write("&"+keys +"\n")
            for keyss, val in values.items():
                conf.write(str(keyss).ljust(20, " ") + "  =  ")
                if '1.0d' in str(val).lower() or 'd0' in str(val).lower():
                    conf.write(str(val).lower().ljust(len(str(val))))
                elif 'NATOMS' in str(val):
                    conf.write(str(val))
                else:
                    json.dump(val, conf)
                conf.write(", \n")
            conf.write("/" + "\n \n")
        conf.write("\n")
        conf.write("ATOMIC_SPECIES \n")
        pseudo_allele = get_pseudopotential()                                       ############### configure pseudo file
        for ii in range(len(element)):
            conf.write(str(element[ii]).ljust(4, " ") + '  ')  
            conf.write(str(parameters.atom_mass[element[ii]]).rjust(10, " ") + "  ")         
            conf.write(pseudo_allele[ii])
            conf.write("\n")
        conf.write("\n")
        conf.write('##COORDINATES##')
        ############### KPOINTS #########################
        conf.write("\n")
        conf.write("K_POINTS {automatic} \n")
        kpoints1 = inpp.kpt_shift
        for i in range(len(inpp.kpt_supercell)):
            conf.write(str(inpp.kpt_supercell[i]).ljust(4, " "))
        for kpp in range(len(kpoints1)):
            conf.write(str(kpoints1[kpp]).ljust(4, " "))
        conf.write("\n \n")
        conf.write('##CELL##')
        conf.write("\n \n")
    conf.close()
def write_input_shengBTE(para, filepath):
    conf = open(osp.join(filepath, 'CONTROL'), 'w')
    conf.write("&allocations \n")
    for key_allo, val_allo in para['allocations'].items():
        conf.write("   " + key_allo + "=")
        if 'ngrid' in key_allo:
            conf.write(str(val_allo))
        else:
            json.dump(val_allo, conf)
        conf.write(', \n')
    conf.write('&end \n')
    conf.write("\n")
    conf.write('&crystal \n')
    for key_cry, val_cry in para['crystals'].items():
        if  'elements' in key_cry:
            conf.write("   " + key_cry + '=')
            for iee in val_cry:
                conf.write('\"')  
                # conf.write(iee + " ")
                conf.write(iee)
                # json.dump(iee , conf)
                conf.write('\" ')
            conf.write(', \n')
        elif 'types' in key_cry:
            conf.write("   " + key_cry + '=')
            for itt in val_cry.split():
                conf.write(itt + " ")
            # json.dump(val_cry, conf)
            conf.write(', \n')
        elif 'born' in key_cry:
            conf.write(str(val_cry) + '\n')
        else:
            conf.write("   " + key_cry + '=' + str(val_cry) + ', \n')
    conf.write('&end \n')
    conf.write("\n")
    conf.write('&parameters \n')
    for key_pars, val_pars in para['parameter'].items():
        conf.write("   " + key_pars + '=')
        json.dump(val_pars, conf)
        conf.write(', \n')
    conf.write('&end \n')
    conf.write("\n")
    conf.write('&flags \n')
    for key_flag, val_flag in para['flags'].items():
        conf.write("   " + key_flag + '=')
        json.dump(val_flag, conf)
        conf.write(', \n')
    conf.write('&end \n')
    conf.write("/" + "\n")
    conf.write("\n")
    conf.close()
def write_qe2pert(par, filepath, filename):
    qp = open(osp.join(filepath, filename), 'w')
    qp.write("&qe2pert \n")
    for key_pert, val_pert in par['qe2pert'].items():
        qp.write("   " + key_pert + "=")
        json.dump(val_pert, qp)
        qp.write(', \n')  
    qp.write('/ \n') 
    qp.close()
def write_perturbo(par, filepath, filename):
    qp = open(osp.join(filepath, filename), 'w')
    qp.write("&perturbo \n")
    for key_pert, val_pert in par['perturbo'].items():
        qp.write("   " + key_pert + "=")
        json.dump(val_pert, qp)
        qp.write(', \n')  
    qp.write('/ \n') 
    qp.close()

class WriCls():
    def __init__(self, filepath=os.getcwd()):
        self.filepath = filepath

    def QE_input_RelaxScfNscf(self, type, par, filepath=os.getcwd()):
        conf = open(osp.join(filepath, type+'.in'), 'w')
        conf = write_QEinput_CONTROL(conf, par['CONTROL'])
        conf = write_QEinput_SYSTEM(conf, par['SYSTEM'])
        conf = write_QEinput_ELECTRONS(conf, par['ELECTRONS'])
        conf = write_QEinput_IONS(conf, par['IONS'])
        conf = write_QEinput_CELL(conf, par['CELL'])
        conf = write_atomP(conf)
        if type == 'bands':
            write_kpoints_band(conf)
        else:
            conf = write_kpoint_scf(conf)
        conf.close()
    def QE_input_RelaxScfNscf_kpt(self, type, par, filepath=os.getcwd()):
        conf = open(osp.join(filepath, type+'.in'), 'w')
        conf = write_QEinput_CONTROL(conf, par['CONTROL'])
        conf = write_QEinput_SYSTEM(conf, par['SYSTEM'])
        conf = write_QEinput_ELECTRONS(conf, par['ELECTRONS'])
        conf = write_QEinput_IONS(conf, par['IONS'])
        conf = write_QEinput_CELL(conf, par['CELL'])
        conf = write_atomP(conf)
        conf = write_kpoint_scf_kpt(conf, inpp.kpt)
        conf.close()
    def QE_input_RelaxScfNscf_ibrav0(self, type, par, filepath=os.getcwd()):
        conf = open(osp.join(filepath, type+'.in'), 'w')
        conf = write_QEinput_CONTROL(conf, par['CONTROL'])
        conf = write_QEinput_SYSTEM(conf, par['SYSTEM'])
        conf = write_QEinput_ELECTRONS(conf, par['ELECTRONS'])
        conf = write_QEinput_IONS(conf, par['IONS'])
        conf = write_QEinput_CELL(conf, par['CELL'])
        conf = write_atomP_ibrav0(conf)
        if type == 'bands':
            write_kpoints_band(conf)
        else:
            conf = write_kpoint_scf(conf)
        conf.close()
    def QE_input_RelaxScfNscf2(self, par, filepath):
        conf = open(osp.join(filepath, 'nscf.in'), 'w')
        conf = write_QEinput_CONTROL(conf, par['CONTROL'])
        conf = write_QEinput_SYSTEM(conf, par['SYSTEM'])
        conf = write_QEinput_ELECTRONS(conf, par['ELECTRONS'])
        conf = write_QEinput_IONS(conf, par['IONS'])
        conf = write_QEinput_CELL(conf, par['CELL'])
        conf = write_atomP(conf)
        conf = write_kmesh_pl(conf, filepath)
        conf.close()
    def QE_input_RelaxScfNscf3(self, par, filepath):
        conf = open(osp.join(filepath, 'nscf.in'), 'w')
        conf = write_QEinput_CONTROL(conf, par['CONTROL'])
        conf = write_QEinput_SYSTEM(conf, par['SYSTEM'])
        conf = write_QEinput_ELECTRONS(conf, par['ELECTRONS'])
        conf = write_QEinput_IONS(conf, par['IONS'])
        conf = write_QEinput_CELL(conf, par['CELL'])
        conf = write_atomP(conf)
        conf = write_kmesh_pl2(conf, filepath, inpp.qpoints)
        conf.close()
    def QE_input_postband(self, par=parameters, filepath=os.getcwd()):
        write_band_post_input(par, filepath)
    def QE_dos(self, par, filepath):
        write_dos(par, filepath)
    def QE_phonon(self, par, filepath):
        write_phonon(par, filepath)
    def QE_phononTEST(self, par, filepath):
        write_phononTEST(par, filepath)
    def QE_q2r(self, par, filepath):
        write_q2r(par, filepath)
    def QE_matdyn(self, par, filepath):
        write_matdyn(par, filepath)
    def QE_phdos(self, par, filepath):
        write_phdos(par, filepath)
    def w90(self, par, filepath):
        write_w90(par, filepath)
    def w90_win(self, par, filepath):
        write_w90_win(par, filepath)
    def epw(self, par, filepath, filename):
        write_epw(par, filepath, filename)
    def phonopy_2nd(self, par, filepath=os.getcwd()):
        write_phonopy_2nd(par, filepath)
    def ShengBTE_3rd(self, par, filepath=os.getcwd()):
        write_3rd(par, filepath)
    def ShengBTE_4th(self, par, filepath=os.getcwd()):
        write_4th(par, filepath)
    def ShengBTE(self, par, filepath=os.getcwd()):
        write_input_shengBTE(par, filepath)
    def qe2pert(self, par, filepath, filename):
        write_qe2pert(par, filepath, filename)
    def perturbo(self, par, filepath, filename):
        write_perturbo(par, filepath, filename)

