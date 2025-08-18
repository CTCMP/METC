import os, sys
import os.path as osp 
from configfile.inputpara import path_QE
from configfile.inputpara import python_path
import configfile.inputpara as parameters
import AutoinputQE as inpp
import configfile.config as configg
import module.funcs as bas 
curPath = os.path.abspath(os.path.dirname(__file__))
scriptpath = curPath.replace(curPath.split(r'/')[-1], '')
sys.path.append(os.getcwd())
sys.dont_write_bytecode = True

mag = bas.mag_judge()
#--------------------------- job module ----------------------------
cal_order = configg.cal_order
item = ''
cal_order_w90 = cal_order.replace(cal_order.split()[-1], str(24))
# cal_order = "srun "
NP = int(inpp.num_process*inpp.num_nodes)
NPP = len(inpp.kpoints_band.split("\n")[1:-1])
cal_order_ephmat = ''
for cc in range(len(cal_order.split())-1):
    cal_order_ephmat = cal_order_ephmat + cal_order.split()[cc] + ' '
cal_order_ephmat = cal_order_ephmat + ' ' 
# NP = 4
NP_pert = 4        # parallel process for pertrubo
try:
    if inpp.mpi_par.strip() == '0':
        npool=False 
    else:
        npool = True  
    if npool:
        mpinpool = "  " + inpp.mpi_par + "  "
        mpipo = ' -npool ' + str(NPP)
    else:
        mpinpool = ''

    npool_must = True
    if npool_must:
        mpinpool_must = ' -npool ' + str(NP)
        mpipo_must = ' -npool ' + str(NPP)
except:
    npool=False 
    if npool:
        mpinpool = ' -npool ' + str(NP)
        mpipo = ' -npool ' + str(NPP)
    else:
        mpinpool = ''

    npool_must = True
    if npool_must:
        mpinpool_must = ' -npool ' + str(NP)
        mpipo_must = ' -npool ' + str(NPP)

if mag == 2:
    band_post = ["mpirun -np 20 " + " bands.x " + " < postBand_up.in > postband_up.log 2>&1 \n",
                 "mpirun -np 20 " + " bands.x " + " < postBand_dw.in > postband_dw.log 2>&1 \n"]
    w90_post = [cal_order_w90 + "  " + "pw2wannier90.x < w90_up.in > pw2Wannier90_up.log 2>&1 \n", 
                cal_order_w90 + "  " + "pw2wannier90.x < w90_dw.in > pw2Wannier90_dw.log 2>&1 \n"]
    w90_distangle = [cal_order + "  " + "wannier90.x pwscf_up \n", cal_order + "  " + "wannier90.x pwscf_dw \n"]
    w90_pp = ["wannier90.x -pp pwscf_up \n", "wannier90.x -pp pwscf_dw \n"]
else:
    band_post = ["mpirun -np 20 " + " bands.x " + " < postBand.in > postband.log 2>&1 \n"]    
    w90_post = [cal_order_w90 + "  " + "pw2wannier90.x  < w90.in > pw2Wannier90.log 2>&1 \n"]
    w90_distangle = [cal_order + "  wannier90.x pwscf \n"]
    w90_pp = ["wannier90.x -pp pwscf \n"]

def write_checkd0(conf, par):
    import json
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

item = ''
ph_xTEST = configg.QE_path.replace('/bin', '')
submitscripts = configg.submit_scripts
submit_order = {
    'pw_vc-relax'    : cal_order + "  " + "pw.x" + mpinpool + " < vc-relax.in > vc-relax.log 2>&1 \n",
    'pw_relax'       : cal_order + "  " + "pw.x" + mpinpool + " < relax.in > relax.log 2>&1 \n",
    'pw_scf'         : cal_order + "  " + "pw.x" + mpinpool + " < scf.in > scf.log 2>&1 \n",
    'pw_nscf'        : cal_order + "  " + "pw.x" + mpinpool + " < nscf.in > nscf.log 2>&1 \n",
    'pw_bands'       : cal_order + "  " + "pw.x" + mpinpool + " < bands.in > bands.log 2>&1 \n",
    'pro_bands'      : cal_order + "  " + "projwfc.x" + mpinpool + " < proBand.in > proBand.log 2>&1 \n",
    'band_plot'      : str(python_path) + ' -u ' + str(osp.join(scriptpath, 'scripts', 'band.py')) + ' \n',
    'pw_dos'         : cal_order + "  " + "projwfc.x" + mpinpool + " < dos.in > dos.log 2>&1 \n",
    'dos_plot'       : str(python_path) + ' -u ' + str(osp.join(scriptpath, 'scripts', 'dos.py')) + ' \n',
    'phonon'         : cal_order + "  " + "ph.x" + mpinpool + " < phonon.in > ph.log 2>&1 \n",
    'phonon_TEST'    : cal_order + "  " + osp.join(ph_xTEST, 'test-suite','not_epw_comp', 'ph.x') + mpinpool + " < phononEPW.in > phEPW.log 2>&1 \n",
    'q2r'            : cal_order + "  " + "q2r.x" + mpinpool + " < q2r.in > q2r.log 2>&1 \n",
    'matdyn'         : cal_order + "  " + "matdyn.x" + mpinpool + " < matdyn.in > matdyn.log 2>&1 \n",
    'phondos'        : cal_order + "  " + "matdyn.x" + mpinpool + " < phdos.in > phdos.log 2>&1 \n",
    'phband_plot'    : str(python_path) + ' -u ' + str(osp.join(scriptpath, 'scripts', 'phonon_band.py')) + ' \n',
    'pw2wannier90'   : w90_post,
    'wannier90_pp'   : w90_pp,
    'wannier90'      : w90_distangle,
    'nihe_plot'      : str(python_path) + ' -u ' + str(osp.join(scriptpath, 'scripts', 'nihe.py')) + ' MLWFs \n',
    'epw'            : cal_order + " epw.x -npool " + str(inpp.num_process) + " < epw.in > epw.out 2>&1 \n",
    'phcollect_epw'  : str(python_path) + ' -u ' + str(osp.join(scriptpath, 'scripts', 'phCollect.py')) + ' \n',
    'last'           : 'date >> time \n',
    'ShengBTE'       : cal_order + " ShengBTE 2>BTE.err >BTE.out \n",
'perturbo_QE2pert'   : cal_order + " qe2pert.x " + mpinpool_must + "P -i QE2pert.in > QE2pert.out 2>&1 \n",
}

def start_submit(filepath = os.getcwd(), filename="submit_QE.sh"):
    subscp = open(osp.join(filepath, filename), 'w')
    for keys, vals in submitscripts.items():
        # json.dump(vals, subscp)
        subscp.write(vals)
        subscp.write("\n")
    subscp.write("export PATH=" + path_QE + ":$PATH \n")
    subscp.write("\n")
    subscp.write("\n")
    return subscp 

def wannier90_submit(subscp, filename='pwscf.win'):
    kpoints = inpp.kpt_nscf
    str_kpt = str(kpoints[0]) + "  " + str(kpoints[1]) + "  " + str(kpoints[2])
    # subscp.write('echo "mp_grid      =  '+ str_kpt +'" >>  \n')
    subscp.write('echo "mp_grid      =  '+ str_kpt +'" >> '+ filename +' \n')
    subscp.write('echo "begin kpoints" >> '+ filename +' \n')
    subscp.write('kmesh.pl ' + str_kpt +' wannier >> '+ filename +' \n')
    subscp.write('echo "end kpoints" >> '+ filename +' \n')
    # subscp.write(submit_order['wannier90_pp'])
    for i in range(len(submit_order['pw2wannier90'])):
        subscp.write(submit_order['wannier90_pp'][i])
        subscp.write(submit_order['pw2wannier90'][i])
        subscp.write(submit_order['wannier90'][i])
    subscp.write(submit_order['nihe_plot'])
    return subscp
def wannier90_submit0(subscp, filename='pwscf.win'):
    kpoints = inpp.kpt_nscf
    str_kpt = str(kpoints[0]) + "  " + str(kpoints[1]) + "  " + str(kpoints[2])
    # subscp.write('echo "mp_grid      =  '+ str_kpt +'" >>  \n')
    subscp.write('echo "mp_grid      =  '+ str_kpt +'" >> '+ filename +' \n')
    subscp.write('echo "begin kpoints" >> '+ filename +' \n')
    subscp.write('kmesh.pl ' + str_kpt +' wannier >> '+ filename +' \n')
    subscp.write('echo "end kpoints" >> '+ filename +' \n')
    # subscp.write(submit_order['wannier90'])
    # for i in range(len(submit_order['pw2wannier90'])):
    #     subscp.write(submit_order['pw2wannier90'][i])
    #     subscp.write(submit_order['wannier90'][i])
    # subscp.write(submit_order['nihe_plot'])
    return subscp

def band_submit(subscp):
    subscp.write(submit_order['pw_bands'])
    for ppb in range(len(band_post)):
        subscp.write(band_post[ppb])
    subscp.write(submit_order['pro_bands'] + '\n')
    subscp.write(submit_order['band_plot'])
    return subscp
def phonon_submit_single(subscp):
    subscp.write(submit_order['q2r'])
    subscp.write(submit_order['matdyn'])
    subscp.write(submit_order['phondos'])
    subscp.write(submit_order['phband_plot'])
    return subscp
def phonon_single(subscp, num):
    latt, element, num_element, tot_ele, cor = bas.get_poscar(os.getcwd(), "POSCAR")
    paras_phon = parameters.PHONON
    for keys, vals in inpp.phonon_input.items():
        paras_phon[keys] = vals
    for i in range(len(element)):
        paras_phon['amass('+ str(i+1) +')'] = parameters.atom_mass[element[i]]
    subscp.write('######### phonon ########## \n')
    subscp.write('for qp in $(seq 1 1 '+str(int(num))+') \n')
    subscp.write('do \n')
    subscp.write('if [ ! -f ./tmp/phonon_${qp} ] \n')
    subscp.write('then \n')
    subscp.write('cat > ./tmp/phonon_${qp} << EOF \n')
    subscp.write(str('phonons of systems') + "\n")    
    subscp.write("&inputph \n")  
    del paras_phon['start_q']    
    paras_phon['outdir'] = './tmp'    
    subscp = write_checkd0(subscp, paras_phon)
    subscp.write("start_q".ljust(20, ' ') + "  = ${qp}, \n")
    subscp.write("last_q".ljust(20, ' ') + "  = ${qp}, \n")
    subscp.write("/" + "\n")
    subscp.write("EOF")
    subscp.write("\n")
    subscp.write(cal_order + "  " + "ph.x" + mpinpool + " < ./tmp/phonon_${qp} > ./tmp/ph_${qp}.log 2>&1 \n")
    subscp.write('fi \n')
    subscp.write("done \n")
    subscp.write('############################ \n')
    return subscp
def phonon_single_pert(subscp, num):
    latt, element, num_element, tot_ele, cor = bas.get_poscar(os.getcwd(), "POSCAR")
    paras_phon = parameters.PHONON
    for keys, vals in inpp.phonon_input.items():
        paras_phon[keys] = vals
    for i in range(len(element)):
        paras_phon['amass('+ str(i+1) +')'] = parameters.atom_mass[element[i]]
    subscp.write('######### phonon ########## \n')
    subscp.write('for qp in $(seq 1 1 '+str(int(num))+') \n')
    subscp.write('do \n')
    subscp.write('if [ ! -f ./tmp/phonon_${qp} ] \n')
    subscp.write('then \n')
    subscp.write('cat > ./tmp/phonon_${qp} << EOF \n')
    subscp.write(str('phonons of systems') + "\n")    
    subscp.write("&inputph \n")  
    del paras_phon['start_q']    
    paras_phon['outdir'] = './tmp'    
    paras_phon['fildyn'] = "pwscf.dyn.xml"  
    paras_phon['fildvscf'] = "dvscf"
    subscp = write_checkd0(subscp, paras_phon)
    subscp.write("start_q".ljust(20, ' ') + "  = ${qp}, \n")
    subscp.write("last_q".ljust(20, ' ') + "  = ${qp}, \n")
    subscp.write("/" + "\n")
    subscp.write("EOF")
    subscp.write("\n")
    subscp.write(cal_order + "  " + "ph.x" + mpinpool + " < ./tmp/phonon_${qp} > ./tmp/ph_${qp}.log 2>&1 \n")
    subscp.write('fi \n')
    subscp.write("done \n")
    subscp.write('############################ \n')
    return subscp
def phonon_submit(subscp):
    try:
        if inpp.phonon_input['start_q'] == 0:
            num_q = int(input("Please input number of q points in phonon calculations: "))
            subscp = phonon_single(subscp, num_q)
        else:
            subscp.write(submit_order['phonon'])
            subscp.write(submit_order['q2r'])
            subscp.write(submit_order['matdyn'])
            subscp.write(submit_order['phondos'])
            subscp.write(submit_order['phband_plot'])
    except:
        subscp.write(submit_order['phonon'])
        subscp.write(submit_order['q2r'])
        subscp.write(submit_order['matdyn'])
        subscp.write(submit_order['phondos'])
        subscp.write(submit_order['phband_plot'])
    return subscp
def phonon_submitTEST(subscp):
    subscp.write(submit_order['phonon'])
    subscp.write(submit_order['phonon_TEST'])
    subscp.write(submit_order['q2r'])
    subscp.write(submit_order['matdyn'])
    # subscp.write(submit_order['phonon_TEST'])
    # subscp.write(submit_order['phondos'])
    subscp.write(submit_order['phband_plot'])
    print('PH_PATH:', ph_xTEST)
    return subscp
def phonon_submit_perturbo_single(subscp):
    subscp.write('cp pwscf.dyn0 pwscf.dyn0.xml \n')
    subscp.write(submit_order['q2r'])
    subscp.write(submit_order['matdyn'])
    subscp.write(submit_order['phondos'])
    subscp.write(submit_order['phband_plot'])
    return subscp
def phonon_submit_perturbo(subscp):
    try:
        if inpp.phonon_input['start_q'] == 0:  
            num_q = int(input("Please input number of q points in phonon calculations: ")) 
            subscp = phonon_single_pert(subscp, num_q)
        else:
            subscp.write(submit_order['phonon'])
            subscp.write('cp pwscf.dyn0 pwscf.dyn0.xml \n')
            subscp.write(submit_order['q2r'])
            subscp.write(submit_order['matdyn'])
            subscp.write(submit_order['phondos'])
            subscp.write(submit_order['phband_plot'])
    except:
        subscp.write(submit_order['phonon'])
        subscp.write('cp pwscf.dyn0 pwscf.dyn0.xml \n')
        subscp.write(submit_order['q2r'])
        subscp.write(submit_order['matdyn'])
        subscp.write(submit_order['phondos'])
        subscp.write(submit_order['phband_plot'])
    return subscp
def phonopy_2nd():
    filepath = os.getcwd()
    subscp = open(osp.join(filepath, 'submit_QE.sh'), 'w')
    for keys, vals in submitscripts.items():
        subscp.write(vals)
        subscp.write("\n")
    subscp.write("export PATH=" + path_QE + ":$PATH \n \n \n")
    subscp.write("random_integer=$((RANDOM % 10000))\n")
    subscp.write('random_number=$(echo "scale=4; $random_integer / 1000.0" | bc) \n')
    subscp.write('sleep $random_number \n \n')
    subscp.write("for ff in `ls supercell-*.in` \n")
    subscp.write('do \n')
    subscp.write('  nums=`echo $ff | tr -cd "[0-9]"` \n')
    subscp.write('  echo $nums test.log \n')
    subscp.write('  if [ ! -e job-$nums ]; then \n')
    subscp.write('    mkdir job-$nums \n')
    subscp.write('    cat 2nd.in $ff >> ./job-$nums/$ff  \n')
    subscp.write('    rm $ff  \n')
    subscp.write('    cd  job-$nums  \n ')
    subscp.write("   " + cal_order + " pw.x " + mpinpool + " < $ff > supercell.log \n")
    subscp.write('    cd  ../  \n')
    subscp.write('  fi \n')
    subscp.write('done  \n')
    subscp.write('phonopy -f job-*/supercell.log \n')
    subscp.write('phonopy --qe -c scf.in -p KPATH.phonopy \n')
    return subscp

def shengbte_3rd():
    filepath = os.getcwd()
    subscp = open(osp.join(filepath, 'submit_QE.sh'), 'w')
    for keys, vals in submitscripts.items():
        subscp.write(vals)
        subscp.write("\n")
    subscp.write("export PATH=" + path_QE + ":$PATH \n \n \n")
    subscp.write("random_integer=$((RANDOM % 10000))\n")
    subscp.write('random_number=$(echo "scale=4; $random_integer / 1000.0" | bc) \n')
    subscp.write('sleep $random_number \n \n')
    subscp.write("for ff in `ls -a | grep 'DISP.scf_sc.in.*'` \n")
    subscp.write('do \n')
    subscp.write('nums=`echo $ff | tr -cd "[0-9]"` \n')
    subscp.write('if [ ! -e job-$nums ]; then \n')
    subscp.write('mkdir job-$nums \n')
    subscp.write('mv $ff job-$nums  \n')
    subscp.write('cd  job-$nums  \n \n')
    subscp.write(cal_order + " pw.x " + mpinpool + " < $ff > $ff.log \n")
    subscp.write('cd  ../  \n')
    subscp.write('fi \n')
    subscp.write('done  \n')
    subscp.write( inpp.third_order_collect + "\n") 
    return subscp
def shengbte_4th():
    filepath = os.getcwd()
    subscp = open(osp.join(filepath, 'submit_QE.sh'), 'w')
    for keys, vals in submitscripts.items():
        subscp.write(vals)
        subscp.write("\n")
    subscp.write("export PATH=" + path_QE + ":$PATH \n \n \n")
    subscp.write("random_integer=$((RANDOM % 10000))\n")
    subscp.write('random_number=$(echo "scale=4; $random_integer / 1000.0" | bc) \n')
    subscp.write('sleep $random_number \n \n')
    subscp.write("for ff in `ls -a | grep 'DISP.scf_sc.in.*'` \n")
    subscp.write('do \n')
    subscp.write('nums=`echo $ff | tr -cd "[0-9]"` \n')
    subscp.write('if [ ! -e job-$nums ]; then \n')
    subscp.write('mkdir job-$nums \n')
    subscp.write('mv $ff job-$nums  \n')
    subscp.write('cd  job-$nums  \n \n')
    subscp.write(cal_order + " pw.x " + mpinpool + " < $ff > $ff.log \n")
    subscp.write('cd  ../  \n')
    subscp.write('fi \n')
    subscp.write('done  \n')
    subscp.write( inpp.fourth_order_collect + "\n") 
    return subscp
def end_submit(subscp):
    subscp.write("\n")
    subscp.write("\n")
    subscp.write(submit_order['last'])
    subscp.write("\n")
    subscp.close()
    os.system('chmod +x ' + osp.join(os.getcwd(), 'submit_QE.sh'))
    os.system('chmod +x ' + osp.join(os.getcwd(), 'AutoinputQE.py'))
    subscp.close()
def end_submit_jb(subscp, filename):
    subscp.write("\n")
    subscp.write("\n")
    subscp.write(submit_order['last'])
    subscp.write("\n")
    subscp.close()
    os.system('chmod +x ' + osp.join(os.getcwd(), filename))
    os.system('chmod +x ' + osp.join(os.getcwd(), 'AutoinputQE.py'))
    subscp.close()

###################### funcs
def job_epw_single(filepath=os.getcwd(), filename='epw.in'):
    subscp = start_submit(filepath)
    filename = filename.split('.')[0]
    subscp.write(cal_order + " epw.x -npool " + str(inpp.num_process) + "  < " + filename + ".in > " + filename + ".out 2>&1 \n")
    end_submit(subscp)

class job_submit():
    def __init__(self, filepath=os.getcwd()):
        self.filepath = filepath
    @staticmethod
    def job_notAuto(type):
        subscp = start_submit()
        if type == 'wannier90':
            wannier90_submit(subscp)
        else:
            subscp.write(submit_order['pw_' + type])
        end_submit(subscp)
    @staticmethod
    def job_AutoW90():
        subscp = start_submit()
        subscp.write(submit_order['pw_scf'])
        subscp = band_submit(subscp)
        subscp.write(submit_order['pw_nscf'])
        if mag == 2:
            subscp = wannier90_submit0(subscp, 'pwscf_up.win')
            subscp.write('\n')
            subscp = wannier90_submit(subscp, 'pwscf_dw.win')
        else:
            subscp = wannier90_submit(subscp, 'pwscf.win')
        end_submit(subscp)
    @staticmethod
    def job_distangle():
        subscp = start_submit()
        subscp.write("\n")
        subscp.write("\n")
        subscp.write("rm *.werr \n ")
        subscp.write("rm *.png \n ")
        subscp.write("rm *.err \n ")
        subscp.write("rm *.out \n ")
        subscp.write("ln -s ../pwscf*.amn ./ \n")
        subscp.write("ln -s ../pwscf*.mmn ./ \n")
        subscp.write("ln -s ../pwscf*.eig ./ \n")
        subscp.write("ln -s ../band*.dat* ./ \n")
        subscp.write("ln -s ../scf.log ./ \n \n")
        for i in range(len(submit_order['pw2wannier90'])):
            subscp.write(submit_order['wannier90'][i])
        subscp.write(submit_order['nihe_plot'])
        end_submit(subscp)    
    @staticmethod
    def job_AutoPhon():
        try:
            if inpp.phonon_input['start_q'] == 0:
                subscp = start_submit(os.getcwd(), "sub1.sh")
                # subscp.write(submit_order['pw_scf'])
                # subscp.write(submit_order['pw_nscf'])
                subscp = phonon_submit(subscp)
                end_submit_jb(subscp, 'sub1.sh')   
                # last 
                subscp2 = start_submit(os.getcwd(), "last.sh")
                subscp2 = phonon_submit_single(subscp2)
                end_submit_jb(subscp2, "last.sh") 
                os.system('mv submit_QE.sh scf.sh')
            else:
                subscp = start_submit()
                subscp.write(submit_order['pw_scf'])
                # subscp.write(submit_order['pw_nscf'])
                subscp = phonon_submit(subscp)
                end_submit(subscp)
        except:
            subscp = start_submit()
            subscp.write(submit_order['pw_scf'])
            # subscp.write(submit_order['pw_nscf'])
            subscp = phonon_submit(subscp)
            end_submit(subscp)
    @staticmethod
    def job_AutoPhonTest():
        subscp = start_submit()
        subscp.write(submit_order['pw_scf'])
        # subscp.write(submit_order['pw_nscf'])
        subscp = phonon_submitTEST(subscp)
        end_submit(subscp)
    @staticmethod
    def job_AutoPhon_perturbo():
        try:
            if inpp.phonon_input['start_q'] == 0:
                subscp = start_submit(os.getcwd(), "sub1.sh")
                # subscp.write(submit_order['pw_scf'])
                # subscp.write(submit_order['pw_nscf'])
                subscp = phonon_submit_perturbo(subscp)
                end_submit_jb(subscp, 'sub1.sh')    
                # last 
                subscp2 = start_submit(os.getcwd(), "last.sh")
                subscp2 = phonon_submit_perturbo_single(subscp2)
                end_submit_jb(subscp2, 'last.sh') 
                os.system('mv submit_QE.sh scf.sh')
            else:
                subscp = start_submit()
                subscp.write(submit_order['pw_scf'])
                # subscp.write(submit_order['pw_nscf'])
                subscp = phonon_submit_perturbo(subscp)
                end_submit(subscp)
        except:
            subscp = start_submit()
            subscp.write(submit_order['pw_scf'])
            # subscp.write(submit_order['pw_nscf'])
            subscp = phonon_submit_perturbo(subscp)
            end_submit(subscp)
    @staticmethod
    def job_AutoBand():
        subscp = start_submit()
        subscp.write(submit_order['pw_scf'])
        subscp = band_submit(subscp)
        end_submit(subscp)
    @staticmethod
    def job_AutoNscf():
        subscp = start_submit()
        subscp.write(submit_order['pw_scf'])
        subscp.write(submit_order['pw_nscf'])
        end_submit(subscp)
    @staticmethod
    def job_AutoBandsNscf():
        subscp = start_submit()
        subscp.write(submit_order['pw_scf'])
        subscp = band_submit(subscp)
        subscp.write(submit_order['pw_nscf'])
        end_submit(subscp)
    @staticmethod
    def job_AutoDos():
        subscp = start_submit()
        subscp.write(submit_order['pw_scf'])
        subscp.write(submit_order['pw_nscf'])
        subscp.write(submit_order['pw_dos'])
        subscp.write(submit_order['dos_plot'])
        end_submit(subscp)
    @staticmethod
    def job_epw(filepath=os.getcwd(), filename='epw.in'):
        subscp = start_submit(filepath)
        filename = filename.split('.')[0]
        subscp.write(cal_order + " epw.x -npool " + str(inpp.num_process) + "  < " + filename + ".in > " + filename + ".out 2>&1 \n")
        end_submit(subscp)
    @staticmethod
    def job_AutoEpw(filepath=os.getcwd(), filenames=[]):
        subscp = start_submit(filepath)
        for filename in filenames:
            filename1 = filename.split('.')[0]
            subscp.write(cal_order + " epw.x -npool " + str(inpp.num_process) + "  < " + filename1 + ".in > " + filename1 + ".out 2>&1 \n")
        end_submit(subscp)

    @staticmethod
    def job_AutoPhonEpw():
        subscp = start_submit()
        subscp.write(submit_order['pw_scf'])
        subscp.write(submit_order['pw_nscf'])
        subscp = phonon_submit(subscp)
        subscp.write(submit_order['phcollect_epw'])
        end_submit(subscp)
    @staticmethod
    def job_ShengBTE():
        subscp = start_submit()
        subscp.write("export PATH=" + configg.shengbte_path + ":$PATH \n")
        subscp.write(submit_order['ShengBTE'])
        end_submit(subscp)
    @staticmethod
    def job_2nd():
        subscp = phonopy_2nd()
        end_submit(subscp)
    @staticmethod
    def job_3rd():
        subscp = shengbte_3rd()
        end_submit(subscp)
    @staticmethod
    def job_4th():
        subscp = shengbte_4th()
        end_submit(subscp)
    @staticmethod
    def job_qe2pert(filepath, filename):
        subscp = start_submit(filepath)
        subscp.write("export PATH=" + configg.perturbo_path + ":$PATH \n")
        subscp.write("\n")
        subscp.write("\n")
        subscp.write("mkdir tmp \n")
        subscp.write("cd tmp \n")
        subscp.write("ln -sf ../../nscf/tmp/pwscf.save . \n")
        subscp.write("cd ../ \n")
        subscp.write("ln -sf ../w90/pwscf_u.mat \n")
        subscp.write("ln -sf ../w90/pwscf_u_dis.mat \n")
        subscp.write("ln -sf ../w90/pwscf_centres.xyz \n")
        subscp.write(cal_order + " qe2pert.x " + mpinpool_must + " -i "+ filename + " > " + filename.split('.')[0] + ".out 2>&1 \n")
        end_submit(subscp)
    @staticmethod
    def job_pert(filepath, filename):
        subscp = start_submit(filepath)
        subscp.write("export PATH=" + configg.perturbo_path + ":$PATH \n")
        if not osp.exists(osp.join(os.getcwd(), "../../qe2pert/pwscf_epwan.h5")):
            subscp.write("ln -sf ../../qe2pert/pwscf_epr.h5 . \n")
        else:
            subscp.write("ln -sf ../../qe2pert/pwscf_epwan.h5 . \n")
        # subscp.write('export OMP_NUM_THREADS=4 \n')
        subscp.write(cal_order + " perturbo.x " + mpinpool_must + " -i "+ filename + " > " + filename.split('.')[0] + ".out 2>&1 \n")
        end_submit(subscp)
    @staticmethod
    def job_pert_ephmat(filepath, filename):
        subscp = start_submit(filepath)
        subscp.write("export PATH=" + configg.perturbo_path + ":$PATH \n")
        subscp.write('export OMP_NUM_THREADS=' +str(int(NP/NP_pert))+ ' \n')
        if not osp.exists(osp.join(os.getcwd(), "../../qe2pert/pwscf_epwan.h5")):
            subscp.write("ln -sf ../../qe2pert/pwscf_epr.h5 . \n")
        else:
            subscp.write("ln -sf ../../qe2pert/pwscf_epwan.h5 . \n")
        subscp.write(cal_order_ephmat + str(NP_pert) + " perturbo.x -npools " + str(NP_pert) + " -i "+ filename + " > " + filename.split('.')[0] + ".out 2>&1 \n")
        end_submit(subscp)
    @staticmethod
    def job_pert_trans(filepath, filename):
        subscp = start_submit(filepath)
        subscp.write("export PATH=" + configg.perturbo_path + ":$PATH \n")
        subscp.write('export OMP_NUM_THREADS=' +str(int(NP/NP_pert))+ ' \n')
        if not osp.exists(osp.join(os.getcwd(), "../../qe2pert/pwscf_epwan.h5")):
            subscp.write("ln -sf ../../qe2pert/pwscf_epr.h5 . \n")
        else:
            subscp.write("ln -sf ../../qe2pert/pwscf_epwan.h5 . \n")
        subscp.write(cal_order_ephmat + str(NP_pert) +  " perturbo.x -npools " + str(NP_pert) + " -i "+ filename + " > " + filename.split('.')[0] + ".out 2>&1 \n")
        end_submit(subscp)


