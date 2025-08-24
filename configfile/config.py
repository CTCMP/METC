import os 
import AutoinputQE as inpp
import os.path as osp
############ intel environment ##################
intel = {
'i0'       : "source /share/intel/2022.2/setvars.sh",
#'i1'       : "module load cmake/3.23.1",
'i2'       : "module load gcc/10.2.0",
#'i3'       : 'cd $PBS_O_WORKDIR',
}

#################################################
# cal and submit order
cal_order = "mpirun -np " + str(int(inpp.num_process*inpp.num_nodes))
submit_order = 'qsub  '

############### software path ###################
QE_path = '/share/home/jinxin/soft/QE/q-e-qe-7.2/bin'
#QE_path = '/share/home/jinxin/soft/QE/qe-7.0/bin'

python_path = '/share/home/jinxin/dxy/soft/miniconda39/bin/python'

path_presudo = '/share/home/jinxin/soft/QE/POT'
path_presudo_My_defined = '/share/home/jinxin/soft/QE/POT/QE-pot-master' # you can download it in : git@github.com:XVDing/QE-pot.git

############## other modules ###################
third_order_path = "/share/home/jinxin/soft/shengbte/thirdorder"
######## kp + kc
#shengbte_path = "/share/home/jinxin/soft/shengbte/diff_four_phonon"
shengbte_path = "/share/home/jinxin/soft/shengbte/FourPhonon_ruan/FourPhonon-main"

#shengbte_path = "/share/home/jinxin/soft/shengbte/shengbtev3v4/ShengBTE"
shengbte_path = "/share/home/jinxin/dxy/soft/FourPhonon-main"
perturbo_path = "/share/home/jinxin/soft/QE/q-e-qe-7.2/perturbo-2.2.0/bin"
#perturbo_path = '/share/home/jinxin/soft/QE/qe-7.0/perturbo-2.0.2/bin'
################ submit scripts configure #######################################
submit_scripts = {
'1'              : '#!/bin/sh',
'2'              : '#PBS -N ' + os.getcwd().split('/')[-1],
'3'              : '#PBS -j oe',           ## sr850 sr860-768 sr860-1536 sr650
'4'              : '#PBS -q workq',
'5'              : '#PBS -l select=1:ncpus='+str(inpp.num_process)+':host='+inpp.cal_nodes+':mpiprocs='+str(inpp.num_process),
'8'              : '#PBS -V',
'i5'             : 'cd $PBS_O_WORKDIR',
'i6'             : '',
'17'             : 'echo -n "start time  " > time',
'18'             : 'export I_MPI_HYDRA_BOOTSTRAP=ssh',
'last'           : 'date >> time',
}
for keys, vals in intel.items():
    submit_scripts[keys] = intel[keys]
###################################################################################

