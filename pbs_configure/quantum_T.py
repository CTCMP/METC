import os 
import AutoinputQE as inpp
import os.path as osp
############ intel environment ##################
intel = {
# 'i0'       : "source /public1/soft/modules/module.sh",
# 'i1'       : "source /work/software/intel/bin/compilervars.sh -arch intel64 -platform linux",
'i2'       : "module load intel-2018",
# 'i3'       : 'module load mpi/intel/17.0.7-thc',
}
def conf_intel():
    for keys, vals in intel.items():
        os.system(vals)
#################################################
# cal and submit order
cal_order = "mpirun -np " + str(int(inpp.num_process*inpp.num_nodes))
submit_order = 'sbatch  '

############### software path ###################
QE_path = '/home/xlyang/jx/soft/QE/perturbo/q-e-qe-7.0/bin'
#QE_path = '/home/xlyang/jx/soft/QE/qe-7.2/bin'
#QE_path = '/home/xlyang/jx/soft/QE/QE-EPW-Modified/qe-7.2/bin'
#QE_path = '/home/xlyang/jx/soft/QE/qe-7.2-DFPT_EPW/qe-7.2/bin'

python_path = '~/jx/soft/miniconda38/bin/python'
path_presudo = '/home/xlyang/jx/soft/QE/pot'
path_presudo_My_defined = '/home/xlyang/jx/soft/QE/pot' # you can download it in : git@github.com:XVDing/QE-pot.git

############## other modules ###################
third_order_path = "/home/xlyang/jx/soft/shengbte/thirdorder"
# shengbte_path = "/work/wangr/dxy/software/shengbte/shengbte_epw/ShengBTE"
shengbte_path = "//home/xlyang/jx/soft/shengbte/shengbtev3v4/ShengBTE"
perturbo_path = "/home/xlyang/jx/soft/QE/perturbo/q-e-qe-7.0/perturbo-2.0.2/bin"

################ submit scripts configure #######################################
submit_scripts = {
'1'              : '#!/bin/bash',
'a'              : '#SBATCH --partition='+str(inpp.cal_nodes),
'2'              : '#SBATCH --job-name='+os.getcwd(),
'3'              : '#SBATCH --nodes=' + str(inpp.num_nodes),           ## sr850 sr860-768 sr860-1536 sr650
'4'              : '#SBATCH --ntasks-per-node=' + str(inpp.num_process),
'5'              : '#SBATCH --error=%j.err',
'6'              : '#SBATCH --output=%j.out',
'8'              : 'CURDIR=`pwd`',
'9'              : 'rm -rf $CURDIR/nodelist.$SLURM_JOB_ID',
'10'             : 'NODES=`scontrol show hostnames $SLURM_JOB_NODELIST`',
'11'             : 'for i in $NODES',
'12'             : 'do',
'13'             : 'echo "$i:$SLURM_NTASKS_PER_NODE" >> $CURDIR/nodelist.$SLURM_JOB_ID',
'14'             : 'done',
'15'             : 'echo $SLURM_NPROCS',
'16'             : 'echo "process will start at : "',
# '12'             : 'source /work/software/intel/bin/compilervars.sh -arch intel64 -platform linux',
'17'             : 'echo -n "start time  " > time',
'last'           : 'date >> time',
}
for keys, vals in intel.items():
    submit_scripts[keys] = intel[keys]
###################################################################################

