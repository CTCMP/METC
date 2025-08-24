import os 
import AutoinputQE as inpp
import os.path as osp
############ intel environment ##################
intel = {
# 'i0'       : "source /public1/soft/modules/module.sh",
'i1'       : "source /work/software/intel/bin/compilervars.sh -arch intel64 -platform linux",
# 'i2'       : "module load mpi/intel/2021.1",
# 'i3'       : 'module load mpi/intel/17.0.7-thc',
}
#################################################
# cal and submit order
cal_order = "mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np " + str(inpp.num_process)
submit_order = 'bsub < '

############### software path ###################
#QE_path = '/work/wangr/dxy/software/QE/epw-qe6.5/bin'
#QE_path = '/work/wangr/dxy/software/QE/perturbo-QE7.0/bin'
#QE_path = '/work/wangr/dxy/software/QE/QE-6.8-EPW-Patch/qe-6.8/bin'
# QE_path = '/work/wangr/dxy/software/QE/epw-QE7.2-use/bin'
#QE_path = '/work/wangr/dxy/software/QE/QE7.2-perturbo/q-e-qe-7.2/bin'
#QE_path = '/work/wangr/dxy/software/QE/QE7.2_modified_perturbo2.2/qe-7.2/bin'
QE_path = '/work/wangr/dxy/software/QE/QE7.2-all/qe-7.2/bin'

python_path = '~/data/dxy/soft/python/miniconda/bin/python'
path_presudo = '/work/wangr/data/dxy/soft/QE/pot'
path_presudo_My_defined = '/work/wangr/dxy/scripts/QEkit/pseudo' # you can download it in : git@github.com:XVDing/QE-pot.git

############## other modules ###################
third_order_path = "/work/wangr/dxy/software/shengbte/third_order/thirdorder"
fourthorder_path = "/work/wangr/jx/soft/fourphonon/Fourthorder-main/"
python2_path = "/work/wangr/tool/anaconda2/anaconda2/bin/python"
# shengbte_path = "/work/wangr/dxy/software/shengbte/shengbte_epw/ShengBTE"
#shengbte_path = "/work/wangr/dxy/software/shengbte/epw/ShengBTE"
shengbte_path = "/work/wangr/dxy/software/shengbte/FourPhonon_ruan/FourPhonon-1.1"
#shengbte_path = "/work/wangr/dxy/software/shengbte/shengbtev3v4/ShengBTE"
#shengbte_path = "/work/wangr/dxy/software/shengbte/shengbte_fourphonon_iteration"
#perturbo_path = "/work/wangr/dxy/software/QE/perturbo-QE7.0/perturbo-2.0.2/bin"
#perturbo_path = "/work/wangr/dxy/software/QE/QE7.2-perturbo/q-e-qe-7.2/perturbo-2.2.0/bin"
#perturbo_path = '/work/wangr/dxy/software/QE/QE7.2_modified_perturbo2.2/qe-7.2/perturbo-2.2.0/bin'
perturbo_path = '/work/wangr/dxy/software/QE/QE7.2-all/qe-7.2/perturbo-2.2.0/bin'

################ submit scripts configure #######################################
submit_scripts = {
'1'              : '#!/bin/bash',
'2'              : '#BSUB -q ' + inpp.cal_nodes,           ## sr850 sr860-768 sr860-1536 sr650
'3'              : '#BSUB -n ' + str(inpp.num_process),
'4'              : '#BSUB -e %J.err',
'5'              : '#BSUB -o %J.out',
'6'              : '#BSUB -J ' + os.getcwd(),
'7'              : '#BSUB -R "span[ptile=' + str(inpp.num_process) + ']"',
'8'              : '#BSUB -R "select[hname!=\'r13n18\']"',
'9'              : 'hostfile=`echo $LSB_DJOB_HOSTFILE`',
'10'             : 'NP=`cat $hostfile | wc -l`',
'11'             : 'cd $LS_SUBCWD',
# '12'             : 'source /work/software/intel/bin/compilervars.sh -arch intel64 -platform linux',
'13'             : 'echo -n "start time  " > time',
'last'           : 'date >> time',
}
for keys, vals in intel.items():
    submit_scripts[keys] = intel[keys]
###################################################################################

