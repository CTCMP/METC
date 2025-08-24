import os 
import AutoinputQE as inpp
import os.path as osp
############ intel environment ##################
intel = {
'i0'       : "source /public1/soft/modules/module.sh",
# 'i1'       : "source /work/software/intel/bin/compilervars.sh -arch intel64 -platform linux",
# 'i2'       : "module load mpi/intel/2021.1",
# 'i3'       : 'module load mpi/intel/17.0.7-thc',
'i3'       : 'module load mpi/intel/18.0.2-thc',
# 'i4'       : 'module load mpi/oneAPI/2022.1',
}
#################################################
# cal and submit order
cal_order = "mpirun -np " + str(inpp.num_process)
submit_order = 'sbatch '

############### software path ###################
# path_presudo = "/work/wangr/data/dxy/soft/perturbo/pseudopotential"
# QE_path = '/work/wangr/data/dxy/soft/perturbo/QE/qe-7.0/bin'
# python_path = '~/data/dxy/soft/python/miniconda/bin/python'


# QE_path = '/public1/soft/qe/6.7/bin'
# QE_path = '/public1/home/scb1717/work/xiang-CeAlSi/dxy/soft/QE/QE/q-e-qe-7.2/bin'
QE_path = '/public1/home/scb1717/work/xiang-CeAlSi/dxy/soft/QE/QE/qe-6.8/bin'

python_path = '~/work/xiang-CeAlSi/dxy/soft/minconda/minconda/bin/python'
path_presudo = '/public1/home/scb1717/work/xiang-CeAlSi/dxy/soft/QE/pot'
path_presudo_My_defined = '/public1/home/scb1717/work/xiang-CeAlSi/dxy/soft/QE/QE-pot-master' # you can download it in : git@github.com:XVDing/QE-pot.git

############## other modules ###################
# third_order_path = "/public1/home/scb1717/work/xiang-CeAlSi/dxy/soft/ShengBTE/thirdorder"
# shengbte_path = "/public1/home/scb1717/work/xiang-CeAlSi/dxy/soft/ShengBTE/ShengBTE"
# perturbo_path = "/work/wangr/data/dxy/soft/perturbo/QE/qe-7.0/perturbo-2.0.2/bin"
third_order_path = '/public1/home/scb1717/work/xiang-CeAlSi/dxy/soft/ShengBTE/thirdorder'
shengbte_path = '/public1/home/scb1717/work/xiang-CeAlSi/dxy/soft/ShengBTE/ShengBTE'
perturbo_path = ""

################ submit scripts configure #######################################
submit_scripts = {
'1'              : '#!/bin/bash',
'2'              : '#SBATCH -p ' + inpp.cal_nodes,           ## amd_256, amd_512, amd_2T
'3'              : '#SBATCH -N ' + str(inpp.num_nodes),
# '4'              : '#SBATCH -J ' + os.getcwd(),
'5'              : '#SBATCH -n ' + str(inpp.num_process),
'6'              : '#SBATCH --mem-per-cpu=' + str('2048MB'),

'13'             : 'echo -n "start time  " > time',
'last'           : 'date >> time',
}

for keys, vals in intel.items():
    submit_scripts[keys] = intel[keys]
cal_order = "mpirun -np " + str(inpp.num_process)
#################################################################################
