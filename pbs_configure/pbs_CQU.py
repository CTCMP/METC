import os 
import AutoinputQE as inpp
import os.path as osp
############ intel environment ##################
intel = {
'i0'       : "module load intel/2022.2.0",
'i1'       : "module load hdf5/1.14.3",
# 'i2'       : "module load mpi/intel/2021.1",
# 'i3'       : 'module load mpi/intel/17.0.7-thc',
}
#################################################
# cal and submit order
cal_order = "mpirun -n " + str(inpp.num_process)
submit_order = 'sbatch '

############### software path ###################
#QE_path = '/work/wangr/dxy/software/QE/epw-qe6.5/bin'
#QE_path = '/work/wangr/dxy/software/QE/perturbo-QE7.0/bin'
#QE_path = '/work/wangr/dxy/software/QE/QE-6.8-EPW-Patch/qe-6.8/bin'
#QE_path = '/work/wangr/dxy/software/QE/epw-QE7.2-use/bin'
#QE_path = '/work/wangr/dxy/software/QE/QE7.2-perturbo/q-e-qe-7.2/bin'
#QE_path = '/work/wangr/dxy/software/QE/QE7.2_modified_perturbo2.2/qe-7.2/bin'
#QE_path = '/work/wangr/dxy/software/QE/QE7.2-all/qe-7.2/bin'
#QE_path = '/work/wangr/dxy/software/QE/QE7.2-perturbo-TDEP/q-e-qe-7.2/bin'
#QE_path = '/share/home/yangxiaolong/software/qe-7.2/bin'
QE_path = '/share/home/yangxiaolong/dxy/soft/qe-7.2/bin'

python_path = '/share/home/yangxiaolong/software/anaconda3/bin/python'
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
#perturbo_path = '/work/wangr/dxy/software/QE/QE7.2-all/qe-7.2/perturbo-2.2.0/bin'
#perturbo_path = '/share/home/yangxiaolong/software/qe-7.2/perturbo-2.2.0/bin'
perturbo_path = '/share/home/yangxiaolong/dxy/soft/qe-7.2/perturbo-2.2.0/bin'

################ submit scripts configure #######################################
submit_scripts = {
'1'              : '#!/bin/bash',
'2'              : '#SBATCH --job-name='+"/".join(str(os.getcwd()).split("/")[-2:]),
'3'              : '#SBATCH --partition=' + inpp.cal_nodes,           ## sr850 sr860-768 sr860-1536 sr650
'4'              : '#SBATCH --nodes=' + str(inpp.num_nodes),
'5'              : '#SBATCH --ntasks-per-node=72',
'6'              : '#SBATCH --output=%j.out',
'7'              : '#SBATCH --error=%j.err',
'8'              : "",
'9'              : "",
'10'             : "export OMPI_MCA_btl='^uct,ofi'",
'11'             : "export OMPI_MCA_pml='ucx'",
'12'             : "export OMPI_MCA_mtl='^ofi'",
'13'             : "",
'14'             : "",
'15'             : 'echo -n "start time  " > time',
'last'           : 'date >> time',
}
for keys, vals in intel.items():
    submit_scripts[keys] = intel[keys]
###################################################################################

