import AutoinputQE as inpp
# ======================================================================================
# python path 
python_path = '/share/home/jinxin/dxy/soft/miniconda39/bin/python'

# ======================================================================================
# cal and submit order
cal_order = "mpirun -np " + str(int(inpp.num_process*inpp.num_nodes))
submit_order = 'qsub  '

# ======================================================================================
# software path 
QE_path = '/share/home/jinxin/soft/QE/q-e-qe-7.2/bin'

# ======================================================================================
# persudo path: you can download it in git@github.com:XVDing/QE-pot.git
persudo_path = '/share/home/jinxin/soft/QE/POT'

# ======================================================================================
# other modules 
third_order_path = "/share/home/jinxin/soft/shengbte/thirdorder"
fourth_order_path = "/share/home/jinxin/soft/shengbte/thirdorder"
shengbte_path = "/share/home/jinxin/dxy/soft/FourPhonon-main"
perturbo_path = "/share/home/jinxin/soft/QE/q-e-qe-7.2/perturbo-2.2.0/bin"

# ======================================================================================
# submit scripts configure
# inpp.num_nodes * inpp.num_process means total cpu cores
submit_script = f'''
#!/bin/sh
#PBS -N calc
#PBS -j oe
#PBS -q workq
#PBS -l select=1:ncpus={inpp.num_nodes * inpp.num_process}:host={inpp.node}:mpiprocs={inpp.num_nodes * inpp.num_process}
#PBS -V
cd $PBS_O_WORKDIR

echo -n "start time  " > time
export I_MPI_HYDRA_BOOTSTRAP=ssh
date >> time
source /share/intel/2022.2/setvars.sh
module load gcc/10.2.0

'''
###################################################################################

