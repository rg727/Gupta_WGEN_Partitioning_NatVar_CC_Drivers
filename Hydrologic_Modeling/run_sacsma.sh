for T in 0 
do
SLURM="#!/bin/bash\n\

#SBATCH -t 100:00:00\n\
#SBATCH --job-name=matlab_test\n\
#SBATCH -p normal\n\
#SBATCH --export=ALL\n\
#SBATCH --nodes=1\n\       
#SBATCH --exclusive\n\          
#SBATCH --output=matlab_test.txt\n\
#SBATCH --error=matlab_error.txt\n\
#SBATCH --ntasks-per-node=16\n\

module load matlab/R2021a\n\

mkdir /home/fs02/pmr82_0001/rg727/check_parallel_matlab/parallel_job_info_${T}\n\

matlab -nodesktop -nosplash -r \"T=${T}; C=0;\" < run_sacsma.m\n\

rm -rf /home/fs02/pmr82_0001/rg727/check_parallel_matlab/parallel_job_info_${T}"

echo -e $SLURM | sbatch 
sleep 0.5
done
  
