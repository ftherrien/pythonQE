#!/bin/bash -x
#SBATCH --account=171206131221
#SBATCH --job-name='name'
#SBATCH --nodes={n_nodes}
#SBATCH --ntasks={n_proc}
##SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time={time}
#SBATCH -o output.txt
#SBATCH -e errlog.txt