#!/bin/bash
#SBATCH -J experiment_5     			# Job name
#SBATCH -o experiment_5.o%j 			# Name of stdout output file (%j expands to jobId)
#SBATCH -e experiment_5.o%j 			# Name of stderr output file(%j expands to jobId)
#SBATCH -n 10         					# Total number of mpi tasks requested
#SBATCH -p normal     					# Submit to the 'normal' or 'development' queue
#SBATCH -t 48:00:00 	  	   			# Run time (hh:mm:ss)
#SBATCH -A ParallelBNP   	            # Allocation name to charge job against
#SBATCH --mail-user=ochavez@utexas.edu	# specify email address for notifications
#SBATCH --mail-type=begin               # email when job begins
#SBATCH --mail-type=end                 # email when job ends

# load R module
module load Rstats           

# call R code from RMPISNOW
ibrun RMPISNOW < experiment_template_cmc_4gs_bcwts.R 
