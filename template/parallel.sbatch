#!/bin/bash
#SBATCH -t 70:00:00
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -e %j.err
#SBATCH -o %j.out
#SBATCH --priority=TOP
#SBATCH --mem=300gb

#burst
#SBATCH --account=chunlidai
#SBATCH --qos=chunlidai-b

# User variables
joblistfile='/home/chunli/chunliwork/work/landslide/testparallel/joblist'
          # Joblist file (a list of subdirectories seperated by lines)
#ppn=20 #if use setsm for coregistration (parallel)
ppn=38

# Load GNU Parallel
#module load parallel/20190222/intel-19.0.5 #loni
module load parallel/20150122
#module load gcc/5.2.0
#module load parallel/20170322

# Change work directory
cd $SLURM_SUBMIT_DIR

# Prepare joblist
joblist=(`cat $joblistfile`)

tmpdir=/scratch/local/$SLURM_JOB_ID/
echo Local tmp directory: $tmpdir

# Run parallelly

parallel --workdir $SLURM_SUBMIT_DIR -j$ppn --sshdelay=10\
         'echo Work on subdirectory: {}; subdir1=`basename {//}`; subdir2=`basename {}`; echo mkdir /scratch/local/$SLURM_JOB_ID/$subdir1/$subdir2; mkdir -p /scratch/local/$SLURM_JOB_ID/$subdir1/$subdir2; mv {//}/$subdir2 /scratch/local/$SLURM_JOB_ID/$subdir1/; rm -rf {//}/$subdir2; ln -sf /scratch/local/$SLURM_JOB_ID/$subdir1/$subdir2 {//}/$subdir2; ln -sf $PWD/job_group.pbs {}job_group.pbs; cd {}; ./job_group.pbs > out1; rm job_group.pbs; cd /scratch/local/$SLURM_JOB_ID/$subdir1; tar -cf $subdir2.tar.gz $subdir2; mv $subdir2.tar.gz {//}; cd {//}; rm $subdir2; tar -xf $subdir2.tar.gz; rm $subdir2.tar.gz ' ::: ${joblist[@]}

# Clear /scratch/local/$SLURM_JOB_ID/
cd /scratch/local/$SLURM_JOB_ID/
find . -type d -user $USER -exec rm -rf "{}" +

