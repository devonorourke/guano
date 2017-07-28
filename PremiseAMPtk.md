# Using Premise - a few notes
One of the complications in using this compute cluster is that you have to use a job submission program called *Slurm*. The following code and instructions may differ slightly than a typical user operating all tasks using their own personal machine or even a single virtual machine like an Amazon Instance. The idea here is to leverage the power of many available pieces of hardware to speed up the process.  

## Step 1 - log in and transfer files if necessary:  
First, log in:  
`ssh devon@pinky.sr.unh.edu` and enter your password.  

Next, transfer the data. We typically receive Illumina **fq.gz** files in an @cobb directory, and move them from there into Premise for work. The following command assumes you are in the directory with the files you want to move.  
```
rsync -avzh Sample_*/*.gz devon@premise.sr.unh.edu:/mnt/lustre/macmaneslab/devon/p5data/lane1
```
Continue transferring files into appropriate directories as needed for all projects and lanes you may be analyzing.  

## Step 2 - Create slurm scripts
To create the scripts you'll want to run there are several potential arguments you can pass on to Slurm through a single script; these have a few caveats relative to a typical *bash* style script. Information provided by the ever helpful Toni Westbrook:  
- Slurm doesn't know Bash environmental variables.
> When the Slurm script runs it’s doing so on a different server than the head node, and therefore in a different login session.  This means that nothing you do in your terminal session (load modules, alter your PATH statement, whatever else) is going to have any effect on anything from the script’s point of view.  The script’s point of view is akin to you logging in fresh with nothing loaded yet.  For this reason, you need to put any module loads, alterations to your path, and anything else, in the actual body of the Slurm script  
- Always specify full paths. Again, Slurm doesn't know variables.  
> SBATCH commands are not interpreted by Bash, they’re interpreted by the sbatch.  Essentially, sbatch scans your script for anything that starts with #SBATCH, makes not of the option you told it, and then passes the script over to the interpreter specified in the shebang line (e.g. bash).  This is an important detail, as you can’t saying something like “#SBATCH -D $HOME/slurmtest”, as $HOME is an environmental variable that’s interpreted by bash – but sbatch processes it before bash does.  So to sbatch, you’re literally asking it to run in a directory called “$HOME” – it has no concept of environment variables. There are ways of getting around this, but the easiest is just to put the full path in – the workaround is almost not worth the effort.  
- The number of tasks is not the number of cores or nodes.  
> “#SBATCH --ntasks” is not telling Slurm how many threads you need, it’s telling Slurm how many tasks (processes) you need.  If you have one program that spins up 24 threads, you actually only want to use 1 ntask.  Using 24 ntasks is going to spin up 24 copies of amptk (which is not what you want).  
- Use **srun** in your scripts to help with debugging.  
> You always want to prepend “srun” before your actual long-running command in your slurm script.  This isn’t strictly necessary, but it makes the Slurm recordkeeping and resource monitoring work better.  You can check out the threaded.slurm to see an example of how this looks too.  
- Finally, use the MacManes parition because you can!  
> Assign the partition “macmanes,shared”, so that way if Matt’s partition is currently in use, but other servers free up, your job will take them  


So what does this code look like?  Here are some slurm scripts for code blocks that could be run individually or combined into a longer single script if you don't need to delve into contamination checks with your reads. In addition, these may be modified as needed, but generally are sufficient with the exception that you may want to rename file names, directories, etc.  

### amptk Illumina SLURM script
This script would work on a single lane of data, present in the directory `/mnt/lustre/macmaneslab/devon/slurmtest/p5data/lane1` and is labeled as **amptk-IlluminaL1.sh**. Notably, the equvalent script for Lane 2 data follows the same structure but substitutes L1 with L2.  

```
#!/bin/bash

## command executed is: sbatch amptk-illumina-test.sh
## see /mnt/lustre/hcgs/shared/slurm-templates for help on creating these scripts

#SBATCH -D /mnt/lustre/macmaneslab/devon/slurmtest/p5data/lane1
#SBATCH -p macmanes,shared
#SBATCH --job-name="orourke-amptkIlluminaL1"
#SBATCH --ntasks=1
#SBATCH --output=oroAmpIllL1.output

module purge
module load linuxbrew/colsa

PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun amptk illumina \
-i /mnt/lustre/macmaneslab/devon/slurmtest/p5data/lane1 \
-o p5L1 \
--rescue_forward on \
--require_primer off \
--min_len 160 \
--full_length \
--cpus 24 \
--read_length 250 \
-f GGTCAACAAATCATAAAGATATTGG \
-r GGWACTAATCAATTTCCAAATCC
```
