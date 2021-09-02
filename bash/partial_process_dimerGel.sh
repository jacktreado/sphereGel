#!/bin/bash

# directories with code
maindir=~/flowers/phasesep

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/jdt45/flowers/dimerGel

# directory for simulations specific to sphereGel
simtypedir=$outputdir/dimerGel

# directory to save matfiles
savedir=$simtypedir/matfiles

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p $savedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
N=$1
dl0=$2
dg=$3
del=$4
l2=$5
partition=$6
time=$7

# name strings
basestr=dgel_N"$N"_dl0"$dl0"_dg"$dg"_del"$del"_l2"$l2"
runstr="$basestr"_PRT_PROCESS
searchstr="$basestr"_seed

# access directory specific for this simulation
simdatadir=$simtypedir/$basestr
if [[ ! -d $simdatadir ]]
then
    echo -- sim directory "$simdatadir" does not exist, ending.
    exit 1
fi

# get mafile string to save data
savestr="$savedir"/partial_"$basestr".mat

# create matlab command
MCODE="addpath ~/sphereGel/viz/; partialProcessDimerGel('$simdatadir','$searchstr','$savestr'); quit"

# setup slurm files
slurmf=slurm/"$runstr".slurm
job_name="$runstr"
runout=out/"$runstr".out
rm -f $slurmf

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo \#SBATCH --mem-per-cpu=50G >> $slurmf
echo module load MATLAB >> $slurmf
echo matlab -nodisplay -r \""$MCODE"\" >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf


# ====================
#       INPUTS
# ====================
# 1. N
# 2. dl0
# 3. dg
# 4. del
# 5. l2
# 6. partition
# 7. time




