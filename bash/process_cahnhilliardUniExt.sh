#!/bin/bash

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/jdt45/flowers

# directory for simulations specific to cahnHilliardUniExt
simtypedir=$outputdir/chUniExt

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
NT=$1
NSKIPSTRAIN=$2
Lx=$3
Ly=$4
Lz0=$5
phi0=$6
startSeed=$7
numSeeds=$8
partition=$9
time="${10}"

let endSeed=$startSeed+$numSeeds-1

basestr=chUniExt_NT"$NT"_NSS"$NSKIPSTRAIN"_Lx"$Lx"_Ly"$Ly"_Lz0"$Lz0"_phi0"$phi0"
runstr="$basestr"_PROCESS
searchstr="$basestr"_seed

# access directory specific for this simulation
simdatadir=$simtypedir/$basestr
if [[ ! -d $simdatadir ]]
then
    echo -- sim directory "$simdatadir" does not exist, ending.
    exit 1
fi

# get mafile string to save data
savestr="$savedir"/"$basestr".mat
mvstr="$savedir"/"$basestr".mp4

# create matlab command
MCODE="addpath ~/sphereGel/viz/; processCahnHilliardUniExt('$simdatadir','$searchstr','$savestr','$mvstr'); quit"

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
# 1. NT
# 3. NSKIPSTRAIN
# 4. Lx
# 5. Ly
# 6. Lz
# 7. phi0
# 8. start seed (startSeed)
# 9. # of seeds (numSeeds)
# 10. partition
# 11. time














