#!/bin/bash

# output directories
outputdir=/gpfs/loomis/project/fas/ohern/jdt45/flowers
simtypedir=$outputdir/chUniExtQD

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
NT=$1
NPRINT=$2
NSKIPSTRAIN=$3
A=$4
Lx=$5
Ly=$6
Lz0=$7
phi0=$8
startSeed=$9
numSeeds="${10}"
partition="${11}"
time="${12}"

let endSeed=$startSeed+$numSeeds-1

basestr=chUniExtQD_NT"$NT"_NSS"$NSKIPSTRAIN"_A"$A"_Lx"$Lx"_Ly"$Ly"_Lz0"$Lz0"_phi0"$phi0"
runstr="$basestr"_s"$startSeed"_e"$endSeed"

# create specific directory for sim
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# create task file
taskf=tasks/"$runstr".task
rm -f $taskf

# add commands to task file
let fcount=0
for seed in `seq $startSeed $endSeed`; do
 	# count files
    let fcount=$fcount+1

    # string for task file
    fpattern="$basestr"_seed"$seed"
    seeddir="$simdatadir"/"$fpattern"
    mkdir -p $seeddir
    ftype="$seeddir"/"$fpattern"
	MCODE="addpath ~/sphereGel/src; cahnHilliardUniExtQD($NT,$NPRINT,$NSKIPSTRAIN,$A,$Lx,$Ly,$Lz0,$phi0,$seed,'$ftype'); quit"

	# add to task file
	echo matlab -nodisplay -r \""$MCODE"\" >> $taskf
done

# test if task file was created
if [[ ! -f "$taskf" ]]
then
    echo task file not created, ending before job submission
    exit 1
fi

# setup slurm files
slurmf=slurm/"$runstr".slurm
job_name="$runstr"
runout=out/"$runstr"-%a.out
rm -f $slurmf

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH --array=1-$fcount >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo \#SBATCH -t $time >> $slurmf
echo module load MATLAB >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch $slurmf

# ====================
#       INPUTS
# ====================
# 1. NT
# 2. NPRINT
# 3. NSKIPSTRAIN
# 4. A
# 5. Lx
# 6. Ly
# 7. Lz
# 8. phi0
# 9. start seed (startSeed)
# 10. # of seeds (numSeeds)
# 11. partition
# 12. time


