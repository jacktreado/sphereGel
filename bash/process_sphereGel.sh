#!/bin/bash

# directories with code
maindir=~/flowers/phasesep

# directory for all output for cell simulations
outputdir=/gpfs/project/fas/ohern/jdt45/flowers

# directory for simulations specific to sphereGel
simtypedir=$outputdir/sphereGel

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
dr=$2
dphi=$3
dlz=$4
l2=$5
partition=$6
time=$7
startSeed="${8}"
endSeed="${9}"

# name strings
basestr=sgel_N"$N"_dr"$dr"_dphi"$dphi"_dlz"$dlz"_l2"$l2"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"_PROCESS

# access directory specific for this simulation
simdatadir=$simtypedir/$basestr
if [[ ! -d $simdatadir ]]
then
    echo -- sim directory "$simdatadir" does not exist, ending.
    exit 1
fi

# create task file
taskf=tasks/"$runstr".task
rm -f $taskf

# loop over files
flist="$simdatadir"/*.xyz
let fcount=0

# LOOP OVER FILES
for f in $flist; do
    # parse file name
    file=${f##*/}
    baseid=${file%%.xyz}
    seed=${baseid#*seed*}

    # check seeds
    if [[ $seed -gt $endSeed ]]; then
        echo - - FILE $file past seed search, skipping...
        continue
    elif [[ $seed -lt $startSeed ]]; then
        echo - - FILE $file before seed search, skipping...
        continue
    else
        echo - - Processing sphereGel sim. with filename: "$file"
        let fcount=$fcount+1
    fi

    # get mafile string
    savestr="$savedir"/"$baseid".mat

    # create matlab command
    MCODE="addpath ~/sphereGel/viz/; processSphereGel('$f','$savestr'); quit"

    # append to runString
    runString="matlab -nodisplay -r \""$MCODE"\" "

    # echo to task file
    echo "$runString" >> $taskf
done

# test if task file was created
if [[ ! -f "$taskf" ]]
then
    echo task file not created, ending before job submission
    exit 1
fi

# get number of jobs to submit to each array
let arraynum=$fcount
echo -- total number of array runs = $arraynum

# setup slurm files
slurmf=slurm/"$runstr".slurm
job_name="$runstr"
runout=out/"$runstr"-%a.out
rm -f $slurmf

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH --array=1-$arraynum >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
echo sbatch -t $time $slurmf


# ====================
#       INPUTS
# ====================
# 1. N
# 2. dr
# 3. dphi
# 4. dlz
# 5. l2
# 6. partition
# 7. time
# 8. startSeed
# 9. endSeed




