
# ===============================================================================
# run_local.sh
# Purpose: Iterate over molecular conformers to run Psi4 calculations one at a time.
# Usage: edit HDIR varialble, then call "sh ./run_local.sh"
# Comment out the psi4 line and do a dry run before actually starting calculations!
# ===============================================================================

HDIR=/home/limvt/Documents/work_get-drugbank/alkanes

cd $HDIR

# collect all subdirectories of conformer input files
ALLJOBS=$(find ./ -mindepth 2 -maxdepth 2 -type d)
ALLJOBSARRAY=($ALLJOBS)

# print the number of calculations found
N=${#ALLJOBSARRAY[@]}
echo "DIRECTORY SIZE: $N"

# iterate over directories and run
for i in $(seq 0 $((N-1))); do

  WDIR=${ALLJOBSARRAY[i]}
  cd $WDIR
  
  echo "STARTING JOB $i:" `pwd`
  #psi4 -i input.dat -o output.dat
  cd $HDIR

done
