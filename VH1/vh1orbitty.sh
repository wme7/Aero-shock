#$ -S /bin/bash
#$ -pe mpich 16
#$ -cwd
#$ -N TEST16

source /etc/profile.d/apps.sh

mpirun -np $NSLOTS -machinefile $TMPDIR/machines ./vh1-mpi

