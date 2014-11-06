#!/bin/sh

#PBS -N pore10r 
#PBS -l nodes=1:ppn=2:del_int_20_128
#PBS -l walltime=1:00:00
#PBS -M kgsteen@ku.edu
#PBS -m ae
#PBS -d /data/kgs32/cass1/forEd/
#PBS -e /data/kgs32/cass1/forEd/pore10r.err
#PBS -o /data/kgs32/cass1/forEd/pore10r.out

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

echo Starting Time is `date`

cd ${PBS_O_WORKDIR}
export OMP_NUM_THREADS=2
/users/kgs32/Cassandra/Src/cassandra.parallel Si_meOH.inp

echo Ending Time is `date`
exit 0

