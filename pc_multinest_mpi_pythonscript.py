# -- script for automation of speed testing --

#variable initialisation
num_nodes = 1               #number of phoenix nodes
num_taskspernode = 0        #number of MPI processes (tasks) per node
num_threadspertask = 0      #number of OpenMP threads (CPUs) per task

str(num_nodes)              #convert to string to use in filename
str(num_taskspernode)
str(num_threadspertask)

from subprocess import call

#creates job files with range of parameter values
for num_taskspernode in range(1,17):
    for num_threadspertask in range(1,17):
        if num_taskspernode*num_threadspertask == 32:

           #num_nodes*num_threadspernode = number of processes (eg. mpirun -np 16)
           num_processes = num_nodes*num_taskspernode
           filename = "pc_multinest_mpi_" + str(num_nodes) + '_' + str(num_taskspernode) + '_' + str(num_threadspertask)
           
           template = """#!/bin/bash

           #SBATCH --export=ALL
           #SBATCH -p batch
           #SBATCH -N {num_nodes}                           # number of nodes
           #SBATCH --ntasks-per-node={num_taskspernode}     # MPI tasks per node
           #SBATCH -c {num_threadspertask}                  # OMP threads (CPUs) per task

           export OMP_NUM_THREADS={num_threadspertask}

           #SBATCH --time=3-00:00:00      # time allocation

           #SBATCH --mem=32GB             # memory for all nodes

           #SBATCH -J filename    #job name - set to the same as the output

           #SBATCH --mail-type=FAIL
           #SBATCH --mail-type=END
           #SBATCH --mail-user=a1686947@student.adelaide.edu.au

           # Run the job from directory in which sbatch command was run
           cd $SLURM_SUBMIT_DIR

           mpirun -np {num_processes} ppr:{num_taskspernode}:node ./pc_multinest/pc_multinest_mpi output/{filename}
  
                  # np stands for number of processes (should be num_nodes*num_taskspernode)
                  # ppr stands for processes per resource
                  # pr stands for processing elements - binds number of porcessing elements to each process
                  # this should be the name of the file to sbatch (including path from home) and then the ouput directory and filename
           """
           context = {
           "num_nodes":num_nodes,
           "num_taskspernode":num_taskspernode,
           "num_threadspertask":num_threadspertask,
           "num_processes":num_processes,
           "filename":filename
           }
           
           with open(filename, 'w') as myfile:
                    myfile.write(template.format(**context))

           print "created file ", filename
           #call(["sbatch", filename])
           #print "submitted    ", filename

        else:
                pass

#call (["squeue", "-u", "a1686947"])

#unknown what modifications to do here - have commented out all loops and kept initialisation that I know works
