#!/bin/sh
#SBATCH --account=civeng
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=16
#SBATCH --time=170:00:00
#SBATCH --job-name="Test-D-Horizontal-T120"

#SBATCH --mail-user=rutger.marquart@uct.ac.za
#SBATCH --mail-type=fail,end

# Set this environment if your software is being called from a container
SINGULARITY_PATH=/opt/exp_soft/singularity-containers/openfoam
CONTAINER_ID=OpenFoam-2306
FOAM_BASHRC=/openfoam/bash.rc

# try to compile solver as part of a SLURM job
cd applications/solvers/my_seaIce_06112024
singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && wclean"
singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && wmake"

# change directory to your calculation directory
cd /scratch/mrqrut001/simulation/postdoc/Github/Test case/Dynamics/Horizontal/T120

# Execute bashrc and OpenFoam Commands
singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && blockMesh"
singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && decomposePar"
singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && mpirun -np 16 my_seaIce_06112024 -parallel"
# singularity exec $SINGULARITY_PATH/openfoam/openfoam-v2106.sif reconstructPar
