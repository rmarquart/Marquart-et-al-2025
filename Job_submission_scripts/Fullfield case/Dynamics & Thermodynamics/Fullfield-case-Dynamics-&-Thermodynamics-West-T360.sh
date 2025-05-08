#!/bin/sh
#SBATCH --account=civeng
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=36
#SBATCH --time=170:00:00
#SBATCH --job-name="Fullfield-D&T-West-T360"

#SBATCH --mail-user=mrqrut001@myuct.ac.za
#SBATCH --mail-type=fail,end

# Number of simulation T{num}:
num=360

# Specify the number of categories for both ice types and define the index value
numCategories_h1=6
numCategories_h2=6

index=2500

# Export the variable
export num numCategories_h1 numCategories_h2 index

# Set Python script path:
cd "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_real/T${num}_cont"

# Set this environment if your software is being called from a container
SINGULARITY_PATH=/opt/exp_soft/singularity-containers/openfoam
CONTAINER_ID=OpenFoam-2306
FOAM_BASHRC=/openfoam/bash.rc

# try to compile solver as part of a SLURM job:
cd /home/mrqrut001/applications/solvers/my_seaIce_06112024
singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && wclean"
singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && wmake"

# change directory to your calculation directory:
cd "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_real/T${num}"

# Execute bashrc and OpenFoam commands:
singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && blockMesh"

# Assign a value to the variables deltat and ntimesteps
deltat=$(echo "15.2" | bc)
ntimesteps1=$(echo "76" | bc)

formatted_deltat=$(printf "%.1f" $deltat)
formatted_ntimesteps1=$(printf "%.1f" $ntimesteps1)

# Export the variable
export formatted_ntimesteps1 formatted_deltat

# Loop to execute OpenFOAM and PYTHON several times:
ii=$deltat

while (( $(awk -v ii="$ii" -v ntimesteps1="$ntimesteps1" 'BEGIN {print (ii <= ntimesteps1)}') ))
do
	ii=$(printf "%.1f" $ii)

	# Export the variable
    export ii

	# OpenFOAM commands:
	echo "Start OpenFOAM..."

	# change directory to your calculation directory:
	cd "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_real/T${num}"

	# Function to check if a string ends with ".0"
	ends_with_zero() {
		[[ $1 =~ \.0$ ]]
	}

	# Iterate over folders
	for folder in */; do
		folder_name="${folder%/}"  # Remove trailing slash
		# Check if the folder name is not "0" and ends with ".0"
    	if [ "$folder_name" != "0" ] && ends_with_zero "$folder_name"; then
    		# Remove ".0" suffix from the folder name
    		new_folder_name="${folder_name%.*}"
    		# Rename the folder
    		mv "$folder_name" "$new_folder_name"
    		# echo "Folder '$folder_name' renamed to '$new_folder_name'"
    	fi
	done

	# Execute bashrc and OpenFoam commands:
	singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && decomposePar"
	singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && mpirun -np 36 my_seaIce_06112024 -parallel"
	singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && reconstructPar"  

	rm -r processor*

	# Function to check if a string is an integer
	is_integer() {
	    [[ $1 =~ ^[0-9]+$ && $1 != 0 ]]
	}

	# Iterate over folders
	for folder in */; do
	    folder_name="${folder%/}"  # Remove trailing slash
	    # Check if the folder name represents an integer (except for "0")
	    if is_integer "$folder_name"; then
	        # Add .0 to the folder name
	        new_folder_name="${folder_name}.0"
	        # Rename the folder
	        mv "$folder_name" "$new_folder_name"
	        # echo "Folder '$folder_name' renamed to '$new_folder_name'"
	    fi
	done
	echo ""

	cp -r "$ii" "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_storage_real/T${num}_storage_real"
	
	cd "$ii"/
	rm -r uniform 

	# Set Python binary path:
	module load python/miniconda3-py3.12-usr

	# Set Python script path:
	cd "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_real/T${num}_cont"

	# Run Python script:
	echo "Start updating controlDict..."
	python create_controlDict.py
	echo "End updating controlDict..."
	echo ""

	ii=$(echo "$ii + $deltat" | bc)

	echo "End OpenFOAM..."
	echo ""	
done

# Set Python binary path:
module load python/miniconda3-py3.12-usr

# Set Python script path:
cd "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_real/T${num}_cont"

echo "Start importing files..."
python importalpha_tplus.py
python importh_tplus.py
python importeta_tplus.py
echo "End importing files..."
echo ""

echo "Start creating thickness categories..."
python categories_part1.py
echo "End creating thickness categories..."
echo ""

echo "Start creating initial previous time step files..."
python initial_previous_index_values.py
echo "End creating initial previous time step files..."
echo ""

# Set path to scripts thermodynamics model:
cd "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_real/T${num}_cont"

# Assign a value to the variable ntimesteps
ntimesteps2=$(echo "3009.6" | bc)

end_condition=$(echo "$ntimesteps2" | bc)

# Loop to execute OpenFOAM and PYTHON several times:
while (( $(echo "$ii <= $end_condition" | bc -l) )); do

	# Export the variable
    export ii

    echo "Start copying and removing previous time step files"
    python copy_or_remove_previous_index_values_h1.py
    python copy_or_remove_previous_index_values_h2.py
    echo "End copying and removing previous time step files"
    echo ""
	echo "Start thermodynamics model..." 
	python RM_ESIM3_h1_FixedInitial_v2.py
	echo ""
	python RM_ESIM3_h2_FixedInitial_v2.py
	echo "End thermodynamics model..."
	echo ""  

	echo "Start calculating offsets..."
	python categories_part2.py
	echo "End calculating offsets..."
	echo ""

	# OpenFOAM commands:
	echo "Start OpenFOAM..."

	# change directory to your calculation directory:
	cd "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_real/T${num}"

	# Function to check if a string ends with ".0"
	ends_with_zero() {
		[[ $1 =~ \.0$ ]]
	}

	# Iterate over folders
	for folder in */; do
		folder_name="${folder%/}"  # Remove trailing slash
		# Check if the folder name is not "0" and ends with ".0"
    	if [ "$folder_name" != "0" ] && ends_with_zero "$folder_name"; then
    		# Remove ".0" suffix from the folder name
    		new_folder_name="${folder_name%.*}"
    		# Rename the folder
    		mv "$folder_name" "$new_folder_name"
    		# echo "Folder '$folder_name' renamed to '$new_folder_name'"
    	fi
	done

	# Execute bashrc and OpenFoam commands:
	singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && decomposePar"
	singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && mpirun -np 36 my_seaIce_06112024 -parallel"
	singularity exec $SINGULARITY_PATH/$CONTAINER_ID.sif bash -c "source $FOAM_BASHRC && reconstructPar"

	rm -r processor*

	# Function to check if a string is an integer
	is_integer() {
	    [[ $1 =~ ^[0-9]+$ && $1 != 0 ]]
	}

	# Iterate over folders
	for folder in */; do
	    folder_name="${folder%/}"  # Remove trailing slash
	    # Check if the folder name represents an integer (except for "0")
	    if is_integer "$folder_name"; then
	        # Add .0 to the folder name
	        new_folder_name="${folder_name}.0"
	        # Rename the folder
	        mv "$folder_name" "$new_folder_name"
	        echo "Folder '$folder_name' renamed to '$new_folder_name'"
	    fi
	done

	cp -r "$ii" "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_storage_real/T${num}_storage_real"
	
	# Iterate through all subdirectories
	for dir in [0-9]*/; do
		dir=${dir%/}  # Remove trailing slash
		if [ -d "$dir/uniform" ]; then
			# If 'uniform' folder exists, delete it
			# echo "Deleting $dir/uniform"
			rm -r "$dir/uniform"
		fi
	done

	echo "End OpenFOAM..."
	echo ""
	echo ""

	# Set Python binary path:
	module load python/miniconda3-py3.12-usr

	# Set Python script path:
	cd "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_real/T${num}_cont"

	# Run Python scripts:
	echo "Start importing files..."
	python importalpha_tplus.py
	python importh_tplus.py
	python importeta_tplus.py
	echo "End importing files..."
	echo ""

	echo "Start updating thickness categories and include offsets..."
	python categories_part3.py
	echo "End updating thickness categories and include offsets..."
	echo ""	

    # Set Python script path:
	cd "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_storage_real/T${num}"

	# Run Python scripts:
	echo "Start creating new thickness files..."
	python create_h.py
	echo "End creating new thickness files..."
	echo ""
	
	cp h "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_storage_real/T${num}_storage_real/$ii/" &&
	cp h "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_real/T${num}/$ii/" &&
	rm h h_dummy.txt

	# Set Python script path:
	cd "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_real/T${num}_cont"

	# Run Python script:
	echo "Start updating controlDict..."
	python update_controlDict.py
	echo "End updating controlDict..."
	echo ""

	# Set Python script path:
	cd "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_real/T${num}_cont"

	# Set Python script path:
	cd "/scratch/mrqrut001/simulation/postdoc/Github/Fullfield case/Dynamics & Thermodynamics/West/Python_real/T${num}_cont"

	# Run Python scripts:
	echo "Start importing files..."
	python importalpha_tplus.py
	python importh_tplus.py
	python importeta_tplus.py
	echo "End importing files..."
	echo ""

	echo "Start creating new thickness categories..."
	python categories_part4.py
	echo "End creating new thickness categories..."
	echo ""

	ii=$(echo "$ii + $deltat" | bc)

	echo "End Python..."
	echo ""
	echo ""
done
