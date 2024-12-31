#!/bin/bash -l
#PBS -l select=2:system=polaris
#PBS -l place=scatter
#PBS -l walltime=01:00:00
#PBS -l filesystems=home:eagle
#PBS -N af3_dc2_debug
#PBS -q debug
#PBS -A [PROJECT]
#PBS -M mwoodc2@uic.edu
#PBS -m bae


# Set number of AF3 runs per node
N_PER_NODE=2

# Define base AF3 macro-directory
BASE_DIR="/lus/eagle/projects/DirectContacts2/miles/af3"

# Define local scratch directory to store MSA database
LOCAL_SCRATCH="/local/scratch/"

# Define AF3 root directory
AF3_DIR="$BASE_DIR/alphafold3"
AF3_LOC="$LOCAL_SCRATCH/af3_db"

# Define input, output directories
INPUT_DIR="$BASE_DIR/input_af3/json_input"
OUTPUT_DIR="$BASE_DIR/output_af3"

# Define AF3 parameter directory
PARAM_DIR="$BASE_DIR/af3_params"

# Define AF3 database directory
AF3_DATA="$BASE_DIR/af3_db"
AF3_DATA_COMPRESS="$BASE_DIR/../af3_db"



# Enable GPU-MPI
export MPICH_GPU_SUPPORT_ENABLED=1

cd $BASE_DIR

# Load necessary modules for running containers
module use /soft/modulefiles
module load spack-pe-base/0.8.1
module load apptainer

# Set proxy for internet access
export HTTP_PROXY=http://proxy.alcf.anl.gov:3128
export HTTPS_PROXY=http://proxy.alcf.anl.gov:3128
export http_proxy=http://proxy.alcf.anl.gov:3128
export https_proxy=http://proxy.alcf.anl.gov:3128

#ADDITIONAL_PATH=/opt/cray/pe/pals/1.2.12/lib
#module load cray-mpich-abi

# Get number of files in input directory
NUM_FILES=$(ls -1 | wc -l)
NNODES=`wc -l < $PBS_NODEFILE`

#FILE_START=101
#FILE_END=$((FILE_START + NNODES - 1))

# Make array of files, based on number of nodes
#file_array=($(ls ${INPUT_DIR} | sort | sed -n "$((FILE_START + 1)),$((FILE_END + 1))p"))
file_array=($(ls ${INPUT_DIR} | sort))
num_files=${#file_array[@]}


# Split PBS node file into separate host files containing addressable node names
split --lines=1 --numeric-suffixes=1 --suffix-length=3 $PBS_NODEFILE local_hostfile.


# Pre-stage AF3 MSA database to each compute node being used
echo "Copying AF3 database files to local scratch on each node..."

for lh in local_hostfile*; do
	
	mpiexec -n 1 --ppn 1 --hostfile ${lh} \
		time bash -c "rsync -aP $AF3_DATA_COMPRESS $LOCAL_SCRATCH && tar -xzf $AF3_LOC/mmcif_files.tar.gz -C $AF3_LOC/" &

	sleep 5s

done

wait

echo "Successfully copied"



# Loop over distinct nodes, start an AF3 instance on each
i=0
for lh in local_hostfile*; do
	(
		start_idx=$((i * num_files / NNODES))
		end_idx=$(((i + 1) * num_files / NNODES - 1))

		for ((j=start_idx ; j<=end_idx ; j++)); do
	
			INPUT_FILE=${file_array[$j]}
	
			echo "Processing file: $INPUT_FILE on node $lh"
	
			time mpiexec -n 1 --ppn 1 --hostfile ${lh} \
				apptainer exec \
				--nv \
				--bind $INPUT_DIR:/root/input_af3 \
				--bind $OUTPUT_DIR:/root/output_af3 \
				--bind $PARAM_DIR:/root/models \
				--bind $AF3_LOC:/root/databases \
				$AF3_DIR/alphafold3.sif \
				python3 $AF3_DIR/run_alphafold.py \
				--json_path=/root/input_af3/${INPUT_FILE} \
				--model_dir=/root/models \
				--db_dir=/root/databases \
				--output_dir=/root/output_af3
		done
	) &

	i=$((i + 1))

done

wait

rm local_hostfile*
