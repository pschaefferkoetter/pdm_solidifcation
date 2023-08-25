#!/bin/bash
#SBATCH --nodes=1                       # Number of requested nodes
#SBATCH --time=72:00:00                 # Max walltime of 2 minutes
#SBATCH --qos=blanca-ceae               # Specify your Blanca QoS
#SBATCH --partition=blanca-ceae         # Specify your Blanca nodes
#SBATCH --output=matlab_%j.out          # Output file name

# Load Matlab module
module load matlab
  

# Run matlab without a GUI
matlab -nodisplay -nodesktop -r "case_1"
