#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=mbarnet
#SBATCH --time=72:00:00
#SBATCH -o eq_%A.out
#SBATCH -e eq_%A.err

# sourcing
export BSS="/home/anna/anaconda3/bin/activate /home/anna/mambaforge/envs/pipeline_annamherz" # to use the conda env to make sure sire works correctly - sourced in each sh script
#export BSS="/home/anna/anaconda3/bin/activate /home/anna/anaconda3/envs/biosimspace-dev"
export amber="/home/anna/amber22" # sourced in each script
export gromacs="/usr/local/gromacs-23.1/bin/GMXRC" # sourced in each script

# sourcing - as needed in the othe sh scripts
source $BSS
source $amber/amber.sh
source $gromacs

#export PYTHONPATH="$HOME/BioSimSpace/python:$PYTHONPATH"
export PYTHONPATH="$HOME/biosimspace_annamherz/python:$PYTHONPATH"
export PYTHONPATH="$HOME/Documents/code/python:$PYTHONPATH"

date
start=`date +%s`

echo "Folder for these runs is : $MAINDIRECTORY"
echo "Analysis for the transformation $1, $2."

echo running eq detection ....
python run_convergence.py tyk2
python run_convergence.py mcl1
python run_convergence.py p38
python run_convergence.py hif2a
python run_convergence.py syk
python run_convergence.py cmet

end=`date +%s`
runtime=$((end-start))
echo "runtime was $runtime seconds, which is $((runtime/60)) minutes"
