#!/bin/bash
#SBATCH -A SNIC2018-3-18 # group to which you belong
#SBATCH -p hebbe  # partition (queue)
#SBATCH -N 1  # number of nodes
#SBATCH -n 20  # number of cores
#SBATCH -t 1-00:00 # runtime limit (D-HH:MM)
#SBATCH -o logs/blast_run.sh.%N.%j.out  # STDOUT
#SBATCH -J blast_run.sh  # partition (queue)
#SBATCH --mail-type=END,FAIL    # notifications for job done & fail
#SBATCH --mail-user=alezel@chalmers.se # send-to address

mkdir -p ./logs
echo "logged to: logs/blast_run.sh.%N.%j.out"
data_dir="data"
output_dir="r_blast"
module load Anaconda3
source activate my_root
mkdir -p $output_dir

command="blastp -query $data_dir/Enzymes_2018_11_27  -db blast_db/ensembl_proteins_all -out $output_dir/Enzymes_2018_11_27_VS_ensembl_proteins_all  -outfmt 6 -num_threads $SLURM_NPROCS"

eval $command
