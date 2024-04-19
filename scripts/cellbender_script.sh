#!/bin/bash
#SBATCH --job-name=cellbender_raw_matrix                            # Job name
#SBATCH --nodes=2 #if you do 2 nodes then you can do 96 cores of 5300                                       # Number of cores job will run on
#SBATCH --time=7-00:00         #max is 7 days                         # Time limit hrs:min:sec
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=v.nikhilesh23k@gmail.com   # Email to which notifications will be sent
#SBATCH -A lu2023-2-6
#SBATCH --mem-per-cpu=5300 #5.3MB is the most you can do if you specify 48 cores 
#SBATCH --tasks-per-node=48 #here you set the cores and can be up to 48
pwd; hostname; date

cellbender remove-background --input /home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nick/Pleuro/outs/20230704_splenocytes_count_r2_10x_annoV2/outs/raw_feature_bc_matrix --output /home/nikvaku/snic2022-6-312/LabMemberScratchDir/Nikhilesh/Raw_data/Anno/matrix_filtered_clear/second_spleen/output.h5 --learning-rate 0.00005 --fpr 0.01 --epochs 150
