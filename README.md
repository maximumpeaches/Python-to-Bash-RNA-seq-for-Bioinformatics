# README #

### Purpose ###
This project was designed to accelerate the RNA-seq pipeline at Dr. Jiang's lab. 

### How it works ###
It consists of a Python script which, when given some user-defined parameters, will output a Bash script and submit it to SLURM on USF's HPC cluster. The input file to the Python script contains the locations of a directory containing RNA-seq reads, a GTF file, and a genome index in FASTA format. To run the program the user would edit the input file and then type `python slurm_generator.py` to create a few Bash files and submit them to the USF Research Computing cluster. The input file and Bash script are for illustration purposes and will look different after the user modifies the input file.

### Author info ###
Written by Maxwell Pietsch for Dr. Rays Jiang's bioinformatics group at USF

maxwell.pietsch@gmail.com

http://maxwell9.myweb.usf.edu
