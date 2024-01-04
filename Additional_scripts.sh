#### Additional shell scripts ###

#LEFSE

singularity shell lefse.sif
lefse_format_input.py input.txt input.in -c 1 -u 2 -o 1000000
lefse_run.py input.in input.res

#VIBRANT to predict lifestyles (lysogenic/lytic) for contigs with minimum sequence length of 1000bp and containing at least 4 ORFs (open readings frames) 

singularity shell vibrant.sif
VIBRANT_run.py -i viral_contigs.fasta -t 114 -folder VIBRANT_results -virome