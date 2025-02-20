#MAFFT for alignments
mafft --auto "$INPUT_FILE" > "$OUTPUT_MAFFT"

#TrimAl for alignemnt trimming
trimal -in $OUTPUT_MAFFT.fasta -out $OUTPUT_TrimAL.fasta -gappyout

#IQ-TREE for phylogeny inference
iqtree2 -s OUTPUT_TrimAL.fasta -bb 1000 -nt AUTO -pre ID