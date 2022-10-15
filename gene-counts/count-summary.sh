# curate data
python3 get-stats-gb.py NewData/mito.genomic.gbff NewData/mt Prelims/stats-codon.csv Prelims/stats-residue.csv > NewOutputs/mt-process.txt
python3 get-stats-gb.py NewData/plastid.genomic.gbff NewData/pt Prelims/stats-codon.csv Prelims/stats-residue.csv > NewOutputs/pt-process.txt
python3 get-raw-counts.py NewData/mt-dna.fasta NewData/mt-gene-raw-occurrence.csv
python3 get-raw-counts.py NewData/pt-dna.fasta NewData/pt-gene-raw-occurrence.csv

# extract gene barcodes
python3 get-feature-labels.py Prelims/stats-residue.csv Prelims/stats-codon.csv  Species,Compartment,GeneLabel, > NewData/mt-stats-manual.csv
python3 process-labels.py 10 NewData/mt-dna.fasta Prelims/mt-manual-replace.csv NewData/mt-stats-manual.csv NewData/mt-barcodes-manual.csv NewData/mt-species-manual.txt NewData/mt-gene-occurrence-manual.csv
python3 get-feature-labels.py Prelims/stats-residue.csv Prelims/stats-codon.csv  Species,Compartment,GeneLabel, > NewData/pt-stats-manual.csv
python3 process-labels.py 10 NewData/pt-dna.fasta Prelims/pt-manual-replace.csv NewData/pt-stats-manual.csv NewData/pt-barcodes-manual.csv NewData/pt-species-manual.txt NewData/pt-gene-occurrence-manual.csv 

# make plots
Rscript belen-paper.R
