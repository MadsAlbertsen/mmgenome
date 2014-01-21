#!/bin/bash

##### Description
# This small shell script wraps the data generation for R in Albertsen et. al., 2013

##### Needed input files
# Fasta file with all assembled scaffolds (keep the naming as >1, >2 etc): assembly.fa

##### Needed software
# Prodigal
# HMMER 3.0
# BLAST
# MEGAN
# Perl scripts from : git clone https://github.com/MadsAlbertsen/mmgenome.git        

clear
echo "---Metagenomics workflow script v.2.0.0---"

echo ""
echo "Calculating tetranucleotide frequency"
perl mmgenome/scripts/calc.kmerfreq.pl -i assembly.fa -o assembly.kmer.tab

echo ""
echo "Calculating gc content"
perl mmgenome/scripts/calc.gc.pl -i assembly.fa -o assembly.gc.tab

echo ""
echo "Finding essential genes - Predicting proteins (Prodigal)"
prodigal -a temp.orfs.faa -i assembly.fa -m -o temp.txt -p meta -q
cut -f1 -d " " temp.orfs.faa > assembly.orfs.faa

echo ""
echo "Finding essential genes - running HMM search"
hmmsearch --tblout assembly.hmm.orfs.txt --cut_tc --notextw mmgenome/scripts/essential.hmm assembly.orfs.faa > hmm.temp.txt
tail -n+4  assembly.hmm.orfs.txt | sed 's/ * / /g' | cut -f1,4 -d " " | sed 's/_/ /' > assembly.orfs.hmm.id.txt
grep -v "#" assembly.hmm.orfs.txt | cut -f1 -d " " > list.of.positive.orfs.txt
perl mmgenome/scripts/extract.using.header.list.pl -l list.of.positive.orfs.txt -s assembly.orfs.faa -o assembly.orfs.hmm.faa

echo ""
echo "Finding essential genes - Blasting positive hits"
blastp -query assembly.orfs.hmm.faa -db refseq_protein -evalue 1e-5 -num_threads 60 -max_target_seqs 5 -outfmt 5 -out assembly.orfs.hmm.blast.xml

echo ""
echo "Finding essential genes - Extracting consensus taxonomic assignment"
MEGAN +g -x "import blastfile= assembly.orfs.hmm.blast.xml meganfile=temp.rma;recompute toppercent=5;recompute minsupport=1;update;collapse rank=Species;update;select nodes=all;export what=CSV format=readname_taxonpath separator=tab file=assembly.orfs.hmm.blast.tax.txt;update;close"
perl mmgenome/scripts/hmm.majority.vote.pl -i assembly.orfs.hmm.blast.tax.txt -o tax.txt

echo ""
echo "Generating connection network"
perl mmgenome/scripts/network.pl -i mapping.sam -f 2

echo ""
echo "Removing temp files"
rm hmm.temp.txt
rm list.of.positive.orfs.txt
rm assembly.orfs.hmm.blast.tax.txt
rm temp.orfs.faa
rm temp.txt
rm temp.rma
rm assembly.orfs.hmm.blast.xml
rm assembly.orfs.hmm.faa
rm assembly.hmm.orfs.txt

echo ""
echo "done"

