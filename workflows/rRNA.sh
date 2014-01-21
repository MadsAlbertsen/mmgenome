echo ""
echo "Extracting 16S sequences"
makeblastdb -in assembly.fa -dbtype nucl > temp.blast2
blastn -query /space/users/malber06/Desktop/Databases/gg_90id/gg_90.fasta -db assembly.fa -num_threads 60 -max_target_seqs 5 -outfmt 6 -evalue 1e-10 -out temp.blast.txt
perl /space/users/malber06/mmgenome/scripts/extract.long.hits.from.blast.pl -b temp.blast.txt -d assembly.fa -m 500 -o 16S.fa

echo ""
echo "Extracting 23S sequences"
blastn -query /space/users/malber06/Desktop/Databases/silva_lsu_90id/silva.lsu.111.90id.fa -db assembly.fa -num_threads 60 -max_target_seqs 5 -outfmt 6 -evalue 1e-10 -out temp.blast.txt
perl /space/users/malber06/mmgenome/scripts/extract.long.hits.from.blast.pl -b temp.blast.txt -d assembly.fa -m 500 -o 16S.fa

rm temp.blast.txt
rm assembly.fa.nhr
rm assembly.fa.nin
rm assembly.fa.nsq
