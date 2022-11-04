#!/bin/bash
for i in {1..23}
	do 
	wget "https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.$i.fa.gz"	
	done
gunzip Homo*
cat Homo* >> ref.fa
