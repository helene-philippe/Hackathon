#!/bin/bash
for i in {1..22}
	do 
	wget "https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.$i.fa.gz"	
	done
wget "https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz"
wget "https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz"
gunzip Homo*
cat Homo* >> ref.fa
wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.chr.gtf.gz
gunzip *.gtf.gz
mv *gtf ref.gtf
