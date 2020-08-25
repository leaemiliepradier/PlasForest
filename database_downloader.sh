#!/usr/bin/env bash

PFdirectory=$pwd
cd $HOME
git clone https://github.com/StuntsPT/NCBI_Mass_Downloader
cd NCBI_Mass_Downloader
python NCBI_downloader.py -q '(bacteria[filter] AND biomol_genomic[PROP] AND refseq[filter] AND is_nuccore[filter] AND plasmid[filter]) AND ("0001/01/01"[PDAT] : "2019/09/01"[PDAT])' -d "nucleotide" -o ${PFdirectory}/plasmid_refseq.fasta
cd $PFdirectory

makeblastdb -in plasmid_refseq.fasta -dbtype nucl -parse_seqids

