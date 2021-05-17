#!/usr/bin/env bash

#######################################################################################
###                                                                                 ###
###     PlasForest 1.1                                                              ###
###     Copyright (C) 2020  Léa Pradier, Tazzio Tissot, Anna-Sophie Fiston-Lavier,  ###
###                         Stéphanie Bedhomme. (leaemiliepradier@gmail.com)        ###
###                                                                                 ###
###     This program is free software: you can redistribute it and/or modify        ###
###     it under the terms of the GNU General Public License as published by        ###
###     the Free Software Foundation, either version 3 of the License, or           ###
###     (at your option) any later version.                                         ###
###                                                                                 ###
###     This program is distributed in the hope that it will be useful,             ###
###     but WITHOUT ANY WARRANTY; without even the implied warranty of              ###
###     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               ###
###     GNU General Public License for more details.                                ###
###                                                                                 ###
###     You should have received a copy of the GNU General Public License           ###
###     along with this program. If not, see <http://www.gnu.org/licenses/>.        ###
###                                                                                 ###
#######################################################################################

PFdirectory=$PWD
cd $HOME
git clone https://github.com/StuntsPT/NCBI_Mass_Downloader
cd NCBI_Mass_Downloader
echo "The database is about to be downloaded."
echo "This may take a while."
while true; do
    python3 NCBI_downloader.py -q '(bacteria[filter] AND biomol_genomic[PROP] AND refseq[filter] AND is_nuccore[filter] AND plasmid[filter]) AND ("0001/01/01"[PDAT] : "2019/09/01"[PDAT])' -d "nucleotide" -o ${PFdirectory}/plasmid_refseq.fasta
    if test -f ${PFdirectory}/plasmid_refseq.fasta; then
        if [ $(grep -o ">" ${PFdirectory}/plasmid_refseq.fasta | wc -l) -ge 35000 ]; then
            break
        fi
    fi
    rm ${PFdirectory}/plasmid_refseq.fasta
done

cd $PFdirectory
while true; do
    ntheoseq=$(wc -l < list_ids.txt)
    nseq=$(grep -o ">" plasmid_refseq.fasta | wc -l)
    if [ $nseq -lt $ntheoseq ]; then
        python3 check_and_download_database.py check
    else
        break
    fi
done
makeblastdb -in plasmid_refseq.fasta -dbtype nucl -parse_seqids

