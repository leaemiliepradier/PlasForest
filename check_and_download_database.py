#!/usr/bin/env python
# -*- coding: utf-8 -*-

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


from Bio import SeqIO
import sys
import os
import re
from Bio import Entrez


def check_missing():
    with open("list_ids.txt","r+") as f:
        list_ids = [x.split("\n")[0] for x in f.readlines()]
    list_downloaded = []
    if os.path.isfile("plasmid_refseq.fasta"):
        with open("plasmid_refseq.fasta","r+") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                list_downloaded.append(record.id)
        list_missing = [x for x in list_ids if x not in list_downloaded]
    else:
        list_missing = list_ids
    return list_missing


def get_response(nmissing):
    print("There are",nmissing,"missing sequences.")
    print("PlasForest may have reduced performances without these sequences.")
    response=input("Do you want to download them? (y/n) ")
    while response != "y" and response != "n":
        response=input("Do you want to download them? (y/n) ")
    if response=="y":
        print("Good choice!")
    else:
        print("Got it. PlasForest will not access the whole database.")
        print("If you want to do it later, you can call this function:")
        print("")
        print("python check_and_download_database.py download")
        print("")
    return response

def get_email():
    regex = '^(\w|\.|\_|\-)+[@](\w|\_|\-|\.)+[.]\w{2,3}$'
    print("The missing sequences will be downloaded from NCBI Entrez.")
    print("Entrez requires your email address to download a large batch of sequences.")
    print("It will only be communicated to NCBI and will not be collected for any purpose.")
    email=input("Please enter an email address: ")
    while not re.search(regex,email):
        email=input("Please enter a valid email address: ")
    return email


def download_missing(list_ids, email):
    Entrez.email = email
    request = Entrez.epost("nucleotide",id=",".join(map(str,list_ids)))
    result = Entrez.read(request)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    handle = Entrez.efetch(db="nucleotide",rettype="fasta", retmode="text",webenv=webEnv, query_key=queryKey)
    for record in SeqIO.parse(handle,"fasta"):
        with open("plasmid_refseq.fasta","a+") as out:
            SeqIO.write(record,out,"fasta")
            
def make_blast_database():
    os.system("makeblastdb -in plasmid_refseq.fasta -dbtype nucl -parse_seqids")

if sys.argv[1] == "check":
    list_missing = check_missing()
    nmissing = len(list_missing)
    if nmissing > 0:
        response = get_response(nmissing)
        if response=="y":
            email=get_email()
            download_missing(list_missing, email)
elif sys.argv[1] == "download":
    list_missing = check_missing()
    nmissing = len(list_missing)
    if nmissing > 0:
        email=get_email()
        download_missing(list_missing,email)
        make_blast_database()


