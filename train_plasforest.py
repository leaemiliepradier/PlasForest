#!/usr/bin/env python
# -*- coding: utf-8 -*-

#######################################################################################
###                                                                                 ###
###     PlasForest 1.2                                                              ###
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

import os
import sys, getopt
import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import pandas as pd
import numpy as np
import pickle
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool

def main(argv):
   inputfile = ""
   outputfile = ""
   labelfile = ""
   verbose = False
   nthreads = 1
   database="plasmid_refseq.fasta"
   print("PlasForest: a homology-based random forest classifier for plasmid identification.")
   print("(C) Lea Pradier, Tazzio Tissot, Anna-Sophie Fiston-Lavier, Stephanie Bedhomme. 2020.")
   try:
      opts, args = getopt.getopt(argv,"hi:o:l:v", ["help","input=","output=","label=","threads=","database="])
   except getopt.GetoptError:
      print('./train_plasforest.py -i <inputfile> -o <outputfile> -l <labelfile>')
      print('\t -i, --input <inputfile>: a FASTA input file')
      print('\t -o, --output <outputfile>: a SAV file')
      print('\t -l, --labels <labelfile>: a CSV file')
      print('List of options:')
      print('\t --threads <number of threads>: number of threads (default: 1)')
      print('\t --database <BLAST+ db>: new plasmid database (default: plasmid_refseq.fasta)')
      print('\t -v: verbose mode')
      print('\t -h, --help: show this message and quit')
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('-h','--help'):
         print('./train_plasforest.py -i <inputfile> -o <outputfile> -l <labelfile>')
         print('\t -i, --input <inputfile>: a FASTA input file')
         print('\t -o, --output <outputfile>: a SAV file')
         print('\t -l, --labels <labelfile>: a CSV file')
         print('List of options:')
         print('\t --threads <number of threads>: number of threads (default: 1)')
         print('\t --database <BLAST+ db>: new plasmid database (default: plasmid_refseq.fasta)')
         print('\t -v: verbose mode')
         print('\t -h, --help: show this message and quit')
         sys.exit()
      elif opt in ("-i", "--input"):
         inputfile = arg
         if (not arg.endswith(".fasta")) and (not arg.endswith(".fna")) and (not arg.endswith(".fa")):
            print('./train_plasforest.py -i <inputfile> -o <outputfile> -l <labelfile>')
            print('Error: input file must be in FASTA format.')
            sys.exit()
      elif opt in ("-o", "--output"):
         outputfile=arg
         if not arg.endswith(".sav"):
            print('./train_plasforest.py -i <inputfile> -o <outputfile> -l <labelfile>')
            print('Error: output file is a SAV file.')
            sys.exit()
      elif opt in ("-l", "--labels"):
         labelfile=arg
         if not arg.endswith(".csv"):
            print('./train_plasforest.py -i <inputfile> -o <outputfile> -l <labelfile>')
            print('Error: label file is a CSV file.')
            sys.exit()
      elif opt=="--threads":
         nthreads = arg
      elif opt=="-v":
         verbose = True
      elif opt=="--database":
         database=arg
   if verbose: print("Training PlasForest on "+inputfile+".")
   blast_table = inputfile+"_blast.out"
   features_table = inputfile+"_features.csv"
   list_records = seq_checker(inputfile, verbose)
   blast_launcher(inputfile, blast_table, verbose, nthreads, database)
   features = get_features(list_records, blast_table, verbose, nthreads)
   make_training(features, labelfile, outputfile, verbose, nthreads)
   if(os.path.isfile(inputfile)):
      os.remove(blast_table)
      if verbose: print("Temporary files are deleted.")
   print('The trained classifier is exported to '+outputfile)

### CHECK WHICH SEQUENCES NEED TO BE CHECKED ###
def seq_checker(inputfile, verbose):
   nb_seqs_analyze = 0
   list_records = []
   with open(inputfile, "r+") as handle:
      for record in SeqIO.parse(handle, "fasta"):
         nb_seqs_analyze+=1
         list_records.append(record)
   if verbose: print(str(nb_seqs_analyze)+" contigs will be used to train PlasForest.")
   return list_records

### BLAST LAUNCHER ###
def blast_launcher(inputfile, blast_table, verbose, nthreads,database):
    blastn_cline = NcbiblastnCommandline(query = inputfile, db = "plasmid_refseq.fasta", evalue = 0.001, outfmt = 6, out = blast_table, num_threads=nthreads)
    if verbose: print("Starting BLASTn...")
    stdout, stderr = blastn_cline()
    if verbose: print("BLASTn is over!")

### GET FEATURES ###
def get_features(list_records, blast_table, verbose, nthreads):
    if verbose: print("Computing the features")
    pool = Pool(int(nthreads), read_blast, [blast_table])
    seqFeatures = pd.concat(pool.map(get_seq_feature, list_records))
    blastFeatures = pd.concat(pool.map(get_blast_feature, seqFeatures["ID"].unique().tolist()))
    pool.close()
    pool.join()
    features = pd.merge(seqFeatures, blastFeatures, on="ID")
    features[["Maximum coverage", "Median coverage","Average coverage"]] = features[["Maximum coverage", "Median coverage","Average coverage"]].div(features["Contig size"], axis="index")
    features["Variance of coverage"] = features["Variance of coverage"].div(features["Contig size"], axis="index").div(features["Contig size"], axis="index")
    return(features)

def read_blast(blast_table):
    global BLAST
    BLAST = pd.read_table(blast_table, sep='\t', header=None, names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend', 'sstart','send','evalue','bitscore'])

def get_seq_feature(record):
    return(pd.DataFrame({"ID":[str(record.id)], "G+C content":[float(SeqUtils.GC(record.seq))], "Contig size":[len(str(record.seq))]}))

def get_blast_feature(ID):
    idBLAST = BLAST[BLAST["qseqid"].astype(str)==ID]
    maxCover=0; averCover=0; medianCover=0; varCover=0
    nHits=len(idBLAST.index)
    if nHits > 0:
        maxCover=max(idBLAST["length"])
        averCover = np.mean(idBLAST["length"])
        medianCover = np.median(idBLAST["length"])
        varCover = np.var(idBLAST["length"])
    return(pd.DataFrame({"ID":[ID],"Number of hits":[nHits],"Maximum coverage":[maxCover],"Median coverage":[medianCover],"Average coverage":[averCover],"Variance of coverage":[varCover]}))

### TRAINING ###
def make_training(features, labelfile, outputfile, verbose, nthreads):
	labels=pd.read_csv(labelfile, sep=",")
	training_set = pd.merge(features, labels, on="ID")
	training_features = training_set[["G+C content","Contig size","Number of hits","Maximum coverage","Median coverage","Average coverage","Variance of coverage"]]
	training_labels = training_set[["Plasmid"]]
	if verbose: print("Start the training")
	plasforest = RandomForestClassifier(n_estimators=500, n_jobs=int(nthreads), oob_score=True)
	plasforest.fit(training_features, training_labels)
	if verbose: print("Training is over.")
	pickle.dump(plasforest, open(outputfile,"wb"))

### EXECUTE MAIN ###
if __name__ == "__main__":
   main(sys.argv[1:])




