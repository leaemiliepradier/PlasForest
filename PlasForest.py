#!/usr/bin/env python

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
from sklearn.ensemble import RandomForestClassifier
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
import time

start_time = time.time()

### GLOBALS ###
global attributed_IDs, attributed_identities
attributed_IDs = []; attributed_identities = [];
plasforest = pickle.load(open("plasforest.sav","rb"))

### MAIN ###
def main(argv):
    inputfile = ""
    outputfile = ""
    besthits = False
    showFeatures = False
    verbose = False
    reattribute = False
    nthreads = 1
    batch=0
    # INPUT OPTIONS #
    print("PlasForest: a homology-based random forest classifier for plasmid identification.")
    print("(C) Lea Pradier, Tazzio Tissot, Anna-Sophie Fiston-Lavier, Stephanie Bedhomme. 2020.")
    if not os.path.isfile("plasforest.sav"):
        print("Error: plasforest.sav is not in the directory.")
        if os.path.isfile("plasforest.sav.tar.gz"):
            print("Please decompress plasforest.sav.tar.gz before launching PlasForest.")
        sys.exit(2)
    if not os.path.isfile("plasmid_refseq.fasta"):
        print("Error: plasmid_refseq.fasta is not in the directory.")
        print("Please execute ./database_downloader.sh before launching PlasForest.")
        sys.exit(2)
    try:
        opts, args = getopt.getopt(argv,"hi:o:bfvr", ["help","input=","output=","threads=","size_of_batch="])
    except getopt.GetoptError:
        print('./PlasForest.py -i <inputfile> -o <outputfile>')
        print('\t -i, --input <inputfile>: a FASTA input file')
        print('\t -o, --output <outputfile>: a CSV file')
        print('List of options:')
        print('\t -b: show best hit from the plasmid database for each contig')
        print('\t -f: keep the features used by the classifier in the output')
        print('\t --threads <int>: number of threads (default: 1)')
        print('\t --size_of_batch <int>: number of sequences per batch')
        print('\t -r: reattribute contigs which are already described as plasmid or chromosome')
        print('\t -v: verbose mode')
        print('\t -h, --help: show this message and quit')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h','--help'):
            print('./PlasForest.py -i <inputfile> -o <outputfile>')
            print('\t -i, --input <inputfile>: a FASTA input file')
            print('\t -o, --output <outputfile>: a CSV file')
            print('List of options:')
            print('\t -b: show best hit from the plasmid database for each contig')
            print('\t -f: keep the features used by the classifier in the output')
            print('\t --threads <int>: number of threads (default: 1)')
            print('\t --size_of_batch <int>: number of sequences per batch')
            print('\t -r: reassign contigs which are already described as plasmid or chromosome')
            print('\t -v: verbose mode')
            print('\t -h, --help: show this message and quit')
            sys.exit()
        elif opt in ("-i", "--input"):
            inputfile = arg
            if (not arg.endswith(".fasta")) and (not arg.endswith(".fna")) and (not arg.endswith(".fa")):
                print('./PlasForest.py -i <inputfile> -o <outputfile>')
                print('Error: input file must be in FASTA format.')
                sys.exit()
            outputfile = inputfile+".csv"
        elif opt in ("-o", "--output"):
            outputfile=arg
            if not arg.endswith(".csv"):
                print('./PlasForest.py -i <inputfile> -o <outputfile>')
                print('Error: output file is a column-separated file.')
                sys.exit()
        elif opt=="--threads":
            nthreads = arg
        elif opt=="-b":
            besthits = True
        elif opt=="-f":
            showFeatures = True
        elif opt=="-r":
            reattribute = True
        elif opt=="-v":
            verbose = True
        elif opt=="--size_of_batch":
            if not arg.isnumeric():
                print("Error: size_of_batch must be integer.")
                sys.exit()
            batch = int(arg)
    if verbose: print("Applying PlasForest on "+inputfile+".")
    tmp_fasta = inputfile+"_tmp.fasta"
    blast_table = inputfile+"_blast.out"
    # THE WHOLE DATASET WILL BE ANALYZED AT ONCE #
    if batch==0:
        list_records = seq_checker(inputfile, tmp_fasta, verbose, reattribute)
        if(os.path.isfile(tmp_fasta)):
            blast_launcher(tmp_fasta, blast_table, verbose, nthreads)
            features = get_features(list_records, blast_table, verbose, nthreads)
            finalfile = plasforest_predict(features, showFeatures, besthits, verbose, attributed_IDs, attributed_identities, nthreads)
            os.remove(tmp_fasta)
            os.remove(blast_table)
        if verbose: print("Temporary files are deleted.")
        else:
            if verbose: print("Contig descriptions already mentioned their identities.")
            finalfile = pd.DataFrame({"ID":attributed_IDs, "Prediction":attributed_identities})
    # THE INPUT DATASET WILL BE ANALYZED IN SEPARATE BATCHES #
    else:
        nb_seqs_to_analyze = seq_counter(inputfile, verbose, batch)
        nb_seqs_analyzed = 0
        finalfile = pd.DataFrame()
        while nb_seqs_analyzed < nb_seqs_to_analyze:
            list_records = seq_checker_batch(inputfile, tmp_fasta, verbose, reattribute, batch, nb_seqs_analyzed)
            nb_seqs_analyzed = nb_seqs_analyzed + len(list_records)
            if(os.path.isfile(tmp_fasta)):
                blast_launcher(tmp_fasta, blast_table, verbose, nthreads)
                features = get_features(list_records, blast_table, verbose, nthreads)
                tmp_finalfile = plasforest_predict(features, showFeatures, besthits, verbose, attributed_IDs, attributed_identities, nthreads)
                os.remove(tmp_fasta)
                os.remove(blast_table)
                finalfile = finalfile.append(tmp_finalfile)
    # THE OUTPUT IS PRINTED TO A CSV FILE #
    finalfile.to_csv(os.path.abspath(outputfile), sep=',', encoding='utf-8', index=False)
    print('Predictions are printed in '+outputfile)
    with open("time_plasforest_"+str(nthreads)+"_threads.dat","a+") as f:
        f.write(str(time.time() - start_time)+"\n")

### COUNT SEQUENCES ###
def seq_counter(inputfile, verbose, batch):
    nb_seqs_analyze = 0
    with open(inputfile, "r+") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            nb_seqs_analyze+=1
    if verbose:
        ntotbatch = nb_seqs_analyze // batch
        if nb_seqs_analyze % batch != 0:
            ntotbatch+=1
        print("PlasForest will analyze the input in ",str(ntotbatch)," different batches.")
    return nb_seqs_analyze

### CHECK WHICH SEQUENCES NEED TO BE CHECKED ###
def seq_checker(inputfile, tmp_fasta, verbose, reattribute):
    nb_seqs_analyze = 0
    list_records = []
    words_to_eliminate = ["plasmid", "chromosome", "Plasmid", "Chromosome"]
    words_plasmid = ["plasmid","Plasmid"]
    with open(inputfile, "r+") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if all(word not in str(record.id) for word in words_to_eliminate) and all(word not in str(record.description) for word in words_to_eliminate):
                nb_seqs_analyze+=1
                list_records.append(record)
                with open(tmp_fasta,"a+") as tmp:
                    SeqIO.write(record, tmp, "fasta")
            else:
                if reattribute:
                    nb_seqs_analyze+=1
                    list_records.append(record)
                    with open(tmp_fasta,"a+") as tmp:
                        SeqIO.write(record, tmp, "fasta")
                else:
                    attributed_IDs.append(str(record.id))
                    if any(word in str(record.id) for word in words_plasmid) or any(word in str(record.description) for word in words_plasmid):
                        attributed_identities.append("Plasmid")
                    else:
                        attributed_identities.append("Chromosome")
                
    if verbose: print(str(nb_seqs_analyze)+" contigs will be analyzed by PlasForest.")
    return list_records

def seq_checker_batch(inputfile, tmp_fasta, verbose, reattribute, batch, nb_already_analyzed):
    nb_seqs_analyze = 0
    list_records = []
    words_to_eliminate = ["plasmid", "chromosome", "Plasmid", "Chromosome"]
    words_plasmid = ["plasmid","Plasmid"]
    with open(inputfile, "r+") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if all(word not in str(record.id) for word in words_to_eliminate) and all(word not in str(record.description) for word in words_to_eliminate):
                nb_seqs_analyze+=1
                if nb_seqs_analyze > nb_already_analyzed and nb_seqs_analyze <= nb_already_analyzed+batch:
                    list_records.append(record)
                    with open(tmp_fasta,"a+") as tmp:
                        SeqIO.write(record, tmp, "fasta")
            else:
                if reattribute:
                    nb_seqs_analyze+=1
                    if nb_seqs_analyze > nb_already_analyzed and nb_seqs_analyze <= nb_already_analyzed+batch:
                        list_records.append(record)
                        with open(tmp_fasta,"a+") as tmp:
                            SeqIO.write(record, tmp, "fasta")
                else:
                    attributed_IDs.append(str(record.id))
                    if any(word in str(record.id) for word in words_plasmid) or any(word in str(record.description) for word in words_plasmid):
                        attributed_identities.append("Plasmid")
                    else:
                        attributed_identities.append("Chromosome")
    if verbose: print("Contigs #"+str(nb_already_analyzed+1)+" to #"+str(min(nb_seqs_analyze,nb_already_analyzed+batch))+" will be analyzed by PlasForest.")
    return list_records

### BLAST LAUNCHER ###
def blast_launcher(tmp_fasta, blast_table, verbose, nthreads):
    blastn_cline = NcbiblastnCommandline(query = tmp_fasta, db = "plasmid_refseq.fasta", evalue = 0.001, outfmt = 6, out = blast_table, num_threads=nthreads)
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
    maxCover=0; averCover=0; medianCover=0; varCover=0; bestHit="NA"
    nHits=len(idBLAST.index)
    if nHits > 0:
        maxCover=max(idBLAST["length"])
        averCover = np.mean(idBLAST["length"])
        medianCover = np.median(idBLAST["length"])
        varCover = np.var(idBLAST["length"])
        theBestHit = idBLAST[idBLAST["length"]==max(idBLAST["length"])]
        bestHit = theBestHit["sseqid"].unique().tolist()[0]
    return(pd.DataFrame({"ID":[ID],"Number of hits":[nHits],"Maximum coverage":[maxCover],"Median coverage":[medianCover],"Average coverage":[averCover],"Variance of coverage":[varCover],"Best hit":[bestHit]}))

### PREDICT WITH PLASFOREST ###
def plasforest_predict(features, showFeatures, besthits, verbose, attributed_IDs, attributed_identities, nthreads):
    if verbose: print("Starting predictions with the random forest classifier...")
    plasforest.n_jobs = int(nthreads)
    features["num_predict"] = plasforest.predict(features[["G+C content","Contig size","Number of hits","Maximum coverage","Median coverage","Average coverage","Variance of coverage"]])
    wordy_prediction = []
    for index, row in features.iterrows():
        if row["num_predict"]==0: wordy_prediction.append("Chromosome")
        else: wordy_prediction.append("Plasmid")
    features["Prediction"] = wordy_prediction
    features = features.drop(columns=["num_predict"])
    if showFeatures:
        finalfile=features
    elif besthits:
        finalfile=features[["ID","Prediction", "Best hit"]]
    else:
        finalfile=features[["ID","Prediction"]]
    if verbose: print("Predictions made!")
    if len(attributed_IDs)>0:
        attributed_df = pd.DataFrame({"ID":attributed_IDs, "Prediction":attributed_identities})
        finalfile = pd.concat([finalfile, attributed_df])
    return(finalfile)

### EXECUTE MAIN ###
if __name__ == "__main__":
   main(sys.argv[1:])

