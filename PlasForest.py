#!/usr/bin/env python3.6
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
from joblib import Parallel, delayed
import multiprocessing

### GLOBALS ###
nHits = []; maxCover = []; averCover = []; medianCover = []; varCover = []; firstHit = []
plasforest = pickle.load(open("plasforest.sav","rb"))
prediction = [];

### MAIN ###
def main(argv):
    inputfile = ""
    outputfile = ""
    firsthits = False
    showFeatures = False
    print("PlasForest: a homology-based random forest classifier for plasmid identification.")
    print("(C) Lea Pradier, Tazzio Tissot, Anna-Sophie Fiston-Lavier, Stephanie Bedhomme. 2020.")
    try:
        opts, args = getopt.getopt(argv,"hi:o:bf", ["help","input=","output="])
    except getopt.GetoptError:
        print('./PlasForest.py -i <inputfile> -o <outputfile>')
        print('\t -i, --input=<inputfile>: a FASTA input file')
        print('\t -o, --output=<outputfile>: a CSV file')
        print('List of options:')
        print('\t -b: show best hit from the plasmid database for each contig')
        print('\t -f: keep the features used by the classifier in the output')
        print('\t -h: show this message and quit')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h','--help'):
            print('./PlasForest.py -i <inputfile> -o <outputfile>')
            print('\t -i, --input=<inputfile>: a FASTA input file')
            print('\t -o, --output=<outputfile>: a CSV file')
            print('List of options:')
            print('\t -b, --besthits: show best hit from the plasmid database for each contig')
            print('\t -f, --features: keep the features used by the classifier in the output')
            print('\t -h: show this message and quit')
            sys.exit()
        elif opt in ("-i", "--input"):
            inputfile = arg
            if (not arg.endswith(".fasta")) and (not arg.endswith(".fna")):
                print('./PlasForest.py -i <inputfile> -o <outputfile>')
                print('Error: input file must be in FASTA format.')
                sys.exit()
        elif opt in ("-o", "--output"):
            outputfile = arg
            if not arg.endswith(".csv"):
                print('./PlasForest.py -i <inputfile> -o <outputfile>')
                print('Error: output file is a column-separated file.')
                sys.exit()
        elif opt=="-b":
            firsthits = True
        elif opt=="-f":
            showFeatures = True
    print("Applying PlasForest on "+inputfile+".")
    tmp_fasta = inputfile+"_tmp.fasta"
    blast_table = inputfile+"_blast.out"
    seq_checker(inputfile, tmp_fasta)
    blast_launcher(tmp_fasta, blast_table)
    features = get_features(tmp_fasta, blast_table, firsthits)
    finalfile = plasforest_predict(features, showFeatures, firsthits)
    finalfile.to_csv(os.path.abspath(outputfile), sep=',', encoding='utf-8', index=False)
    os.remove(tmp_fasta)
    os.remove(blast_table)
    print("Temporary files are deleted.")
    print('Predictions are printed in '+outputfile)

### CHECK WHICH SEQUENCES NEED TO BE CHECKED ###
def seq_checker(inputfile, tmp_fasta):
    nb_seqs_analyze = 0
    words_to_eliminate = ["plasmid", "chromosome", "Plasmid", "Chromosome"]
    with open(inputfile, "r+") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if all(word not in str(record.id) for word in words_to_eliminate) and all(word not in str(record.description) for word in words_to_eliminate):
                nb_seqs_analyze+=1
                with open(tmp_fasta,"a+") as tmp:
                    SeqIO.write(record, tmp, "fasta")
            else:
                nb_seqs_analyze+=1
                with open(tmp_fasta,"a+") as tmp:
                    SeqIO.write(record, tmp, "fasta")
    print(str(nb_seqs_analyze)+" contigs will be analyzed by PlasForest.")

### BLAST LAUNCHER ###
def blast_launcher(tmp_fasta, blast_table):
    nthreads=os.cpu_count()
    blastn_cline = NcbiblastnCommandline(query = tmp_fasta, db = "plasmid_refseq.fasta", evalue = 0.001, outfmt = 6, out = blast_table, num_threads=nthreads)
    print("Starting BLASTn...")
    stdout, stderr = blastn_cline()
    print("BLASTn is over!")

### GET FEATURES ###
def get_features(tmp_fasta, blast_table, firsthits):
    print("Computing the features")
    IDs = []
    GC = []
    cSize = []
    with open(tmp_fasta, "r+") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            IDs.append(str(record.id))
            GC.append(float(SeqUtils.GC(record.seq)))
            cSize.append(len(str(record.seq)))
    nthreads=os.cpu_count()
    Parallel(n_jobs = nthreads, require='sharedmem')(delayed(get_feature)(ID, blast_table, firsthits) for ID in IDs)
    maxCoverage = [float(x)/float(y) for x,y in zip(maxCover, cSize)]
    medianCoverage = [float(x)/float(y) for x,y in zip(medianCover, cSize)]
    averCoverage = [float(x)/float(y) for x,y in zip(averCover, cSize)]
    varCoverage = [float(x)/(float(y)*float(y)) for x,y in zip(varCover, cSize)]
    if firsthits:
        features = pd.DataFrame({"ID":IDs,"G+C content":GC,"Contig size":cSize,"Number of hits":nHits,"Maximum coverage":maxCoverage,"Median coverage":medianCoverage,"Average coverage":averCoverage,"Variance of coverage":varCoverage, "First hit":firstHit})
    else:
        features = pd.DataFrame({"ID":IDs,"G+C content":GC,"Contig size":cSize,"Number of hits":nHits,"Maximum coverage":maxCoverage,"Median coverage":medianCoverage,"Average coverage":averCoverage,"Variance of coverage":varCoverage})
    return(features)

def get_feature(ID,blast_table, firsthits):
    BLAST = pd.read_table(blast_table, sep='\t', header=None, names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend', 'sstart','send','evalue','bitscore'])
    idBLAST = BLAST[BLAST["qseqid"].astype(str)==ID]
    nHits.append(len(idBLAST.index))
    if len(idBLAST.index) > 0:
        maxCover.append(max(idBLAST["length"]))
        averCover.append(np.mean(idBLAST["length"]))
        medianCover.append(np.median(idBLAST["length"]))
        varCover.append(np.var(idBLAST["length"]))
        if firsthits:
            theFirstHit = idBLAST[idBLAST["length"]==max(idBLAST["length"])]
            firstHit.append(theFirstHit["sseqid"].unique().tolist()[0])
    else:
        maxCover.append(0)
        averCover.append(0)
        medianCover.append(0)
        varCover.append(0)
        if firsthits:
            firstHit.append("NA")

### PREDICT WITH PLASFOREST ###
def plasforest_predict(features, showFeatures, firsthits):
    print("Starting predictions with the random forest classifier...")
    temp_table = features[["G+C content","Contig size","Number of hits","Maximum coverage","Median coverage","Average coverage","Variance of coverage"]]
    nthreads=os.cpu_count()
    Parallel(n_jobs = nthreads, require='sharedmem')(delayed(make_prediction)(row) for index, row in temp_table.iterrows())
    wordy_prediction = []
    for p in prediction:
        if p==0: wordy_prediction.append("Chromosome")
        else: wordy_prediction.append("Plasmid")
    if showFeatures:
        features["Prediction"] = pd.Series(wordy_prediction, index=features.index)
        finalfile=features
    elif firsthits:
        finalfile=pd.DataFrame({"ID":features["ID"], "Prediction":wordy_prediction, "Best hit": features["First hit"]})
    else:
        finalfile=pd.DataFrame({"ID":features["ID"], "Prediction":wordy_prediction})
    print("Predictions made!")
    return(finalfile)

def make_prediction(temp_table):
    topredict = temp_table.to_frame().T
    pre = plasforest.predict(topredict)
    prediction.append(pre)


### EXECUTE MAIN ###
if __name__ == "__main__":
   main(sys.argv[1:])


