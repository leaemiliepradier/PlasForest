#!/usr/bin/env Rscript
rm(list=ls())
require(optparse)
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  return(script.dir)
}
setwd(getScriptPath())

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="FASTA file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="output.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (FASTA file name).\n", call.=FALSE)
}
if ((tail(strsplit(opt$file, split="[.]")[[1]],1) != "fasta") & (tail(strsplit(opt$file, split="[.]")[[1]],1) != "fna")){
  print_help(opt_parser)
  stop("Input file must be a FASTA file.\n", call.=FALSE)
}

#require(caret)
require(seqinr)

launchBlastn <- function(file, blast.output){
  blast.command = sprintf("./tools/blastn -query %s -out %s -db ./db/plasmids_DB -evalue 0.001 -outfmt '6 qseqid sseqid qstart qend qlen slen evalue length'",
                          file,
                          blast.output)
  print(sprintf("Launch blastn on FASTA file %s", file))
  system(command = blast.command)
}

parseBlastnResults <- function(file, blast.output){
  fasta.file <- read.fasta(file, seqtype = "DNA")
  blast.output.file <- read.table(blast.output, sep="\t")
  seq.IDs = names(fasta.file)
  seq.lengths = c(); for(s in fasta.file){seq.lengths <- c(seq.lengths, length(s))}
  blast.results <- data.frame(ID = seq.IDs, queryNS = NA, ratMaxS = NA, ratMedS = NA, ratMeanS = NA, len = seq.lengths)
  for (id in seq.IDs){
    temp.blast <- blast.output.file[blast.output.file$qseqid == id, ]
    blast.results$queryNS[blast.results$ID==id] <- dim(temp.blast)[1]
    blast.results$ratMaxS[blast.results$ID==id] <- ifelse(dim(temp.blast)[1]>0, max(temp.blast$length)/unique(temp.blast$qlen), 0)
    blast.results$ratMedS[blast.results$ID==id] <- ifelse(dim(temp.blast)[1]>0, median(temp.blast$length)/unique(temp.blast$qlen), 0)
    blast.results$ratMeanS[blast.results$ID==id] <- ifelse(dim(temp.blast)[1]>0, mean(temp.blast$length)/unique(temp.blast$qlen), 0)
  }
  print("Blastn results successfully parsed.")
  return(blast.results)
}

predictPlasmids <- function(blast.results){
  print("Starting prediction of contig identity...")
  PlasForest <- readRDS("tools/plasmids_full_balanced_dataset.rds")
  prediction.plasforest.tmp <- predict(PlasForest,newdata=blast.results)
  prediction.plasforest <- as.numeric(prediction.plasforest.tmp)
  prediction.plasforest[prediction.plasforest == 1] <- "Plasmid"
  prediction.plasforest[prediction.plasforest == 0] <- "Chromosome"
  output.PlasForest <- data.frame(ID = blast.results$ID, prediction = prediction.plasforest)
  return(output.PlasForest)
}

blast.output = paste(tools::file_path_sans_ext(opt$file),".out",sep="")
launchBlastn(opt$file, blast.output)
blast.results <- parseBlastnResults(opt$file, blast.output)
output <- predictPlasmids(blast.results)

write.table(output, opt$out, sep="\t")

print(sprintf("Identity of contigs in file %s has been successfully predicted. Results are displayed in file %s.", opt$file, opt$out))




