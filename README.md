# PlasForest
A random forest classifier of contigs to identify plasmids in draft genomes and metagenomes

## Dependencies
- R 3.5+
and the following R packages
  - getopt
  - seqinr

## Install
Clone this git on your computer:
````
$ git clone https://github.com/leaemiliepradier/PlasForest
````

## Minimal use
PlasForest requires at least an input FASTA file, and generates an output file named `output.txt`.
````
$ Rscript --slave PlasForest.R --input /path/to/your/file.fasta
````
You can also specify the path to your desired output file.
````
$ Rscript --slave PlasForest.R --input /path/to/your/file.fasta --output /path/to/your/output.txt
````
