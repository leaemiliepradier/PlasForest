# PlasForest
A random forest classifier of contigs to identify contigs of plasmid origin in draft genomes

## Dependencies
- Python 3.6+
- Python libraries:
  - BioPython
  - scikit-learn 0.22.2.post1
  - Numpy
  - Pandas
  - Joblib
- NCBI Blast+

## Install
Clone this git on your computer:
````
$ git clone https://github.com/leaemiliepradier/PlasForest
````
Untar the random forest classifier:
````
$ tar -xzvf plasforest.sav.tar.gz
````
The use of PlasForest requires to download a database of plasmid sequences (2.5GB). This script should be launched from the directory in which PlasForest is stored:
````
$ ./database_downloader.sh
````

## Minimal use
PlasForest requires at least an input FASTA file, and generates an column-separated output file.
````
$ python PlasForest.py -i /path/to/your/inputfile.fasta
````
You can also specify the name of the output file.
````
$ python PlasForest.py -i /path/to/your/inputfile.fasta -o /path/to/your/outputfile.csv
````

## Additional options
Option ```-b``` will add a third column to the output file, with the identification number of the best hit in the plasmid database.

Option ```-f``` will add seven columns to the output file, with the features measured for each contig.
