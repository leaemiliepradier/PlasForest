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
Python dependencies can be manually installed using ```pip```:
````
pip install --user biopython numpy pandas joblib scikit-learn==0.22.2.post1
````
on Ubuntu/Debian, NCBI Blast+ can be installed through the repositories with the command ```sudo apt-get install ncbi-blast+```. For other systems, it is advised to download the binaries from the NCBI FTP repository:
````
# Download latest BLAST version
# For systems other than Linux, check the FTP for a compatible version
$ wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.1+-x64-linux.tar.gz
# Decompress
$ tar -zxvf ncbi-blast-2.10.1+-x64-linux.tar.gz
# Add BLAST location to system PATH
$ export PATH=$HOME/ncbi-blast-2.10.1+/bin:$PATH
````

Untar the random forest classifier:
````
$ cd PlasForest/
$ tar -zxvf plasforest.sav.tar.gz
````
The use of PlasForest requires to download a database of plasmid sequences (2.5GB). This script should be launched from the directory in which ```PlasForest.py``` and ```plasforest.sav``` are stored:
````
$ ./database_downloader.sh
````
## Manual
### Minimal use
PlasForest requires at least an input FASTA file, and generates an column-separated output file.
````
$ python PlasForest.py -i /path/to/your/inputfile.fasta
````
You can also specify the name of the output file.
````
$ python PlasForest.py -i /path/to/your/inputfile.fasta -o /path/to/your/outputfile.csv
````

### Additional options
Option ```-b``` will add a third column to the output file, with the identification number of the best hit in the plasmid database.

Option ```-f``` will add seven columns to the output file, with the features measured for each contig.

Option ```-v```can activate verbose mode.
