# PlasForest
A random forest classifier of contigs to identify contigs of plasmid origin in contig and scaffold genomes

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
$ cd PlasForest-master/
$ tar -zxvf plasforest.sav.tar.gz
````
The use of PlasForest requires to download a database of plasmid sequences (2.5GB). This script should be launched from the directory in which ```PlasForest.py``` and ```plasforest.sav``` are stored. Be careful this step may take a while.
````
$ chmod 755 database_downloader.sh
$ ./database_downloader.sh
````
### Testing the installation
To test that PlasForest has been correctly installed, you can run the pipeline on the test dataset:
````
$ python3 PlasForest.py -i test.fasta
````
which should output ```test.fasta.csv``` without error.

## Manual
The PlasForest pipeline is able to process 2000 contigs (~100 genomes) in about 20 minutes on a computer with 6 CPUs and 16GB of RAM.
### Minimal use
PlasForest requires at least an input FASTA file, and generates an column-separated output file.
````
$ python3 PlasForest.py -i /path/to/your/inputfile.fasta
````
You can also specify the name of the output file.
````
$ python3 PlasForest.py -i /path/to/your/inputfile.fasta -o /path/to/your/outputfile.csv
````

### Additional options
Option ```-b``` will add a third column to the output file, with the identification number of the best hit in the plasmid database.

Option ```-f``` will add seven columns to the output file, with the features measured for each contig.

Option ```-v``` can activate verbose mode.

Option ```-r``` will allow to re-assign contigs that are already described as plasmid or chromosome.

Option ```--threads``` allows to define the number of CPUs on which PlasForest will be run.

# Citation
If you use PlasForest for your research, please cite the following papers:
- Pradier L, Tissot T, Fiston-Lavier AS, Bedhomme S (2020). PlasForest: a homology-based random forest classifier for plasmid detection in genomic datasets. Available from https://github.com/leaemiliepradier/PlasForest
- Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL (2009) BLAST+: architecture and applications. BMC Bioinformatics 10:421
- Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Freidberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25:1422-1423

## Support
Any issues connected with PlasForest should be addressed to Léa Pradier (lea.pradier (at) cefe.cnrs.fr).
