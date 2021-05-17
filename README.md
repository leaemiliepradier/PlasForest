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
$ cd PlasForest/
$ tar -zxvf plasforest.sav.tar.gz
````
The use of PlasForest requires to download a database of plasmid sequences (2.5GB). This script should be launched from the directory in which ```PlasForest.py``` and ```plasforest.sav``` are stored. Be careful this step may take a while.
````
$ chmod 755 database_downloader.sh
$ ./database_downloader.sh
# or alternatively
$ python check_and_download_database.py download
````
### Testing the installation
To test that PlasForest has been correctly installed, you can run the following script:
````
$ ./test_plasforest.sh
````
which will test if all required files are here and will run PlasForest on a test dataset.

## Manual
The PlasForest pipeline is able to process 1000 contigs (~50 genomes) in about 8 minutes using 1 CPU and 2.3 GB of memory. Using 16 CPUs, it will take only 2 minutes, but will require 20 GB of memory.
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

Option ```--threads <int>``` allows to define the number of CPUs on which PlasForest will be run.

Option ```--size_of_batch <int>``` alows to define ow many sequences can be used in the single batch.

### Training PlasForest on custom data
The classifier is already provided for PlasForest and there is no need to train it again. However, if you want to train PlasForest on a custom training dataset and/or using a custom plasmid database, you can do it by using ```train_plasforest.py```.

This pipeline requires you to provide an input FASTA file for your training set, a CSV file containing the labels for each sequence of the training set, and the output name for your custom classifier:
````
$ python3 train_plasforest.py -i /path/to/your/trainingset.fasta -l /path/to/your/labelfile.csv -o /path/to/your/output_classifier.sav
````
Make sure that your label file has two columns: one named "ID" containing all the identifiers of your training set sequences; and the other names "Plasmid" taking the values 0 if a sequence comes from a chromosome and 1 if it comes from a plasmid.

Optionally, you can also
- run in verbose mode with option ```-v```
- specify a custom plasmid database with option ```--database <BLAST+ db>```
- define the number of CPUs on which the training will be run with option ```--threads <int>```

# Citation
If you use PlasForest for your research, please cite the following papers:
- Pradier L, Tissot T, Fiston-Lavier AS, Bedhomme S (2020). PlasForest: a homology-based random forest classifier for plasmid detection in genomic datasets. bioRxiv 2020.10.05.326819; doi: https://doi.org/10.1101/2020.10.05.326819
- Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL (2009) BLAST+: architecture and applications. BMC Bioinformatics 10:421
- Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Freidberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25:1422-1423

## Support
Any issues connected with PlasForest should be addressed to Léa Pradier (lea.pradier (at) cefe.cnrs.fr).
