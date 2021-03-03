#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -en "${BLUE}Starting to test your PlasForest install...\n"
echo -en "\r${BLUE}Checking if all the files are here."
sleep 1
if [ ! -f PlasForest.py ]; then
    echo -en "\r${RED}ERROR: You did not download the Python script."
    exit
fi

echo -en "\r${BLUE}Checking if all the files are here.."
sleep 0.2
if [ ! -f plasmid_refseq.fasta ]; then
    echo -en "\r${RED}ERROR: You must first download the plasmid database by using database_downloader.sh"
    exit
fi

echo -en "\r${BLUE}Checking if all the files are here..."
sleep 0.2

if [ ! -f plasforest.sav ]; then
    echo -en "\r${RED}ERROR: You did not download or decompress the trained classifier."
    exit
fi

echo -en "\r${BLUE}Checking if all the files are here...."
sleep 0.2

if [ ! -f plasmid_refseq.fasta.nsq ]; then
    echo -en "\r${RED}ERROR: The plasmid database must be compiled to BLAST+."
    exit
fi

echo -en "\r${BLUE}Checking if all the files are here.... ${GREEN}OK\n"

sleep 1
echo -en "\r${BLUE}We now run PlasForest on the test dataset..."
python3.6 PlasForest.py -i test.fasta > silent_plasforest.tmp
rm silent_plasforest.tmp
STATUS="$(cmp --silent test.fasta.csv test_results.csv; echo $?)"
if [[ $STATUS -ne 0 ]]; then
    echo -en "\r${BLUE}We now run PlasForest on the test dataset... ${RED}ERROR"
    exit
fi
echo -en "\r${BLUE}We now run PlasForest on the test dataset... ${GREEN}OK"
sleep 0.5
echo -en "\n${GREEN}The test is successful, you can now use PlasForest.\n${NC}"


