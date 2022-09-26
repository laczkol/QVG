#!/bin/bash

mkdir test/output

apptainer exec qvg.sif QVG.sh -r ./test/ref.fa -samples-list ./test/slist -s ./test -o ./test/output -annot yes -g ./test/ref.gff -np 1

cons=$(find ./test/output/ -name "samples*masked*fa" | sort | head -1)

if [[ -f $cons ]]; then
	echo "#######################"
	echo "Consensus genome found"
	echo -e "#######################\n"
else
	echo "#######################"
	echo "Consensus genome was not found. Please check or download the Apptainer container."
	echo -e "#######################\n"
	exit 1
fi	

if [[ -f ./test/output/non_haploid_sites.txt ]]; then
	echo "#######################"
	echo "Within host variable sites found"
	echo -e "#######################\n"
else
	echo "#######################"
	echo "Within host variable sites were not found. Please check or download the Apptainer container."
	echo -e "#######################\n"
	exit 1
fi

if [[ -f ./test/output/S1/S1.gff ]]; then
	echo "#######################"
	echo "Transferred annotations found"
	echo -e "#######################\n"
else
	echo "#######################"
	echo "Transferred annotations were not found. Please check or download the Apptainer container."
	echo -e "#######################\n"
	exit 1
fi

echo "TEST RUN IS SUCCESSFUL"
