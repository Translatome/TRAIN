#!/bin/bash

# Author: Eric Rivals
# Address: LIRMM, 161, rue Ada, F-34095 Montpellier Cedex 5
# Email: rivals@lirmm.fr
# creation date: 27.02.2023
# License: CeCILL

################################################## DOCUMENTATION
# go_format.sh
# reformat the GO hierarchy file to make an easy to parse file in TSV format
# 2 parameters required input and output filenames

# input: hierarchy in a text formatted file (in GO format)
# output: a TSV formatted file with same information

# example command:
# ./go_format.sh human_BP_MF_CC_concise.txt new_human_BP_MF_CC_concise.txt

##################################################

DEBUG=0

if [ $# -ge 1 ]; then
    gofile=$1
fi

if [ $# -ge 2 ]; then
    outfile=$2
fi

if [ $# -lt 2 ]; then
    echo "Usage:" $0 "GO_filename output_filename"
    echo "        2 compulsory arguments required; not enough try again"
    exit 1
fi

tmpdir="mytmp"
mkdir $tmpdir

# fetch the ids and isolate them
cut -c 1-31 ${gofile} | awk '{print $NF}' > ${tmpdir}/"id.txt"

# fetch intermediate fields
cut -c 32-71 ${gofile} | awk '{for (i=1; i< NF; i++){printf("%s\t", $i);} printf("%s\n", $NF)}' >  ${tmpdir}/"inter.txt"

# fetch description field and surround it by double quotes
cut -c 72- ${gofile} | awk '{printf("\"%s\"\n", $0)}' >  ${tmpdir}/"desc.txt"


# paste all together separated by tabulation
pushd ${tmpdir} >/dev/null

paste -d '\t' id.txt inter.txt  desc.txt > ../${outfile}

popd >/dev/null

if [ $DEBUG -ge 1 ]; then
    exit 0;
else
    rm -r ${tmpdir}
fi
