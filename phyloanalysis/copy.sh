#!/bin/sh
if [ $# -lt 3 ]
then
	echo "USAGE: $0 <PATIENT> <INPUT> <OUTPUT>"
	exit 1
fi
DONOR=$1
OUT=$3
STEM=$2 

cp $STEM/${DONOR}.cfg $OUT/.

for x in `ls $STEM/snp/*.ANNOT.tsv`
do
	echo "copying $x to ${DONOR}_snp_pileup.tsv"
	cp $x $OUT/${DONOR}_snp_pileup.tsv
done
for x in `ls $STEM/indel/*.ANNOT.tsv`
do
	echo "copying $x to ${DONOR}_indel_pileup.tsv"
	cp $x $OUT/${DONOR}_indel_pileup.tsv
done
cp $STEM/drivers.bed $OUT/${DONOR}_drivers.bed
