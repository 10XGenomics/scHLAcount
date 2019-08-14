#!/bin/bash

HLANUC="hla_nuc.fasta"
HLAGEN="hla_gen.fasta"

GTFILE=$1

while read -r line; do grep -F -m 1 $line $HLANUC >> tmpallele.txt; done < $GTFILE

samtools faidx  $HLANUC $(cut -f1 -d' ' tmpallele.txt | tr '>' ' ' | tr '\n' ' ') > cds.fasta
samtools faidx  $HLAGEN $(cut -f1 -d' ' tmpallele.txt | tr '>' ' ' | tr '\n' ' ') > gen.fasta

while read -r line; do IFS=' ' read -r f1 f2 <<<"$line"; sed -i "" "s/$f1/$f1 $f2/g" cds.fasta; done < tmpallele.txt
while read -r line; do IFS=' ' read -r f1 f2 f3 <<<"$line"; sed -i "" "s/$f1/$f1 $f2/g" gen.fasta; done < tmpallele.txt

rm tmpallele.txt

