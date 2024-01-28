#!/usr/bin/bash
set -e

N_ID=$( basename $(pwd) | cut -d '_' -f3 )
T_ID=$( basename $(pwd) | cut -d '_' -f1 )

grep -v '^Chr' *.hg38_multianno.txt | awk -F "\t" '{print $0"\t""'$T_ID'""\t""'$N_ID'"}' > tmp-mutatanno.txt

head -n 1 *hg38_multianno.txt | awk -F "\t" '{print $0"\t""Tumor_Sample_Barcode""\t""Matched_Norm_Sample_Barcode"}' > table-head.txt

cat table-head.txt tmp-mutatanno.txt > ${T_ID}_vs_${N_ID}_mutatanno.txt