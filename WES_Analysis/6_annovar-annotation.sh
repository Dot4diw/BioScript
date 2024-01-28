#!/usr/bin/bash
set -e

ANNOVAR="/software/annovar"
ANNOVAR_DB="/database/ANNOVAR/humandb20211202/"

VARIANT_VCF=$(pwd)"/"$(ls *_somatic-mutations_filter.vcf)

echo "START ANNOVAR ANNOTATION......"  `date`
echo "VARIANT VCF FILE : ${VARIANT_VCF}"
${ANNOVAR}/table_annovar.pl \
    ${VARIANT_VCF} \
    ${ANNOVAR_DB} \
    -buildver hg38 \
    -out  $( basename ${VARIANT_VCF} .vcf )_anno \
    -remove \
    -protocol refGene,knownGene,dbnsfp42a,avsnp150,dbscsnv11,intervar_20180118,cosmic95,gnomad211_exome,clinvar_20210501 \
    -operation g,g,f,f,f,f,f,f,f \
    -nastring . \
    -vcfinput

echo "FINISHED ANNOVAR ANNOTATION!!!"  `date`