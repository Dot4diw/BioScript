#!/usr/bin/bash
set -e
### wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa
### wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz
### wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_indels.vcf.gz

gatk="/software/gatk-4.5.0.0/gatk"
# gzip -d mgp_REL2021_snps.vcf.gz
# gzip -d mgp_REL2021_indels.vcf.gz

echo -n "[`date` - Filter the snp database file --- ] -> "
# get the header of vcf
zcat mgp_REL2021_snps.vcf.gz | head -n 1000 | grep "^#" | cut -f 1-8 \
    > mgp_REL2021_snps.pass.vcf
# keep onely passing site.
zcat mgp_REL2021_snps.vcf.gz | grep -v "^#" | cut -f 1-8 | grep -w "PASS" \
    | awk '{print "chr"$0}'  >> mgp_REL2021_snps.pass.vcf
echo "OK"

echo -n "[`date` - Filter the indels database file --- ] -> "
# get the header of vcf
zcat mgp_REL2021_indels.vcf.gz | head -1000  | grep "^#" | cut -f 1-8 \
    > mgp_REL2021_indels.pass.vcf
# keep onely passing site.
cat mgp_REL2021_indels.vcf.gz | grep -v "^#" | cut -f 1-8 \
    | awk '{print "chr"$0}'  >> mgp_REL2021_indels.pass.vcf
echo "OK"

echo -n "[`date` - Create reference sequence dictionary ... ] -> "
${gatk} CreateSequenceDictionary -R mm39.fa -O mm39.fa.dict
echo "OK"

echo -n "[`date` - Update vcf sequence dictionary --- ] -> "
${gatk} UpdateVcfSequenceDictionary \
    -I mgp_REL2021_snps.pass.vcf \
    -O mgp_REL2021_snps.pass.dict.vcf \
    --SEQUENCE_DICTIONARY mm39.fa.dict
echo "+++++++++++++++++++++done+++++++++++++++++++++"
${gatk} UpdateVcfSequenceDictionary \
    -I mgp_REL2021_indels.pass.vcf \
    -O mgp_REL2021_indels.pass.dict.vcf \
    --SEQUENCE_DICTIONARY mm39.fa.dict
echo "OK"

echo -n "[`date` - Sort vcf --- ] -> "
${gatk} SortVcf \
    -I mgp_REL2021_snps.pass.dict.vcf \
    -O mgp_REL2021_snps.pass.sorted.vcf \
    --SEQUENCE_DICTIONARY mm39.fa.dict

${gatk} SortVcf \
    -I mgp_REL2021_indels.pass.dict.vcf \
    -O mgp_REL2021_indels.pass.sorted.vcf \
    --SEQUENCE_DICTIONARY mm39.fa.dict
echo "OK"

echo -n "[`date` - bgzip files --- ] -> "
bgzip mgp_REL2021_snps.pass.sorted.vcf
bgzip mgp_REL2021_indels.pass.sorted.vcf
echo "OK"

rm *.vcf.idx
rm *.pass.vcf *.dict.vcf

echo "[`date` - Build index of vcf file. --- ] -> "
${gatk} IndexFeatureFile -I mgp_REL2021_snps.pass.sorted.vcf.gz
${gatk} IndexFeatureFile -I mgp_REL2021_indels.pass.sorted.vcf.gz
echo "OK"
echo "[`date` All Done!]"
