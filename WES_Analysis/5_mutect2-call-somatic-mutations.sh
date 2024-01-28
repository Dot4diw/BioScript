#!/usr/bin/bash
set -e

GATK="/software/gatk-4.2.6.1/gatk"

REF="/database/gatk-hg38-ref/Homo_sapiens_assembly38.fasta"
BED="/database/gatk-hg38-ref/hg38.exon.bed"

N_ID=$( basename $(pwd) | cut -d '_' -f3 )
T_ID=$( basename $(pwd) | cut -d '_' -f1 )

N_SAMPLE="/data/0_analysis/"${N_ID}"/5_bqsr/"${N_ID}"_bqsr.bam"
T_SAMPLE="/data/0_analysis/"${T_ID}"/5_bqsr/"${T_ID}"_bqsr.bam"


echo "Start Mutect2 : "  `date`

${GATK}  --java-options "-Xmx26G -Djava.io.tmpdir=./"  Mutect2 -R ${REF} \
    -I ${T_SAMPLE} -tumor ${T_ID} \
    -I ${N_SAMPLE} -normal ${N_ID} \
    -L ${BED} \
    -O ./${T_ID}_vs_${N_ID}_mutect2.vcf

${GATK}  FilterMutectCalls \
    -R ${REF} \
    -V ./${T_ID}_vs_${N_ID}_mutect2.vcf \
    -O ./${T_ID}_vs_${N_ID}_somatic-mutations.vcf

echo "End Mutect2 : " `date`

cat ./${T_ID}_vs_${N_ID}_somatic-mutations.vcf \
    | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }' \
    > ./${T_ID}_vs_${N_ID}_somatic-mutations_filter.vcf

echo "End of All : " `date`
echo "All Done!"