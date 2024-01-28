#!/usr/bin/bash
set -e

#####Software#####
FASTQC="/software/fastqc/0.11.9/fastqc"
FASTP="/software/fastp/fastp"
BWA="/share/app/bwa/0.7.17/bwa"
SAMTOOLS="/software/samtools-1.15.1/samtools"
GATK="/software/gatk-4.2.6.1/gatk"

#####Database#####
REF="/database/gatk-hg38-ref/Homo_sapiens_assembly38.fasta"
REF_INDEX="/database/gatk-hg38-ref/Homo_sapiens_assembly38.fasta"
SNP="/database/gatk-hg38-ref/dbsnp_146.hg38.vcf.gz"
INDELS="/database/gatk-hg38-ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

WORK_DIR=$( pwd )

#Make directory to save the output files.
mkdir ${WORK_DIR}/1_fastqc
mkdir ${WORK_DIR}/2_fastp
mkdir ${WORK_DIR}/3_bwa-align
mkdir ${WORK_DIR}/4_markduplicates
mkdir ${WORK_DIR}/5_bqsr

for SAMPLE in $( ls *1.fq.gz )
do
    SAMPLE_ID=$( basename ${SAMPLE} _1.fq.gz )
    echo $(date) ': Start of fastqc for current sample!'
    cd ${WORK_DIR}/1_fastqc/
    $FASTQC -o ${WORK_DIR}/1_fastqc/ -t 8 ${WORK_DIR}/${SAMPLE_ID}_1.fq.gz ${WORK_DIR}/${SAMPLE_ID}_2.fq.gz \
        1> ${WORK_DIR}/1_fastqc/fastqc.log 2>&1
    echo $(date) ': End of fastqc for current sample --- SUCCESS'

    echo $(date) ': Start of fastp for current sample!'
    cd ${WORK_DIR}/2_fastp/
    $FASTP -i ${WORK_DIR}/${SAMPLE_ID}_1.fq.gz -o ${WORK_DIR}/2_fastp/${SAMPLE_ID}_1_clean.fq.gz \
        -I ${WORK_DIR}/${SAMPLE_ID}_2.fq.gz -O ${WORK_DIR}/2_fastp/${SAMPLE_ID}_2_clean.fq.gz \
        1> ${WORK_DIR}/2_fastp/fastp.log 2>&1

    echo $(date) ': End of fastp for current sample --- SUCCESS'

    echo $(date) ': Start of [BWA-Align] for current sample!'
    echo "Start Align ......"
    cd ${WORK_DIR}/3_bwa-align/
    PLATFORM=$( sed 's/tr>/tr>\n/g' ${WORK_DIR}/1_fastqc/${SAMPLE_ID}_1_fastqc.html | grep Encoding | cut -d " " -f  3 )
    ${BWA} mem -t 8 -M \
        -R "@RG\tID:${SAMPLE_ID}\tPL:${PLATFORM}\tLB:WES_LIB\tSM:${SAMPLE_ID}" ${REF_INDEX} \
        ${WORK_DIR}/2_fastp/${SAMPLE_ID}_1_clean.fq.gz \
        ${WORK_DIR}/2_fastp/${SAMPLE_ID}_2_clean.fq.gz \
        > ${WORK_DIR}/3_bwa-align/${SAMPLE_ID}.sam

    echo $(date) ': Done of [BWA-Align] for current sample --- SUCCESS'
    echo '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

    echo $(date) ': Convert the sam to bam and sorted it......'
    ${SAMTOOLS} view -b -S ${WORK_DIR}/3_bwa-align/${SAMPLE_ID}.sam > ${WORK_DIR}/3_bwa-align/${SAMPLE_ID}.bam
    ${SAMTOOLS} sort ${WORK_DIR}/3_bwa-align/${SAMPLE_ID}.bam -o ${WORK_DIR}/3_bwa-align/${SAMPLE_ID}.sorted.bam
    ${SAMTOOLS} flagstat ${WORK_DIR}/3_bwa-align/${SAMPLE_ID}.sorted.bam > ${WORK_DIR}/3_bwa-align/${SAMPLE_ID}.sorted.bam.flagstat
    echo $(date) ': Samtools Done!'

    echo $(date) ': Mark and remove the PCR duplicates......'
    cd ${WORK_DIR}/4_markduplicates/
    ${GATK} --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates \
        -I ${WORK_DIR}/3_bwa-align/${SAMPLE_ID}.sorted.bam \
        --REMOVE_DUPLICATES true \
        -O ${WORK_DIR}/4_markduplicates/${SAMPLE_ID}.sorted.MarkAndRemoveDuplicates.bam \
        -M ${WORK_DIR}/4_markduplicates/${SAMPLE_ID}.sorted.bam.metrics \
        1> ${WORK_DIR}/4_markduplicates/mark_and_remove_duplicates.log 2>&1

    echo $(date) ': [Samtools index] After Mark and remove the PCR duplicates......'
    ${SAMTOOLS} index -@ 6 -m 4G -b ${WORK_DIR}/4_markduplicates/${SAMPLE_ID}.sorted.MarkAndRemoveDuplicates.bam
    echo $(date) ': Done the mark and remove pcr duplicates. --- SUCCESS'

    echo $(date) ': Start GATK ApplyBQSR......'
    cd ${WORK_DIR}/5_bqsr/
    $GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" BaseRecalibrator \
        -R ${REF} \
        -I ${WORK_DIR}/4_markduplicates/${SAMPLE_ID}.sorted.MarkAndRemoveDuplicates.bam \
        --known-sites ${SNP} \
        --known-sites ${INDELS} \
        -O ${WORK_DIR}/5_bqsr/${SAMPLE_ID}_recal_data.table \
        1> ${WORK_DIR}/5_bqsr/BaseRecalibrator.log 2>&1

    $GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" ApplyBQSR \
        -R ${REF} \
        -I ${WORK_DIR}/4_markduplicates/${SAMPLE_ID}.sorted.MarkAndRemoveDuplicates.bam \
        -bqsr ${WORK_DIR}/5_bqsr/${SAMPLE_ID}_recal_data.table \
        -O ${WORK_DIR}/5_bqsr/${SAMPLE_ID}_bqsr.bam \
        1> ${WORK_DIR}/5_bqsr/ApplyBQSR.log 2>&1
    echo $(date) ': End GATK ApplyBQSR --- SUCCESS'
    cd  ${WORK_DIR}/
    echo $(date) ': All Done!'
done