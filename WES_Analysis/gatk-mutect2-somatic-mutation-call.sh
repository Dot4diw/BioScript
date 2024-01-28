#!/usr/bin/bash
set -e

#####Software#####
FASTQC="/share/fastqc/0.11.9/fastqc"
FASTP="/software/fastp/fastp"
BWA="/share/app/bwa/0.7.17/bwa"
SAMTOOLS="/software/samtools-1.15.1/samtools"
GATK="/software/gatk-4.2.6.1/gatk"
ANNOVAR="/software/annovar"

#####Database#####
REF="/database/gatk-hg38-ref/Homo_sapiens_assembly38.fasta"
BED="/database/gatk-hg38-ref/hg38.exon.bed"
REF_INDEX="/database/gatk-hg38-ref/Homo_sapiens_assembly38.fasta"
SNP="/database/gatk-hg38-ref/dbsnp_146.hg38.vcf.gz"
INDELS="/database/gatk-hg38-ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
ANNOVAR_DB="database/ANNOVAR/humandb20211202/"

#######INPUT FQ FILE PATH CONTENT##########
#EXAMPLE FILE NAME fastq_file_path.txt#####
######     /path1/S1_1.fq.gz      #########
######     /path1/S1_2.fq.gz      #########
######     /path1/S2_1.fq.gz      #########
######     /path1/S2_2.fq.gz      #########
###########################################

######INPUT TN SAMPLE ID FILE CONTENT######
######EXAMPLE FILE NAME tn_pariid.txt######
######      SSSSS1_vs_SSS-N        ########
######      SSSSS2_vs_SSS-N        ########
######      SSSSS3_vs_SSS-N        ########
###########################################

# FASTQC QUAILTY CONTORL, BWA ALIGN, MARKDUPLICATE and APPLYBQSR
FASTQ2BQSRBAM(){
    WORK_DIR=`pwd`
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

        # Clean the big output files of BWA.
        echo 'Clean the big output files of BWA-Align.'
        rm ${WORK_DIR}/3_bwa-align/${SAMPLE_ID}.sam ${WORK_DIR}/3_bwa-align/${SAMPLE_ID}.bam

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
}

# MUTECT2 CALL THE MUTATION AND ANNOTATION BY ANNOVAR
MUTAT_CALL_AND_ANNOTATION(){
    N_ID=$( basename $(pwd) | cut -d '_' -f3 )
    T_ID=$( basename $(pwd) | cut -d '_' -f1 )

    N_SAMPLE="/data/0_analysis/"${N_ID}"/5_bqsr/"${N_ID}"_bqsr.bam"
    T_SAMPLE="/data/0_analysis/"${T_ID}"/5_bqsr/"${T_ID}"_bqsr.bam"


    echo "Start Mutect2 : "  `date`
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "Normal : ${N_SAMPLE}"
    echo "Tumor  : ${T_SAMPLE}"
    echo "Starting mutation call : ${T_ID} vs ${N_ID}"
    echo

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
    #!/usr/bin/bash
    set -e

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

    # ADD THE SAMPLE BARCODES TO ANNOTATION RESULT FILES
    grep -v '^Chr' *.hg38_multianno.txt | awk -F "\t" '{print $0"\t""'$T_ID'""\t""'$N_ID'"}' > tmp-mutatanno.txt

    head -n 1 *hg38_multianno.txt | awk -F "\t" '{print $0"\t""Tumor_Sample_Barcode""\t""Matched_Norm_Sample_Barcode"}' > table-head.txt

    cat table-head.txt tmp-mutatanno.txt > ${T_ID}_vs_${N_ID}_mutatanno.txt
}

echo "GATK MUTECT2 PIPELINE STARTING......"  `date`
PROJECTED_DIR=`pwd`

for SAMPEPATH in $( cat fastq_file_path.txt )
do
    SAMPEID=$( basename $SAMPLEPATH | cut -d "_" -f1 )

    if [ ! -d $SAMPLEID ]; then
        mkdir $N_ID
    else
        :
    fi

    cd $PROJECTED_DIR/$SAMPLEID
    FASTQ2BQSRBAM &

    cd $PROJECTED_DIR
done
wait

cd $PROJECTED_DIR
for VSID in $( cat tn_pariid.txt )
do
    if [ ! -d $VSID ]; then
        mkdir $VSID
    else
        echo "${VSID} EXISTED"
    fi

    cd $PROJECTED_DIR/$VSID
    MUTAT_CALL_AND_ANNOTATION &

    cd $PROJECTED_DIR
done
echo "GATK MUTECT2 PIPELINE ALL DONE!!!"  `date`
