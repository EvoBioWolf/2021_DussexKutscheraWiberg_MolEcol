#! /bin/bash -l
#SBATCH -A b2013182
#SBATCH -p core -n 2
#SBATCH -t 0-18:00:00
#SBATCH -C usage_mail

#job id: $SLURM_JOB_ID
#job name (-J): $SLURM_JOB_NAME
#tmp directory: $SNIC_TMP

#java heap size: 2GB less than max. memory, i.e. 6GB per core
mem=12g

bams="/proj/b2013182/nobackup/PHYLO/Ccornix_v2.5/processed_bams"
ref="/proj/b2013182/private/assemblies/Ccornix"

#####
##Copy files from home to tmp directory
#####
cp ${bams}/${1}* ${SNIC_TMP}/
cp ${ref}/genome_HC_allpaths41687_v2.5.fasta ${SNIC_TMP}/
cp ${ref}/genome_HC_allpaths41687_v2.5.fasta.fai ${SNIC_TMP}/
cp ${ref}/genome_HC_allpaths41687_v2.5.dict ${SNIC_TMP}/

# Mark duplicates and index it
cd ${SNIC_TMP}
java -Xmx${mem} -jar /sw/apps/bioinfo/picard/1.141/milou/picard.jar MarkDuplicates INPUT=${1}.bam OUTPUT=${1}_mrkdup.bam METRICS_FILE=metrics.txt CREATE_INDEX=true

# Check if bam files are intact
java -Xmx${mem} -jar /sw/apps/bioinfo/picard/1.141/milou/picard.jar ValidateSamFile INPUT=${1}_mrkdup.bam OUTPUT=${1}_mrkdup_validateSamFile.txt MODE=SUMMARY
#
#module load bioinfo-tools samtools/1.3
#samtools flagstat ${1}_mrkdup.bam > ${1}_mrkdup_stats.txt
#java -Xmx${mem} -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T CountLoci -R genome_HC_allpaths41687_v2.5.fasta -I ${1}_mrkdup.bam -o ${1}_mrkdup_CountLoci.txt
#java -Xmx${mem} -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T DepthOfCoverage -R genome_HC_allpaths41687_v2.5.fasta -I ${1}_mrkdup.bam -o ${1}_mrkdup_DepthOfCoverage.txt

#####
##Copy files from tmp directory to project directory
#####
cp ${SNIC_TMP}/*mrkdup* ${bams}/
