#! /bin/bash -l
#SBATCH -A b2013182
#SBATCH -p core -n 2
#SBATCH -t 2-00:00:00
#SBATCH -C usage_mail

#job id: $SLURM_JOB_ID
#job name (-J): $SLURM_JOB_NAME
#tmp directory: $SNIC_TMP

#java heap size: 2GB less than max. memory, i.e. 6GB per core
mem=6g
bams="/proj/b2013182/nobackup/PHYLO/Cmoneduloides_v1/processed_bams"
ref="/proj/b2013182/private/assemblies/Cmoneduloides/otago-temp"

######### Usage #########
#bash /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/gatk/step1.3_start_slurm_jobs_Cfrugilegus_Csplendens.sh

#already done
#bash /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/gatk/step1.3_start_slurm_jobs_Cmoneduloides_Cwoodfordi.sh

#####
##Copy files from home to tmp directory
#####
cp $bams/${1}*_sorted_rg_dedup.ba* ${SNIC_TMP}/
cp $ref/nc_crow-final.assembly.fasta.gapcloser-1.fa ${SNIC_TMP}/
cp $ref/nc_crow-final.assembly.fasta.gapcloser-1.fa.fai ${SNIC_TMP}/
cp $ref/nc_crow-final.assembly.fasta.gapcloser-1.dict ${SNIC_TMP}/

#########Input files#########
###Bam files after first round of sorting, indexing and marking duplicates. 16-04-19: shell scripts to start job so far only includes C. moneduloides and C. woodfordi data.

#########merge bams per individual#########
cd ${SNIC_TMP}
java -Xmx${mem} -jar /sw/apps/bioinfo/picard/1.141/milou/picard.jar MergeSamFiles $(printf ' INPUT= %s' ${1}_[A-Z][0-9]*_sorted_rg_dedup.bam) OUTPUT=${1}_merged.bam CREATE_INDEX=true

#########check bam files#########
java -Xmx${mem} -jar /sw/apps/bioinfo/picard/1.141/milou/picard.jar ValidateSamFile INPUT=${1}_merged.bam OUTPUT=${1}_merged.txt MODE=SUMMARY
module load bioinfo-tools samtools/1.3
samtools flagstat ${1}_merged.bam > ${1}_merged_stats.txt
java -Xmx${mem} -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T CountLoci -R $ref/nc_crow-final.assembly.fasta.gapcloser-1.fa -I ${1}_merged.bam -o ${1}_merged_CountLoci.txt
java -Xmx${mem} -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T DepthOfCoverage -R $ref/nc_crow-final.assembly.fasta.gapcloser-1.fa -I ${1}_merged.bam -o ${1}_merged_DepthOfCoverage.txt

#####
##Copy files from tmp directory to project directory
#####
cp ${SNIC_TMP}/*merged* $bams
