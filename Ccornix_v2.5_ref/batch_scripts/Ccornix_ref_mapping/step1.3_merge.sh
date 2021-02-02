#! /bin/bash -l
#SBATCH -A b2013182
#SBATCH -p core -n 2
#SBATCH -t 2-00:00:00
#SBATCH -C usage_mail

#job id: $SLURM_JOB_ID
#job name (-J): $SLURM_JOB_NAME
#tmp directory: $SNIC_TMP
#module load bioinfo-tools samtools/1.3

#java heap size: 2GB less than max. memory, i.e. 6GB per core
mem=12g
bams="/proj/b2013182/nobackup/PHYLO/Ccornix_v2.5/processed_bams"
ref="/proj/b2013182/private/assemblies/Ccornix/"

#####
##Copy files from home to tmp directory
#####
cp ${bams}/${1}*mrkdup.ba* ${SNIC_TMP}/
cp ${ref}/genome_HC_allpaths41687_v2.5.fasta ${SNIC_TMP}/
cp ${ref}/genome_HC_allpaths41687_v2.5.fasta.fai ${SNIC_TMP}/
cp ${ref}/genome_HC_allpaths41687_v2.5.dict ${SNIC_TMP}/

# Input files:
# Print input files to log
echo ${1}
# Bam files are mapped, sorted, readgroups added, duplicates marked
# Merge bams for the individual given as input
cd ${SNIC_TMP}
java -Xmx${mem} -jar /sw/apps/bioinfo/picard/1.141/milou/picard.jar MergeSamFiles $(printf ' INPUT=%s' ${1}_[A-Z][0-9]*_sorted_rg_mrkdup.bam) OUTPUT=${1}_merged.bam CREATE_INDEX=true

# Check bam files
java -Xmx${mem} -jar /sw/apps/bioinfo/picard/1.141/milou/picard.jar ValidateSamFile INPUT=${1}_merged.bam OUTPUT=${1}_merged.txt MODE=SUMMARY
#samtools flagstat ${1}_merged.bam > ${1}_merged_stats.txt
#java -Xmx${mem} -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T CountLoci -R $ref/nc_crow-final.assembly.fasta.gapcloser-1.fa -I ${1}_merged.bam -o ${1}_merged_CountLoci.txt
#java -Xmx${mem} -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T DepthOfCoverage -R $ref/nc_crow-final.assembly.fasta.gapcloser-1.fa -I ${1}_merged.bam -o ${1}_merged_DepthOfCoverage.txt

#####
##Copy files from tmp directory to project directory
#####
cp ${SNIC_TMP}/*merged* ${bams}/
