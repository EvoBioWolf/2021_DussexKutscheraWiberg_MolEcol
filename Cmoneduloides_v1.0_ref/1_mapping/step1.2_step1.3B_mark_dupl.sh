#! /bin/bash -l
#SBATCH -A b2013182
#SBATCH -p core -n 2
#SBATCH -t 1-00:00:00
#SBATCH -C usage_mail

#job id: $SLURM_JOB_ID
#job name (-J): $SLURM_JOB_NAME
#tmp directory: $SNIC_TMP

#java heap size: 2GB less than max. memory, i.e. 6GB per core
mem=12g

bams="/proj/b2013182/nobackup/PHYLO/Cmoneduloides_v1/processed_bams"
ref="/proj/b2013182/private/assemblies/Cmoneduloides/otago-temp"

#########Input files#########
#Bam files incl. read groups, sorted plus indexed

######### Usage #########
#cd /proj/b2013182/nobackup/PHYLO/Cmoneduloides_v1/processed_bams

###step 1.2###
#for f in $(ls *sorted_rg.bam | grep Cfru | sed 's/.bam//g'); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.2_step1.3B_mark_dupl.sh $f; done
#for f in $(ls *sorted_rg.bam | grep Csplendens | sed 's/.bam//g'); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.2_step1.3B_mark_dupl.sh $f; done

#already done
#for f in $(ls *sorted_rg.bam | grep Cmon_ | sed 's/.bam//g'); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/gatk/step1.2_step1.3B_mark_dupl.sh $f; done
#for f in $(ls *sorted_rg.bam | grep Cwoo | sed 's/.bam//g'); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/gatk/step1.2_step1.3B_mark_dupl.sh $f; done


#####
##Copy files from home to tmp directory
#####
cp $bams/${1}.ba* ${SNIC_TMP}/
cp $ref/nc_crow-final.assembly.fasta.gapcloser-1.fa ${SNIC_TMP}/
cp $ref/nc_crow-final.assembly.fasta.gapcloser-1.fa.fai ${SNIC_TMP}/
cp $ref/nc_crow-final.assembly.fasta.gapcloser-1.dict ${SNIC_TMP}/

#########mark duplicates and index it#########
cd ${SNIC_TMP}
java -Xmx${mem} -jar /sw/apps/bioinfo/picard/1.141/milou/picard.jar MarkDuplicates INPUT=${1}.bam OUTPUT=${1}_dedup.bam METRICS_FILE=metrics.txt CREATE_INDEX=true

#########check if bam files are intact#########
java -Xmx${mem} -jar /sw/apps/bioinfo/picard/1.141/milou/picard.jar ValidateSamFile INPUT=${1}_dedup.bam OUTPUT=${1}_dedup_validateSamFile.txt MODE=SUMMARY

module load bioinfo-tools samtools/1.3
samtools flagstat ${1}_dedup.bam > ${1}_dedup_stats.txt

java -Xmx${mem} -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T CountLoci -R $ref/nc_crow-final.assembly.fasta.gapcloser-1.fa -I ${1}_dedup.bam -o ${1}_dedup_CountLoci.txt
java -Xmx${mem} -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T DepthOfCoverage -R $ref/nc_crow-final.assembly.fasta.gapcloser-1.fa -I ${1}_dedup.bam -o ${1}_dedup_DepthOfCoverage.txt

#####
##Copy files from tmp directory to project directory
#####
cp ${SNIC_TMP}/*dedup* $bams
