#! /bin/bash -l
#SBATCH -A b2013182
#SBATCH -p core -n 6
#SBATCH -t 2-00:00:00
#SBATCH -C usage_mail

#job id: $SLURM_JOB_ID
#job name (-J): $SLURM_JOB_NAME
#tmp directory: $SNIC_TMP

#java heap size: 2GB less than max. memory, i.e. 6GB per core
mem=36g
bams="/proj/b2013182/nobackup/PHYLO/Ccornix_v2.5/processed_bams"
ref="/proj/b2013182/private/assemblies/Ccornix"

#####
##Copy files from home to tmp directory
#####
cp ${bams}/${1}* ${SNIC_TMP}/
cp ${ref}/genome_HC_allpaths41687_v2.5.fasta ${SNIC_TMP}/
cp ${ref}/genome_HC_allpaths41687_v2.5.fasta.fai ${SNIC_TMP}/
cp ${ref}/genome_HC_allpaths41687_v2.5.dict ${SNIC_TMP}/

# Realign indels merged bam files from one species combined
cd ${SNIC_TMP}
java -Xmx${mem} -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -R genome_HC_allpaths41687_v2.5.fasta $(printf ' -I %s' ${1}.bam) -o ${1}_raln_targets.list -nt 6
java -Xmx${mem} -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T IndelRealigner -R genome_HC_allpaths41687_v2.5.fasta $(printf ' -I %s' ${1}.bam) -targetIntervals ${1}_raln_targets.list -nWayOut _indraln.bam

#####
##Copy files from tmp directory to project directory
#####
cp ${SNIC_TMP}/*indraln* ${bams}/

