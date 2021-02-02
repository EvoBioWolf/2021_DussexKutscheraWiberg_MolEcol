#! /bin/bash -l
#SBATCH -A b2013182
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -C usage_mail
#SBATCH -J bamstats

#job id: $SLURM_JOB_ID
#job name (-J): $SLURM_JOB_NAME
#tmp directory: $SNIC_TMP

module load bioinfo-tools 
module load samtools/1.3

# Java heap size: 2GB less than max. memory, i.e. 6GB per core
mem=6g

# Variables from command line: directory of bam files
ending=${1}
handle=${2}
dir=${3}

# Copy files from home to tmp directory
cp ${dir}/*${ending} ${SNIC_TMP}/
# Move to tmp directory
cd ${SNIC_TMP}

# Run samtools flagstat for each bamfile.
for bam in $(ls *.bam); do echo ${bam} && samtools flagstat ${bam}; done >> ${handle}_stats.txt

# Run picard validate same for each bamfile.
for bam in $(ls *.bam | sed 's/.bam//g'); do java -Xmx${mem} -jar /sw/apps/bioinfo/picard/1.141/milou/picard.jar ValidateSamFile INPUT=${bam}.bam OUTPUT=${bam}_validateSamFile.txt MODE=SUMMARY; done

##Copy files from tmp directory to project directory
#####
cp ${SNIC_TMP}/${handle}_stats.txt ~/
cp ${SNIC_TMP}/*validateSamFile.txt ${dir}/
