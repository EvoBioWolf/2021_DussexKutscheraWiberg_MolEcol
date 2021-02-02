#! /bin/bash -l
#SBATCH -A b2013182
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0-04:00:00
#SBATCH -C usage_mail
#SBATCH -J readgroup_fix

#job id: $SLURM_JOB_ID
#job name (-J): $SLURM_JOB_NAME
#tmp directory: $SNIC_TMP

module load bioinfo-tools 
module load samtools/1.3

# Variables from command line: filename without the .bam suffix
file=${1}

# Java heap size: 2GB less than max. memory, i.e. 6GB per core
mem=6g

# Directories
#raw_bams="/proj/b2013182/private/bams/Ccornix_ref/bams"
raw_bams="/proj/b2013182/nobackup/POPseq/vcf/raw_bams/Ccornix_v2.5"
proc_bams="/proj/b2013182/nobackup/PHYLO/Ccornix_v2.5/processed_bams"

# Lookup table for read group ids. Bam file name without file extension. 
# 16-04-18: includes 5 C. moneduloides and 5 C. woodfordi individuals. Renamed lanes according to flow cell ids used.
table="/proj/b2013182/nobackup/PHYLO/Ccornix_v2.5/reads_table_readgroups_info.csv"

# Copy files from home to tmp directory
cp $raw_bams/${file}.bam ${SNIC_TMP}/
# Move to tmp directory
cd ${SNIC_TMP}
#------------------------------------------------------------#
#  Read group format                                         #
# @RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1 #
#------------------------------------------------------------#

# Sort (Picard), add read group information and index the file

# RGID (String)	Read Group ID Default value: 1. This option can be set to 'null' to clear the default value. Be sure to change from default of 1! 
# Assigned consecutive number for each sample (one number per pair) starting at 1001. Lookup table: column 2.
rgid=`awk '{FS=OFS=","}{if ($1 == LOOKUP) print $2 }' LOOKUP=${file} ${table}`

#RGSM (String)	Read Group sample name Required. Lookup table: column 4.
rgsm=`awk '{FS=OFS=","}{if ($1 == LOOKUP) print $4 }' LOOKUP=${file} ${table}`

#RGLB (String)	Read Group Library Required. Lookup table: column 6.
rglb=`awk '{FS=OFS=","}{if ($1 == LOOKUP) print $6 }' LOOKUP=${file} ${table}`

#RGPU (String)	Read Group platform unit (eg. run barcode) Required. RGPU = RGSM.RGID_date_RGLB. Lookup table: column 10.
rgpu=`awk '{FS=OFS=","}{if ($1 == LOOKUP) print $10 }' LOOKUP=${file} ${table}`

java -Xmx${mem} -jar /sw/apps/bioinfo/picard/1.141/milou/picard.jar AddOrReplaceReadGroups INPUT=${file}.bam OUTPUT=${file}_rg.bam RGID=$rgid RGSM=$rgsm RGPL=illumina RGLB=$rglb RGPU=$rgpu SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

# Check bam files
#java -Xmx${mem} -jar /sw/apps/bioinfo/picard/1.141/milou/picard.jar ValidateSamFile INPUT=${file}_rg.bam OUTPUT=${file}_rg_validateSamFile.txt MODE=SUMMARY
#samtools flagstat ${file}_rg.bam > ${file}_rg_stats.txt

# Check if read groups have been added successfully
echo ${file}
awk '{FS=","}{OFS=" "}{if ($1 == LOOKUPVAL) print $1,$2,$4,$6,$9,$10 }' LOOKUPVAL=${file} ${table}
samtools view -H ${file}_rg.bam | grep '^@RG'

#####
##Copy files from tmp directory to project directory
#####
cp ${SNIC_TMP}/*rg* $proc_bams
