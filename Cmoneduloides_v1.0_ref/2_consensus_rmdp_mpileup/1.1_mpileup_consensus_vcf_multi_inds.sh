#!/bin/bash -l
#SBATCH -A b2013182
#SBATCH -p core -n 3
#SBATCH -t 10-00:00:00
#SBATCH -C usage_mail

#job id: $SLURM_JOB_ID
#job name (-J): $SLURM_JOB_NAME
#tmp directory: $SNIC_TMP

#Usage: loop through bam lists for each set of bams per consensus (/proj/b2013182/nobackup/PHYLO/Cmoneduloides_v1/processed_bams)
#cd /proj/b2013182/nobackup/PHYLO/Cmoneduloides_v1/processed_bams
#for f in $(ls Cmone_inds_bam.list); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/consensus_rmdp_mpileup/1.1_mpileup_consensus_vcf_multi_inds.sh $f; done
# for f in $(ls Cbrach_inds_bam.list); do sbatch  /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/consensus_rmdp_mpileup/1.1_mpileup_consensus_vcf_multi_inds.sh $f; done
# for f in $(ls Ccorone_3sp_inds_bam.list); do sbatch  /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/consensus_rmdp_mpileup/1.1_mpileup_consensus_vcf_multi_inds.sh $f; done
# for f in $(ls Cdau_inds_bam.list); do sbatch  /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/consensus_rmdp_mpileup/1.1_mpileup_consensus_vcf_multi_inds.sh $f; done
# for f in $(ls Cfru_inds_bam.list); do sbatch  /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/consensus_rmdp_mpileup/1.1_mpileup_consensus_vcf_multi_inds.sh $f; done
# for f in $(ls Cmon_inds_bam.list); do sbatch  /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/consensus_rmdp_mpileup/1.1_mpileup_consensus_vcf_multi_inds.sh $f; done
# for f in $(ls Csplendens_inds_bam.list); do sbatch  /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/consensus_rmdp_mpileup/1.1_mpileup_consensus_vcf_multi_inds.sh $f; done
# for f in $(ls Cwoo_inds_bam.list); do sbatch  /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/consensus_rmdp_mpileup/1.1_mpileup_consensus_vcf_multi_inds.sh $f; done
# for f in $(ls Cancestral_inds_bam.list); do sbatch  /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/consensus_rmdp_mpileup/1.1_mpileup_consensus_vcf_multi_inds.sh $f; done
# for f in $(ls Ckub_inds_bam.list); do sbatch  /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/consensus_rmdp_mpileup/1.1_mpileup_consensus_vcf_multi_inds.sh $f; done



#bam list files to loop through
#/proj/b2013182/nobackup/PHYLO/Cmoneduloides_v1/processed_bams/Cmone_inds_bam.list


# Copy files from home to tmp directory
for f in `cat /proj/b2013182/nobackup/PHYLO/Cmoneduloides_v1/processed_bams/${1}`; do cp ${f} ${SNIC_TMP}/; done
cp /proj/b2013182/private/assemblies/Cmoneduloides/otago-temp/nc_crow-final.assembly.fasta.gapcloser-1* ${SNIC_TMP}/



# Get species ID from bam list files as variable for later steps
species=`echo ${1} | cut -d'_' -f1`

# Move to tmp/ directory
cd ${SNIC_TMP}

# Generate a new bam file list in the tmp directory, containing only the bam file names
for f in `cat /proj/b2013182/nobackup/PHYLO/Cmoneduloides_v1/processed_bams/${1}`; do basename ${f} >> short_${1}; done

module load bioinfo-tools samtools/1.3 bcftools/1.3 tabix/0.2.6 vcftools/0.1.14

# Make consensus (mpileup -g = bcf output; -t DP = print depth in the output bcf file; -u = uncompressed; -s = output mapping quality; bcftools call -c = consensus caller; -M = keep masked ref; -O v = output type vcf)
samtools mpileup -f nc_crow-final.assembly.fasta.gapcloser-1.fa -b short_${1} -g -t DP -u -s | bcftools call -c -M -O v - > consensus_${species}.vcf

# Zip and index the vcf files
bgzip consensus_${species}.vcf
tabix -p vcf consensus_${species}.vcf.gz

# Remove INDELS from consensus.vcf
java -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T SelectVariants -R nc_crow-final.assembly.fasta.gapcloser-1.fa --variant:VCF consensus_${species}.vcf.gz -o consensus_${species}_noIndels.vcf.gz --selectTypeToExclude INDEL

# Index the vcf files
tabix -p vcf consensus_${species}_noIndels.vcf.gz

# NEW: Generate VCF file of sites that will be masked by “N” in the consensus fasta: 1) sites of DP < 3 and 2) truly polymorphic sites  
java -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T SelectVariants -R nc_crow-final.assembly.fasta.gapcloser-1.fa --variant:VCF consensus_${species}.vcf.gz -select "DP < 5" -o consensus_${species}_lowDP.vcf.gz

vcftools --gzvcf consensus_${species}.vcf.gz --maf 0.00001 --max-maf 0.99999 --min-alleles 2 --remove-indels --recode --recode-INFO-all --out consensus_${species}_trueSNPs

# NEW: Convert VCF files of sites to be masked into tab delimited files of the site positions, concatenate and sort them.
bgzip -cd consensus_${species}_lowDP.vcf.gz | grep -v '^#' | awk '{print $1 "\t" $2}' > consensus_${species}_lowDP.list
grep -v '^#' consensus_${species}_trueSNPs.recode.vcf | awk '{print $1 "\t" $2}' > consensus_${species}_trueSNPs.list
cat consensus_${species}_lowDP.list consensus_${species}_trueSNPs.list | sort -k1,1 -k2,2n -u > consensus_${species}_lowDP_trueSNPs.list

# Zip and index the vcf file of true SNPs to safe space
bgzip consensus_${species}_trueSNPs.recode.vcf
tabix -p vcf consensus_${species}_trueSNPs.recode.vcf.gz

# Convert consensus .vcf (without indels) to .fasta, masking sites from the list file
bcftools consensus --mask consensus_${species}_lowDP_trueSNPs.list --fasta-ref nc_crow-final.assembly.fasta.gapcloser-1.fa consensus_${species}_noIndels.vcf.gz > consensus_${species}_unfixed.fa

# Rename headers  so that they have the 'scaffold_10' format instead of '1 scaffold_10:1' one
perl -pne 's/>\d+\s+(.*):.*/>\1/' consensus_${species}_unfixed.fa > consensus_${species}.fa

# Move files to output directory
cp ${SNIC_TMP}/consensus_${species}*.vcf* /proj/b2013182/private/bams/Cmon_Ref/Consensus_per_species_nodup
cp ${SNIC_TMP}/consensus_${species}*.list /proj/b2013182/private/bams/Cmon_Ref/Consensus_per_species_nodup
cp ${SNIC_TMP}/consensus_${species}.fa /proj/b2013182/private/bams/Cmon_Ref/Consensus_per_species_nodup


