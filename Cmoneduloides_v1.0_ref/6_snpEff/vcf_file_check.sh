#! /bin/bash -l
#SBATCH -A b2010060
#SBATCH -p core -n 1
#SBATCH -t 48:00:00
#SBATCH -C usage_mail

cd /proj/b2013182/nobackup/POPseq/vcf/Cmoneduloides_v1/Cmoneduloides/5inds_gatk_HC_ERC/snpEff

java -Xmx4g -jar /home/verena/glob/software/snpEff/SnpSift.jar vcfCheck Cmon_gatkHC_allSites_vqsr99_Q6_pass.ann.vcf > Cmon_gatkHC_allSites_vqsr99_Q6_pass.ann.vcfCheck.txt

