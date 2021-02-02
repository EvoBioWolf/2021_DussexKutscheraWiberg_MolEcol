#! /bin/bash -l
#SBATCH -A b2013182
#SBATCH -p core -n 1
#SBATCH -t 48:00:00
#SBATCH -C usage_mail

vcf="/proj/b2013182/nobackup/POPseq/vcf/Cmoneduloides_v1"

#for f in $(ls *vqsr99.9_Q6_pass.recode.vcf.gz | grep Cmon | sed 's/_gatkHC_allSites_vqsr99.9_Q6_pass.recode.vcf.gz//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/snpEff/annotate_vcf.sh $f; done
#for f in $(ls *vqsr99.9_Q6_pass.recode.vcf.gz | grep Cwoo | sed 's/_gatkHC_allSites_vqsr99.9_Q6_pass.recode.vcf.gz//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/snpEff/annotate_vcf.sh $f; done

cd $vcf/${1}*/*gatk_HC_ERC
java -Xmx4g -jar /home/verena/glob/software/snpEff/snpEff.jar NCcrow.v1 ${1}_gatkHC_allSites_vqsr99.9_Q6_pass.recode.vcf.gz > ${1}_gatkHC_allSites_vqsr99.9_Q6_pass.recode.ann.vcf

# compress and index vcf file
module load bioinfo-tools tabix/0.2.6
bgzip -c ${1}_gatkHC_allSites_vqsr99.9_Q6_pass.recode.ann.vcf > ${1}_gatkHC_allSites_vqsr99.9_Q6_pass.recode.ann.vcf.gz
tabix -p vcf ${1}_gatkHC_allSites_vqsr99.9_Q6_pass.recode.ann.vcf.gz
