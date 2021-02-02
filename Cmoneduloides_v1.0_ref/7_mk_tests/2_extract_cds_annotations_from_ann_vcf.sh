#! /bin/bash -l
#SBATCH -A b2013182
#SBATCH -p core -n 1
#SBATCH -t 5-00:00:00
#SBATCH -C usage_mail

home="/proj/b2013182/nobackup/POPseq/vcf/Cmoneduloides_v1/Cmoneduloides/5inds_gatk_HC_ERC/snpEff"

cd $home

module load bioinfo-tools tabix/0.2.6
bgzip -cd Cmon_gatkHC_allSites_vqsr99_Q6_pass.ann.vcf.gz | /home/verena/glob/software/snpEff/scripts/vcfEffOnePerLine.pl | java -Xmx4g -jar /home/verena/glob/software/snpEff/SnpSift.jar extractFields - CHROM POS REF ALT "GEN[*].GT" "GEN[*].DP" "ANN[*].GENE:" "ANN[*].FEATUREID" "ANN[*].EFFECT" "ANN[*].ERRORS" > Cmon_allSites_annotations_160818.txt

# input file moved to /proj/b2013182/private/POPseq/vcf/Cmoneduloides_v1/Cmoneduloides/5inds_gatk_HC_ERC/snpEff
# output file moved to /proj/b2013182/nobackup/POPseq/mk_tests/keep