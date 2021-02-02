#! /bin/bash -l

# bash script to start slurm jobs to merge all bam files per individual using the GATK v3.4.0 pipeline
# Merging script is located in
#scripts=/proj/b2013182/nobackup/batch_scripts/Ccornix_ref/Ccornix_ref_mapping
# OR
scripts=/pica/h1/axelw/batch_scripts/Ccornix_ref_mapping
# bam files are located in /proj/b2013182/nobackup/PHYLO/Ccornix_v2.5/processed_bams/
bams=/proj/b2013182/nobackup/PHYLO/Ccornix_v2.5/processed_bams
# Move to the bams directory
cd ${bams}

# Remember to move all the slurm logs back to the batch_scripts folder

# C. cornix
echo "C. cornix"
for f in $(ls *_mrkdup.bam | grep IRL_Lm_R01 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep B_So_H04 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep ISR_TA_H01 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep ITA_Ro_H11 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep ITA_Ro_H12 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep PL_Wa_H22 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep RUS_No_H03 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done #Only 1 sample
for f in $(ls *_mrkdup.bam | grep S_Up_H03 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep S_Up_H09 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done

# C. corone
echo "C. corone"
for f in $(ls *_mrkdup.bam | grep D_Ko_C04 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep D_Ko_C13 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep D_Ko_C15 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep D_Ra_C05 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep D_Ra_C16 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep E_Vi_C01 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep E_Vi_C57 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep E_Vi_C58 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done

# C. orientalis
echo "C. orientalis"
for f in $(ls *_mrkdup.bam | grep RUS_Pr_O01 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep RUS_Pr_O02 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep RUS_Pr_O03 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep RUS_Pr_O04 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep RUS_Pr_O05 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep RUS_Tv_O01 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done #Only 1 sample

# C. brach
echo "C. brachyrhynchos"
for f in $(ls *_mrkdup.bam | grep USA_CA_B01 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep USA_CA_B03 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep USA_CA_B03 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep USA_CA_B08 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep USA_NJ_B02 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep USA_NY_B03 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep USA_NY_B04 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep USA_NY_B05 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done

# C. dauricus
echo "C. dauricus"
for f in $(ls *_mrkdup.bam | grep DDCHN_Gu_D02 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done #Only 1 sample
for f in $(ls *_mrkdup.bam | grep MON_Kh_D01 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep MON_Kh_D02 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep MON_Kh_D03 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done

# C. frugilegus
echo "C. frugilegus"
for f in $(ls *_mrkdup.bam | grep IRL_Lm_R01 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep IRL_Lm_R05 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done #Only 1 sample
for f in $(ls *_mrkdup.bam | grep RUS_Mu_R03 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep RUS_Mu_R04 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done

# C. moneduloides
echo "C. moneduloides"
for f in $(ls *_mrkdup.bam | grep FS66096 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep NC3 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep NC5 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep NC6 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep NC8 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done

# C. monedula
echo "C. monedula"
for f in $(ls *_mrkdup.bam | grep S_Ri_J01 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep S_Ri_J02 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep S_Ri_J03 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep S_Ri_J08 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done

# C. splendens
echo "C. splendens"
for f in $(ls *_mrkdup.bam | grep CKW77 | sed 's/_L00[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep CSW7069 | sed 's/_L00[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep EVL981 | sed 's/_L00[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep LKW131 | sed 's/_L00[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep PJG222 | sed 's/_L00[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done

# C. woodfordii
echo "C. woodfordii"
for f in $(ls *_mrkdup.bam | grep CROW12 | sed 's/_L00[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep CROWR21 | sed 's/_L00[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep MLK164 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep MLK166 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done
for f in $(ls *_mrkdup.bam | grep MLK92 | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done


# C. kubaryi # Only one individual sampled
echo "C. kubaryi"
for f in $(ls *_mrkdup.bam | grep "Rota" | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done

# C. hawaiiensis # Only one individual sampled
echo "C. hawaiiensis"
for f in $(ls *_mrkdup.bam | grep "Oli-Studbook67" | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done

# C. corax # Only one individual sampled
echo "C. corax"
for f in $(ls *_mrkdup.bam | grep "2407-51841" | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done

# C. tasmanicus # Only one individual sampled
echo "C. tasmanicus"
for f in $(ls *_mrkdup.bam | grep "B44919" | sed 's/_S[0-9]_sorted_rg_mrkdup.bam//g' | uniq); do sbatch ${scripts}/step1.3_merge.sh ${f}; done


