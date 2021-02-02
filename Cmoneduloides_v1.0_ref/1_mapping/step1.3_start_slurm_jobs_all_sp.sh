cd /proj/b2013182/nobackup/PHYLO/Cmoneduloides_v1/processed_bams

### Cfru, Csplendens
for f in $(ls *_dedup.bam | grep IRL_Lm_R01 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep IRL_Lm_R05 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep RUS_Mu_R03 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep RUS_Mu_R04 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep CKW77 | sed 's/_L00[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep CSW7069 | sed 's/_L00[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep EVL981 | sed 's/_L00[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep LKW131 | sed 's/_L00[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep PJG222 | sed 's/_L00[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done

### Cmon, Cwoo
for f in $(ls *_dedup.bam | grep FS66096 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep NC3 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep NC5 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep NC6 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep NC8 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep MLK92 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep MLK164 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep MLK166 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep CROW12 | sed 's/_L00[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep CROWR21 | sed 's/_L00[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done

### Cmone

for f in $(ls *_dedup.bam | grep S_Ri_J01 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep S_Ri_J02 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep S_Ri_J03 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep S_Ri_J08 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done


### Cdau
for f in $(ls *_dedup.bam | grep FS66096 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep NC3 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep NC5 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep NC6 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep NC8 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep MLK92 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep MLK164 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep MLK166 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep CROW12 | sed 's/_L00[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep CROWR21 | sed 's/_L00[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done

### C.corx C.haw
for f in $(ls *_dedup.bam | grep 2407-51841 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep Oli-Studbook67 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done

### C.tas
for f in $(ls *_dedup.bam | grep B44919 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done


### C.corone cornix orientalis
for f in $(ls *_dedup.bam | grep D_Ko_C04 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep D_Ko_C13 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep D_Ko_C15 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep D_Ra_C05 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep D_Ra_C16 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep E_Vi_C01 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep E_Vi_C57 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep E_Vi_C58 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done

for f in $(ls *_dedup.bam | grep B_So_H04 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep ISR_TA_H01 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep ITA_Ro_H11 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep ITA_Ro_H12 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep PL_Wa_H22 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep RUS_No_H03| sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep S_Up_H03 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep S_Up_H09 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done

for f in $(ls *_dedup.bam | grep RUS_Pr_O01 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep RUS_Pr_O02 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep RUS_Pr_O03 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep RUS_Pr_O04 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep RUS_Pr_O05 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep RUS_Tv_O01 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done

### Cbrach

for f in $(ls *_dedup.bam | grep USA_CA_B01 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep USA_CA_B03 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep USA_CA_B08 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep USA_NJ_B02 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep USA_NY_B04 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
for f in $(ls *_dedup.bam | grep USA_NY_B05 | sed 's/_S[0-9]_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done

### Ckub 
for f in $(ls bwa_Ckub*_dedup.bam | grep Ckub | sed 's/_S1_sorted_rg_dedup.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step1.3_merge.sh $f; done
