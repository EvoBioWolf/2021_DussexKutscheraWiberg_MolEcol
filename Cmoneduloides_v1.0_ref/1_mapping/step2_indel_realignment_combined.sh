#! /bin/bash -l
#SBATCH -A b2013182
#SBATCH -p core -n 2
#SBATCH -t 2-00:00:00
#SBATCH -C usage_mail

#job id: $SLURM_JOB_ID
#job name (-J): $SLURM_JOB_NAME
#tmp directory: $SNIC_TMP

#java heap size: 2GB less than max. memory, i.e. 6GB per core
mem=36g
bams="/proj/b2013182/nobackup/PHYLO/Cmoneduloides_v1/processed_bams"
ref="/proj/b2013182/private/assemblies/Cmoneduloides/otago-temp"

#Usage; doesn't seem to work
#for f in $(ls *merged.bam | grep Cfru | sed 's/_[A-Z]*_[A-Z][a-z]_[A-Z]*[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done
#for f in $(ls *merged.bam | grep Csplendens | sed 's/_[A-Z]*[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done
#for f in $(ls *merged.bam | grep Ccorx | sed 's/_[0-9]*-[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done
#for f in $(ls *merged.bam | grep Ctas | sed 's/_[A-Z]*[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done
#for f in $( ls *merged.bam | grep Chaw | sed 's/_[A-Z]*[a-z]*-[A-Z]*[a-z]*[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done
#for f in $( ls *merged.bam | grep C.corone | sed 's/_[A-Z]*_[A-Z]*[a-z]_[A-Z]*[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done
#for f in $(ls *_merged.bam | grep C.cornix | sed 's/_[A-Z]*_[A-Z]*[A-Z]*[a-z]*_[A-Z]*[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done
#for f in $( ls *merged.bam | grep C.orientalis | sed 's/_[A-Z]*_[A-Z]*[a-z]_[A-Z]*[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done
#for f in $( ls *merged.bam | grep Cdau | sed 's/_[A-Z]*_[A-Z]*[a-z]_[A-Z]*[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done
#for f in $( ls *merged.bam | grep Cfru | sed 's/_[A-Z]*_[A-Z]*[a-z]_[A-Z]*[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done
#for f in $(ls *_merged.bam | grep Cmon_ | sed 's/_[A-Z]*[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done
#for f in $(ls *_merged.bam | grep Cmone | sed 's/_[A-Z]*_[A-Z]*[a-z]_[A-Z]*[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done
#for f in $(ls *_merged.bam | grep Ckub | sed 's/_[A-Z]*[a-z]*_[0-9]*-[0-9]*-[0-9]*_merged.bam//g' | uniq); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/mapping/step2_indel_realignment_combined.sh $f; done


#####
##Copy files from home to tmp directory
#####
cp $bams/${1}*_merged.ba* ${SNIC_TMP}/
cp $ref/nc_crow-final.assembly.fasta.gapcloser-1.fa ${SNIC_TMP}/
cp $ref/nc_crow-final.assembly.fasta.gapcloser-1.fa.fai ${SNIC_TMP}/
cp $ref/nc_crow-final.assembly.fasta.gapcloser-1.dict ${SNIC_TMP}/

#########realign indels for all bam files from one species combined#########
cd ${SNIC_TMP}
java -Xmx${mem} -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -R nc_crow-final.assembly.fasta.gapcloser-1.fa $(printf ' -I %s' ${1}*_merged.bam) -o ${1}_realignment_targets.list -nt 6
java -Xmx${mem} -jar /sw/apps/bioinfo/GATK/3.4.0/GenomeAnalysisTK.jar -T IndelRealigner -R nc_crow-final.assembly.fasta.gapcloser-1.fa $(printf ' -I %s' ${1}*_merged.bam) -targetIntervals ${1}_realignment_targets.list -nWayOut _realigned.bam

#####
##Copy files from tmp directory to project directory
#####
cp ${SNIC_TMP}/*realign* $bams


#start batch scripts to check the quality of the new bam files
cd $bams

          
for f in $(ls *realigned.bam | sed 's/.bam//g'); do sbatch /proj/b2013182/nobackup/batch_scripts/Cmon_Ref/gatk/bam_quality_control.sh $f; done
