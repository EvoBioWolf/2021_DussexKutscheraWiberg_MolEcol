#! /bin/bash -l
#SBATCH -A b2010060
#SBATCH -p core -n 2
#SBATCH -t 48:00:00
#SBATCH -C usage_mail

cd /home/verena/glob/software/snpEff
java -Xmx8g -jar snpEff.jar build -gff3 -v NCcrow.v1

