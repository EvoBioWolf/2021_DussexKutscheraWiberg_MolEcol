#!/bin/bash
# This script starts slurm jobs for all species in the 7species set

echo "3sp_Ccorone"
sbatch Ccornix_ref_extract_genes.bsh 3sp_Ccorone_mskVars 7species
echo "Cdau"
sbatch Ccornix_ref_extract_genes.bsh Cdau_mskVars 7species
echo "Cfru"
sbatch Ccornix_ref_extract_genes.bsh Cfru_mskVars 7species
echo "Csple"
sbatch Ccornix_ref_extract_genes.bsh Csple_mskVars 7species
echo "Cmon"
sbatch Ccornix_ref_extract_genes.bsh Cmon_mskVars 7species
echo "Ctas"
sbatch Ccornix_ref_extract_genes.bsh Ctas_mskVars 7species
echo "Ccorx"
sbatch Ccornix_ref_extract_genes.bsh Ccorx_mskVars 7species
echo "All Jobs Submitted"
