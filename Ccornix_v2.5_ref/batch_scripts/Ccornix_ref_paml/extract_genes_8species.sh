#!/bin/bash
# This script starts slurm jobs for all species in the 8species set

echo "3sp_Ccorone"
sbatch Ccornix_ref_extract_genes.bsh 3sp_Ccorone_mskVars 8species
echo "Cdau"
sbatch Ccornix_ref_extract_genes.bsh Cdau_mskVars 8species
echo "Cfru"
sbatch Ccornix_ref_extract_genes.bsh Cfru_mskVars 8species
echo "Csple"
sbatch Ccornix_ref_extract_genes.bsh Csple_mskVars 8species
echo "Cmon"
sbatch Ccornix_ref_extract_genes.bsh Cmon_mskVars 8species
echo "Ctas"
sbatch Ccornix_ref_extract_genes.bsh Ctas_mskVars 8species
echo "Ccorx"
sbatch Ccornix_ref_extract_genes.bsh Ccorx_mskVars 8species
echo "Chaw"
sbatch Ccornix_ref_extract_genes.bsh Chaw_mskVars 8species
echo "All Jobs Submitted"
