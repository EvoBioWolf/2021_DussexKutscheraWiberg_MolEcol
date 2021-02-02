#!/bin/bash
# This script starts slurm jobs for all species in the 5species set

echo "3sp_Ccornix"
sbatch Ccornix_ref_extract_genes.bsh 3sp_Ccornix 5species
echo "Cdau"
sbatch Ccornix_ref_extract_genes.bsh Cdau 5species
echo "Cfru"
sbatch Ccornix_ref_extract_genes.bsh Cfru 5species
echo "Csple"
sbatch Ccornix_ref_extract_genes.bsh Csple 5species
echo "Cmon"
sbatch Ccornix_ref_extract_genes.bsh Cmon 5species
echo "All Jobs Submitted"
