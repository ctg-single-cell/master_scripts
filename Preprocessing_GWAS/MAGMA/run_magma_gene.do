#!/bin/bash
#SBATCH -N 1 
#SBATCH -t 20:00:00
#SBATCH -p normal

# Preparing MAGMA on the input files
magma=/home/brouwer2/sources/magma_1.10/magma
ref_input=/home/brouwer2/sources/magma_1.10/g1000_eur

phenos=/home/brouwer2/sources/sumstats/input_GWAS
gene_loc_dir=/home/brouwer2/sources/sumstats/magma_gene/gene_locations

gunzip ${phenos}/*gz 

# we assume that the annotation files are already created
# all magma's are run using multi option - new release version 1.10 and are based on entrez IDs

for annot in all_ensemble protein_coding_ensemble all_entrez protein_coding_entrez
do
if [ ${batch} == "first" ]
then
# SZ
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/SZ3_reformatted  pval=PVAL ncol=NEFF snp-id=ID --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/SZ3_${annot}_magma &

# Bipolar disorder
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/BD_reformatted  pval=PVAL ncol=NEFFDIV2 snp-id=ID --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/BP_${annot}_magma &

# ADHD
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.syn onyms --pval ${phenos}/daner_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta  pval=P ncol=Neff snp-id=SNP --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/ADHD_${annot}_magma &

# ASD
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/iPSYCH-PGC_ASD_Nov2017  pval=P N=46350 snp-id=SNP --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/ASD_${annot}_magma &

# MDD 
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/PGC_UKB_depression_genome-wide.txt  pval=P N=500199 snp-id=MarkerName --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/MDD_${annot}_magma &

# BMI
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt  pval=P ncol=N snp-id=SNP --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/BMI_${annot}_magma &

# Addiction - cannabis
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/CUD_EUR_casecontrol_public_11.14.2020 pval=P ncol=N snp-id=SNP --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/AddCannabis_${annot}_magma &

# Addiction - alcohol
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/AUDIT_UKB_2018_AJP.txt  pval=P_T ncol=N snp-id=rsid --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/AddAlc_${annot}_magma &

# Parkinsons

$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/PD_reformatted  pval=p ncol=Neff snp-id=rsID --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/PD_${annot}_magma &

# Alzheimers

$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/AD_reformatted  pval=p ncol=N snp-id=SNP --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/AD_${annot}_magma &

# Anorexia 
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/pgcAN2.2019-07.vcf.tsv  pval=PVAL ncol=NEFFDIV2 snp-id=ID --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/AN_${annot}_magma &

# Anxiety
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/anxiety.meta.full.cc.tbl  pval=P.value ncol=TotalN snp-id=SNPID --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/ANX_${annot}_magma &

# PTSD
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/pts_aam_freeze2_overall.results pval=P ncol=Neff snp-id=SNP --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/PTSD_${annot}_magma &

# intelligence
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/SavageJansen_2018_intelligence_metaanalysis.txt  pval=P ncol=N_analyzed snp-id=SNP --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/IQ_${annot}_magma &

# neuroticism
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/sumstats_neuroticism_ctg_format.txt  pval=P ncol=N snp-id=RSID --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/Neur_${annot}_magma &

# height
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/Meta-analysis_Wood_et_al+UKBiobank_2018.txt pval=P snp-id=SNP ncol=N --gene-model multi --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/Height_${annot}_magma
wait 
fi

if [ ${batch} == "second" ]
then
# heart failure
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/HERMES_Jan2019_HeartFailure_summary_data.txt pval=p snp-id=SNP ncol=N --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/HF_${annot}_magma &

# Addiction - drinks per week
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/GSCAN_DrnkWk_2022_GWAS_SUMMARY_STATS_EUR.txt pval=P snp-id=RSID ncol=N --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/DrnkWk_${annot}_magma &

# Addiction - Cigarettes per day
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt pval=P snp-id=RSID ncol=N --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/CigDay_${annot}_magma &

# Addiction - Smoking initiation (ever smoked regularly)
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt pval=P snp-id=RSID ncol=N --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/SmkInit_${annot}_magma &

# Addiction - Age of smoking initiation
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.txt pval=P snp-id=RSID ncol=N --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/AgeSmk_${annot}_magma &

# Addiction - Smoking cessation (former versus current smokers)
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.txt pval=P snp-id=RSID ncol=N --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/SmkCes_${annot}_magma &

# CAD 
$magma --bfile $ref_input/g1000_eur synonyms=$ref_input/g1000_eur.synonyms --pval ${phenos}/CAD_reformatted pval=p_value snp-id=SNP ncol=n --gene-annot ${gene_loc_dir}/${annot}.genes.annot --out ${annot}/CAD_${annot}_magma &

fi

done
wait

