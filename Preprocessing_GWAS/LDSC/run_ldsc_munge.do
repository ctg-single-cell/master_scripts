#!/bin/bash
#SBATCH -N 1 
#SBATCH -t 2:00:00
#SBATCH -p normal

hapmap=/project/prjsbrouwer2/sumstats/munged_sumstats 
phenos=/project/prjsbrouwer2/sumstats/input_GWAS
ldsc_path=/home/brouwer2/sources/ldsc

for i in 1
do
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/SZ3_reformatted--merge-alleles ${hapmap}/w_hm3.snplist --out SZ3 --snp ID --N-con-col NCON --N-cas-col NCAS"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/BD_reformatted --merge-alleles ${hapmap}/w_hm3.snplist --out BD --snp ID --N-con-col NCON --N-cas-col NCAS"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/daner_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta --merge-alleles ${hapmap}/w_hm3.snplist --out ADHD --daner-n"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/iPSYCH-PGC_ASD_Nov2017 --merge-alleles ${hapmap}/w_hm3.snplist --out ASD --N 46350"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/PGC_UKB_depression_genome-wide.txt --merge-alleles ${hapmap}/w_hm3.snplist --signed-sumstats logOR,0 --out MDD --N 500199"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt --merge-alleles ${hapmap}/w_hm3.snplist --out BMI --a1 Tested_Allele --a2 Other_Allele"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/CUD_EUR_casecontrol_public_11.14.2020 --merge-alleles ${hapmap}/w_hm3.snplist --out AddCannabis "
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/AUDIT_UKB_2018_AJP.txt --merge-alleles ${hapmap}/w_hm3.snplist --out AddAlc --p P_T --a1 a_1 --a2 a_0 --signed-sumstats BETA_T,0"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/PD_reformatted --merge-alleles ${hapmap}/w_hm3.snplist --out PD --N-con-col N_controls --N-cas-col N_cases"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/AD_reformatted --merge-alleles ${hapmap}/w_hm3.snplist --out AD --a1 testedAllele --a2 otherAllele --signed-sumstats z,0"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/AN_reformatted --merge-alleles ${hapmap}/w_hm3.snplist --out AN --snp ID --N-cas-col NCAS --N-con-col NCON --a1 ALT --a2 REF"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/anxiety.meta.full.cc.tbl --merge-alleles ${hapmap}/w_hm3.snplist --out ANX --N-col TotalN"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/pts_aam_freeze2_overall.results --merge-alleles ${hapmap}/w_hm3.snplist --out PTSD --daner-n"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/SavageJansen_2018_intelligence_metaanalysis.txt --merge-alleles ${hapmap}/w_hm3.snplist --out IQ --N-col N_analyzed"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/sumstats_neuroticism_ctg_format.txt --merge-alleles ${hapmap}/w_hm3.snplist --out Neur --snp RSID --ignore SNP"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/Meta-analysis_Wood_et_al+UKBiobank_2018.txt --merge-alleles ${hapmap}/w_hm3.snplist --out height --a1 Tested_Allele --a2 Other_Allele"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.txt --merge-alleles ${hapmap}/w_hm3.snplist --out AgeSmk"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/GSCAN_DrnkWk_2022_GWAS_SUMMARY_STATS_EUR.txt --merge-alleles ${hapmap}/w_hm3.snplist --out DrnkWk"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt --merge-alleles ${hapmap}/w_hm3.snplist --out CigDay"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.txt --merge-alleles ${hapmap}/w_hm3.snplist --out SmkCes"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt --merge-alleles ${hapmap}/w_hm3.snplist --out SmkInit"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/HERMES_Jan2019_HeartFailure_summary_data.txt --merge-alleles ${hapmap}/w_hm3.snplist --out HF"
echo "python /home/brouwer2/sources/ldsc/munge_sumstats.py --sumstats ${phenos}/CAD_reformatted --merge-alleles ${hapmap}/w_hm3.snplist --signed-sumstats odds_ratio,1 --out CAD"
done > list_munge


parallel -N 1 --delay .2 -j 15 --joblog munge_GWAS.log --resume -a list_munge



