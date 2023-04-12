# GWAS preprocessing
 
This directory holds the scripts that were used to preprocess GWAS summary statistics for later use in MAGMA or LDSC

Included in this directory are currently (see below for description of data sources):
 
 phenotype | abbreviation 
 |:-- |:-- | 
 Addiction - alcohol use | AddAlc 
 Addiction - cannabis use | AddCannabis
 ADHD |  
 Alzheimer's disease | AD
 Anorexia Nervosa | AN
 Anxiety | ANX
 ASD | 
 Bipolar disorder | BP
 BMI | 
 Height |
 Heart Failure | HF
 Intelligence | IQ
 MDD |
 Neuroticism | Neur
 Parkinson's | PD
 PTSD |
 Schizophrenia - UPDATED 7/6/22 | SZ
 Addiction - Drinks per week (DrnkWk), added 11/1/23 | DrnkWk 	
 Addiction - Cigarettes per day (CigDay), added 11/1/23 | CigDay
 Addiction - Smoking initiation (SmkInit), added 11/1/23 | SmkInit 
 Addiction - Smoking cessation (SmkCes), added 11/1/23 | SmkCes
 Addiction - Age of initiation (AgeSmk), added 11/1/23 | AgeSmk
 Coronary Artery Disease - added 12/1/23 | CAD



Raw data files are currently stored in /project/prjsbrouwer2/sumstats on lisa.

---
## General prerocessing of summary stats 

Some of the files need a little preprocessing before we can use them:

AD/PD/CAD need the addition of rsIDs (scripts in preprocessing directory, make sure to add the correct paths to the working directory, input file directory and the 1000G phase 3 bimfile : 

```
R --no-save < Add_rs_AD.R
R --no-save < Add_rs_PD.R
R --no-save < Add_rs_CAD.R
```

reformat BD file to include the header (commented out in original file)
```
cat pgc-bip2021-all.vcf.tsv | grep CHROM | awk -F'#' '{ print $2 }' > ${phenos}/BD_reformatted
cat pgc-bip2021-all.vcf.tsv | grep -v "#"  >> ${phenos}/BD_reformatted
```

for AN and SZ remove all the extra information before the actual sumstats, munge_stats does not like that
```
cat pgcAN2.2019-07.vcf.tsv | grep -v '##' > AN_reformatted
cat PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv | grep -v '##' > SZ3_reformatted
```

## Processing for use with MAGMA

See MAGMA directory for scripts that process the GWAS above for use with MAGMA (different choices of gene naming system and genes to include (protein coding only or all genes)).

## Processing for use with LDSC

See LDSC directory for instructions to process the GWAS above for use with LDSC. 

---

# Description of data sources

* Addiction - alcohol use
Sanchez-Roige S et al. Genome-Wide Association Study Meta-Analysis of the Alcohol Use Disorders Identification Test (AUDIT) in Two Population-Based Cohorts. Am J Psychiatry. 2019 Feb 1;176(2):107-118. doi: 10.1176/appi.ajp.2018.18040369. Epub 2018 Oct 19. PMID: 30336701; PMCID: PMC6365681.

File AUDIT_UKB_2018_AJP.txt.gz downloaded from PGC https://www.med.unc.edu/pgc/download-results/ on 14/01/2022

* Addiction - cannabis use
Johnson EC et al. A large-scale genome-wide association study meta-analysis of cannabis use disorder. Lancet Psychiatry. 2020 Dec;7(12):1032-1045. doi: 10.1016/S2215-0366(20)30339-4. Epub 2020 Oct 20. PMID: 33096046; PMCID: PMC7674631.

File CUD_EUR_casecontrol_public_11.14.2020.gz downloaded from PGC https://www.med.unc.edu/pgc/download-results/ on 14/01/2022

* ADHD 
Demontis D, et al. Discovery of the first genome-wide significant risk loci for attention deficit/hyperactivity disorder. Nat. Genet. 51, 63-75 (2019).

File daner_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta downloaded from PGC https://www.med.unc.edu/pgc/download-results/ on 13/01/2022

* Alzheimer's
Wightman D et al. A genome-wide association study with 1,126,563 individuals identifies new risk loci for Alzheimer's disease.  Genet. 2021 Sep;53(9):1276-1282 

File (excluding 23andMe) PGCALZ2sumstatsExcluding23andMe.txt.gz downloaded from https://ctg.cncr.nl/software/summary_statistics on 14/01/2022

* Anorexia Nervosa
Watson HJ, et al. Genome-wide association study identifies eight risk loci and implicates metabo-psychiatric origins for anorexia nervosa. Nat Genet. 2019 Aug;51(8):1207-1214. doi: 10.1038/s41588-019-0439-2. Epub 2019 Jul 15. PMID: 31308545; PMCID: PMC6779477.

File pgcAN2.2019-07.vcf.tsv.gz downloaded from from PGC https://www.med.unc.edu/pgc/download-results/ on 14/01/2022

* Anxiety 
Otowa T, et al. Meta-analysis of genome-wide association studies of anxiety disorders [published correction appears in Mol Psychiatry. 2016 Oct;21(10):1485]. Mol Psychiatry. 2016;21(10):1391-1399. doi:10.1038/mp.2015.197

File anxiety.meta.full.cc.tbl.gzip downloaded from from PGC https://www.med.unc.edu/pgc/download-results/ on 14/01/2022

* ASD
Grove J, et al. Identification of common genetic risk variants for autism spectrum disorder. Nat Genet. 2019 Mar;51(3):431-444. doi: 10.1038/s41588-019-0344-8. Epub 2019 Feb 25. PMID: 30804558; PMCID: PMC6454898.

File iPSYCH-PGC_ASD_Nov2017 downloaded from PGC https://www.med.unc.edu/pgc/download-results/ on 13/01/2022

* BMI
Yengo L Meta-analysis of genome-wide association studies for height and body mass index in ~700000 individuals of European ancestry, Human Molecular Genetics, Volume 27, Issue 20, 15 October 2018, Pages 3641-3649, https://doi.org/10.1093/hmg/ddy271

File Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt downloaded from https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#BMI_and_Height_GIANT_and_UK_BioBank_Meta-analysis_Summary_Statistics on 14/01/2022

* BP
Mullins N et al. Genome-wide association study of more than 40,000 bipolar disorder cases provides new insights into the underlying biology. Nat Genet. 2021 Jun;53(6):817-829. doi: 10.1038/s41588-021-00857-4. Epub 2021 May 17. PMID: 34002096; PMCID: PMC8192451.

Files (3, BP combined, type 1 and type 2) pgc-bip2021-all.vcf.tsv.gz; pgc-bip2021-BDI.vcf.tsv.gz; pgc-bip2021-BDII.vcf.tsv.gz; available downloaded from PGC https://www.med.unc.edu/pgc/download-results/ on 13/01/2022
We use the combined phenotype for these analyses

* Height
Yengo L Meta-analysis of genome-wide association studies for height and body mass index in ~700000 individuals of European ancestry, Human Molecular Genetics, Volume 27, Issue 20, 15 October 2018, Pages 3641-3649, https://doi.org/10.1093/hmg/ddy271

File Meta-analysis_Wood_et_al+UKBiobank_2018.txt downloaded from https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#BMI_and_Height_GIANT_and_UK_BioBank_Meta-analysis_Summary_Statistics on 25/01/2022

* Heart Failure
Shah, S., Henry, A., Roselli, C. et al. Genome-wide association and Mendelian randomisation analysis provide insights into the pathogenesis of heart failure. Nat Commun 11, 163 (2020). https://doi.org/10.1038/s41467-019-13690-5

File HERMES_Jan2019_HeartFailure_summary_data.txt downloaded from https://cvd.hugeamp.org/dinspector.html?dataset=GWAS_HERMES_eu on 21/11/22

* MDD
Howard DM et al. Genome-wide meta-analysis of depression identifies 102 independent variants and highlights the importance of the prefrontal brain regions. Nat Neurosci. 2019 Mar;22(3):343-352. doi: 10.1038/s41593-018-0326-7. Epub 2019 Feb 4. PMID: 30718901; PMCID: PMC6522363.

File PGC_UKB_depression_genome-wide.txt.gz (excluding 23andMe) downloaded from https://datashare.ed.ac.uk/handle/10283/3203 on 13/01/2022

* Neuroticism
Nagel, M., Jansen, P.R., Stringer, S. et al. Meta-analysis of genome-wide association studies for neuroticism in 449,484 individuals identifies novel genetic loci and pathways.Nat Genet 50, 920-927 (2018). https://doi.org/10.1038/s41588-018-0151-7 

File sumstats_neuroticism_ctg_format.txt.gz (excluding 23andMe) downloaded from https://ctg.cncr.nl/software/summary_statistics/ on 25/01/2022


* Intelligence
Savage et al. Genome-wiide association meta-analysis (N=269,867) identifies new genetic and functional links to intelligence. Nature Genetics, 2018 Jul;50(7):912-919

File SavageJansen_2018_intelligence_metaanalysis.txt,gz downloaded from https://ctg.cncr.nl/software/summary_statistics/ on 14/01/22

* Parkinson's
Nalls MA et al. Identification of novel risk loci, causal insights, and heritable risk for Parkinson's disease: a meta-analysis of genome-wide association studies. Lancet Neurol. 2019 Dec;18(12):1091-1102. doi: 10.1016/S1474-4422(19)30320-5. PMID: 31701892; PMCID: PMC8422160.

File nallsEtAl2019_excluding23andMe_allVariants.tab.gz (excluding 23andMe) downloaded from https://drive.google.com/drive/folders/10bGj6HfAXgl-JslpI9ZJIL_JIgZyktxn on 14/01/2022

* PTSD
Nievergelt CM et al. International meta-analysis of PTSD genome-wide association studies identifies sex- and ancestry-specific genetic risk loci. Nat Commun. 2019 Oct 8;10(1):4558. doi: 10.1038/s41467-019-12576-w. PMID: 31594949; PMCID: PMC6783435.

File pts_aam_freeze2_overall.results.gz downloaded from PGC https://www.med.unc.edu/pgc/download-results/ on 14/01/2022

* SZ
Mapping genomic loci prioritises genes and implicates synaptic biology in schizophrenia
The Schizophrenia Working Group of the Psychiatric Genomics Consortium, Stephan Ripke, James TR Walters, Michael C O'Donovan doi: https://doi.org/10.1101/2020.09.12.20192922

Files ~/sources/PGC/summary_stats/autosomes/meta/daner_PGC_SCZ_w3_90_0518d_eur.gz
and ~/sources/PGC/summary_stats/chrX/meta/daner_scz_w3_HRC_chrX_deduped_0618a.gz 
Used with PGC permission before publication of the paper. 

AS OF JUNE 7, 2022: we use the file that is freely available from the PGC website (and ignore the X chromosome, turns out it is of mixed ancestry). 
file PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv  https://www.med.unc.edu/pgc/download-results/ downloaded on June 2, 2022

Different set of addiction related pehnotypes asociated with the paper
Saunders, G.R.B., Wang, X., Chen, F. et al. Genetic diversity fuels gene discovery for tobacco and alcohol use. Nature (2022).
https://doi.org/10.1038/s41586-022-05477-4
Five phenoytypes:
Drinks per week (DrnkWk), Cigarettes per day (CigDay), Smoking initiation (SmkInit), Smoking cessation (SmkCes), and Age of initiation (AgeSmk).
File EUR_stratified.zip containing these sumstats (European ancestry) downloaded from https://conservancy.umn.edu/handle/11299/241912
This contains files GSCAN_XXX_2022_GWAS_SUMMARY_STATS_EUR.txt.gz for XXX in DrnkWk, CigDay, SmkInit, SmkCes, AgeSmk

Coronary Artery Disease (CAD)
Aragam, K.G., Jiang, T., Goel, A. et al. Discovery and systematic characterization of risk variants and genes for coronary artery disease in over a million participants. Nat Genet 54, 1803-1815 (2022). https://doi-org.vu-nl.idm.oclc.org/10.1038/s41588-022-01233-6
File GCST90132314_buildGRCh37.tsv downloaded from http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/ on 12/01/23




