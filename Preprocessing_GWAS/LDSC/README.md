# Munge sumstats

Installation of LDSC: follow instructions on https://github.com/bulik/ldsc; no need to download any reference data.
LDSC makes use of a reference SNP list. Originally obtained from https://data.broadinstitute.org/alkesgroup/LDSCORE/${hapmap}/w_hm3.snplist.bz2. These data are currently behind a paywall, we make use of the originally downloaded file located on the system. Set these paths:

```
hapmap=/project/prjsbrouwer2/sumstats/munged_sumstats 
phenos=/project/prjsbrouwer2/sumstats/input_GWAS
ldsc_path=/home/brouwer2/sources/ldsc
```

The run_ldsc_munge.do script contains the program call munge the sumstats with GWAS specific settings.  
Submit to the cluster using 

```
sbatch run_ldsc_munge.do
```


