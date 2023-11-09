import os

rule all:
    input:
        expand(os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/ldsc-seg_out/baselineLD_v2.1/{trait}_ldsc-seg_sumstats.aggregated_cell_type_results.txt"), method=config["spec_methods"], scrna=config["scrnas"], trait=config["traits"]),
        expand(os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/ldsc-seg_out/baselineLD_v2.1/{trait}_ldsc-seg_sumstats.aggregated_cell_type_results_P.txt"), method=config["spec_methods"], scrna=config["scrnas"], trait=config["traits"])

rule aggregate_cts:
   output:
      file=os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/ldsc-seg_out/baselineLD_v2.1/{trait}_ldsc-seg_sumstats.aggregated_cell_type_results.txt")
   params:
      script="/gpfs/home6/tphung/tmp/20230831_ldsc/extract_results_ldsc-seg.sh",
      inputdir=os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/ldsc-seg_out/baselineLD_v2.1/"),
      gwas = "{trait}"
   shell:
      """
      {params.script} {params.inputdir} {params.gwas} {output.file}
      """

rule compute_p:
    input:
        os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/ldsc-seg_out/baselineLD_v2.1/{trait}_ldsc-seg_sumstats.aggregated_cell_type_results.txt")
    output:
        os.path.join(config["ldsc_dir"],"spec_{method}/{scrna}/ldsc-seg_out/baselineLD_v2.1/{trait}_ldsc-seg_sumstats.aggregated_cell_type_results_P.txt")
    params:
    shell:
        """
        Rscript post_ldsc_processing.R --input {input} --output {output}
        """