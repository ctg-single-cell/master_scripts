# This is a snakemake file for running the initial_setup.py script
import os


# current_dir = "/gpfs/work5/0/vusr0480/Preprocessing_scRNA/"
ids = []
with open(os.path.join("/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/scrnaseq_data_master.csv"), "r") as f:
    for line in f:
        if not line.startswith("Internal_ID"):
            ids.append(line.rstrip("\n").split(",")[0])
# print(ids)

rule all:
    input:
        expand("/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/{id}/readme.md", id=ids)

rule initial_setup:
    output:
        "/gpfs/work5/0/vusr0480/Preprocessing_scRNA/data/{id}/readme.md"
    params:
        id = "{id}"
    shell:
        """
        python initial_setup.py --id {params.id} --work_dir /gpfs/work5/0/vusr0480/Preprocessing_scRNA/data
        """