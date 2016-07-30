import os

configfile: "config.yaml"

SAMPLES = config["SAMPLES"]
ASSEMBLIES = config["ASSEMBLIES"]
CENTRIFUGE_BASE = config["CENTRIFUGE_BASE"]

ENV_ANVI = "source activate anvio2"

shell.prefix(ENV_ANVI + '; ')

# SAMPLES = ["A", "B"]
# ASSEMBLIES = ["X","Y"]
# CENTRIFUGE_BASE = "/home/jgsanders/miniconda/envs/anvio2/centrifuge"


rule all:
    input: 
        expand("data/anvio/{assembly}/{assembly}_SAMPLES-SUMMARY/index.html",
               assembly=ASSEMBLIES,
               sample=SAMPLES)


rule bowtie2_index:
    input:
        "data/assemblies/{assembly}.fa"
    output:
        idx="data/assemblies/{assembly}.rev.2.bt2"
    shell:
        "bowtie2-build {input} data/assemblies/{wildcards.assembly}"


rule bowtie2_map:
    input:
        fa="data/assemblies/{assembly}.fa",
        idx="data/assemblies/{assembly}.rev.2.bt2",
        R1="data/samples/{sample}.R1.fq.gz",
        R2="data/samples/{sample}.R2.fq.gz"
    params:
        threads=12,
        idx_base="data/assemblies/{assembly}"
    output:
        temp("data/mapped_reads/{assembly}.{sample}.bam")
    shell:
        "bowtie2 -x {params.idx_base} -p {params.threads} --no-unal " \
                  "-q -1 {input.R1} -2 {input.R2} | " \
                  "samtools view -bS - > {output}"

rule samtools_sort:
    input:
        "data/mapped_reads/{assembly}.{sample}.bam"
    output:
        "data/sorted_reads/{assembly}.{sample}.bam"
    shell:
        "samtools sort -O bam -o {output} {input}"

rule samtools_index:
    input:
        "data/sorted_reads/{assembly}.{sample}.bam"
    output:
        "data/sorted_reads/{assembly}.{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule anvi_gen_contigs_database:
    input:
        "data/assemblies/{assembly}.fa"
    output:
        "data/anvio/{assembly}/{assembly}.db"
    shell:
        "anvi-gen-contigs-database -f {input} -o {output}"

rule anvi_run_hmms:
    input:
        "data/anvio/{assembly}/{assembly}.db"
    output:
        touch("data/anvio/{assembly}/{assembly}.db.run-hmms.done")
    params:
        threads=30
    shell:
        """
        anvi-run-hmms -c {input} --num-threads {params.threads}
        """

rule anvi_export_gene_calls:
    input:
        "data/anvio/{assembly}/{assembly}.db"
    output:
        "data/anvio/{assembly}/{assembly}.gene-calls.fa"
    shell:
        "anvi-get-dna-sequences-for-gene-calls -c {input} -o {output}"

rule anvi_run_centrifuge:
    input:
        fa="data/anvio/{assembly}/{assembly}.gene-calls.fa",
        db="data/anvio/{assembly}/{assembly}.db"
    output:
        hits="data/anvio/{assembly}/{assembly}.centrifuge_hits.tsv",
        report="data/anvio/{assembly}/{assembly}.centrifuge_report.tsv",
        done=touch("data/anvio/{assembly}/{assembly}.db.added_centrifuge.done")
    params:
        threads=8,
        centrifuge_base=CENTRIFUGE_BASE,
        centrifuge_models=CENTRIFUGE_BASE + '/b+h+v/b+h+v'
    shell:
        """
        export CENTRIFUGE_BASE={params.centrifuge_base}
        centrifuge -f --threads {params.threads} \
        -x {params.centrifuge_models} \
        {input.fa} \
        -S {output.hits} \
        --report-file {output.report}

        ln -s {output.hits} data/anvio/{wildcards.assembly}/centrifuge_hits.tsv
        ln -s {output.report} data/anvio/{wildcards.assembly}/centrifuge_report.tsv

        cd data/anvio/{wildcards.assembly}

        anvi-import-taxonomy -c {wildcards.assembly}.db \
        -i centrifuge_report.tsv centrifuge_hits.tsv \
        -p centrifuge

        rm centrifuge_hits.tsv
        rm centrifuge_report.tsv

        cd ../../../
        """

rule anvi_profile:
    input:
        sorted="data/sorted_reads/{assembly}.{sample}.bam",
        idx="data/sorted_reads/{assembly}.{sample}.bam.bai",
        db="data/anvio/{assembly}/{assembly}.db",
    output:
        aux="data/sorted_reads/{assembly}.{sample}.bam-ANVIO_PROFILE/AUXILIARY-DATA.h5",
        prof="data/sorted_reads/{assembly}.{sample}.bam-ANVIO_PROFILE/PROFILE.db",
        info="data/sorted_reads/{assembly}.{sample}.bam-ANVIO_PROFILE/RUNINFO.cp",
        log="data/sorted_reads/{assembly}.{sample}.bam-ANVIO_PROFILE/RUNLOG.txt"
    shell:
        """
        anvi-profile -i {input.sorted} \
        -c {input.db} \
        --overwrite-output-destinations \
        -o data/sorted_reads/{wildcards.assembly}.{wildcards.sample}.bam-ANVIO_PROFILE
        """

rule anvi_merge:
    input:
        profiles=lambda wildcards: expand("data/sorted_reads/{assembly}.{sample}.bam-ANVIO_PROFILE/RUNINFO.cp",
                        assembly=wildcards.assembly,
                        sample=SAMPLES),
        db="data/anvio/{assembly}/{assembly}.db"
        centrifuge_done="data/anvio/{assembly}/{assembly}.db.added_centrifuge.done",
        hmms_done="data/anvio/{assembly}/{assembly}.db.run-hmms.done"
    output:        
        aux="data/anvio/{assembly}/SAMPLES_MERGED/AUXILIARY-DATA.h5",
        prof="data/anvio/{assembly}/SAMPLES_MERGED/PROFILE.db",
        info="data/anvio/{assembly}/SAMPLES_MERGED/RUNINFO.mcp"
    shell:
        """
        anvi-merge {input.profiles} \
        -o data/anvio/{wildcards.assembly}/SAMPLES_MERGED \
        -c {input.db}
        """

rule anvi_summarize:
    input:
        prof="data/anvio/{assembly}/SAMPLES_MERGED/PROFILE.db",
        db="data/anvio/{assembly}/{assembly}.db"
    output:
        "data/anvio/{assembly}/{assembly}_SAMPLES-SUMMARY/index.html"
    shell:
        """
        anvi-summarize -p {input.prof} \
        -c {input.db} \
        -o data/anvio/{wildcards.assembly}/{wildcards.assembly}_SAMPLES-SUMMARY \
        -C CONCOCT
        """
