# make assembly folder with links
# make samples folder with links

# index assemblies

# map samples
# sort and index samples

# make anvio db

# annotate sc genes

# export gene calls

# centrifuge

# profile

# merge

# summarize



SAMPLES = ["A", "B"]


rule all:
    input:
        "report.html"

# =expand("sorted_reads/{sample}.bam", sample=config["samples"]),

rule link_files:
    input:
        assembly=

rule bowtie2_index:
    input:
        fa=expand("data/assemblies/{assembly}.fa", assembly=config["assemblies"])
    output:
        idx=expand("data/assemblies/{assembly}.rev.2.bt2", assembly=config["assemblies"])
    params:
        assembly=config["assemblies"]
    run:
        shell("bowtie2-build {input.fa} {params.assembly}")


rule bowtie2_map:
    input:
        fa=expand("data/assemblies/{assembly}.fa", assembly=config["assemblies"]),
        idx=expand("data/assemblies/{assembly}.rev.2.bt2", assembly=config["assemblies"]),
        R1=expand("data/samples/{sample}.R1.fa.gz", sample=config["samples"]),
        R2=expand("data/samples/{sample}.R2.fa.gz", sample=config["samples"])
    params:
        assembly=config["assemblies"]
        threads=8
    output:
        bam=expand("data/mapped_reads/{assembly}.{sample}.bam",
                   assembly=config["assemblies"],
                   sample=config["sample"])
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("bowtie2 -x {params.assembly} -p {params.threads} --no-unal "
                  "-q -1 {input.R1} -2 {input.R2} | "
                  "samtools view -bS - > {output.bam}")


rule samtools_sort:
    input:
        bam=expand("data/mapped_reads/{assembly}.{sample}.bam",
                   assembly=config["assemblies"],
                   sample=config["sample"])
    output:
        bam=expand("data/sorted_reads/{assembly}.{sample}.bam",
                   assembly=config["assemblies"],
                   sample=config["sample"])
    shell:
        "samtools sort -O bam -o {output.bam} {input.bam}"


rule samtools_index:
    input:
        bam=expand("data/sorted_reads/{assembly}.{sample}.bam",
                   assembly=config["assemblies"],
                   sample=config["sample"])
    output:
        bai=expand("data/sorted_reads/{assembly}.{sample}.bam.bai",
                   assembly=config["assemblies"],
                   sample=config["sample"])
    shell:
        "samtools index {input.bam}"


rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"


rule report:
    input:
        "calls/all.vcf"
    output:
        "report.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as vcf:
            n_calls = sum(1 for l in vcf if not l.startswith("#"))

        report("""
        An example variant calling workflow
        ===================================

        Reads were mapped to the Yeast
        reference genome and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], T1=input[0])