# snakemake_anvio
snakemake rules for running anvio

This Snakefile will run a basic [anvi'o](http://merenlab.org/software/anvio/) analysis
of a set of assemblies and a set of samples, using the [snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home)
workflow manager. 

## Files

### Snakefile
The file containing the snakemake rules

### config.yaml.template
An example file containing the specific sample and assembly names that will
be used in the workflow, as well as additional variable specifications.

Must be renamed 'config.yaml' for use.

### cluster.json
A file specifying compute cluster job resource parameters for particular
rules in the Snakefile.

### requirements.txt
A file specifying conda requirements for the snakefile conda environment

### launch.sh
Example launch command for launching on a Torque (qsub) job manager

### anvio_install.sh
The commands I used for setting up my Anvi'o conda environment on Centos6

### dag.svg
Visualization of the workflow steps run when executing the example config.yaml

### data [directory]
Fake example data structure. For your own data, replace the files in data/assemblies
and data/samples with the assembled contigs and the per-sample paired-end gzipped
fastqs, respectively. Files must be named [assembly].fa and [sample].R[1,2].fq.gz 
