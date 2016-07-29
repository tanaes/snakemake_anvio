#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "script.sh <out_dir (with config.yaml)> [additional snakemake arguments]*"
    exit
fi

snakemake -j 100 --local-cores 8 --cluster-config cluster.json --cluster "qsub -q {cluster.queue} -k eo -m n -l nodes=1:ppn={cluster.n} -l mem={cluster.mem}gb -l walltime={cluster.time}" --directory "$@"
#snakemake --directory ~/pipeline_test/exp1/ --dag | dot -Tsvg > dag.svg 
#-o /home/snurk/logs/ -e /home/snurk/logs/ 
