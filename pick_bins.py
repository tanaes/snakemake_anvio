#!/usr/bin/env python

import pandas as pd
import os
from shutil import copyfile


def calc_scores(bins_df, cov_df, sample_name, C = 1, R = 1, A = 10):
    join_df = bins_df.join(cov_df)

    join_df['score'] = join_df['percent_complete'] * C - 
                   join_df['percent_redundancy'] * R +
                   join_df[sample_name] / join_df[sample_name].max() * A

    return(join_df)

def copy_fastas(bins_dir, bins, out_dir):
    for bin_name in bins:
        copyfile(os.path.join(bins_dir,bin_name,'{}-contigs.fa'.format(bin_name)), 
                 os.path.join(out_dir,'{}-contigs.fa'.format(bin_name)))
    return()


bins_dir = snakemake.input['bins_dir']
bin_summary_fp = snakemake.input['bins_fp']
coverage_fp = snakemake.input['cov_fp']

out_dir = snakemake.output['bins_outdir']
out_info_fp = snakemake.output['bins_info']

C = snakemake.config["C"]
R = snakemake.config["R"]
A = snakemake.config["A"]
N = snakemake.config["N"]

assembly = snakemake.wildcards['assembly']
sample = snakemake.wildcards['sample']

# get bin summary -- sort by completeness
# bins    taxon   total_length    num_contigs     N50     GC_content      percent_complete        percent_redundancy
# Bin_10  Unknown 1724296 212     8076    43.944531913    25.4583543022   7.08043819317
# Bin_11  Unknown 2818857 203     16305   45.52894129     68.8161950567   1.03008134033

bins_df = pd.read_csv(bin_summary_fp, sep='\t', header=0, index_col = 0)

# get mean coverage in native sample
# bins    s18772_18772_TRUSEQNANO_400BP_S8_L001   s18772_18783_TRUSEQNANO_400BP_S5_L001
# Bin_10  3.20549815941   0.82594390421

cov_df = pd.read_csv(coverage_fp, sep='\t', header=0, index_col = 0)

# find top N bins with less than R redundancy 
# score: C * completion  - (R * redundancy) +  A * (coverage /  max_coverage) 
# C = 1, R = 1, A = 10

sample_name = 's{0}_{1}'.format(assembly.upper(), sample.upper())

bins_df = calc_scores(bins_df, cov_df, sample_name, C = C, R = R, A = A)

topN = bins_df.nlargest(N, 'score')

topN_names = topN.index

copy_fastas(bins_dir, bins, out_dir)

topN.to_csv(out_info_fp, sep='\t')



