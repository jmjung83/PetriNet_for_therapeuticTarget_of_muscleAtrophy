'''
Created on 2016. 12. 8.

@author: jmjung
'''

import gseapy
import glob
import pandas as pd

if __name__ == '__main__':
    target_df = pd.read_table('target_candidates.txt', sep='\t', index_col = 0)
    target_list = list(target_df["gene"])
    gseapy.enrichr(gene_list=target_list, description='BP2015', gene_sets='GO_Biological_Process_2015', outdir='enrichedGOTerms', cutoff=0.001)
    gseapy.enrichr(gene_list=target_list, description='KEGG2016', gene_sets='KEGG_2016', outdir='enrichedGOTerms')

