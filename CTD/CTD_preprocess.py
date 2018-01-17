'''
Created on 2017. 11. 10.

@author: jmjung
'''
import pandas as pd
if __name__ == '__main__':
    header = "GeneSymbol\tGeneID\tDiseaseName\tDiseaseID\tDirectEvidence\tInferenceChemicalName\tInferenceScore\tOmimIDs\tPubMedIDs"
    fout = open("CTD_genes_related_to_atrophy.txt",'w+')
    fout.write(header+'\n')
    

    fin = open('CTD_genes_diseases.tsv','r')
    
    line_cnt = 0
    for line in fin:
        if line.startswith('#'): continue
        
        line_cnt += 1
        if line_cnt % 1000000 == 0: print line_cnt #51,181,644
        
        atts = line.strip().split('\t')
        # D009133: Muscular Atrophy

        if atts[3] in ["MESH:D009133"]:
            fout.write(line)
        
    print("total_line",line_cnt)
    fin.close()
    fout.close()
    