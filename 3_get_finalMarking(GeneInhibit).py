'''
Created on 2016. 4. 21.

@author: jmjung
'''
import glob
import numpy as np
import scipy.stats
import os
import pandas as pd

def get_gene_sym():
    in_file = open("gene(GE)_HomoSapiens.txt")
    
    entz_dic = {}
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
        
        line = line.replace("\n","").replace("\r","")
        [BSML, entz, sym, synon] = line.split("\t")[:4]
        
        entz_dic[entz] = sym
        
    in_file.close()
    
    return entz_dic

def getFinalMark(STATE,time):
    #gene_sym_dic = get_gene_sym()
    DIR = "./%s/Sim_result_%s/"%(STATE,time)
    folder_list = [x[0] for x in os.walk(DIR) if x[0] != DIR]
    
    all_finalMark_ave_df = pd.DataFrame()
    
    print("###### %s state %s time"%(STATE, time))
    print("#### Get 10 final markings of 10 iterations for each gene inhibition state")
    print("#### Get and save their averaged final markings in terms of 5 predefined biomarkers")
    print("## Entrez_GeneSymbol")
    for fDir in folder_list:
        all_files = glob.glob(fDir + "/result_place_R*.txt")
         
        finalMark_stacked_df = pd.DataFrame()
        for in_file in all_files:
            sim_df = pd.read_table(in_file, sep='\t', index_col = 0)
            finalMark_df = sim_df.ix[[100],["4087", "4088", "5970", "2308", "2309"]]
            finalMark_stacked_df = finalMark_stacked_df.append(finalMark_df)
 
        ## averaged final markings of the each gene inhibit state
        finalMark_ave_ds = finalMark_stacked_df.mean(axis=0)
        all_finalMark_ave_df[fDir.split('/')[-1]] = finalMark_ave_ds
        print(fDir.split('/')[-1])
        print(finalMark_ave_ds)
    
    pd.DataFrame.to_csv(all_finalMark_ave_df, path_or_buf='./GeneInhibit/final_marking_bioMarker(TH35)_%stime.txt'%(time),sep='\t')
    
if __name__ == '__main__':
    getFinalMark("GeneInhibit",1)
    getFinalMark("GeneInhibit",2)
    getFinalMark("GeneInhibit",3)
    