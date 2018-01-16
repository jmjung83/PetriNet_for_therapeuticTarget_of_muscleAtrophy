'''
Created on 2016. 4. 21.

@author: jmjung
'''

'''
Created on 2016. 3. 4.

@author: jmjung
'''
import glob
import numpy as np
import scipy.stats
import os
import pandas as pd

def get_gene_sym(in_file_name):
    in_file = open(in_file_name)
    
    entz_dic = {}
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
        
        line = line.replace("\n","").replace("\r","")
        [BSML, entz, sym, synon] = line.split("\t")[:4]
        
        entz_dic[entz] = sym
        
    in_file.close()
    
    return entz_dic

### Main Function ###
if __name__ == '__main__':
    sig_TH=1
    # get saturated tokens in the six proteins for the reference state
    dic_refr = {}
    in_file = open("./REFR/final_marking_bioMarker(TH35).txt")
    for line in in_file.xreadlines():
        line = line.replace("\n","").replace("\r","")
        [entz, sym, val] = line.split("\t")
        dic_refr[entz]=float(val)    
        
    in_file.close()

    # get saturated tokens in the six proteins for the myostatin state
    dic_myst = {}
    in_file = open("./MYST/final_marking_bioMarker(TH35).txt")
    for line in in_file.xreadlines():
        line = line.replace("\n","").replace("\r","")
        [entz, sym, val] = line.split("\t")
        dic_myst[entz]=float(val)    
        
    in_file.close()
    
    # get their difference
    stat_df = pd.DataFrame({"REFR":dic_refr,"MYST":dic_myst})
    stat_df["MYST-REFR"] = stat_df["MYST"] - stat_df["REFR"]
    
    # filtered in the 5 validated markers
    stat_df = stat_df.loc[["4087", "4088", "5970", "2308", "2309"]]
    stat_df.index = stat_df.index.astype(int)
    
    
    opp_dir_dic = {}
    for ii in range(1,4):
        # get saturated tokens in the 5 validated proteins in Gene-Inhibited state
        simGI_df = pd.read_table('./GeneInhibit/final_marking_bioMarker(TH35)_%stime.txt'%(ii), sep='\t', index_col = 0)
        if ii == 1:
            for gene in list(simGI_df):
                opp_dir_dic[gene] = [] 
        for gene in list(simGI_df):
            MR_diff = stat_df["MYST-REFR"]
            GM_diff = simGI_df[gene] - stat_df["MYST"]
            # get the number of the opposite direction markers
            opp_dir_num = sum(MR_diff*GM_diff<0)
            opp_dir_dic[gene].append(opp_dir_num)
    
    target_num = 0
    all_num = 0
    target_df = pd.DataFrame(columns=["entrez","gene","first","second","third"])
    all_df = pd.DataFrame(columns=["entrez","gene","first","second","third"])
    for key,val in opp_dir_dic.items():
        if sum(np.array(val)>=4)==3:
            target_num += 1
            target_df.loc[target_num] = key.split("_")+val
        
        all_num += 1
        all_df.loc[all_num] = key.split("_")+val
    print target_df.to_string()
    pd.DataFrame.to_csv(target_df, path_or_buf='target_candidates.txt',sep='\t')
    pd.DataFrame.to_csv(all_df, path_or_buf='all_proteins.txt',sep='\t')
    
    
