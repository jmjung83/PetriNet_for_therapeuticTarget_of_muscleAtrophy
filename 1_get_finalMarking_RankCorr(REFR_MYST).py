'''
Created on 2016. 3. 4.

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

def getFinalMark_and_RankCor(STATE):
    gene_sym_dic = get_gene_sym()
    
    DIR = "./%s/Sim_result/"%(STATE)
    folder_list = [x[0] for x in os.walk(DIR) if x[0] != DIR]
    
    print("###### %s state"%(STATE))
    print("#### Get 10 final markings of 10 iterations")
    print("#### Calculate 45 rank correlations among them and average the correlations")
    print("#### Get and save their averaged final markings in terms of 6 predefined biomarkers")
    print("## Enabling threshold || Ave.rank correlation among 10 iterations")
    for fDir in folder_list:
        finalMark_Dlist = []
        
        all_files = glob.glob(fDir + "/result_place_R*.txt")
        
        finalMark_stacked_df = pd.DataFrame()
        for in_file in all_files:
            sim_df = pd.read_table(in_file, sep='\t', index_col = 0)
            
            finalMark_list = list(sim_df.loc[100])
            finalMark_Dlist.append(finalMark_list)
            
            finalMark_df = sim_df.ix[[100],["4087", "4088", "5970", "207", "2308", "2309"]]
            finalMark_stacked_df = finalMark_stacked_df.append(finalMark_df)
        
        ## Ave. rank correlation        
        rankCor_2D = scipy.stats.spearmanr(finalMark_Dlist, axis=1)[0]
        rankCor_2D_ave = (sum(sum(rankCor_2D))-len(all_files))/(len(all_files)*(len(all_files)-1))
        print("%s || %s")%(fDir.split('/')[-1],rankCor_2D_ave)
        
        ## averaged final markings of the selected enabling threshold, 35
        finalMark_ave_ds = finalMark_stacked_df.mean(axis=0)
        if fDir.split('/')[-1]=="TH35":
            out_file = open("./%s/final_marking_bioMarker(TH35).txt"%(STATE),'w+') 
            for (entz, ave) in finalMark_ave_ds.iteritems():
                out_file.write("%s\t%s\t%s\n"%(entz,gene_sym_dic[entz],ave))
                
            out_file.close()
if __name__ == '__main__':
    getFinalMark_and_RankCor("REFR")
    getFinalMark_and_RankCor("MYST")
    