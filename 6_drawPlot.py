'''
Created on 2017. 11. 3.

@author: jmjung
'''
import pandas as pd
import matplotlib.pyplot as plt
import glob

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

def get_subAveToken(folder, marker):
    all_files = glob.glob(folder + "/result_place_R*.txt")
     
    sim_added_df = pd.DataFrame()
    for in_file in all_files:
        sim_df = pd.read_table(in_file, sep='\t', index_col = 0)
        sim_df = sim_df[marker]
        sim_added_df = sim_added_df.add(sim_df, fill_value=0)
    
    ## averaged markings
    sim_ave_df = sim_added_df/10.0
    return sim_ave_df

if __name__ == '__main__':
    gene_sym_dic = get_gene_sym()
    markers = ['5802','4088']
    GI1_df = get_subAveToken('./GeneInhibit/Sim_result_1/5802_PTPRS',markers)
    GI2_df = get_subAveToken('./GeneInhibit/Sim_result_2/5802_PTPRS',markers)
    GI3_df = get_subAveToken('./GeneInhibit/Sim_result_3/5802_PTPRS',markers)
    
    GI_df = GI1_df.add(GI2_df, fill_value=0)
    GI_df = GI_df.add(GI3_df, fill_value=0)
    GI_df = GI_df/3.0
    
    RF_df = get_subAveToken('./REFR/Sim_result/TH35',markers)
    
    MT_df = get_subAveToken('./MYST/Sim_result/TH35',markers)
    
    for ii in range(2):
        ax = plt.subplot(1,2,ii+1)
        
        entz = markers[ii]
        plt.plot(range(100),list(RF_df[entz])[1:],'darkslateblue',linewidth=4.0,label='reference')
        plt.plot(range(100),list(MT_df[entz])[1:],'indianred',linewidth=4.0,label='atrophic')
        plt.plot(range(100),list(GI_df[entz])[1:],'sandybrown',linewidth=4.0,label='PTPRS inhibited')
        if ii == 0:
            ax.set_ylim(-1,70)
        elif ii == 1:
            ax.set_ylim(25,60)
        plt.title(gene_sym_dic[entz],fontsize = 20)
        plt.xticks(fontsize=17)
        plt.yticks(fontsize=17)
        plt.xlabel('block', fontsize=20)
        plt.ylabel('status (token)', fontsize=20)
        if ii==0:
            plt.legend(loc='center left', bbox_to_anchor=(0, 0.7),fontsize=20)
    plt.tight_layout()
    plt.show()
    
