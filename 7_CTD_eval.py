'''
Created on 2017. 11. 14.

@author: jmjung
'''
import networkx as nx
import operator
import pandas as pd
import numpy as np
from scipy import stats

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

def get_connectedPlace_to_fiveMarkers(transition_file,markers):
    ## directed network generation
    DG=nx.DiGraph()
    in_file = open(transition_file)
     
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
         
        line = line.replace("\n","").replace("\r","")
        [TR_ID, rel_type, lt, rt] = line.split("\t")
             
        DG.add_edge(lt, rt)
     
    in_file.close()
    
    ## check if a node is connected to the five markers
    conNode_list = []
    for start in DG.nodes():
        pathExistFlg = True
        for end in markers:
            if nx.has_path(DG, start, end) == False:
                pathExistFlg = False
         
        if pathExistFlg == True:
            conNode_list.append(start)
    return conNode_list,DG

def get_target_candidate_SPL():
    markers = ["4087", "4088", "2308", "2309", "5970"]
    INH_PLC_list,DG = get_connectedPlace_to_fiveMarkers("network/transition.txt",markers)
    
    SPL_dic = {}
    for start in INH_PLC_list:
        SPL = 0
        for end in markers:
            SPL += nx.shortest_path_length(DG, source=start, target=end)
        SPL_ave = SPL/5.0
        SPL_dic[start] = SPL_ave
    
    SPL_TG = set()
    for key,val in SPL_dic.items():
        if val <=2.0:
            SPL_TG.add(key)
    
    return SPL_TG, INH_PLC_list

if __name__ == '__main__':
    entz_dic = get_gene_sym()
    ## target candidate from petrinet
    target_df = pd.read_table("target_candidates.txt", sep='\t', index_col = 0)
    PET_TG  = set([str(gene) for gene in target_df['entrez']])

    ## target candidate from shortest path length(SPL)
    SPL_TG, INH_PLC_list = get_target_candidate_SPL()
    
    ## comparison to shortest path analysis
    ALL = set(INH_PLC_list)
    HIT = PET_TG&SPL_TG
    pval = stats.hypergeom.sf(len(HIT),len(ALL),len(PET_TG),len(SPL_TG))
    print("Petrinet & shortest path", pval,[entz_dic[entz] for entz in HIT])

    ## get gold standard from CTD database
    gene_df = pd.read_table('./CTD/CTD_genes_related_to_atrophy.txt', sep='\t')
    
    ## get z-score of inferenece score in CTD
    infScore_ds = gene_df['InferenceScore']
    scoreList = list(infScore_ds.dropna())
    
    scoreDic = {}
    zScoreList = stats.zscore(scoreList)
    for ii in range(len(scoreList)):
        scoreDic[scoreList[ii]] = zScoreList[ii]
        
    GS_dic = {}
    for index, row in gene_df.iterrows():
        if pd.isnull(row['DirectEvidence'])==False:
            if str(row['GeneID']) == '3479': # use IGF1R instead of IGF1
                GS_dic['3480'] = 100
            else:
                GS_dic[str(row['GeneID'])] = 100

    for index, row in gene_df.iterrows():
        if pd.isnull(row['DirectEvidence'])==True:
            if GS_dic.has_key(str(row['GeneID'])):
                continue
            GS_dic[str(row['GeneID'])] = row['InferenceScore']
    
    ## gold standard threshold by its number
    for TH in [20,25,30,35,40]:
        print(TH)
        GS = set()
        for gene, sc in sorted(GS_dic.items(), key=operator.itemgetter(1), reverse=True):
            if gene in INH_PLC_list:
                GS.add(gene)
                
            if len(GS) == TH:
                print(sc, scoreDic[sc])
                break
        
        ALL = set(INH_PLC_list)
        TG = PET_TG
        HIT = TG&GS
        pval = stats.hypergeom.sf(len(HIT),len(ALL),len(GS),len(TG))
        print("Petrinet & gold standard from CTD", len(HIT),len(ALL),len(GS),len(TG), pval,[entz_dic[entz] for entz in HIT])

    
