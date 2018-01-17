'''
Created on 2014. 9. 14.

@author: jmjung
'''
import time
import os
import networkx as nx
import shutil
from v1_Sim_petrinet import petriNet_simulator
from v1_Sim_initialize import place_initialize, transition_initialize, get_characterized_proteins

def measure(init_time):    
    after_time=time.time()
    dif_time=after_time-init_time                                     
    hour=int(dif_time/3660)
    mins=int((dif_time-hour*3660)/60)
    sec=dif_time-hour*3660-mins*60                                 
    print 'Processing Time:'+str(hour) +' hour\t'+str(mins) +' min\t'+str(sec) +' sec'

def get_kinPhos(file_name):
    ##########################
    
    in_file = open(file_name)
    
    kin_set = set()
    pht_set = set() 

    header = True
    for line in in_file.xreadlines():
        if header: header=False; continue
        
        [TNum, rel_type, lt, rt] = line.split("\t")
        
        if rel_type == "phos":  kin_set.add(lt)
        if rel_type == "deph":  pht_set.add(lt)
        
    in_file.close()
    
    return kin_set, pht_set

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

def get_connectedPlace_to_fiveMarkers(transition_file):
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
        for end in ["4087", "4088", "2308", "2309", "5970"]:
            if nx.has_path(DG, start, end) == False:
                pathExistFlg = False
         
        if pathExistFlg == True:
            conNode_list.append(start)
    return conNode_list

### Main Function ###
if __name__ == '__main__':
    ITERATION = 10
    TOTAL_STEP = 100
    
    INIT_TARGETS = ["90","91","92", "93"]
    THRESH = 35
    
    UP_effect = THRESH*2
    mode = "R"
    
    gene_sym = get_gene_sym("../gene(GE)_HomoSapiens.txt")
    rt_phos_inhi_set = get_characterized_proteins()
    INH_PLC_list = get_connectedPlace_to_fiveMarkers("../network/transition.txt")
 
    print("######## Gene-inhibited state")
    for INH_PLC in INH_PLC_list:
        print("###### Gene: %s (%s)"%(INH_PLC,gene_sym[INH_PLC]))
        print("###### enabling threshold %s, %s iteration "%(THRESH,ITERATION))
        cur_path = "./Sim_result/%s_%s"%(INH_PLC, gene_sym[INH_PLC])
        if not os.path.exists(cur_path):
            os.makedirs(cur_path)
            
        for ITER in range(ITERATION):
            print("####  %s iteration"%(ITER))
            init_time=time.time()
            place_dic = place_initialize("../network/place.txt", mode, UP_effect)
            transition_dic =  transition_initialize("../network/transition.txt")
                 
            petriNet_simulator(place_dic,
                               transition_dic,
                               INIT_TARGETS,
                               TOTAL_STEP,
                               ITER,
                               THRESH,
                               UP_effect,
                               rt_phos_inhi_set,
                               mode,
                               cur_path,
                               INH_PLC)
                   
            measure(init_time)


    
    
