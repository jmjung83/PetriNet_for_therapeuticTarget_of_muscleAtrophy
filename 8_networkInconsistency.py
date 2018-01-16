'''
Created on 2017. 11. 28.

@author: jmjung
'''
import networkx as nx

### Main Function ###
if __name__ == '__main__':

    ## directed network generation
    DG=nx.DiGraph()
    in_file = open("./network/transition.txt")
     
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
         
        line = line.replace("\n","").replace("\r","")
        [TR_ID, rel_type, lt, rt] = line.split("\t")
             
        DG.add_edge(lt, rt,rel=rel_type)
     
    in_file.close()
    
    len2_path_cnt = 0
    dir_path_cnt = 0
    for node in DG.nodes():
        path_dic = nx.single_source_shortest_path(DG, node, cutoff=2)
        for key,val in path_dic.items():
            if len(val) == 3:
                print(val[0],val[2])
                len2_path_cnt+=1
                if DG.has_edge(val[0],val[2]):
                    dir_path_cnt+=1
    print(len2_path_cnt)
    print(dir_path_cnt)
    
    
    
    
    
    
