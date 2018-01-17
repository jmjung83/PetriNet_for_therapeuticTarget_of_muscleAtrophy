'''
Created on 2014. 9. 14.

@author: jmjung
'''
import random

# obtain places and their initial values
def place_initialize(place_path, mode, UP_effect):
    print("## place initialization (%s)"%(mode))
    place_dic = {}
    
    place_file = open(place_path)
        
    head_flg = True
    for line in place_file.xreadlines():
        if head_flg: head_flg = False; continue
        
        place = line.replace('\n','').replace('\r','')

        if mode == "Z":
            place_dic[place] = 0
   
        elif mode == "R":
            place_dic[place] = random.randrange(0,UP_effect+1)
                        
    #print ('Place_num:'), (len(place_dic)) 
    place_file.close()
    
    return place_dic

# obtain transitions    
def transition_initialize(transition_path):
    print("## transition initialization")
    transition_dic = {}
    transition_file = open(transition_path)
    
    head_flg = True
    for line in transition_file.xreadlines():
        if head_flg: head_flg = False; continue
        
        line = line.replace('\n','').replace('\r','')
        
        transition = tuple(line.split('\t'))
        transition_dic[transition] = 0
    
    #print ('Transition_num:'), (len(transition_dic.keys())) 
    transition_file.close()
    
    return transition_dic

def get_characterized_proteins():

    file_name = "../network/KEGG/KEGG_relations_expansion_tmp.txt"

    ###########
    phos_set = set(); deph_set = set(); actv_set = set(); inhi_set = set()
    
    in_file = open(file_name)
    
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
        
        line = line.replace('\n','').replace('\r','')
        [pathway, lt, rt, rel_type, rel_name, subtype_val] = line.split("\t")
        
        if rel_type != "PPrel": continue
        if ("/" in lt) or ("/" in rt): continue
        
        lt = lt.replace("hsa:","")
        rt = rt.replace("hsa:","")
        
        if rel_name == "phosphorylation"  : phos_set.add((lt,rt))
        if rel_name == "dephosphorylation": deph_set.add((lt,rt))
        if rel_name == "activation"       : actv_set.add((lt,rt))
        if rel_name == "inhibition"       : inhi_set.add((lt,rt))
                 
    in_file.close()
    
    phos_inhi_set = phos_set&inhi_set
    deph_actv_set = deph_set&actv_set
    
    rt_phos_inhi_set = set()
    
    for (lt,rt) in phos_inhi_set:
        rt_phos_inhi_set.add(rt)

    for (lt,rt) in deph_actv_set:
        rt_phos_inhi_set.add(rt)
    
    return rt_phos_inhi_set



