'''
Created on 2016. 4. 8.

@author: jmjung
'''
def get_entity_mapped_entrez(in_file_name):
    in_file = open(in_file_name)
    
    entz_dic = {}
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
        
        line = line.replace("\n","").replace("\r","")
        [BSML, entz, sym, synon] = line.split("\t")[:4]
        
        entz_dic[entz] = [sym] + synon.split("|")
        
    in_file.close()
    
    return entz_dic

def get_entity_set(in_file_name):
    in_file = open(in_file_name)
    
    entity_set = set()
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
        
        line = line.replace("\n","").replace("\r","")
        [kin, sub] = line.split("\t")[:2]
        
        entity_set.add(kin)
        entity_set.add(sub)
        
    in_file.close()
    
    return entity_set

### Main Function ###
if __name__ == '__main__':
    rel_set = set()
    entz_dic = get_entity_mapped_entrez("gene(GE)_HomoSapiens.txt")
    entity_set = get_entity_set("comKSI (used)_tab.txt")
    
    entity_entz_dic = {}
    
    ### main symbol mapping
    for entity in entity_set:
        for entz, sym_syn in entz_dic.items():
            if entity == sym_syn[0]:
                if not entity_entz_dic.has_key(entity):
                    entity_entz_dic[entity] = [entz]
                else:
                    entity_entz_dic[entity].append(entz)
                


    print len(entity_set)
    print len(entity_entz_dic.keys())
        
    ### synonym mapping for non-mapped entity      
    left_entity_set =  entity_set - set(entity_entz_dic.keys())
    
    for entity in left_entity_set:
        mapped_entz_set = set()
        for entz, sym_syn in entz_dic.items():
            if entity in sym_syn[1:]:
                mapped_entz_set.add(entz)
        
        if len(mapped_entz_set) == 1:
            entity_entz_dic[entity] = list(mapped_entz_set)

    print len(entity_entz_dic.keys())
    
    ########################################

    in_file = open("comKSI (used)_tab.txt")
    
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
         
        line = line.replace("\n","").replace("\r","")
        [kin, sub] = line.split("\t")[:2]
        
        if kin not in entity_entz_dic.keys(): continue 
        if sub not in entity_entz_dic.keys(): continue
        
        for kin_entz in entity_entz_dic[kin]:
            for sub_entz in entity_entz_dic[sub]:
                rel_set.add((kin_entz, "phos", sub_entz))
                   
    in_file.close()

    ########
    out_file = open("phos_network(KSI).txt", "w+")
    out_file.write("\t".join(["kinase (phosphatase)", "phos/deph", "substrate"]) + "\n")
    for rel in rel_set:
        out_file.write("\t".join(rel) + "\n")
    out_file.close()    
