

### Main Function ###
if __name__ == '__main__':
    rel_set = set()
    lt_set = set()
    in_file = open("KEGG_relations_expansion_tmp.txt")
     
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
        
        line = line.replace('\n','').replace('\r','')
        [pathway, lt, rt, rel_type, rel_name, subtype_val] = line.split("\t")
        
        if rel_type != "PPrel": continue
        if ("/" in lt) or ("/" in rt): continue
        
        lt = lt.replace("hsa:","")
        rt = rt.replace("hsa:","")
        
        if rel_name == "phosphorylation"  : rel_set.add((lt, "phos", rt)); lt_set.add(lt) 
        if rel_name == "dephosphorylation": rel_set.add((lt, "deph", rt)); lt_set.add(lt)  
                 
    in_file.close()
    
    print len(lt_set)
    ########
    out_file = open("phosDeph_network(KEGG).txt", "w+")
    out_file.write("\t".join(["kinase (phosphatase)", "phos/deph", "substrate"]) + "\n")
    for rel in rel_set:
        out_file.write("\t".join(rel) + "\n")
    out_file.close()   