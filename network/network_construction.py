'''
Created on 2016. 2. 22.

@author: jmjung
'''


def get_phosDeph_rel(rel_set, dup_cnt, file_name):
    
    in_file = open(file_name)
    
    line_cnt = 0
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg = False; continue
        
        line_cnt += 1
        line = line.replace('\n','').replace('\r','')
        rel = tuple(line.split("\t"))
        
        if rel in rel_set: dup_cnt += 1
        rel_set.add(rel)
    
    print '# of rel in %s: %s'%(file_name.split("/")[1], line_cnt)
        
    in_file.close()
    return dup_cnt

def place_transition_maker(rel_set):
    ##########################
    out_file_trans = open("transition.txt", "w+")
    out_file_trans.write("Transition_ID\tRelation\tLeft\tRight\n")

    place_set = set()
    trans_cnt = 0
    
    phos_cnt = 0
    deph_cnt = 0
    
    kin_set = set()
    pht_set = set() 

    for rel in rel_set:
        [lt, rel_type, rt] = rel

        if rel_type == "phos":  phos_cnt += 1; kin_set.add(lt)
        if rel_type == "deph":  deph_cnt += 1; pht_set.add(lt)
        
        place_set.add(lt)
        place_set.add(rt)
        
        trans_cnt += 1
        out_file_trans.write("\t".join(["T"+repr(trans_cnt), rel_type, lt, rt]) + "\n")

    out_file_trans.close()

    print "# of total relations in myocyte-specific network (transitions): %s"%(trans_cnt)
    print "    # of Phos Rel.: %s"%(phos_cnt)
    print "    # of Deph Rel.: %s"%(deph_cnt)            
    print "    # of proteins(places): %s"%(len(place_set))
    print "        # of kinases: %s"%(len(kin_set))
    print "        # of phosphatases: %s"%(len(pht_set))
    print "        # of both kinases and phosphatases: %s"%(len(kin_set&pht_set))
    print "        # of substrates only: %s"%(len(place_set) - len(kin_set|pht_set))

    ##########################
    out_file_place = open("place.txt", "w+")
    out_file_place.write("Entrez\n")
    
    for entz in place_set:
        out_file_place.write(entz + "\n")
    
    out_file_place.close()
    
    return place_set, kin_set, pht_set

def get_characterized_proteins():

    file_name = "./KEGG/KEGG_relations_expansion_tmp.txt"

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
    

    phos_actv_set = phos_set&actv_set; phos_inhi_set = phos_set&inhi_set
    deph_actv_set = deph_set&actv_set; deph_inhi_set = deph_set&inhi_set
    
    rt_phos_inhi_set = set()
    rt_phos_actv_set = set()
    
    for (lt,rt) in phos_inhi_set:
        rt_phos_inhi_set.add(rt)

    for (lt,rt) in deph_actv_set:
        rt_phos_inhi_set.add(rt)

    for (lt,rt) in phos_actv_set:
        rt_phos_actv_set.add(rt)

    for (lt,rt) in deph_inhi_set:
        rt_phos_actv_set.add(rt)
    
    return rt_phos_inhi_set, rt_phos_actv_set

def remove_conf_rel(rel_set):
    rel_list = list(rel_set)[:]
    
    phos_cnt = 0
    deph_cnt = 0
    for rel in rel_list:
        [lt, rel_type, rt] = rel
        
        if (lt, rel_type) in [('5499','phos'), ('5803','phos'), ('7419','deph'), ('7419','phos'),
                              ('4659','phos'), ('5500','phos'), ('5501','phos')]:
            rel_set.remove(rel)
        
        else:
            if rel_type == "phos":  phos_cnt += 1;
            if rel_type == "deph":  deph_cnt += 1;
    
    print "# of total relations without Dup. and Conf.: %s"%(len(rel_set))
    print "    # of Phos Rel.: %s"%(phos_cnt)
    print "    # of Deph Rel.: %s"%(deph_cnt)

def get_BSML_gene_onto():
    BSML_gene_dic = {}
    ##### gene
    in_file = open("./PhosphoNet/gene(GE)_HomoSapiens.txt")
    
    header = True
    for line in in_file.xreadlines():
        if header: header=False; continue;
        line = line.replace("\n","").replace("\r","")
        [BISL_ID, Entrez_gene_ID] = line.split("\t")[:2]
        
        BSML_gene_dic[BISL_ID] = Entrez_gene_ID 
    
    in_file.close()
    
    return BSML_gene_dic

def get_myocyte_specific_gene(file_name):
    BSML_gene_dic = get_BSML_gene_onto()
    ND_gene_set = set()
    
    in_file = open(file_name)
    head_flg = True
    
    tot_myocyte_gene = 0
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
        
        line = line.replace('\n','').replace('"','')
        line_list = line.split('\t')
        
        gene  = line_list[0]
        organ = line_list[1]
        cell  = line_list[2]
        level = line_list[3]
        conf  = line_list[5]
        
        if gene == "NA": continue
        if conf == "Uncertain": continue
        if (organ == "skeletal muscle") and (cell == "myocytes"):
            tot_myocyte_gene += 1
            if level == "Not detected":
                ND_gene_set.add(BSML_gene_dic[gene])
        
    in_file.close()

    return ND_gene_set

def get_myocyte_specific_rel(rel_set):
    ND_gene_set = get_myocyte_specific_gene("HPA_protein_expression.tsv")
    
    rel_set_list = list(rel_set)[:]
    
    for rel in rel_set_list:
        [lt, rel_type, rt] = rel
        if (lt in ND_gene_set) or (rt in ND_gene_set):
            rel_set.remove(rel)
    
    
### Main Function ###
if __name__ == '__main__':
    targets = ["90","91","92", "93"]
    rel_set = set()
    dup_cnt = 0
    
    dup_cnt = get_phosDeph_rel(rel_set, dup_cnt,"./KEGG/phosDeph_network(KEGG).txt")
    dup_cnt = get_phosDeph_rel(rel_set, dup_cnt,"./PhosphoNet/phos_network(KSI).txt")
    dup_cnt = get_phosDeph_rel(rel_set, dup_cnt,"./DEPOD/deph_network(DEPOD).txt")
    
    print "# of total relations without Dup.: %s"%(len(rel_set))
    print "    # of accumulated duplicates: %s"%(dup_cnt)
    
    remove_conf_rel(rel_set)
    
    get_myocyte_specific_rel(rel_set)
    
    place_set, kin_set, pht_set = place_transition_maker(rel_set)

    rt_phos_inhi_set, rt_phos_actv_set = get_characterized_proteins()

    print "# of kinases: %s"%(len(kin_set))
    print "    # of kinases activated without phos.", len(kin_set&rt_phos_inhi_set)
    print "# of phosphatases: %s"%(len(pht_set))
    print "    # of phosphatases activated without phos.", len(pht_set&rt_phos_inhi_set)
    
    print "# check available marker in networks"
    print "    (predefined bioMarker: 4087, 4088, 5970, 207, 2308, 2309)"
    print "    (myostatin targets: 90, 91, 92, 93)"
    print "    ## available markers in the network :",
    for marker in ["4087", "4088", "5970", "207", "2308", "2309", "90","91","92", "93"]:
        if marker in place_set:
            print marker,
    
    
    
    
    
