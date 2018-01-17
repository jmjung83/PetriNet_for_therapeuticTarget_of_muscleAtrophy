'''
Created on 2016. 4. 8.

@author: jmjung
'''
import urllib2
from bs4 import BeautifulSoup

def crawling(deph):
    fish_url = "http://www.koehn.embl.de/depod/shows.php?name=" + deph.replace(" ","%20")
    
    page = urllib2.urlopen(fish_url)
    html_doc = page.read()
    soup = BeautifulSoup(html_doc)
    
    for ii in range(len(soup.find_all('tr')))[1:]:
        if "EntrezGene ID" in soup.find_all('tr')[ii].get_text():
            entz_string = soup.find_all('tr')[ii].get_text()
            entz_index = entz_string.index("EntrezGene ID")
            entz = entz_string[entz_index + len("EntrezGene ID"):]
            entz = repr(entz)
            entz = entz.replace("'","").replace("u","").replace(" ","")
            return entz.split(",")
    
    return []

def get_substrate_set(in_file_name):
    in_file = open(in_file_name)
    
    sub_set = set()
    
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
        line = line.replace("\n","").replace("\r","")
        [deph, sub, sub_type, org] = line.split("\t")[:4]
        
#         if "Homo sapiens (Human)" not in org:
#             continue
        if sub_type != "protein substrate":
            continue
        if sub == "p210 BCR-ABL":
            continue
        
        sub_set.add(sub)
    
    in_file.close()
    return sub_set

def get_sub_mapped_entrez(sub_set):
    sub_dic = {}
    for sub in sub_set:
        entz_list = crawling(sub)
        if len(entz_list) == 0:
            print "####### no mapped entrez ID for substrate when crawrling !!!! substrate: %s"%sub
        else:
            sub_dic[sub] = entz_list
        
    return sub_dic

def get_deph_mapped_entrez(in_file_name):
    in_file = open(in_file_name)
    
    deph_dic = {}
    
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
        line = line.replace("\n","").replace("\r","")
        [entz, name] = line.split("\t")[3:5]
        
        if not deph_dic.has_key(name):
            deph_dic[name] = [entz]
        else:
            deph_dic[name].append(entz)
    in_file.close()
    
    return deph_dic

#def dephos_net_con(rel_set):
### Main Function ###
if __name__ == '__main__':
    rel_set = set()
    sub_set = get_substrate_set("DEPOD_201408_human_phosphatase-substrate.txt")
    sub_dic = get_sub_mapped_entrez(sub_set)
    
    deph_dic = get_deph_mapped_entrez("DEPOD_201410_human_phosphatases.txt")
    ######
    in_file = open("DEPOD_201408_human_phosphatase-substrate.txt")
    
    head_flg = True
    for line in in_file.xreadlines():
        if head_flg: head_flg=False;continue
        line = line.replace("\n","").replace("\r","")
        [deph, sub, sub_type, org] = line.split("\t")[:4]
         
#         if "Homo sapiens (Human)" not in org:
#             continue
        if sub_type != "protein substrate":
            continue
        if sub == "p210 BCR-ABL":
            continue
        
        if not deph_dic.has_key(deph):
            print "no depho in dic", deph
            continue
            
        if not sub_dic.has_key(sub):
            print "no substrate in dic", sub
            continue
                    
        for deph_entz in deph_dic[deph]:
            for sub_entz in sub_dic[sub]:
                if sub_entz=="N/A": continue
                rel_set.add((deph_entz, "deph", sub_entz))
    
    in_file.close()
    
    ########
    out_file = open("deph_network(DEPOD).txt", "w+")
    out_file.write("\t".join(["kinase (phosphatase)", "phos/deph", "substrate"]) + "\n")
    for rel in rel_set:
        out_file.write("\t".join(rel) + "\n")
    out_file.close()

    
       