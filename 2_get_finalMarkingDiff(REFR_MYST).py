'''
Created on 2016. 4. 21.

@author: jmjung
'''

'''
Created on 2016. 3. 4.

@author: jmjung
'''
import glob
import numpy as np
import scipy.stats
import pandas as pd

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

### Main Function ###
if __name__ == '__main__':
    # get saturated tokens in the six proteins for the reference state
    dic_refr = {}
    in_file = open("./REFR/final_marking_bioMarker(TH35).txt")
    for line in in_file.xreadlines():
        line = line.replace("\n","").replace("\r","")
        [entz, sym, val] = line.split("\t")
        dic_refr["%s(%s)"%(sym,entz)]=float(val)    
        
    in_file.close()

    # get saturated tokens in the six proteins for the myostatin state
    dic_myst = {}
    in_file = open("./MYST/final_marking_bioMarker(TH35).txt")
    for line in in_file.xreadlines():
        line = line.replace("\n","").replace("\r","")
        [entz, sym, val] = line.split("\t")
        dic_myst["%s(%s)"%(sym,entz)]=float(val)    
        
    in_file.close()
    
    # get their difference
    diff_df = pd.DataFrame({"REFR":dic_refr,
                            "MYST":dic_myst})

    diff_df["MYST-REFR"] = diff_df["MYST"] - diff_df["REFR"]
    print(diff_df)
        