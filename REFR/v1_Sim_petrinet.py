'''
Created on 2014. 9. 14.

@author: jmjung
'''

import re
import random

def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def atoi(text):
    return int(text) if text.isdigit() else text

def petriNet_simulator(place_dic,
                         transition_dic,
                         INIT_TARGETS,
                         TOTAL_STEP,
                         ITER,
                         THRESH,
                         UP_effect,
                         rt_phos_inhi_set,
                         mode):
    
    print("## Start petrinet simulation for %s blocks"%(TOTAL_STEP))
    ##########
    # open place result file
    result_place_file = open("./Sim_result/th%s/result_place_%s%s.txt"%(THRESH, mode, ITER),'w+')
    
    # write a header (place ID) of the result place file
    sorted_place = place_dic.keys()
    sorted_place.sort(key=natural_keys)
    
    place_ID_str='step\t'
    for place_ID in sorted_place:
        place_ID_str += place_ID + '\t'

    result_place_file.write(place_ID_str[:-1] + '\n')
    
    # write the initial token for each place ID ("step: 0")
    write_place(0, sorted_place, place_dic, result_place_file)
    ##########
    
    ##########
    # open transition result file
    result_transition_file = open("./Sim_result/th%s/result_transition_%s%s.txt"%(THRESH, mode, ITER),'w+')

    # write a header (tansition ID) of the result transition file
    transition_ID_str='step\t'
    for transition_ID in range(len(transition_dic.keys())+1)[1:]:
        transition_ID_str += "T" + repr(transition_ID) + '\t'

    result_transition_file.write(transition_ID_str[:-1] + '\n')
    ##########
    
    ##########
    # run simulation
    for step in range(TOTAL_STEP+1)[1:]:

        transition_list = transition_dic.keys()
        
        # shuffling transitions
        index_list = range(len(transition_dic.keys()))
        random.shuffle(index_list)
        
        for index in index_list:
            transition = transition_list[index]
            
            # perform the transition
            perform_transition(transition, transition_dic, place_dic, THRESH, rt_phos_inhi_set, UP_effect)

        # normalization places
        normalization(place_dic, 0.1)
      
        # write the results of each step
        write_place(step, sorted_place, place_dic, result_place_file)

        # write transition weight & transition_dic initialization
        write_transition(step, transition_dic, result_transition_file)
    
    result_place_file.close()
    result_transition_file.close()
    print("## End petrinet simulation for %s blocks"%(TOTAL_STEP))

def write_place(step, sorted_place, place_dic, result_file):
    place_token_str = repr(step) + '\t'
    for place_ID in sorted_place:
        place_token_str += repr(place_dic[place_ID]) + '\t'
        
    result_file.write(place_token_str[:-1] + '\n')

def write_transition(step, transition_dic, result_transition_file):
    #transition = (T7, deph, 5788, 3718)
    
    transition_firing_list = [0]*len(transition_dic)
    for transition, marking in transition_dic.items():
        transition_firing_list[int(transition[0][1:])-1] = repr(marking)
        transition_dic[transition] = 0
    
    transition_str = repr(step) + '\t' + "\t".join(transition_firing_list) + "\n"
    result_transition_file.write(transition_str)
    
def perform_transition(transition, transition_dic, place_dic, THRESH, rt_phos_inhi_set, UP_effect):
    #transition = (T7, deph, 5788, 3718)
    [trans_ID, rel_type, lt, rt] = transition

    if lt in rt_phos_inhi_set:
        if place_dic[lt] < THRESH:
            if rel_type == "deph":
                place_dic[rt] -= 1
                if place_dic[rt] < 0: place_dic[rt] = 0
                transition_dic[transition] = -1
                 
            if rel_type == "phos":
                place_dic[rt] += 1
                if place_dic[rt] > UP_effect: place_dic[rt] = UP_effect
                transition_dic[transition] = 1
    else:
        if place_dic[lt] >= THRESH:
            if rel_type == "deph":
                place_dic[rt] -= 1
                if place_dic[rt] < 0: place_dic[rt] = 0
                transition_dic[transition] = -1
                 
            if rel_type == "phos":
                place_dic[rt] += 1
                if place_dic[rt] > UP_effect: place_dic[rt] = UP_effect
                transition_dic[transition] = 1


def normalization(place_dic, PERC):
    for place, val in place_dic.items():
        int_part = int(float(val)*float(PERC))
        
        new_val = val-int_part
        
        place_dic[place] = new_val


