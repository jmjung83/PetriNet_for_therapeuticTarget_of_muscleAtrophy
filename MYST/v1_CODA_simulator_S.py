'''
Created on 2014. 9. 14.

@author: jmjung
'''
import time
import os

from v1_Sim_petrinet import petriNet_simulator
from v1_Sim_initialize import place_initialize, transition_initialize, get_characterized_proteins

def measure(init_time):    
    after_time=time.time()
    dif_time=after_time-init_time                                     
    hour=int(dif_time/3660)
    mins=int((dif_time-hour*3660)/60)
    sec=dif_time-hour*3660-mins*60                                 
    print 'Processing Time:'+str(hour) +' hour\t'+str(mins) +' min\t'+str(sec) +' sec'

### Main Function ###
if __name__ == '__main__':
    
    ITERATION = 10
    TOTAL_STEP = 100
    
    INIT_TARGETS = ["90","91","92","93"]
    rt_phos_inhi_set = get_characterized_proteins()
    
    print("######## Atrophic state")
    for THRESH in [1,5,10,15,20,25,30,35,40,45,50]:
        print("###### enabling threshold %s, %s iteration "%(THRESH,ITERATION))
        UP_effect = THRESH*2
        if not os.path.exists("./Sim_result/" + "TH" + repr(THRESH)):
            os.makedirs("./Sim_result/" + "TH" + repr(THRESH))
             
        for mode in ["R"]:
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
                                   mode)
                 
                measure(init_time)


    
    
