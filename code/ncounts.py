import numpy as np
import os
from sys import argv

if len(argv) > 1:
    s_b = argv[1]
else:
    s_b = 'all'

if s_b == 's':
    mod_directory = '/astro/store/gradscratch/tmp/tzdw/BPASS/singles/'
    out_str = 'single'
elif s_b == 'b':
    mod_directory = '/astro/store/gradscratch/tmp/tzdw/BPASS/binaries/'
    out_str = 'binary'
else:
    mod_directory = '/astro/store/gradscratch/tmp/tzdw/BPASS/all_mods/'
    out_str = 'allmod'

znames = ['zem5','zem4','z001','z002','z003','z004','z006','z008','z010','z014','z020','z030','z040']

for z in znames:
    
    print(z)
    
    models = np.zeros((51,171))
    for i in range(51):
        models[i,0] = 6.0+0.1*i #51 steps of 0.1 dex in logtime
        
    input_filename = '/astro/store/gradscratch/tmp/tzdw/BPASS/BPASSv2.1_imf135_300/input_bpass_{0}'.format(z)
    
    if os.path.isfile(input_filename):
        input_file = open(input_filename)
        totalimfnumber = float(input_file.readline())
        
        while True:
            
            modelname = input_file.readline()
            if modelname == '':
                break
                
            model_imftype = input_file.readline().split()
            modelimf = float(model_imftype[0])
            modeltype = int(model_imftype[1])
            if modeltype >= 1.9:
                mix_imfage = input_file.readline().split()
                mixedimf = float(mix_imfage[0])
                mixedage = float(mix_imfage[1])
                modelimf -= mixedimf
                modelimf = np.clip(modelimf,a_min=0,a_max=None)
                mixedimf = np.clip(mixedimf,a_min=0,a_max=None)
            else:
                mixedage = 0.0
            if modeltype >= 3.9:
                foo = input_file.readline() #stuff about black hole, not used
                    
            modelfile = modelname[26:].rstrip()
            
            if (s_b == 's') & (modelfile.split('/')[0] == 'NEWBINMODS'):
                continue
            elif (s_b == 'b') & (modelfile.split('/')[0] == 'NEWSINMODS'):
                continue
            
            if os.path.isfile(mod_directory+modelfile):
                
                model = open(mod_directory+modelfile)
                
                lasttime = 0.0
                lasttime2 = mixedage
                count = -1
                itlast = 0
                it2last = int(np.clip(np.rint(10.0 * np.log10(mixedage)) - 60.0,0,50))
                
                while True:
                    modelstep = model.readline().split()
                    
                    if modelstep == []:
                        break
                        
                    count += 1
                    time = float(modelstep[1])
                    dt = time - lasttime
                    if dt >= 0.0:
                        it = int(np.clip(np.rint(10.0*np.log10(time))-60,0,50)) #0.1 dex steps from 10^6 yrs, bound to 51 steps (50 + 1 for zero)

                        X_h = float(modelstep[10])
                        T = float(modelstep[3])
                        grav1 = 6.67259e-8 - 1.9891e33*float(modelstep[5])
                        grav2 = np.power((np.power(10.0,float(modelstep[2]))*6.9598e10),2.0)
                        gravity = np.log10(grav1/grav2)
                        C_O_He =(float(modelstep[12])/3.0+float(modelstep[14])/4.0)/float(modelstep[11])
                        logOfG = 3.676*T - 13.253 

                        if (X_h <= 0.4)&(T >= 4.45): 
                            if X_h >= 1e-3: 
                                startype=8 #WNH
                            else:
                                if C_O_He <= 0.03:
                                    startype = 9 #WN
                                else:
                                    startype = 10 #WC
                        else:
                            if T < 3.55:
                                startype = 7 #M
                            elif T < 3.66:
                                startype = 6 #K
                            elif T < 3.9:
                                startype = 5 #F/G
                            elif T < 4.041:
                                startype = 4 #A
                            elif T < 4.48:
                                startype = 3 #B
                            elif (gravity < logOfG) & (T >= 4.519):
                                startype = 2 #Of
                            else: #hotter than 4.48 and not Of
                                startype = 1 #O
                                
                        L = float(modelstep[4])
                        
                        #Now, we want every 0.1 dex below 5.0 give us a shift of 10 to the column number
                        L_chunk = int(np.clip(np.ceil((5.0 - L)/0.1),0,16)) #clipped to 0 for above 5, 16 for below 3.5

                        column_number=startype+(10*L_chunk)
                        
                        if it == itlast: #If model timestep is smaller than 0.1 dex in logtime 
                             models[it,column_number] += dt * modelimf
                        elif (it - itlast) > 1: #if its larger than 1 0.1 dex jump:
                            dt_first_step = np.power(10.0,6.05 + 0.1*itlast) - lasttime #last dt
                            dt_last_step = time - np.power(10.0,5.95 + 0.1*it) #this dt
                            models[itlast,column_number] += dt_first_step * modelimf
                            models[it,column_number] += dt_last_step * modelimf
                            for i in range(itlast+1,it):
                                dt_partial_step = np.power(10.0,6.05 + 0.1*i) - np.power(10.0,5.95 + 0.1*i)
                                models[i,column_number] += dt_partial_step * modelimf
                        else:
                            dt_first_step = np.power(10.0,5.95 + 0.1*it) - lasttime
                            dt_last_step = time - np.power(10.0,5.95 + 0.1*it)
                            models[itlast,column_number] += dt_first_step * modelimf
                            models[it,column_number] += dt_last_step * modelimf

                        if modeltype >= 1.9: #handle the secondary
                            mixedtime = time + mixedage
                            it2 = int(np.clip(np.rint(10.0*np.log10(mixedtime))-60,0,50))
                            if it2 == it2last:
                                models[it2,column_number] += dt * mixedimf
                                
                            elif (it2 - it2last) > 1: #if its larger than 1 0.1 dex jump:
                                dt_first_step = np.power(10.0,6.05 + 0.1*it2last) - lasttime2 
                                dt_last_step = mixedtime - np.power(10.0,5.95 + 0.1*it2) 
                                models[it2last,column_number] += dt_first_step * mixedimf
                                models[it2,column_number] += dt_last_step * mixedimf
                                for i in range(it2last+1,it2):
                                    dt_partial_step = np.power(10.0,6.05 + 0.1*i) - np.power(10.0,5.95 + 0.1*i)
                                    models[i,column_number] += dt_partial_step * mixedimf
                            else:
                                dt_first_step = np.power(10.0,5.95 + 0.1*it2) - lasttime2
                                dt_last_step = mixedtime - np.power(10.0,5.95 + 0.1*it)
                                models[it2last,column_number] += dt_first_step * mixedimf
                                models[it2,column_number] += dt_last_step * mixedimf
                                

               
                        lasttime=time
                        itlast=it
                        if modeltype >= 1.9:
                            it2last=it2
                            lasttime2=mixedtime

                
                model.close()
            
            else:
                
                print('failed on: '+mod_directory+modelfile)
        
        input_file.close()
        
    np.savetxt('ncounts_{}_{}.dat'.format(out_str,z),models)
    
    
