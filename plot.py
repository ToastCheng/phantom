#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 16:13:42 2017

@author: toast
"""

import numpy as np
import matplotlib.pyplot as plt
import sys



#set the parameter
SIZE_MIN = 0.9
SIZE_MAX = 1.1
STEP = 0.01
#if there is argv pass in, then consider argv first
if len(sys.argv)>1:
    SIZE_MIN = float(sys.argv[1])
if len(sys.argv)>2:
    SIZE_MAX = float(sys.argv[2])
if len(sys.argv)>3:
    STEP = float(sys.argv[3])



NUM_STEP = int((SIZE_MAX - SIZE_MIN)/STEP)
NUM_WAVELENGTH = 44

wavelength = [400,410,414,420,430,440,450,460,470,480,
	490,500,510,520,530,540,542,550,560,570,576,580,590,600,610,620,
	630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800]

wavelength32=[400,410,414,420,430,440,450,460,470,480,
	490,500,510,520,530,540,542,550,560,570,576,580,590,600,610,620,
	630,640,650,660,670,680]



# a 3% CV gaussian distribution
def distribution(x):
    return 13.298*np.exp(-555.6*(x-1)**2)

#create an array of size
size_range = [SIZE_MIN + STEP*x for x in range(NUM_STEP+1)]

#weights of each size
weight = np.array([distribution(np.array(size_range))])
weight /= weight.sum()



#for storing the weighted spectrum
#22,42,65 stands for 0.022 0.042 0.065cm SDS
weighted22 = np.zeros(NUM_WAVELENGTH)
weighted42 = np.zeros(NUM_WAVELENGTH)
weighted65 = np.zeros(NUM_WAVELENGTH)


for index in range(NUM_STEP+1):
    
    input_file = "output_DRS/simulation_" + str(index) + ".txt"
    
    data22 = []
    data42 = []
    data65 = []
    
    with open(input_file,'r') as f:
        row = f.read().split('\n')
        del row[NUM_WAVELENGTH]
        for element in row:
            element = element.split('\t')
            data22.append(float(element[1]))
            
            data42.append(float(element[2]))
            data65.append(float(element[3]))
    
    data22 = np.array(data22)
    data42 = np.array(data42)
    data65 = np.array(data65)
    
    weighted22 += data22*weight[0][index]
    weighted42 += data42*weight[0][index]
    weighted65 += data65*weight[0][index]
    
    
    
    size = SIZE_MIN + index*STEP
    size = np.floor(size*100)/100 # for getting rid of some weird little value e.g. 0.940000001micron
    file_name = str(size) + 'micron'
    
    fig, ax = plt.subplots(figsize=(10,6))
    plt.xlabel('wave length[nm]')
    plt.ylabel('reflectance[-]')
    plt.ylim(0,0.0006)
    plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0)) #scientific notation
    
    plt.title(file_name)
    plt.plot(wavelength,data22,'r',label='SDS=0.022')
    plt.plot(wavelength,data42,label='SDS=0.042')
    plt.plot(wavelength,data65,'g',label='SDS=0.065')
    plt.legend()   
    plt.savefig('output_figure/'+file_name+'.png')
    plt.close()

##########################################################################################################
#plot the weighted spectrum and save it into file
    
fig2, ax2 = plt.subplots(figsize=(10,6))
plt.xlabel('wave length[nm]')
plt.ylabel('reflectance[-]')
plt.ylim(0,0.0006)
plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0)) #scientific notation
    
plt.title("weighted spectrum")
plt.scatter(wavelength,weighted22,label='SDS=0.022')
plt.plot(wavelength,weighted22,'r',label='SDS=0.022')
plt.plot(wavelength,weighted42,label='SDS=0.042')
plt.plot(wavelength,weighted65,'g',label='SDS=0.065')
plt.legend()
    
plt.savefig('output_figure/'+'weighted'+'.png')
plt.close()


###########################################################################################################
#save the weighted spectrum data into file

with open("output_DRS/weighted_0.txt","w") as output:
    for i in range(NUM_WAVELENGTH):
        output.write(str(wavelength[i])+"\t"+str(weighted22[i])+"\t"+str(weighted42[i])+"\t"+str(weighted65[i])+"\n")








