#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The wrapper for group jobs of peak finding
Usage:
    wrapper <folder_name> <thld> <min_pix>

"""

import sys,os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('/Users/chufeng/Desktop/BSCCO/NSLS_IND')
import NSLS_IND_util as NSLS


def file_name_gen(img_file_name):
    
    img_file_name=os.path.abspath(img_file_name)
    img_file_name_short=os.path.basename(img_file_name)
    base_name=img_file_name_short[0:img_file_name_short.rfind('_')+1]
    id_num=int(img_file_name_short[img_file_name_short.rfind('_')+1\
                               :img_file_name_short.rfind('.')])
    id_str='%d'%(id_num)
    
    #id_str='%03d'%(id_num)
    txt_file_name_short=base_name+id_str+'.txt'
    txt_file_name=os.path.dirname(img_file_name)+'/'+txt_file_name_short
    
    return txt_file_name_short,id_num
    
    
    

#################################


img_folder_name=sys.argv[1]
thld=int(sys.argv[2])
min_pix=int(sys.argv[3])


img_folder_name=os.path.abspath(img_folder_name)

pwd=os.getcwd()

if os.path.isdir(img_folder_name+'/peak_lists')==True:
    files=os.listdir(img_folder_name+'/peak_lists')
    for file in files:
        os.remove(img_folder_name+'/peak_lists/'+file)
    os.rmdir(img_folder_name+'/peak_lists')
    
peak_lists_dir=img_folder_name+'/peak_lists'
os.mkdir(peak_lists_dir)
print('peak_lists directory is made: ',peak_lists_dir)

print(img_folder_name)

img_file_list=os.listdir(img_folder_name)
img_file_list.sort()

for file in img_file_list:
    if file.endswith('.tif'):
        
        out_file_name_short,id_num=file_name_gen(file)
        out_file_name=os.path.join(peak_lists_dir,out_file_name_short)
        #print(out_file_name)
        file_full=os.path.join(img_folder_name,file)
        label,peak_coord,props=NSLS.peak_finder(file_full,thld,min_pix) #thld=30,min_pix=15
        
        
        NSLS.peak_list_out(label,peak_coord,props,out_file_name)
        print('peakfound:',out_file_name_short)
        print('Approximate completion: %d out of %d'%(id_num,len(img_file_list)))

print('ALL DONE!')

