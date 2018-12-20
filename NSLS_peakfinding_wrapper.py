#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

The wrapper for group jobs of peak finding


Usage:
    NSLS_peakfinding_wrapper.py -i <folder_name> [options]
    
    
Options:
    -i <folder_name>
    -h --help                Show this message
    --list_fmt=<list_fmt>    Output peak list format,'SPIND_py'or 'CFL_NSLS' [default: SPIND_py]
    --photon-energy=<E_ph>   Photon energy of the X-ray [default: 12000]
    --thld=<threshold>       Threshold of the peak finding [default: 30]
    --min-pix=<num_pix>      Minimal number of connected pixels [default: 15]
"""

import sys,os
import numpy as np
import matplotlib.pyplot as plt
import glob
from docopt import docopt


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
if __name__=='__main__':

    
    argv= docopt(__doc__)
    
    print(argv)
    
    
    img_folder_name=argv['-i']
    thld=int(argv['--thld'])
    min_pix=int(argv['--min-pix'])
    E_ph=float(argv['--photon-energy'])
    list_fmt=argv['--list_fmt']


    img_folder_name=os.path.abspath(img_folder_name)

    pwd=os.getcwd()

    if os.path.isdir(img_folder_name+'/peak_list')==True:
        files=os.listdir(img_folder_name+'/peak_list')
        for file in files:
            os.remove(img_folder_name+'/peak_list/'+file)
        os.rmdir(img_folder_name+'/peak_list')
    
    peak_lists_dir=img_folder_name+'/peak_list'
    os.mkdir(peak_lists_dir)
    print('peak_lists directory is made: ',peak_lists_dir)

    print(img_folder_name)

    img_file_list=os.listdir(img_folder_name)

    img_file_list=glob.glob(img_folder_name+'/*.tif')
    img_file_list.sort()

    for file in img_file_list:
        if file.endswith('.tif'):
        
            out_file_name_short,id_num=file_name_gen(file)
            out_file_name_short='event-%d.txt'%(id_num)#This line adapts the file name to the SPIND format. Comment this line
            #if the output peaklist file name needs to match the .tif file name.
            out_file_name=os.path.join(peak_lists_dir,out_file_name_short)
            #print(out_file_name)
            file_full=os.path.join(img_folder_name,file)
            label,peak_coord,props=NSLS.peak_finder(file_full,thld,min_pix) #thld=30,min_pix=15
        
        
            NSLS.peak_list_out(label,peak_coord,props,out_file_name,list_fmt,E_ph)
            print('peakfound:',out_file_name_short)
            print('Approximate completion: %d out of %d'%(id_num,len(img_file_list)))

    print('ALL DONE!')
    if list_fmt=='CFL_NSLS':
        print('The format of the peak list is\n \
              row,col,Intensity,connectivity,pos_x,pos_y,pos_z')
    elif list_fmt=='SPIND_py':
        print('The format of the peak list is\n \
              row,col,Intensity,connectivity,q_x,q_y,q_z,res')
        

