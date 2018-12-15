#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 10:26:25 2018
'NSLS_IND_util' contains the utilities and tools for the indexing prep


@author: chufeng
"""

import sys,os
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure, morphology, feature



def peak_finder(img_file_name,thld,min_pix):
    
    im=plt.imread(img_file_name)
    im=im[:,:,0] # To convert the image to gray scale
    all_labels=measure.label(im>thld)
#    num_blobs=all_labels.max()+1
#    
#    
#    blob_tab=np.zeros((num_blobs,2),dtype=np.int16)#the first column is the blob id, the second
#    blob_tab[:,0]=np.arange(0,num_blobs)   
#
#                           #column is the number of pixels in this blob
#    for blob_id in np.arange(0,num_blobs):
#        
#        ind=np.where(all_labels==blob_id)
#        
#        num_tab[blob_id,1]=np.shape(ind[0])[0]
#        
#        blob_mask=(all_labels==blob_id)
#        blob_seg=blob_mask*im
#        
#        props=measure.regionprops(all_labels,im)
        
        
    props=measure.regionprops(all_labels,im)
    
    
    
    area=np.array([r.area for r in props]).reshape(-1,)
    max_intensity=np.array([r.max_intensity for r in props]).reshape(-1,)
    #coords=np.array([r.coords for r in props]).reshape(-1,)
    label=np.array([r.label for r in props]).reshape(-1,)
    centroid=np.array([np.array(r.centroid).reshape(1,2) for r in props]).reshape((-1,2))
    weighted_centroid=np.array([r.weighted_centroid for r in props]).reshape(-1,)
    
    label_filtered=label[area>min_pix]
    
    area_filtered=area[area>min_pix]
    
    area_sort_ind=np.argsort(area_filtered)[::-1]
    
    label_filtered_sorted=label_filtered[area_sort_ind] 
        
    area_filtered_sorted=area_filtered[area_sort_ind]
    
    weighted_centroid_filtered=np.zeros((len(label_filtered_sorted),2))
    
    for index,value in enumerate(label_filtered_sorted):
        
        weighted_centroid_filtered[index,:]=np.array(props[value-1].weighted_centroid)
    
    print('In image: %s \n %5d peaks are found' %(img_file_name, len(label_filtered_sorted)))
    
    beam_center=np.array([1492.98,2163.41])
    
    plt.figure(figsize=(15,15))
#    
    plt.imshow(im,cmap='gray')
    plt.clim(0,1.5*thld)
    plt.scatter(weighted_centroid_filtered[:,1],weighted_centroid_filtered[:,0],edgecolors='r',facecolors='none')
    plt.scatter(beam_center[1],beam_center[0],marker='*',color='b')
    plt.show()
    
    return label_filtered_sorted,weighted_centroid_filtered,props




def pos2angle(col,row):
    pix = 74.8
    R = 2.315e5
    th1 = 0.7617
    phi1 = 3.0366
    th2 = 0.1796
    phi2 = 2.5335
    phi3 = -0.1246
    alpha = 8.5*np.pi/180

    row = 1536-row
    
    det_orig = R*np.array([np.sin(th1)*np.cos(phi1),np.sin(th1)*np.sin(phi1),np.cos(th1)])
    det_z = np.array([np.sin(th2)*np.cos(phi2), np.sin(th2)*np.sin(phi2),np.cos(th2)])
    th3 = np.arctan(-1.0/(np.cos(phi2-phi3)*np.tan(th2)))
    det_x = np.array([np.sin(th3)*np.cos(phi3),np.sin(th3)*np.sin(phi3),np.cos(th3)])
    det_y = np.cross(det_z,det_x)

    #print(np.inner(det_x,det_x))
    #print(np.inner(det_y,det_y))
    #print(np.inner(det_z,det_z))

#    print('det_orig: ',det_orig)
#    print('det_orig length: ',np.sqrt(np.inner(det_orig,det_orig)))
#    print('det_x: ',det_x)
#    
#    print('det_y: ',det_y)
#    print('det_z: ',det_z)

    pos = det_orig + (col - 1)*pix*det_x + (row -1)*pix*det_y

    A=np.array([det_x[0:2],det_y[0:2]]).T
#    print('A: \n',A)
    b=-det_orig[0:2]

    sol=np.linalg.solve(A,b)
    beam_center_col=1+sol[0]/pix
    beam_center_row=1536-(1+sol[1]/pix)
    beam_center=np.array([beam_center_col,beam_center_row])
#    print('beam_center[col,row]:' , beam_center)
#
#    print('pos: ',pos)
#    print('distance from frame origin: ',np.sqrt(np.inner(pos,pos)))

    M = np.array([[np.cos(alpha),-np.sin(alpha),0],[np.sin(alpha), np.cos(alpha),0],[0,0,1]])

    #pos = np.dot(M,pos)

#    print('after alpha rotation: ',pos)

    tth = np.arccos(pos[2]/np.sqrt(pos[0]**2+pos[1]**2+pos[2]**2))*180.0/np.pi
    delta = np.arcsin(pos[1]/np.sqrt(pos[0]**2+pos[1]**2+pos[2]**2))*180.0/np.pi
    pos_xy = pos*np.array([1,0,1])
    gamma = np.arccos(pos[2]/np.sqrt(pos_xy[0]**2+pos_xy[1]**2+pos_xy[2]**2))*180.0/np.pi
#    print(gamma,delta,tth)
    #return (gamma,delta,tth)
    return (pos)
#if __name__=='__main__':
#
#    col=np.int16(sys.argv[1])
#    row=np.int16(sys.argv[2])
#
#    pos2angle(col,row)

def peak_list_out(label,peak_coord,props,out_file_name):
    
    t_list=np.zeros((label.shape[0],7))
    
    for index,value in enumerate(label):
        t_list[index,0:2]=np.array(props[value-1].weighted_centroid)
        t_list[index,2]=np.array(props[value-1].max_intensity)
        t_list[index,3]=np.array(props[value-1].area)
        pos=pos2angle(props[value-1].weighted_centroid[1],props[value-1].weighted_centroid[0])
        t_list[index,4:]=pos

    np.savetxt(out_file_name,t_list,delimiter=' ',fmt='%.5e')
    print(out_file_name,'saved!')

    return None
    

