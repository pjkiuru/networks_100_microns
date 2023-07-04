# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 18:02:42 2023

Script for generating figures of tomography images and porosity profiles

@author: pkiuru
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in' 

fontname = 'Arial'
fontsize = 8

rcParams['font.size'] = fontsize

rcParams['axes.titlesize'] = fontsize
rcParams['axes.labelsize'] = fontsize
rcParams['axes.titleweight'] = 'normal'

rcParams['xtick.labelsize'] = fontsize
rcParams['ytick.labelsize'] = fontsize

rcParams['legend.handlelength'] = 1.0

subtpos = [0.01,1.05]
subtpos_porosity_t = [0.01,1.03]
subtpos_porosity_b = [0.01,1.06]
subt = ['(a)', '(b)', '(c)', '(d)']

colors = list([(0,0,0), (230/255,159/255,0), (86/255,180/255,233/255),
           (0,158/255,115/255), (240/255,228/255,66/255),
           (0,114/255,178/255), (213/255,94/255,0),(204/255,121/255,167/255)]) 
ls = [':','-.','-','--']
lw = 0.75


draw_porosity_profile_image =False

draw_cross_section_image = True


#Vertical and horizontal porosity profiles
if draw_porosity_profile_image:
    
    df=pd.read_excel('C:\\Users\\pkiuru\\OneDrive - University of Eastern Finland\\scan_analysis_large.xlsx', sheet_name='Distros',
                         header=0,skiprows=[0])
    
    #Vertical porosity profiles
    vpp_1A = df['1A'][2:]
    vpp_3A = df['3A'][2:]
    vpp_5A = df['5A'][2:]
    vpp_6A = df['6A'][2:]
    vpp_7A = df['7A'][2:]
    
    #Horizontal porosity profiles
    hpp_1A = df.iloc[2:,1]
    hpp_3A = df.iloc[2:,9]
    hpp_5A = df.iloc[2:,17]
    hpp_6A = df.iloc[2:,21]
    hpp_7A = df.iloc[2:,25]
    
    
    
    fig = plt.figure(num=100)
    fig.set_size_inches(3.0,4.8)
    plt.clf()
    
    axt = fig.add_subplot(2,1,1)
    axt.plot(0.01*vpp_1A,0.1*np.linspace(0,len(vpp_1A)-1,len(vpp_1A)),ls[0],lw=lw,c=colors[0],label='1A')
    axt.plot(0.01*vpp_3A,0.1*np.linspace(0,len(vpp_3A)-1,len(vpp_3A)),ls[1],lw=lw,c=colors[1],label='3A')
    axt.plot(0.01*vpp_5A,0.1*np.linspace(0,len(vpp_5A)-1,len(vpp_5A)),ls[2],lw=lw,c=colors[2],label='5A')
    axt.plot(0.01*vpp_7A,0.1*np.linspace(0,len(vpp_7A)-1,len(vpp_7A)),ls[3],lw=lw,c=colors[3],label='7A')
    axt.legend(bbox_to_anchor=(0.775, 0.355), bbox_transform=axt.transAxes, handletextpad=0.55)
    axt.text(subtpos_porosity_t[0], subtpos_porosity_t[1], subt[0], transform=axt.transAxes) 
    axt.text(0.805, 0.660, 'Sample', transform=axt.transAxes)  
    axt.set_xlim([0.00, 0.260])
    axt.set_ylim([0.0, 86])
    axt.invert_yaxis()
    
    axt.set_ylabel('Distance from top (mm)')
    axt.set_xlabel('Air-filled porosity')
    
    
    #l = np.sum(~hpp_5A.isna())
    pi_by_four_900 = 0.7853098765432098 # discretized pi/4 for a 900-diameter circle    
    volumes = np.linspace(0,4*pi_by_four_900*450**2,46)
    radii = 0.1*np.sqrt(volumes/pi_by_four_900)/2
    radii = (radii[:-1]+radii[1:])/2
    
    
    axb = fig.add_subplot(2,1,2)
    axb.plot(radii,hpp_1A[~hpp_1A.isna()],ls[0],lw=lw,c=colors[0],label='1A')
    axb.plot(radii,hpp_3A[~hpp_3A.isna()],ls[1],lw=lw,c=colors[1],label='3A')
    axb.plot(radii,hpp_5A[~hpp_5A.isna()],ls[2],lw=lw,c=colors[2],label='5A')
    axb.plot(radii,hpp_7A[~hpp_7A.isna()],ls[3],lw=lw,c=colors[3],label='7A')
    axb.text(subtpos_porosity_b[0], subtpos_porosity_b[1], subt[1], transform=axb.transAxes)
    axb.set_xlim([0, 47])
    axb.set_ylim([-0.02, 0.52])
    
    axb.set_xlabel('Distance from center (mm)')
    axb.set_ylabel('Air-filled porosity')
    
    axb.yaxis.set_label_coords(-0.09, 0.5)
    axb.xaxis.set_label_coords(0.5, -0.14)
    
    axt.xaxis.set_label_coords(0.5, -0.07)
    
    #plt.legend()
    
    axt.set_position([0.125, 0.425, 0.83, 0.53])
    axb.set_position([0.125, 0.065, 0.83, 0.255])  
    
    #plt.savefig('sample_porosity_profiles.pdf')
    

# micro-CT image cross-sections
if draw_cross_section_image: 
    
    # Line width, linestyle (loosely dashed), color
    lw_vlim = 0.75
    ls_vlim = [(0, (5, 10))]
    c_vlim = 'r'
    
    #Intensity scaling limits in tjhe images
    vmin, vmax = 40, 255
    
    im_1A_cs1 = np.load('D:/filtered_1A_cs1.npy')
    im_1A_cs2 = np.load('D:/filtered_1A_cs2.npy')
    
    im_1A_axrot30_cs1 = np.load('D:/filtered_1A_axrot30_cs1.npy')
    im_1A_axrot30_cs2 = np.load('D:/filtered_1A_axrot30_cs2.npy')
    
    im_3A_cs1 = np.load('D:/filtered_3A_cs1.npy')
    im_3A_cs2 = np.load('D:/filtered_3A_cs2.npy')
    
    im_5A_cs1 = np.load('D:/filtered_5A_cs1.npy')
    im_5A_cs2 = np.load('D:/filtered_5A_cs2.npy')
    
    im_6A_cs1 = np.load('D:/filtered_6A_cs1.npy')
    im_6A_cs2 = np.load('D:/filtered_6A_cs2.npy')
    
    im_7A_cs1 = np.load('D:/filtered_7A_cs1.npy')
    im_7A_cs2 = np.load('D:/filtered_7A_cs2.npy')
    
    #Network domain limits
    vlim_1A = [110,970]
    vlim_3A = [140,940]
    vlim_5A = [100,950]
    vlim_6A = [130,970]
    vlim_7A = [160,970]
    
    
    
    fig = plt.figure(num=1)
    fig.set_size_inches(6.5,2.0) #7 x 2.1
    plt.clf()
    
    
    
    ax1 = fig.add_subplot(1,4,1)
    ax1.imshow(im_1A_axrot30_cs2, vmin=vmin, vmax=vmax, cmap='gray', interpolation='nearest')
    #plt.axis('off')
    plt.axhline(y=vlim_1A[0], color=c_vlim, linestyle=ls_vlim[0],lw=lw_vlim)
    plt.axhline(y=vlim_1A[1], color=c_vlim, linestyle=ls_vlim[0],lw=lw_vlim)
    
    #ax1.set_title('a)', loc='left', fontsize = fontsize)
    ax1.text(subtpos[0], subtpos[1], subt[0], transform=ax1.transAxes)
    #ax1.patch.set_edgecolor('black')  
    #ax1.patch.set_linewidth('0.5') 
    ax1.xaxis.set_ticks_position('none')
    ax1.yaxis.set_ticks([]) 
    ax1.xaxis.set_ticklabels([])
    ax1.yaxis.set_ticklabels([])
    
    ax1.spines['left'].set_linewidth(0.5)
    ax1.spines['right'].set_linewidth(0.5)
    ax1.spines['top'].set_linewidth(0.5)
    ax1.spines['bottom'].set_linewidth(0.5)
    
    
    ax2 = fig.add_subplot(1,4,2)
    ax2.imshow(im_3A_cs2, vmin=vmin, vmax=vmax, cmap='gray', interpolation='nearest')
    #plt.axis('off')
    plt.axhline(y=vlim_3A[0], color=c_vlim, linestyle=ls_vlim[0],lw=lw_vlim)
    plt.axhline(y=vlim_3A[1], color=c_vlim, linestyle=ls_vlim[0],lw=lw_vlim)
    
    #ax2.set_title('b)', loc='left', fontsize = fontsize)
    ax2.text(subtpos[0], subtpos[1], subt[1], transform=ax2.transAxes)
    #ax2.patch.set_edgecolor('black')  
    #ax2.patch.set_linewidth('0.5') 
    ax2.xaxis.set_ticks_position('none')
    ax2.yaxis.set_ticks([]) 
    ax2.xaxis.set_ticklabels([])
    ax2.yaxis.set_ticklabels([])
    
    ax2.spines['left'].set_linewidth(0.5)
    ax2.spines['right'].set_linewidth(0.5)
    ax2.spines['top'].set_linewidth(0.5)
    ax2.spines['bottom'].set_linewidth(0.5)
    
    
    ax3 = fig.add_subplot(1,4,3)
    ax3.imshow(im_5A_cs1, vmin=vmin, vmax=vmax, cmap='gray', interpolation='nearest')
    #plt.axis('off')
    plt.axhline(y=vlim_5A[0], color=c_vlim, linestyle=ls_vlim[0],lw=lw_vlim)
    plt.axhline(y=vlim_5A[1], color=c_vlim, linestyle=ls_vlim[0],lw=lw_vlim)
    
    #ax3.set_title('c)', loc='left', fontsize = fontsize)
    ax3.text(subtpos[0], subtpos[1], subt[2], transform=ax3.transAxes)
    #ax3.patch.set_edgecolor('black')  
    #ax3.patch.set_linewidth('0.5') 
    ax3.xaxis.set_ticks_position('none')
    ax3.yaxis.set_ticks([]) 
    ax3.xaxis.set_ticklabels([])
    ax3.yaxis.set_ticklabels([])
    
    ax3.spines['left'].set_linewidth(0.5)
    ax3.spines['right'].set_linewidth(0.5)
    ax3.spines['top'].set_linewidth(0.5)
    ax3.spines['bottom'].set_linewidth(0.5)
    
    
    ax4 = fig.add_subplot(1,4,4)
    ax4.imshow(im_7A_cs1, vmin=vmin, vmax=vmax, cmap='gray', interpolation='nearest')
    #plt.axis('off')
    plt.axhline(y=vlim_7A[0], color=c_vlim, linestyle=ls_vlim[0],lw=lw_vlim)
    plt.axhline(y=vlim_7A[1], color=c_vlim, linestyle=ls_vlim[0],lw=lw_vlim)
    
    #ax4.set_title('d)', loc='left', fontsize = fontsize)
    ax4.text(subtpos[0], subtpos[1], subt[3], transform=ax4.transAxes)
    #ax4.patch.set_edgecolor('black')  
    #ax4.patch.set_linewidth('0.5') 
    ax4.xaxis.set_ticks_position('none')
    ax4.yaxis.set_ticks([]) 
    ax4.xaxis.set_ticklabels([])
    ax4.yaxis.set_ticklabels([])
    
    ax4.spines['left'].set_linewidth(0.5)
    ax4.spines['right'].set_linewidth(0.5)
    ax4.spines['top'].set_linewidth(0.5)
    ax4.spines['bottom'].set_linewidth(0.5)
    
    
    ax1.set_position([0.010, 0.06, 0.24, 0.81])
    ax2.set_position([0.255, 0.06, 0.24, 0.81])
    ax3.set_position([0.500, 0.06, 0.24, 0.81])
    ax4.set_position([0.745, 0.06, 0.24, 0.81])

    #plt.savefig('uCT_cross_sections.pdf',dpi=300)

