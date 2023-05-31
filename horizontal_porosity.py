# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 09:30:37 2023

@author: pkiuru

Script file for
(1) calculation of vertical and horizontal air-filled porosity distributions of
    a binary image (for the network domain)
(2) determination of surface topography, i.e., locations of top and bottom solid
    space voxels
    

Required input:
    Output of the script file net_poros.py

"""

import numpy as np
import matplotlib.pyplot as plt
import porespy.tools as pst
import porespy.metrics as psm

#bin_img_nofloat = np.load('C:\\Users\\pkiuru\\Pictures\\metnet_ct_reco\\20_25_5.npy')
#bin_img_nofloat = np.load('D:/binary_5A.npy')
#img = bin_img_nofloat[110:940,:,:]
#img = bin_img_percolating[110:940,:,:]
#img = np.copy(reg_image[3:-3,:,:])
#img = img > 0
#img_boundary = borders[110:940,:,:]


def horizont_por_prof(img, img_boundary):

    #Otetaan mallikappale sylinterin mallista.
    template = np.copy(~img_boundary[0,:,:])
    
    #Muodostetaan horisontaalihila
    shape = np.shape(img)
    step = 1
    
    grid = []
    section_area = []
    section_porosity = []
    
    for m in range(0,shape[1],step):
        #print(m)
        for n in range(0,shape[2],step):
            if shape[1]-m<step:
                incr1 = shape[1] - m
            else:
                incr1 = step
            if shape[2]-n<step:
                incr2 = shape[2] - n
            else:
                incr2 = step
            grid.append([[m, m + incr1],[n, n + incr2]])
            
    for section in grid:
        dummy = template[section[0][0]:section[0][1],section[1][0]:section[1][1]]
        section_area.append(np.sum(dummy)/np.size(dummy))
        dummy = img[:,section[0][0]:section[0][1],section[1][0]:section[1][1]]
        section_porosity.append(np.sum(dummy)/np.size(dummy))
    
    side_length = int(np.sqrt(len(section_porosity)))
       
    section_area = np.asarray(section_area)
    section_area = np.reshape(section_area,(side_length,side_length))
    
    section_porosity = np.asarray(section_porosity)
    section_porosity = np.reshape(section_porosity,(side_length,side_length))/section_area
    
    return section_porosity


def vertical_boundary_voxels(img):

    #Muodostetaan horisontaalihila
    shape = np.shape(img)
    step = 1
    
    grid = []
    top_voxels = []
    bottom_voxels = []
    
    for m in range(0,shape[1],step):
        #print(m)
        for n in range(0,shape[2],step):
            if shape[1]-m<step:
                incr1 = shape[1] - m
            else:
                incr1 = step
            if shape[2]-n<step:
                incr2 = shape[2] - n
            else:
                incr2 = step
            grid.append([[m, m + incr1],[n, n + incr2]])
            
    for section in grid:

        dummy = img[:,section[0][0]:section[0][1],section[1][0]:section[1][1]]
        solids = np.where(~dummy)[0]
        try:
            top_voxels.append(solids[0])
            bottom_voxels.append(solids[-1])
        except IndexError:
            top_voxels.append(shape[0])
            bottom_voxels.append(0)
    
    side_length = int(np.sqrt(len(bottom_voxels)))
       
    top_voxels = np.asarray(top_voxels)
    top_voxels = np.reshape(top_voxels,(side_length,side_length))
    
    bottom_voxels = np.asarray(bottom_voxels)
    bottom_voxels = np.reshape(bottom_voxels,(side_length,side_length))
    
    return top_voxels, bottom_voxels

#z_top = 50
#z_bottom = 950

borders = None

# Vertical range of the network domain
z_top, z_bottom = domain_params[scan-1,0:2]

# Create binary mask for the region outside the cylindrical sample
if borders is None:
    cylinder = pst.extract_cylinder(bin_img_nofloat.astype(int)+10, axis=0)
    borders = cylinder == 0

# Calculate the horizontal porosity distribution
section_porosity_bin = horizont_por_prof(bin_img_nofloat[z_top:z_bottom,:,:], borders[z_top:z_bottom,:,:])
#section_porosity_reg = horizont_por_prof(reg_image[3:-3,:,:] > 0, borders[z_top:z_bottom,:,:])
#section_porosity_perc = horizont_por_prof(bin_img_percolating[z_top:z_bottom,:,:], borders[z_top:z_bottom,:,:])


#When diameter = 900 voxels, this is the circular area factor ("discretized pi/4").
#pi_by_4 = 0.7853098765432098

# Calculate the vertical porosity profile for the network domain and for the
# whole image domain
#p_prof = psm.porosity_profile(bin_img_nofloat[z_top:z_bottom,:,:], 0) / (np.pi / 4)
p_prof = psm.porosity_profile(bin_img_nofloat[z_top:z_bottom,:,:], 0) / pi_by_4
p_prof_all = psm.porosity_profile(bin_img_nofloat, 0) / pi_by_4


# Vertical air-filled porosity profile, network domain
fig = plt.figure(num=901)
plt.clf()
fig.set_size_inches([4,4])
plt.plot(0.01*p_prof,'b',lw=0.5,label='Vertical afp profile, network domain')
plt.legend()
plt.xlabel('Distance from top (voxels)')
plt.ylabel('Air-filled porosity')
plt.tight_layout()

# Vertical air-filled porosity profile, binary image domain
fig = plt.figure(num=902)
plt.clf()
fig.set_size_inches([4,4])
plt.plot(0.01*p_prof_all,label='Vertical afp profile, whole image')
plt.legend()
plt.xlabel('Distance from top (voxels)')
plt.ylabel('Air-filled porosity')
plt.tight_layout()

#Determine the vertical locations of the top and bottom voxels of the solid matrix
tops, bottoms = vertical_boundary_voxels(bin_img_nofloat)

#Topographic map of the top surface
fig = plt.figure(921)
plt.clf()
fig.set_size_inches([4,4])
plt.imshow(tops, cmap=plt.cm.gist_earth.reversed(), interpolation='nearest',
           vmin=np.min(tops[~borders[0,:,:]]),vmax=np.percentile(tops[~borders[0,:,:]],99.75))
plt.title('Topmost solid voxel')
plt.colorbar()
plt.axis('off')
plt.tight_layout()

#Topographic map of the bottom surface
fig = plt.figure(922)
plt.clf()
fig.set_size_inches([4,4])
plt.imshow(bottoms, cmap=plt.cm.gist_earth, interpolation='nearest',
           vmin=np.percentile(bottoms[~borders[0,:,:]],0.25),vmax=np.max(bottoms[~borders[0,:,:]]))
plt.title('Bottommost solid voxel')
plt.colorbar()
plt.axis('off')
plt.tight_layout()

#Binary mask for the cylindrical image region
inner = np.ravel(~borders[0,:,:])

#Histogram of top and bottom surface topography with 0.1 mm (one-voxel) bins
fig = plt.figure(923)
plt.clf()
fig.set_size_inches([6,4])
ax1=plt.subplot(2,1,1)
dummy = np.ravel(tops)
counts, bins = np.histogram(dummy[inner],bins=np.max(dummy[inner]))
ax1.plot(0.1*bins[1:], counts)
plt.yscale('log')
plt.xlabel('Location of topmost voxel (mm)')
plt.title('Histogram of top surface topography')
ax2=plt.subplot(2,1,2)
dummy = np.ravel(bottoms)
counts, bins = np.histogram(dummy[inner],bins=np.max(dummy[inner])-np.min(dummy[inner]))
ax2.plot(0.1*bins[1:], counts)
plt.yscale('log')
plt.xlabel('Location of bottommost voxel (mm)')
plt.title('Histogram of bottom surface topography')
plt.tight_layout()

#Histogram of surface topography with 1 mm bins
bins = np.linspace(0,100,101)
fig = plt.figure(924)
plt.clf()
fig.set_size_inches([4,4])
ax1=plt.subplot(1,1,1)
dummy = np.ravel(tops)
dumm=ax1.hist(0.1*dummy[inner],bins=bins,label='Top',density=True)
#plt.yscale('log')
plt.xlabel('Location of the voxel closest to the side (mm)')
plt.title('Histogram of surface topography')
dummy = np.ravel(bottoms)
dumm=ax1.hist(0.1*dummy[inner],bins=bins,label='Bottom',density=True)
plt.legend()
plt.tight_layout()

dummy = np.ravel(tops)
mean_top = np.mean(dummy[inner])

# Horizontal air-filled porosity distribution for the network domain; max. value
# of the colorbar is the 99th percentile of the porosity values
fig = plt.figure(911)
plt.clf()
fig.set_size_inches([4,4])
plt.imshow(section_porosity_bin, cmap=plt.cm.plasma, interpolation='nearest',
           vmax=np.percentile(section_porosity_bin[~np.isnan(section_porosity_bin)],99))
plt.title('Horizontal air-filled porosity distribution \n (total macropore space)')
plt.colorbar()
plt.axis('off')
plt.tight_layout()

# Horizontal air-filled porosity distribution for the network domain; max. value
# of the colorbar is the 95th percentile of the porosity values
fig = plt.figure(912)
plt.clf()
fig.set_size_inches([4,4])
plt.imshow(section_porosity_bin, cmap=plt.cm.plasma, interpolation='nearest',
           vmax=np.percentile(section_porosity_bin[~np.isnan(section_porosity_bin)],95))
plt.title('Horizontal air-filled porosity distribution \n (total macropore space)')
plt.colorbar()
plt.axis('off')
plt.tight_layout()

'''
section_porosity_top = horizont_por_prof(bin_img_nofloat[z_top:500,:,:], borders[z_top:500,:,:])
section_porosity_btm = horizont_por_prof(bin_img_nofloat[500:z_bottom,:,:], borders[500:z_bottom,:,:])

fig = plt.figure(913)
plt.clf()
fig.set_size_inches([4,4])
plt.imshow(section_porosity_top, cmap=plt.cm.plasma, interpolation='nearest',
           vmax=np.percentile(section_porosity_bin[~np.isnan(section_porosity_bin)],95))
plt.title('Horizontal air-filled porosity distribution \n (top)')
plt.colorbar()
plt.axis('off')
plt.tight_layout()

fig = plt.figure(914)
plt.clf()
fig.set_size_inches([4,4])
plt.imshow(section_porosity_btm, cmap=plt.cm.plasma, interpolation='nearest',
           vmax=np.percentile(section_porosity_bin[~np.isnan(section_porosity_bin)],95))
plt.title('Horizontal air-filled porosity distribution \n (bottom)')
plt.colorbar()
plt.axis('off')
plt.tight_layout()
'''
'''
fig = plt.figure(912)
plt.clf()
fig.set_size_inches([5,5])
plt.imshow(section_porosity_reg, cmap=plt.cm.Reds, interpolation='nearest')
plt.title('Horizontal porosity distribution, regions')
plt.colorbar()
plt.axis('off')
plt.tight_layout()

fig = plt.figure(913)
plt.clf()
fig.set_size_inches([5,5])
plt.imshow(section_porosity_perc, cmap=plt.cm.Reds, interpolation='nearest')
plt.title('Horizontal porosity distribution, percolated')
plt.colorbar()
plt.axis('off')
plt.tight_layout()


fig = plt.figure(910)
plt.clf()
fig.set_size_inches([5,5])
plt.imshow((section_porosity_bin-section_porosity_reg)/section_porosity_reg, cmap=plt.cm.Reds, interpolation='nearest')
plt.title('Horizontal porosity distribution')
plt.colorbar()
plt.axis('off')
plt.tight_layout()

'''
'''
if False:
    
    radius = np.sqrt((pn['pore.coords'][pn['pore.internal'],1]-np.shape(bin_img_nofloat)[1]*voxel_size/2)**2 +
                     (pn['pore.coords'][pn['pore.internal'],2]-np.shape(bin_img_nofloat)[1]*voxel_size/2)**2)
    vsum = np.sum(pn['pore.volume'][pn['pore.internal']])
    nbins = 22
    z_bin = np.linspace(0.0, np.max(radius),nbins)
    pore_volume_hor_distr = np.zeros((nbins-1,3))
    for i in range(len(z_bin)-1):
        A=np.logical_and(radius>=z_bin[i], radius<z_bin[i+1])
        pore_volume_hor_distr[i,0] = z_bin[i]
        pore_volume_hor_distr[i,1] = z_bin[i+1]
        pore_volume_hor_distr[i,2] = np.sum(pn['pore.volume'][pn['pore.internal']][A]) / vsum
    
    
    relative_gradient = np.diff(np.concatenate([np.array([100]),p_prof_all]))/np.median(p_prof_all)

    fig = plt.figure(num=903)
    plt.clf()
    fig.set_size_inches([5,5])
    plt.plot(relative_gradient,label='Vertical air-filled porosity profile gradient, whole image')
    plt.legend()
    plt.xlabel('Distance from top (voxels)')
    plt.ylabel('Gradient')
    plt.tight_layout()
    
    
    inner_space = np.where(np.logical_and(np.abs(relative_gradient) < 0.1, p_prof_all < 3*np.median(p_prof_all)))[0]
'''    
#porosity_bin_img = np.sum(img)/np.size(img)/caf_real
#porosity_sections_mean = np.mean(section_porosity[~np.isnan(section_porosity)])


