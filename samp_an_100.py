# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:17:27 2020

@author: pkiuru

Script file for preprocessing the 3D microtomography images wit 100 micron voxel
size and extracting pore network

Required input:
    Folder containing the microtomography image files (topfolder)
    File containing the image crop and rotation parameters, image grayscale
    intensity ranges and network domain ranges (paramfile)

Output:
    3D binary image (bin_img_nofloat) with void space = 1 & solid space = 0
    OpenPNM network object with all pores (pn_all) and the largest connected
    network (pn)
    
    These can be saved with the commented commands at the end of this file

"""

import os
import time

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from skimage import filters
from skimage import exposure
from skimage import io, util
from skimage.morphology import ball

import numpy as np
import my_models as mm
from scipy import ndimage#, stats
import porespy as ps
import openpnm as op

#Folder of the image stack folders
topfolder = 'D:/peat_030223'

#File that contains the image rotation parameters and grayscale intensity ranges etc.
paramfile = 'C:\\Users\\pkiuru\\OneDrive - University of Eastern Finland\\METNET_scripts\\Scan_parameters_2.csv'

# Name of the tomography image folder in the folder list
sample_name = '7B'

#sample_name = None
scan = 10

# Switch for reading the image files and stacking a 3D image.
# The first array axis contains the axial coordinate of the sample cylinder.
# Creates also images of central cross sections in all directions (Figures 110-112)
# and two horizontal cross sections vith a vertical distance of 800 voxels (Figure 113).
imageread = True

# Switch for performing straightening and cylinder extraction of the 3D image.
# Creates also voxel grayscale intensity histogram (Figure 22) and horizontal
# cross sections near the top and the bottom of the original and
# straightened images (Figure 532).
cylinder_analysis = True

# Performs conversion from 16-bit to 8-bit grayscale image
bitscale = True

# Performs noise filtering.
# Creates also voxel grayscale intensity histograms for the original and filtered
# 8-bit images (Figure 29).
filtering = True

# Performs solid-void thresholding and creates the binary image (with void
# space having value 1)
# Creates also horizontal cross-section images of the original, filtered, and
# binary images (Figure 501)
binary = True

# Calculates image porosity and vertical porosity profile
# Creates also a plot of vertical porosity profile (Figure 94)
poros = False

# Switch for network extraction with PoreSpy and pore network object generation
# with OpenPNM
netgen = True

# Simple test code for percolation and effective diffusion simulation in
# the pore network - NOT validated
percolation = False

#Image voxel size (m)
voxel_size = 100e-6

if imageread:
    # Generate  list of image stack folders
    image_folders = []
    filenames = []
    samples = []
    
    for root, dirs, files in os.walk(topfolder, topdown=False):
        for name in dirs:
            image_folders.append(os.path.join(root, name))
            samples.append(name.split("_"))

    #Read the image in 16 bit grayscale
    stack_list = []
    
    if sample_name is not None:
        
        for i, item in enumerate(samples):
            if len(item)==2:
                if item[0] == sample_name:
                    dummy = i
                    break
        
        print('Loading sample ' + sample_name + ' from folder ' + os.path.basename(image_folders[dummy]))
        for root, dirs, files in os.walk(image_folders[dummy], topdown=False):
            for filename in files:
                filenames.append(os.path.join(root, filename))
                if(filename[-3:]=='tif'):
                    stack_list.append(io.imread(os.path.join(root, filename)))
        
    else:
        print('Loading ' + os.path.basename(image_folders[scan-1]))
        for root, dirs, files in os.walk(image_folders[scan-1], topdown=False):
            for filename in files:
                filenames.append(os.path.join(root, filename))
                stack_list.append(io.imread(os.path.join(root, filename)))
    
    whole_image = np.asarray(stack_list)
    del stack_list
    
    borders = None  
    
    #Middle z
    mm.osakuva(whole_image[576,:,:], 110, cmap='gray')
    #Middle x
    mm.osakuva(whole_image[:,576,:], 111, cmap='gray') 
    #Middle y
    mm.osakuva(whole_image[:,:,576], 112, cmap='gray')
    
    fig, ax = plt.subplots(1, 2, num=113, figsize=(12, 7), sharex=True, sharey=True)
    ax[0].imshow(whole_image[190, :, :], cmap='gray', interpolation='nearest')
    ax[1].imshow(whole_image[990, :, :], cmap='gray', interpolation='nearest')
    
    
if cylinder_analysis:
    
    if sample_name is not None:
        parameters = (mm.read_params_name(paramfile))
        names = parameters.filenames
        crop_coords = parameters.crop_coords
        rot_params = parameters.rot_params
        scale_params = parameters.scale_params
        domain_params = parameters.domain_params
        
    try:
        scan = names.index(sample_name) + 1
    except ValueError:
        print("Sample does not exist")
        
    else:
        crop_coords, rot_params, scale_params = (mm.read_params(paramfile))
    
    
    #cropped = whole_image[crop_coords[scan-1,0]:crop_coords[scan-1,1],
    #                      crop_coords[scan-1,2]:crop_coords[scan-1,3], 
    #                      crop_coords[scan-1,4]:crop_coords[scan-1,5]]
    
    rotated = np.copy(whole_image)
    
    if rot_params[scan-1,0] != 0:
        print('Rotating around x axis')
        rotated = ndimage.rotate(rotated, rot_params[scan-1,0],
                             axes=(0, 1),
                             mode='reflect', reshape=False)
    
    if rot_params[scan-1,1] != 0:
        print('Rotating around y axis')
        rotated = ndimage.rotate(rotated, rot_params[scan-1,1],
                             axes=(0, 2),
                             mode='reflect', reshape=False)
    
    cropped_rot = np.copy(rotated[crop_coords[scan-1,0]:crop_coords[scan-1,1],
                          crop_coords[scan-1,2]:crop_coords[scan-1,3], 
                          crop_coords[scan-1,4]:crop_coords[scan-1,5]])
    
    # Finds the voxels with intensity = 0 inside the cylinder and changes these
    # values to 1 (value 0 is the mask for the cylindrical cropping algorithm)
    numzeros = np.sum(cropped_rot==0)
    if numzeros > 0:
        
        checker = np.copy(cropped_rot)
        print(numzeros, 'voxels with zero value in the image; set their value to 1.')
        cropped_rot[cropped_rot==0] = 1
    
    print('Extracting cylinder')
    syl_rot = ps.tools.extract_cylinder(cropped_rot, axis=0)
    
    #Binary mask for the region outside the cylindrical peat matrix
    borders = syl_rot == 0
    
    #Fraction of the area of a sphere (with diameter = d) to the area
    #of a square (with side length = d) in a binary image
    caf_real = np.sum(~borders[0,:,:])/np.size(borders[0,:,:])

    if numzeros > 0:
        numzeros_2 = np.sum(checker[borders]==0)
        print(numzeros_2, 'voxels with zero value are outside the core')
        del checker
    
    hist_s, bins_center_s = exposure.histogram(syl_rot[~borders])
    
    plt.figure(num=22, figsize=(5,5))
    plt.clf()
    ax = plt.subplot(111)   
    ax.plot(bins_center_s, hist_s, linewidth=2,label='Original image')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(visible=True, which='both')
    plt.legend()
    
    plt.figure(num=532, figsize=(8, 8))
    plt.clf()
    ax1 = plt.subplot(221)    
    plt.imshow(whole_image[150,:,:], cmap='gray', interpolation='nearest')#, vmin=0, vmax=25000)
    plt.title('Orig 100')
    ax2 = plt.subplot(222, sharex=ax1, sharey=ax1)    
    plt.imshow(whole_image[-250,:,:], cmap='gray', interpolation='nearest')#, vmin=0, vmax=25000)
    plt.title('Orig 1000')
    ax3 = plt.subplot(223, sharex=ax1, sharey=ax1)    
    plt.imshow(rotated[150,:,:], cmap='gray', interpolation='nearest')#, vmin=0, vmax=25000)
    plt.title('Rot 100')
    ax4 = plt.subplot(224, sharex=ax1, sharey=ax1)    
    plt.imshow(rotated[-250,:,:], cmap='gray', interpolation='nearest')#, vmin=0, vmax=25000)
    plt.title('Rot 1000')
    plt.tight_layout()
    
    #Lokaalien maksimien etsintäfunktio piirtää myös kuvan
    #bins_av_s, avg_dir, loclower, lochigher, lower, higher = mm.local_maxima(hist_s, bins_center_s, num=80, rang=100)

    # 0.5 and 99.5 percentiles (P005 and P995) used in the intendity range rescaling
    P005, P995 = np.percentile(syl_rot[~borders], (0.5,99.5))
'''    
if histo:
    
    kokokuva = whole_image[200:800,210:790,210:790]
    
    hist, bins_center = exposure.histogram(kokokuva)
    
    mm.osakuva(kokokuva[0,:,:], 1)
    mm.osakuva(kokokuva[-1,:,:], 11)
    
    plt.figure(num=2, figsize=(7,7))
    plt.clf()
    plt.plot(bins_center, hist, lw=2)
    
    # First quartile (Q1) 
    Q1 = np.percentile(kokokuva, 25, interpolation = 'midpoint') 
    # Third quartile (Q3) 
    Q3 = np.percentile(kokokuva, 75, interpolation = 'midpoint') 
      
    # Interquartile range (IQR) 
    IQR = Q3 - Q1 
    
    ol_range = 1.5
    
    I_min = Q1 - ol_range * IQR
    I_max = Q3 + ol_range * IQR
    
    #G16min = 9500#5000#I_min
    #G16max = 25000#22000#I_max
'''
    
if bitscale:

    G8min = 0
    G8max = 255
    
    h = np.copy(syl_rot)
    
    G16min = P005
    G16max = P995

    print('Scaling to 8 bit with limits', int(G16min), 'and', int(G16max))
    h[h>G16max] = G16max
    h[h<G16min] = G16min           
    image_8bit = (G8min + (G8max-G8min)/(G16max-G16min) * (h-G16min)).astype(np.uint8)
    del h
    print('Scaling done')
    
    print('Calculating histogram')
    try:
        hists, bins_centers = exposure.histogram(image_8bit[~borders])
    except:
        hists, bins_centers = exposure.histogram(image_8bit)
        
    print('Histogram done')

if filtering:    
    
    filter_radius = 2
    
    print('Filtering...')
    print('Filter radius = ' + str(filter_radius))

    #filtered = ndimage.median_filter(image_8bit, footprint=ball(filter_radius))  
    #print('Filtered')
    
    tic = time.perf_counter()
    filtered = util.apply_parallel(ndimage.median_filter,image_8bit,
                                     depth=int(np.ceil(filter_radius))+1,
                                     extra_keywords={'footprint': ball(filter_radius)})
    toc = time.perf_counter()
    print('Parallelized filtering done: time was', toc-tic, 'seconds')
    
    if borders is not None:
        histf, bins_centerf = exposure.histogram(filtered[~borders])
    else:
        histf, bins_centerf = exposure.histogram(filtered)
    
    
    plt.figure(num=29, figsize=(6,6))
    plt.clf()
    plt.plot(bins_centerf, 1e-6*histf, lw=2,label='Filtered')
    try:
        plt.plot(bins_centers, 1e-6*hists, lw=2,label='Rescaled')
    except Exception:
        pass
    plt.ylabel('Number of voxels (millions)')
    plt.xlabel('Dark                     Intensity value                        Light')
    plt.legend()

if binary:
   
    print('Thresholding...')
    
     #Otsu thresholding
    
    if borders is not None:
        val  = filters.threshold_otsu(filtered[~borders])
    else:
        val  = filters.threshold_otsu(filtered)
    print('Thresholded, threshold value is ' + str(val))
    
    #Binary image (void voxels = 1)
    
    bin_img = filtered <= val
    
    #The region outside the cylinder is set to be solid
    
    if borders is not None:
        bin_img[borders] = False
    
    #Identification of floating solid voxels
    
    print('Finding disconnected solids...')
    bin_img2 = np.copy(bin_img)
    disconn = ps.filters.find_disconnected_voxels(~bin_img2)
    floating_solid = np.sum(disconn)
    del bin_img2
    print('Done')
    
    #Removal of floating solid voxels from the final binary image
    
    bin_img_nofloat = np.copy(bin_img)
    bin_img_nofloat[disconn] = True
    del bin_img, disconn
    
    #del whole_image, rotated, syl_rot, cropped_rot
    
    slic = 350

    plt.figure(num=501, figsize=(11, 4))
    plt.clf()
    ax1 = plt.subplot(131)
    plt.imshow(image_8bit[slic,:,:], cmap='gray', interpolation='nearest')#, vmin=0, vmax=255)
    plt.title('Original 8 bit')
    ax2 = plt.subplot(132, sharex=ax1, sharey=ax1)    
    plt.imshow(filtered[slic,:,:], cmap='gray', interpolation='nearest', vmin=0, vmax=255)
    plt.title('Filtered')
    ax3 = plt.subplot(133, sharex=ax1, sharey=ax1)    
    plt.imshow(bin_img_nofloat[slic,:,:], cmap='gray', interpolation='nearest')
    plt.title('Otsu thresholding')
    plt.tight_layout()

if poros:
    
    print('Starting porosity calculations...')
    
    try:    
        p = ps.metrics.porosity(bin_img_nofloat[~borders])    
        p_prof = ps.metrics.porosity_profile(bin_img_nofloat, 0) / caf_real    
    except:
        p = ps.metrics.porosity(bin_img_nofloat)    
        p_prof = ps.metrics.porosity_profile(bin_img_nofloat, 0)
    
    fig = plt.figure(num=94)
    plt.clf()
    fig.set_size_inches([4,4])
    plt.plot(p_prof,label='Air-filled porosity')
    plt.legend()
    plt.xlabel('Distance from top (voxels)')
    plt.ylabel('Air-filled porosity (%)')
    plt.tight_layout()

if netgen:
    
    # 4A 170:710 # 4B 110:940 # 5A 50:950
    try:
        if sample_name is not None:
            parameters = mm.read_params_name(paramfile)
            names = parameters.filenames
            domain_params = parameters.domain_params    
            try:
                scan = names.index(sample_name) + 1
            except ValueError:
                print("Sample does not exist")
        z_top = domain_params[scan-1,0]
        z_bottom = domain_params[scan-1,1]
    except Exception:
        z_top = 200
        z_bottom = 950   
    hmin = 0
    hmax = np.shape(bin_img_nofloat)[1]
    #z_top = 200; z_bottom = 950; hmin = 150; hmax = 750
    
    print('Vertical network domain from', str(z_top), 'to', str(z_bottom) )
    
    ws = op.Workspace()
    ws.clear()
    
    proj1 = ws.new_project(name = 'full_image')
    
    np.random.seed(42)
    binaariverkko = (ps.networks.snow(bin_img_nofloat[z_top:z_bottom,hmin:hmax,hmin:hmax],
                     voxel_size = voxel_size, boundary_faces=['left','right']))
    
    pn = op.network.GenericNetwork(project=proj1)
    pn.update(binaariverkko)
    
    areas = ps.metrics.region_surface_areas(regions=binaariverkko.regions)
    pn['pore.surface_area'] = areas * voxel_size**2
    
    geom = op.geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
    
    proj2 = ws.new_project(name = 'full_image_all')
    pn_all = op.network.GenericNetwork(project=proj2)
    pn_all.update(binaariverkko)

    pn_all['pore.surface_area'] = areas * voxel_size**2
    geom_all = op.geometry.GenericGeometry(network=pn_all, pores=pn_all.Ps, throats=pn_all.Ts)
    
    hech = pn.check_network_health()
    op.topotools.trim(network=pn, pores=hech['trim_pores']) 
    
    reg_image = np.copy(binaariverkko.regions)
    
    vol_tot = (z_bottom-z_top)*np.shape(reg_image)[1]*np.shape(reg_image)[2]
    
    #Draw a voxel image of a pore
    #mm.single_pore(duu,33)
    
    print('Porosity, total image',
          ps.metrics.porosity(bin_img_nofloat[~borders]))
    
    print('Porosity, image section',
          ps.metrics.porosity(bin_img_nofloat[z_top:z_bottom,:,:])/caf_real)
    
    print('Porosity, total pore space',
          np.sum(pn_all['pore.volume'][~pn_all['pore.boundary']])/ vol_tot /caf_real / voxel_size**3)
    
    print('Porosity, trimmed network',
          np.sum(pn['pore.volume'][~pn['pore.boundary']]) / vol_tot /caf_real / voxel_size**3)

    
if percolation:

    
    air = op.phases.Air(network=pn)
    phys_air = op.physics.Standard(network=pn, phase=air, geometry=geom)
    
    air['pore.contact_angle'] = 180.0
    #air['pore.entry_pressure'] = 0.0
    
    qw = np.copy(pn['pore.volume'])
    
    pn['pore.volume_percolation'] = pn['pore.volume']
    pn['pore.volume_percolation'][pn['pore.boundary']] = 0.
    
    print('Air intrusion simulation')
    
    #Ordinary percolation
    perc = op.algorithms.OrdinaryPercolation(network=pn)
    perc.setup(phase=air)
    perc.set_inlets(pores=pn['pore.left'])
    perc.setup(pore_volume='pore.volume_percolation', throat_volume='throat.volume')
    perc.run(points = np.unique(phys_air['throat.entry_pressure']))

    perc_data = perc.get_intrusion_data()
    perc_data_new = np.copy(perc_data).T
    
    plt.figure(num=842, figsize=[7, 7])
    plt.clf()
    plt.semilogy(1-np.asarray(perc_data.Snwp), perc_data.Pcap, 'b-', label='Ordinary Percolation - Bond')
    plt.legend()
    
    #Relative diffusion
    
    print('Effective diffusion simulation')
    
    air.update(perc.results(Pc=1000))

    #phys_air.add_model(propname='throat.entry_pressure', model=op.models.physics.capillary_pressure.washburn)
    phys_air.add_model(model=op.models.physics.multiphase.conduit_conductance,
                   propname='throat.conduit_diffusive_conductance',
                   throat_conductance='throat.diffusive_conductance')    
    
    diff_air = []
    sat= []
    tot_vol = np.sum(pn["pore.volume"]) + np.sum(pn["throat.volume"])
    
    press_diffusion = np.unique(perc['pore.invasion_pressure'])
    
    for Pc in press_diffusion:
        air.update(perc.results(Pc=Pc))
        phys_air.regenerate_models()
        this_sat = 0
        this_sat += np.sum(pn["pore.volume"][air["pore.occupancy"] == 1])
        sat.append(this_sat)
        BC1_pores = pn['pore.left']
        BC2_pores = pn['pore.right']
        FD_1 = op.algorithms.FickianDiffusion(network=pn)
        FD_1.setup(phase=air, conductance='throat.conduit_diffusive_conductance')
        FD_1.set_value_BC(values=0.6, pores=BC1_pores)
        FD_1.set_value_BC(values=0.2, pores=BC2_pores)
        FD_1.run()
        eff_diff = FD_1.calc_effective_diffusivity(domain_area=caf_real*((hmax-hmin)*voxel_size)**2,
                                                   domain_length= (z_bottom-z_top)*voxel_size)
        diff_air.append(eff_diff)
        pn.project.purge_object(FD_1)
        #pdb.set_trace()
    
    sat = np.asarray(sat)
    sat /= tot_vol
    
    rel_diff_air = np.asarray(diff_air)
    rel_diff_air /= rel_diff_air[-1]
    
    print('Network effective diffusivity is', diff_air[-1], '(Relative', diff_air[-1]/air['pore.diffusivity'][0], ')')
    
    plt.figure(188)
    plt.clf()
    plt.plot(sat, np.asarray(rel_diff_air), '^-r')
    plt.xlabel('Air saturation (Air volume / total void space volume)')
    plt.ylabel('Relative diffusivity (relative to totally dry pore conditions)')
    
    '''
    
    fig = plt.figure(num=40,figsize=(9,9))
    plt.clf()
    ax = fig.gca(projection='3d')
    ax.scatter(1000*pn['pore.coords'][:,0],
           1000*pn['pore.coords'][:,1],
           1000*pn['pore.coords'][:,2],
           c='#000080', s = 5e6*(pn['pore.equivalent_diameter'])**2,
           cmap=plt.cm.jet)
    
    fig = plt.figure(num=41,figsize=(9,9))
    plt.clf()
    ax = fig.gca(projection='3d')
    ax.scatter(1000*pn_all['pore.coords'][:,0],
           1000*pn_all['pore.coords'][:,1],
           1000*pn_all['pore.coords'][:,2],
           c='#000080', s = 5e6*(pn_all['pore.equivalent_diameter'])**2,
           cmap=plt.cm.jet)
    
    fig = op.topotools.plot_coordinates(network=pn, pores=pn.pores('internal'), c='b')
    fig.set_size_inches([5, 5])
    Ts = pn.find_neighbor_throats(pores=pn.pores('internal'), mode='xnor')
    fig = op.topotools.plot_connections(network=pn, throats=Ts, fig=fig, c='r')
    fig.axes[0].set_xlim((0.0,np.max(pn['pore.coords'][:,0])))
    fig.axes[0].set_ylim((0.0,np.max(pn['pore.coords'][:,1])))
    fig.axes[0].set_zlim((0.0,np.max(pn['pore.coords'][:,2])))
    plt.tight_layout()   
    
    fig = op.topotools.plot_coordinates(network=pn_all, pores=pn_all.pores('internal'), c='b')
    fig.set_size_inches([5, 5])
    Ts = pn_all.find_neighbor_throats(pores=pn_all.pores('internal'), mode='xnor')
    fig = op.topotools.plot_connections(network=pn_all, throats=Ts, fig=fig, c='r')
    fig.axes[0].set_xlim((0.0,np.max(pn_all['pore.coords'][:,0])))
    fig.axes[0].set_ylim((0.0,np.max(pn_all['pore.coords'][:,1])))
    fig.axes[0].set_zlim((0.0,np.max(pn_all['pore.coords'][:,2])))
    plt.tight_layout()  
    
    
    fig = plt.figure(num=411,figsize=(9,9))
    plt.clf()
    ax = fig.gca(projection='3d')
    ax.scatter(1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][0],0],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][0],1],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][0],2],
           c='#000080', s = 5000e6*(pn_5w['pore.equivalent_diameter'][hech_5['disconnected_clusters'][0]])**2,
           cmap=plt.cm.jet)
    ax.scatter(1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][1],0],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][1],1],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][1],2],
           c='c', s = 5000e6*(pn_5w['pore.equivalent_diameter'][hech_5['disconnected_clusters'][1]])**2,
           cmap=plt.cm.jet)
    ax.scatter(1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][2],0],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][2],1],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][2],2],
           c='m', s = 5000e6*(pn_5w['pore.equivalent_diameter'][hech_5['disconnected_clusters'][2]])**2,
           cmap=plt.cm.jet)
    ax.scatter(1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][3],0],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][3],1],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][3],2],
           c='r', s = 5000e6*(pn_5w['pore.equivalent_diameter'][hech_5['disconnected_clusters'][3]])**2,
           cmap=plt.cm.jet)
    ax.scatter(1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][4],0],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][4],1],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][4],2],
           c='k', s = 5000e6*(pn_5w['pore.equivalent_diameter'][hech_5['disconnected_clusters'][4]])**2,
           cmap=plt.cm.jet)
    ax.scatter(1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][5],0],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][5],1],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][5],2],
           c='g', s = 5000e6*(pn_5w['pore.equivalent_diameter'][hech_5['disconnected_clusters'][5]])**2,
           cmap=plt.cm.jet)
    ax.scatter(1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][6],0],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][6],1],
           1000*pn_5w['pore.coords'][hech_5['disconnected_clusters'][6],2],
           c='y', s = 5000e6*(pn_5w['pore.equivalent_diameter'][hech_5['disconnected_clusters'][6]])**2,
           cmap=plt.cm.jet)
  
    '''
if False:
    
    AAAAA = []
    for i in range(np.shape(syl_rot)[0]):
        Aq = syl_rot[i,:,:]
        Bq = borders[i,:,:]
        Aq[Bq] = 1e6
        B = np.sum(Aq<=3800)
        if B > -10:
            AAAAA.append([i,B])
    AAAAA = np.asarray(AAAAA)
    plt.figure(num=600)
    plt.clf()
    ax1=plt.subplot(1,1,1)
    plt.plot(AAAAA[:,0],AAAAA[:,1])
    ax1.set_yscale('log')
    
    #Calculate more accurate surface areas
    import omatesti
    volume_tests, surface_areas, sphericities = omatesti.regionprops_3D(binaariverkko.regions)
    pn['pore.surface_area_correct'] = surface_areas*voxel_size**2
    
    fig = plt.figure(num=3643)
    fig.set_size_inches([8,4])
    plt.clf()
    plt.subplot(1,2,1)
    plt.hist(sphericities[pn_all['pore.internal']],50)
    plt.subplot(1,2,2)
    plt.scatter(volume_tests[pn_all['pore.internal']],sphericities[pn_all['pore.internal']])
    
    vim1 = ps.networks.generate_voxel_image(pn,max_dim=836)
    mm.osakuva(vim1[3,:,:],7890,cmap=plt.cm.viridis)

if False:
        # Karsitaan ei-kytketty huokostila pois ja tehdään siitä toinen verkko
        bin_img_percolating = ps.filters.trim_nonpercolating_paths(bin_img_nofloat)
        np.random.seed(42)
        binaariverkko_percolating = (ps.networks.snow(bin_img_percolating[z_top:z_bottom,hmin:hmax,hmin:hmax],
                     voxel_size = voxel_size, boundary_faces=['left','right']))
        
        proj3 = ws.new_project(name = 'full_image_percolated')
        pn_perc = op.network.GenericNetwork(project=proj3)
        pn_perc.update(binaariverkko_percolating)
    
        geom_perc = op.geometry.GenericGeometry(network=pn_perc, pores=pn_perc.Ps, throats=pn_perc.Ts)
        
        
        proj4 = ws.new_project(name = 'full_image_percolated_trimmed')
        pn_perc_tr = op.network.GenericNetwork(project=proj4)
        pn_perc_tr.update(binaariverkko_percolating)
    
        geom_perc_tr = op.geometry.GenericGeometry(network=pn_perc_tr, pores=pn_perc_tr.Ps, throats=pn_perc_tr.Ts)
        
        hech_perc = pn_perc_tr.check_network_health()
        
        op.topotools.trim(network=pn_perc_tr, pores=hech_perc['trim_pores']) 
        
        print('Porosity, total percolated pore space',
          np.sum(pn_perc['pore.volume'][~pn_perc['pore.boundary']])/ vol_tot /caf_real / voxel_size**3)
        print('Porosity, trimmed percolated pore space',
          np.sum(pn_perc_tr['pore.volume'][~pn_perc_tr['pore.boundary']])/ vol_tot /caf_real / voxel_size**3)

    #----------------------
    
    
#ws.save_workspace(filename = 'D:/network_4Anew_330_860')
#np.save('D:/binary_4Anew.npy', bin_img_nofloat)