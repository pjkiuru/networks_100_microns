# -*- coding: utf-8 -*-
"""
Created on Wed May  3 10:26:35 2023

@author: pkiuru

Script file for
(1) loading the binary image and the network object
(2) calculating the radial porosity distribution of a binary
    solid-void image
(3) calculating the number fraction and the volume fraction of pores near the
    cylinder walls

Required input:
    Binary image (.npy file)
    Network object (.pnm file)
    Parameter file (Scan_parameters_2.csv)

"""

import openpnm as op
import porespy.tools as pst
import numpy as np
import matplotlib.pyplot as plt
import my_models as mm
import pdb
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

from matplotlib import rcParams

rcParams['figure.max_open_warning'] = 40
rcParams['font.size'] = 8


def pi_by_four(radius):
    square = np.ones((1,radius,radius))
    cylinder = pst.extract_cylinder(square, axis=0)
    pi_by_four = np.sum(cylinder)/np.sum(square)
    
    return pi_by_four


def radial_porosity(img,pi_by_four,nslices=None):
    
    #The image domain is divided into hollow cylinders with equal volumes
    
    a = list(img.shape)
    a.pop(0)
    #r = np.floor(np.amin(a) / 2)
    dim = [range(int(-s / 2), int(s / 2) + s % 2) for s in img.shape]
    inds = np.meshgrid(*dim, indexing='ij')
    inds[0] = inds[0] * 0
    d = np.sqrt(np.sum(np.square(inds), axis=0))
    
    if nslices is None:
        nslices = 2
    
    #rs = np.linspace(0,r,nslices + 1)
    
    n_total = pi_by_four * np.size(img)
    
    Vs = np.linspace(0,n_total,nslices + 1)
    
    rs_by_volume = np.sqrt(Vs/pi_by_four/img.shape[0])/2
    rs = rs_by_volume
    
    por_slices = np.zeros(nslices) 
    n_void = np.zeros(nslices) 
    n_tot = np.zeros(nslices) 
    
    for i in range(nslices):
        mask = np.logical_and(d>=rs[i],d<rs[i+1])
        n_void[i] = np.sum(img[mask])
        n_tot[i] = np.sum(mask)
        por_slices[i] = n_void[i]/n_tot[i]
    
    return por_slices, n_void, n_tot

def radial_porosity_image(por_slices,diam):
    
    square = np.ones((1,diam,diam))
    cylinder = pst.extract_cylinder(square, axis=0)
    
    a = list(cylinder.shape)
    a.pop(0)
    r = np.floor(np.amin(a) / 2)
    dim = [range(int(-s / 2), int(s / 2) + s % 2) for s in cylinder.shape]
    inds = np.meshgrid(*dim, indexing='ij')
    inds[0] = inds[0] * 0
    d = np.sqrt(np.sum(np.square(inds), axis=0))[0]
    
    nslices = np.size(por_slices)
    
    #rs = np.linspace(0,r,nslices + 1)
    
    pi_by_four = np.sum(cylinder)/np.sum(square)
    
    Vs = np.linspace(0,4*pi_by_four*r**2,nslices + 1)
    
    rs_by_volume = np.sqrt(Vs/pi_by_four)/2
    rs = rs_by_volume
    pcol = np.zeros((diam, diam))
    
    for i in range(nslices):
        mask = np.logical_and(d>=rs[i],d<rs[i+1])
        pcol[mask] = por_slices[i]
        
        
    fig = plt.figure(930)
    plt.clf()
    fig.set_size_inches([4,4])
    plt.imshow(pcol/cylinder[0], cmap=plt.cm.plasma, interpolation='nearest',vmin=0)
    plt.title('Radial air-filled porosity distribution \n (network domain, total macropore space)')
    plt.colorbar()
    plt.axis('off')
    plt.tight_layout()
    
    fig = plt.figure(931)
    plt.clf()
    fig.set_size_inches([4,4])
    plt.plot(0.1*rs[1:],por_slices,'b',lw=0.5)
    plt.plot(0.1*rs[1:],por_slices,'bo',markersize=2)
    plt.xlabel('Radius (mm)')
    plt.ylabel('Air-filled porosity')
    plt.title('Radial air-filled porosity distribution \n (network domain, total macropore space)')


# Load the image and network data and calculate total neteork domain volume

if True:
    
    topfolder = 'D:/'
    
    sample_name = '7B'
    
    paramfile = 'C:\\Users\\pkiuru\\OneDrive - University of Eastern Finland\\METNET_scripts\\Scan_parameters_2.csv'
    
    
    bin_img_nofloat = np.load(topfolder + 'binary_' + sample_name + '.npy')
    
    voxel_size = 100e-6
    
    if sample_name is not None:
        parameters = (mm.read_params_name(paramfile))
        names = parameters.filenames
        crop_coords = parameters.crop_coords
        domain_params = parameters.domain_params
        scan = names.index(sample_name) + 1
    
    ws = op.Workspace()
    ws.clear()
    #ws.load_workspace(filename='D:/network_1A_110_970.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_1B_120_970.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_2A_110_970.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_2B_220_970.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_3A_140_940.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_3B_50_950.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_4Anew_328_868.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_4B_110_940.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_5A_100_950.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_5B_200_950.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_6A_130_970.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_6B_130_970.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_7A_160_970.pnm', overwrite=True)    
    ws.load_workspace(filename='D:/network_7B_100_960.pnm', overwrite=True)    
    #ws.load_workspace(filename='D:/network_8A_200_970.pnm', overwrite=True)
    #ws.load_workspace(filename='D:/network_8B_70_930.pnm', overwrite=True)
    
    #Trimmed, i.e., connected network
    proj = ws['full_image']
    pn = proj['net_01']
    geom = proj['geo_01']
    
    #All networks (incl. isolated pores in the domain)
    proj_all = ws['full_image_all']
    pn_all = proj_all['net_01']
    geom_all = proj_all['geo_01']
    
    hech = pn_all.check_network_health()
    
    Ps_front = pn['pore.left']
    
    #Set the volumes of the artificial boundary pores to zero
    pn['pore.volume'][pn['pore.boundary']] = 0.
    pn_all['pore.volume'][pn_all['pore.boundary']] = 0.
    
    #Total sample volume (m3)
    vol_sample_appr = (np.max(pn['pore.coords'][:,0]) * np.max(pn['pore.coords'][:,1])
                  * np.max(pn['pore.coords'][:,2]) * np.pi/4)
    
    height = np.diff(domain_params[scan-1,:])[0]
    width = np.diff(crop_coords[scan-1,2:4])
    
    #square = np.ones((1,int(width),int(width)))
    #cylinder = pst.extract_cylinder(square, axis=0)
    #caf_real = np.sum(cylinder)/np.sum(square)
    caf_real = pi_by_four(int(width))
    
    vol_sample = caf_real*height*width**2*voxel_size**3


porosity_total = np.sum(pn_all['pore.volume']) / vol_sample
print('Porosity, total pore space', porosity_total[0])

porosity_trimmed = np.sum(pn['pore.volume']) / vol_sample
print('Porosity, trimmed network', porosity_trimmed[0])


#Select the network domain from the binary image
z_top, z_bottom = domain_params[scan-1,0:2]
img = bin_img_nofloat[z_top:z_bottom,:,:]
diameter = np.shape(img)[1]

#Calculate the radial porosity distribution. The image domain is divided into
#hollow cylinders with equal volumes.
if True:

    pi_by_4 = pi_by_four(int(width))
    
    por_slices, n_void, n_tot = radial_porosity(img,pi_by_4,45)
    
    n_void_total = np.sum(img)
    n_total = pi_by_4 * np.size(img)
    
    print('Domain porosity',n_void_total/n_total)
    
    #Check for the equality of the volumes of the hollow cylinders
    volfracs = n_tot/n_total
    
    radial_porosity_image(por_slices,diameter)
    
    #radial_porosity_image(np.linspace(0,45,45),diameter)

#Calculate the number fraction and the volume fraction of pores near the cylinder walls
if True:
    
    center = 0.5*diameter*voxel_size
    pore_dist_radial = (pn['pore.coords'][pn['pore.internal'],1]-center)**2 + (pn['pore.coords'][pn['pore.internal'],2]-center)**2
    
    plt.figure(num=1)
    plt.clf()
    plt.hist(np.sqrt(pore_dist_radial),45)
    
    near_wall_1mm = pore_dist_radial > (0.5*diameter*voxel_size-0.001)**2#(0.475*diameter*voxel_size)**2
    near_wall_2mm = pore_dist_radial > (0.5*diameter*voxel_size-0.002)**2#(0.475*diameter*voxel_size)**2
    
    volumes_internal = pn['pore.volume'][pn['pore.internal']]
    near_wall_1mm_volume_fraction = np.sum(volumes_internal[near_wall_1mm]/np.sum(volumes_internal))
    near_wall_2mm_volume_fraction = np.sum(volumes_internal[near_wall_2mm]/np.sum(volumes_internal))
    
    print(np.sum(near_wall_1mm)/np.sum(pn['pore.internal']), 'of pores are < 1 mm the wall ')
    print(near_wall_1mm_volume_fraction, 'of pore volume is < 1 mm from the wall ')
    
    print(np.sum(near_wall_2mm)/np.sum(pn['pore.internal']), 'of pores are < 2 mm the wall ')
    print(near_wall_2mm_volume_fraction, 'of pore volume is < 2 mm from the wall ')


'''
Q = np.zeros([10,10,10]).astype('bool')
Qcub = np.ones([10,10,10])
Q[3,4,4] = 1
Q[5,4,4] = 1
Q[7,4,4] = 1
Q[5:8,1,4] = 1
print(radial_porosity(Q))
Qcyl = pst.extract_cylinder(Qcub, axis=0)
'''