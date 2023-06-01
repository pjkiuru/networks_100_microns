# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 10:48:57 2023

@author: pkiuru

Auxiliary functions for manipulating an OpenPNM network with cubic pores.

"""

import numpy as np

def extract_cyl(target,nx,space):
    xy = np.copy(target['pore.coords'][:,1:3])
    r = nx / 2
    dr = r * space
    xy -= dr
    d = np.sqrt(np.sum(np.square(xy), axis=1))
    mask = d > dr
    #Simport pdb; pdb.set_trace()
    return mask

def find_center(target,nx,space):
    xy = np.copy(target['pore.coords'])
    r = nx / 2
    dr = r * space
    xy -= dr
    d = np.sqrt(np.sum(np.square(xy), axis=1))
    #import pdb; pdb.set_trace()
    min_value = np.where(d == np.min(d))
    if type(min_value)=='int':
        A = min_value
    else:
        A = min_value[0][0]
    return A