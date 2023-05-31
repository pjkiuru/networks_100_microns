# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 13:59:29 2023

@author: pkiuru

Script for drawing 3D images of all pores in the network domain and the largest
connected pore clusters divided into four vertical sections

Required input:
    OpenPNM network object (= output of the script file net_poros.py or output 
    of the script file samp_an_100.py)

"""
# If True, the boundary boundary pores are illustrated in red colour
draw_border = True

a = pn['pore.coords'][:,0]
amax = np.max(pn['pore.coords'][:,0])
internal_pores = pn['pore.internal']

AA = a < 0.25*amax
AB = np.logical_and(a>=0.25*amax, a<=0.5*amax)
AC = np.logical_and(a>=0.5*amax, a<=0.75*amax)
AD = a > 0.75*amax

angles = [20,0]

# The 4 vertical sections of the largest connected pore cluster, viewed from bottom to top

fig = plt.figure(num=500,figsize=(19,5))
plt.clf()

ax = plt.subplot(1,4,1,projection='3d')
ax.scatter(1000*pn['pore.coords'][AA,0],
       1000*pn['pore.coords'][AA,1],
       1000*pn['pore.coords'][AA,2],
       c=a[AA], s = 5e6*(pn['pore.equivalent_diameter'][AA])**2,
       cmap=plt.cm.Blues)
if draw_border:
    ax.scatter(1000*pn['pore.coords'][pn['pore.left'],0],
           1000*pn['pore.coords'][pn['pore.left'],1],
           1000*pn['pore.coords'][pn['pore.left'],2],
           c='r', s = 5e6*(pn['pore.equivalent_diameter'][pn['pore.left']])**2,
       )
ax.w_xaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_yaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_zaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.view_init(elev=angles[0], azim=angles[1])

ax = plt.subplot(1,4,2,projection='3d')
ax.scatter(1000*pn['pore.coords'][AB,0],
       1000*pn['pore.coords'][AB,1],
       1000*pn['pore.coords'][AB,2],
       c=a[AB], s = 5e6*(pn['pore.equivalent_diameter'][AB])**2,
       cmap=plt.cm.Blues)
ax.w_xaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_yaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_zaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.view_init(elev=angles[0], azim=angles[1])

ax = plt.subplot(1,4,3,projection='3d')
ax.scatter(1000*pn['pore.coords'][AC,0],
       1000*pn['pore.coords'][AC,1],
       1000*pn['pore.coords'][AC,2],
       c=a[AC], s = 5e6*(pn['pore.equivalent_diameter'][AC])**2,
       cmap=plt.cm.Blues)
ax.w_xaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_yaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_zaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.view_init(elev=angles[0], azim=angles[1])

ax = plt.subplot(1,4,4,projection='3d')
ax.scatter(1000*pn['pore.coords'][AD,0],
       1000*pn['pore.coords'][AD,1],
       1000*pn['pore.coords'][AD,2],
       c=a[AD], s = 5e6*(pn['pore.equivalent_diameter'][AD])**2,
       cmap=plt.cm.Blues)
if draw_border:
    ax.scatter(1000*pn['pore.coords'][pn['pore.right'],0],
           1000*pn['pore.coords'][pn['pore.right'],1],
           1000*pn['pore.coords'][pn['pore.right'],2],
           c='r', s = 5e6*(pn['pore.equivalent_diameter'][pn['pore.right']])**2,
       )
ax.w_xaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_yaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_zaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.view_init(elev=angles[0], azim=angles[1])

plt.tight_layout()


a = pn_all['pore.coords'][:,0]
amax = np.max(pn_all['pore.coords'][:,0])

AA = a < 0.25*amax
AB = np.logical_and(a>=0.25*amax, a<=0.5*amax)
AC = np.logical_and(a>=0.5*amax, a<=0.75*amax)
AD = a > 0.75*amax

# The 4 vertical sections of the pore space, viewed from bottom to top

fig = plt.figure(num=502,figsize=(19,5))
plt.clf()

ax = plt.subplot(1,4,1,projection='3d')
ax.scatter(1000*pn_all['pore.coords'][AA,0],
       1000*pn_all['pore.coords'][AA,1],
       1000*pn_all['pore.coords'][AA,2],
       c=a[AA], s = 5e6*(pn_all['pore.equivalent_diameter'][AA])**2,
       cmap=plt.cm.Blues)
ax.w_xaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_yaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_zaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.view_init(elev=angles[0], azim=angles[1])

ax = plt.subplot(1,4,2,projection='3d')
ax.scatter(1000*pn_all['pore.coords'][AB,0],
       1000*pn_all['pore.coords'][AB,1],
       1000*pn_all['pore.coords'][AB,2],
       c=a[AB], s = 5e6*(pn_all['pore.equivalent_diameter'][AB])**2,
       cmap=plt.cm.Blues)
ax.w_xaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_yaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_zaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.view_init(elev=angles[0], azim=angles[1])

ax = plt.subplot(1,4,3,projection='3d')
ax.scatter(1000*pn_all['pore.coords'][AC,0],
       1000*pn_all['pore.coords'][AC,1],
       1000*pn_all['pore.coords'][AC,2],
       c=a[AC], s = 5e6*(pn_all['pore.equivalent_diameter'][AC])**2,
       cmap=plt.cm.Blues)
ax.w_xaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_yaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_zaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.view_init(elev=angles[0], azim=angles[1])

ax = plt.subplot(1,4,4,projection='3d')
ax.scatter(1000*pn_all['pore.coords'][AD,0],
       1000*pn_all['pore.coords'][AD,1],
       1000*pn_all['pore.coords'][AD,2],
       c=a[AD], s = 5e6*(pn_all['pore.equivalent_diameter'][AD])**2,
       cmap=plt.cm.Blues)
ax.w_xaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_yaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_zaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.view_init(elev=angles[0], azim=angles[1])

plt.tight_layout()


a = pn['pore.coords'][:,0]
amax = np.max(pn['pore.coords'][:,0])

angles = [30,210]

# LArgest connected pore cluster, viewed from top to bottom

fig = plt.figure(num=504,figsize=(4,4))
plt.clf()

ax = plt.subplot(1,1,1,projection='3d')
ax.scatter(1000*pn['pore.coords'][internal_pores,0],
       1000*pn['pore.coords'][internal_pores,1],
       1000*pn['pore.coords'][internal_pores,2],
       c=a[internal_pores], s = 5e5*(pn['pore.equivalent_diameter'][internal_pores])**2,
       cmap=plt.cm.Blues)
ax.w_xaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_yaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_zaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.view_init(elev=angles[0], azim=angles[1])

plt.tight_layout()