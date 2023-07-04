"""
Script for drawing 3D images of all pores in the network domain and 11 largest
connected pore clusters. Also, an outline of images to Supplementary material.

Required input:
    OpenPNM network object (= output of the script file net_poros.py or output 
    of the script file samp_an_100.py)

"""
# All pores

fig = plt.figure(num=41,figsize=(7,7))
plt.clf()
ax = fig.gca(projection='3d')
ax.scatter(1000*pn_all['pore.coords'][:,0],
       1000*pn_all['pore.coords'][:,1],
       1000*pn_all['pore.coords'][:,2],
       c='#000080', s = 5e6*(pn_all['pore.equivalent_diameter'])**2,
       cmap=plt.cm.jet)

# 11 largest networks

fig = plt.figure(num=411,figsize=(7,7))
plt.clf()
ax = fig.gca(projection='3d')
ax.scatter(1000*pn_all['pore.coords'][hech['disconnected_clusters'][0],0],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][0],1],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][0],2],
       c='#000080', s = 5e6*(pn_all['pore.equivalent_diameter'][hech['disconnected_clusters'][0]])**2,
       cmap=plt.cm.jet)
ax.scatter(1000*pn_all['pore.coords'][hech['disconnected_clusters'][1],0],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][1],1],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][1],2],
       c='c', s = 5e6*(pn_all['pore.equivalent_diameter'][hech['disconnected_clusters'][1]])**2,
       cmap=plt.cm.jet)
ax.scatter(1000*pn_all['pore.coords'][hech['disconnected_clusters'][2],0],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][2],1],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][2],2],
       c='m', s = 5e6*(pn_all['pore.equivalent_diameter'][hech['disconnected_clusters'][2]])**2,
       cmap=plt.cm.jet)
ax.scatter(1000*pn_all['pore.coords'][hech['disconnected_clusters'][3],0],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][3],1],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][3],2],
       c='r', s = 5e6*(pn_all['pore.equivalent_diameter'][hech['disconnected_clusters'][3]])**2,
       cmap=plt.cm.jet)
ax.scatter(1000*pn_all['pore.coords'][hech['disconnected_clusters'][4],0],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][4],1],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][4],2],
       c='k', s = 5e6*(pn_all['pore.equivalent_diameter'][hech['disconnected_clusters'][4]])**2,
       cmap=plt.cm.jet)
ax.scatter(1000*pn_all['pore.coords'][hech['disconnected_clusters'][5],0],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][5],1],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][5],2],
       c='g', s = 5e6*(pn_all['pore.equivalent_diameter'][hech['disconnected_clusters'][5]])**2,
       cmap=plt.cm.jet)
ax.scatter(1000*pn_all['pore.coords'][hech['disconnected_clusters'][6],0],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][6],1],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][6],2],
       c='y', s = 5e6*(pn_all['pore.equivalent_diameter'][hech['disconnected_clusters'][6]])**2,
       cmap=plt.cm.jet)
ax.scatter(1000*pn_all['pore.coords'][hech['disconnected_clusters'][7],0],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][7],1],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][7],2],
       c='gray', s = 5e6*(pn_all['pore.equivalent_diameter'][hech['disconnected_clusters'][7]])**2,
       cmap=plt.cm.jet)
ax.scatter(1000*pn_all['pore.coords'][hech['disconnected_clusters'][8],0],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][8],1],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][8],2],
       c='c', s = 5e6*(pn_all['pore.equivalent_diameter'][hech['disconnected_clusters'][8]])**2,
       cmap=plt.cm.jet)
ax.scatter(1000*pn_all['pore.coords'][hech['disconnected_clusters'][9],0],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][9],1],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][9],2],
       c='k', s = 5e6*(pn_all['pore.equivalent_diameter'][hech['disconnected_clusters'][9]])**2,
       cmap=plt.cm.jet)
ax.scatter(1000*pn_all['pore.coords'][hech['disconnected_clusters'][10],0],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][10],1],
       1000*pn_all['pore.coords'][hech['disconnected_clusters'][10],2],
       c='y', s = 5e6*(pn_all['pore.equivalent_diameter'][hech['disconnected_clusters'][10]])**2,
       cmap=plt.cm.jet)


# Side views of all pores and the largest connected network

a = pn['pore.coords'][:,0]
amax = np.max(pn['pore.coords'][:,0])
angles = [10,270]
internal_pores = pn['pore.internal']

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
#ax.set_xlim((0.0,1000*0.105))
ax.set_ylim((0.0,1000*0.095))
ax.set_zlim((0.0,1000*0.095))
ax.view_init(elev=angles[0], azim=angles[1])
plt.tight_layout()

a = pn_all['pore.coords'][:,0]
amax = np.max(pn_all['pore.coords'][:,0])
internal_pores = pn_all['pore.internal']

fig = plt.figure(num=506,figsize=(4,4))
plt.clf()
ax = plt.subplot(1,1,1,projection='3d')
ax.scatter(1000*pn_all['pore.coords'][internal_pores,0],
       1000*pn_all['pore.coords'][internal_pores,1],
       1000*pn_all['pore.coords'][internal_pores,2],
       c=a[internal_pores], s = 5e5*(pn_all['pore.equivalent_diameter'][internal_pores])**2,
       cmap=plt.cm.Blues)
ax.w_xaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_yaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
ax.w_zaxis.set_pane_color((0.6, 0.6, 0.6, 0.6))
#ax.set_xlim((0.0,1000*0.105))
ax.set_ylim((0.0,1000*0.095))
ax.set_zlim((0.0,1000*0.095))
ax.view_init(elev=angles[0], azim=angles[1])
plt.tight_layout()




# Largest cluster, second largest cluster, all pores; the sample is upright

#In the beginning, pn_all contains all pores in the pore space

hech_all = pn_all.check_network_health()
op.topotools.trim(network=pn_all, pores=hech_all['disconnected_clusters'][0]) 

#Now pn_all contains all pores but the largest connected cluster

colors = list([(0,0,0), (230/255,159/255,0), (86/255,180/255,233/255),
           (0,158/255,115/255), (240/255,228/255,66/255),
           (0,114/255,178/255), (213/255,94/255,0),(204/255,121/255,167/255)]) 

hech_tr = pn_all.check_network_health()

angles = [-150,135]

#----------------------------------------------------------------------
#Auxiliary function for drawing the outline of a cylinder
def data_for_cylinder_along_z(center_x,center_y,radius,height_z):
    z = np.linspace(0, height_z, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

#Cylinder with radius 45 mm and height z (1A 86; 3A 80; 5A 85;7A 81) mm
Xc,Yc,Zc = data_for_cylinder_along_z(45,45,45,86)
#----------------------------------------------------------------------

fig = plt.figure(num=507)
#fig.set_size_inches([4.5,6.0]); scaling_factor = 1e6
fig.set_size_inches([3.0,4.0]); scaling_factor = 5e5
plt.clf()
ax = plt.subplot(1,1,1,projection='3d')

first = hech_tr['disconnected_clusters'][0]

dummy = np.array(colors[0]).reshape(1,-1)
ax.scatter(1000*pn_all['pore.coords'][first,1],
       1000*pn_all['pore.coords'][first,2],
       1000*pn_all['pore.coords'][first,0],
       c=dummy, s =scaling_factor*(pn_all['pore.equivalent_diameter'][first])**2,
       cmap=plt.cm.Blues)

op.topotools.trim(network=pn_all, pores=first) 

internal_pores = pn_all['pore.internal']
#internal_pores = np.logical_and(pn_all['pore.coords'][:,0] < 0.003,pn_all['pore.coords'][:,0] > 0.00011)
dummy = np.array(colors[5]).reshape(1,-1)
ax.scatter(1000*pn_all['pore.coords'][internal_pores,1],
       1000*pn_all['pore.coords'][internal_pores,2],
       1000*pn_all['pore.coords'][internal_pores,0],
       c=dummy, s = scaling_factor*(pn_all['pore.equivalent_diameter'][internal_pores])**2,
       cmap=plt.cm.Blues)
internal_pores = pn['pore.internal']
#internal_pores = np.logical_and(pn['pore.coords'][:,0] < 0.003,pn['pore.coords'][:,0] > 0.00011)
dummy = np.array(colors[1]).reshape(1,-1)
ax.scatter(1000*pn['pore.coords'][internal_pores,1],
       1000*pn['pore.coords'][internal_pores,2],
       1000*pn['pore.coords'][internal_pores,0],
       c=dummy, s =scaling_factor*(pn['pore.equivalent_diameter'][internal_pores])**2,
       cmap=plt.cm.Blues)

ax.plot_surface(Xc, Yc, Zc, alpha=0.2,color='gray')

#ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
#ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
#ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
#ax.set_xlim((0.0,1000*0.065))
#ax.set_ylim((0.0,1000*0.065))
#ax.set_zlim((0.0,1000*0.085))
ax.view_init(elev=angles[0], azim=angles[1])
plt.axis('off')
ax.grid(b=None)
plt.axis('scaled')
ax.set_position([-0.425, -0.425, 1.8, 1.8])


#plt.savefig('pores_5A_small_200.png',dpi=200)