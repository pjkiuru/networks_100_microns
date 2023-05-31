"""
Script for drawing 3D images of all pores in the network domain and 11 largest
connected pore clusters

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