# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 14:03:20 2022

@author: pkiuru

Script for calculating network metrics for the cylindrical domain networks

Required input:
    OpenPNM network object (= output of the script file net_poros.py)

Output:
    Network metrics calculation results (printed in the console)

"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import my_models as mm
import openpnm as op
import networkx as nx

# Switch for creating a NetworkX object from the total largest connected pore
# network (including boundary pores) and from the internal pores of the network
make_networkx_object = True

# Calculation of vertical geometrical tortuosity
tortuosity_switch = True

#The following network metrics use the internal network only (boundary pores excluded)

#Calculation of pore coordination number
degree_switch = True

#Calculation of clustering coefficient
clust_switch = True

#Calculation of closeness centrality and betweenness centrality
centrality_switch = True

#Reduction of nodes included in tortuosity / centrality calculations by
#a factor 'space' (1 if the number of internal pores < 6500, otherwise 5)
space = 1

print('Internal pores in the domain:', np.sum(pn_all['pore.internal']))
print('Internal pores in the largest network:', np.sum(pn['pore.internal']))
print('Top boundary pores in the largest network:', np.sum(pn['pore.left']))
print('Bottom boundary pores in the largest network:', np.sum(pn['pore.right']))

if np.sum(pn['pore.internal']) > 6500:
    space = 5

if make_networkx_object:
    
    # Create a NetworkX object from the total pore network and from the internal
    # pores of the network
    
    pn['throat.conduitlength'] = pn['throat.total_length']
    pnx = op.io.NetworkX.to_networkx(pn)
    
    #Remove the boundary pores from the NetworkX object
    internal_nodes = [x for x,y in pnx.nodes(data=True) if y['internal']==True]        
    pnx_int = pnx.subgraph(nodes=internal_nodes)

if tortuosity_switch:

    path_lengths = mm.tortuosity(pnx,space)
    
    his = mm.parse_histogram(np.histogram(path_lengths, bins=25, density=True))
    
    fit_alpha, fit_loc, fit_beta = stats.gamma.fit(path_lengths)
        
       
    xax = np.linspace(stats.gamma.ppf(0.005, fit_alpha, fit_loc, fit_beta),
                          stats.gamma.ppf(0.995, fit_alpha, fit_loc, fit_beta), 100)
    
    print('')
    print('Vertical tortuosity:')
    print(np.mean(path_lengths))
    #print(stats.gamma.stats(fit_alpha, fit_loc, fit_beta,'mvsk'))        
            
    plt.figure(num=572,figsize=(5,5))
    plt.clf()
    plt.subplot(1,1,1)
    plt.bar(his.bin_centers, his.pdf, label='Calculated', width = his.bin_centers[1]-his.bin_centers[0])
    #plt.plot(xax, stats.gamma.pdf(xax, fit_alpha, fit_loc, fit_beta), c='sandybrown', label='Fitted gamma distribution')
    plt.xlabel('Tortuosity')
    plt.ylabel('Probability density')
    plt.legend()
          
    
if degree_switch:
    #Pore degree or coordination number
    
    #Pore coordination number calculated from the OpenPNM network

    connections = np.zeros(pn.Np)
    
    for j, k in enumerate(pn['throat.conns']):
        if pn['pore.internal'][pn['throat.conns'][j,1]]:
            connections[pn['throat.conns'][j,0]] += 1
            connections[pn['throat.conns'][j,1]] += 1
    
    connections_all = np.zeros(pn.Np)
    
    for j, k in enumerate(pn['throat.conns']):
        if pn['pore.all'][pn['throat.conns'][j,1]]:
            connections_all[pn['throat.conns'][j,0]] += 1
            connections_all[pn['throat.conns'][j,1]] += 1
            

    
    nnodes = pnx_int.number_of_nodes()
    
    degrees = np.asarray([val for (node, val) in pnx_int.degree()])
    av_degree = np.sum(degrees)/nx.classes.function.number_of_nodes(pnx_int)
    
    print('Average pore coordination number')
    print(av_degree)
    
    #This way isolated pores can be excluded 
    av_degree_corr = np.sum(degrees)/(nx.classes.function.number_of_nodes(pnx_int)-np.sum(degrees==0))
    print('Corrected average pore coordination number')
    print(av_degree_corr)
    
    degreehist = np.histogram(degrees, bins=np.arange(np.min(degrees),np.max(degrees)+2),density=False)
    
    
    fig1 = plt.figure(num=75, figsize=(5,5))
    plt.clf()
    ax2 = fig1.add_subplot(1,1,1)
    ax2.bar(degreehist[1][0:-1], degreehist[0], width = np.diff(degreehist[1]),align='edge')
    #ax2.set_yscale('log')
    ax2.set_xlabel('Pore coordination number')
    ax2.set_ylabel('Number of pores')
    
    fig1 = plt.figure(num=85, figsize=(5,5))
    plt.clf()
    ax1 = fig1.add_subplot(1,1,1)
    ax1.plot(pn['pore.coords'][pn['pore.internal'],0], degrees, 'b.')
    ax1.set_xlabel('Pore vertical coordinate')
    ax1.set_ylabel('Coordination number')
    
if clust_switch:

    #Clustering coefficient
    clust = np.asarray(list(nx.algorithms.cluster.clustering(pnx_int).values()))
    
    #Network average clustering coefficient
    clust_av = np.mean(clust)
    
    print('Network average clustering coefficient, single network')
    print(clust_av)
    
    clusthist = np.histogram(clust, bins=10,density=False)

    fig1 = plt.figure(num=76, figsize=(5,5))
    plt.clf()
    ax2 = fig1.add_subplot(1,1,1)
    ax2.bar(clusthist[1][0:-1], clusthist[0], width = np.diff(clusthist[1]),align='edge')
    ax2.set_yscale('log')
    ax2.set_xlabel('Pore clustering coefficient')
    ax2.set_ylabel('Number of pores')


    fig1 = plt.figure(num=6)
    fig1.set_size_inches(9,5)
    plt.clf()
    cnvscc_hist = []
    for i in range(2,11):
        try:
            ax1 = fig1.add_subplot(3,3,i-1)
            qwer = np.unique(clust[degrees==i])
            qwert = np.asarray(qwer[-1] + np.diff(qwer)[0])
            np.concatenate([qwer,np.array([qwert])])
            dummy1 =  ax1.hist(clust[degrees==i], bins=np.concatenate([qwer,np.array([qwert])]), density=False, align='left')
            cnvscc_hist.append(dummy1)
            ax1.set_ylim(0.9, 1.1*np.max(dummy1[0]))
            ax1.set_yscale('log')
            if i-1 == 1 or i-1 == 6 or i-1 == 11:
                ax1.set_ylabel('Pore clustering coefficient')
            if i-1 > 10:
                ax1.set_xlabel('Pore coordination number')
            ax1.set_title('CN = '+str(i))
        except Exception:
            pass
        
    fig1 = plt.figure(num=86, figsize=(5,5))
    plt.clf()
    ax1 = fig1.add_subplot(1,1,1)
    ax1.plot(pn['pore.coords'][pn['pore.internal'],0], clust, 'b.')
    ax1.set_xlabel('Pore vertical coordinate')
    ax1.set_ylabel('Clustering coefficient')
    
    '''
    # Network transitivity (global clustering coefficient)
    transitiv = nx.algorithms.cluster.transitivity(pnx_int)
    
    print('Network transitivity')
    print(transitiv)
    '''
    
if centrality_switch:    
    
    print('Calculating closeness centrality')
    
    ccentr = []
    for i in range(0, pnx_int.number_of_nodes(), space):
        ccentr.append(nx.algorithms.centrality.closeness_centrality(pnx_int, u=i, distance='conduitlength'))
        if i%1000==0:
            print(i)
    
    #Closeness centrality
    ccentr = np.asarray(ccentr)
    ccentrhist = np.histogram(ccentr, bins=10,density=False)
    
    try:
        testing = pn['pore.equivalent_diameter'][0]
        del testing
    except Exception:
        pn['pore.equivalent_diameter']=pn['pore.diameter']
        
    
    fig1 = plt.figure(num=777)
    fig1.set_size_inches(7,7)
    plt.clf()
    ax1 = fig1.add_subplot(2,2,1)
    ax1.hist(ccentr,bins=20)
    ax1.set_xlabel('Pore closeness centrality')
    ax1.set_ylabel('Number of pores')
    if space == 1: 
        ax2 = fig1.add_subplot(2,2,2)
        ax2.plot(pn['pore.coords'][pn['pore.internal'],0], ccentr, 'b.')
        ax2.set_xlabel('Vertical coordinate')
        ax2.set_ylabel('Closeness centrality')
        ax2.set_ylim([0,1.05*np.max(ccentr)])
        
        ax3 = fig1.add_subplot(2,2,3)
        ax3.scatter(pn['pore.equivalent_diameter'][pn['pore.internal']], ccentr, c='b')
        ax3.set_xlabel('Pore equivalent diameter')
        ax3.set_ylabel('Closeness centrality')
        ax3.set_ylim([0,1.05*np.max(ccentr)])
        ax3.set_xlim([-0.0005,0.005])
        
        ax4 = fig1.add_subplot(2,2,4)
        ax4.scatter(pn['pore.coords'][pn['pore.internal'],0],
                    pn['pore.equivalent_diameter'][pn['pore.internal']], c='b')
        ax4.set_xlabel('Vertical coordinate')
        ax4.set_ylabel('Pore equivalent diameter')
        ax4.set_xlim([-0.001,0.031])
        ax4.set_ylim([-0.0005,0.005])
    plt.tight_layout()
    
    print('Average closeness centrality')
    print(np.mean(ccentr))
    
    #Betweenness centrality
    bc_all = np.asarray(list(nx.algorithms.centrality.betweenness_centrality(pnx_int, k=np.int(pnx_int.number_of_nodes()/space), weight='conduitlength').values()))
    
    print('Average betweenness centrality for the whole network')
    print(np.mean(bc_all))        
    
    bcentrhist = np.histogram(bc_all, bins=20,density=False)
    
    fig78 = plt.figure(num=788)
    fig78.set_size_inches(7,7)
    plt.clf()
    ax1 = fig78.add_subplot(2,2,1)
    ax1.bar(bcentrhist[1][0:-1], bcentrhist[0], width = np.diff(bcentrhist[1]),align='edge')
    ax1.set_yscale('log')
    ax1.set_xlabel('Pore betweenness centrality')
    ax1.set_ylabel('Number of pores')
    if space == 1: 
        ax2 = fig78.add_subplot(2,2,2)
        ax2.plot(pn['pore.coords'][pn['pore.internal'],0], bc_all, 'b.')
        ax2.set_xlabel('Vertical coordinate')
        ax2.set_ylabel('Betweenness centrality')
        ax2.set_ylim([0,1.05*np.max(bc_all)])
        
        ax3 = fig78.add_subplot(2,2,3)
        ax3.scatter(pn['pore.equivalent_diameter'][pn['pore.internal']], bc_all, c='b')
        ax3.set_xlabel('Pore equivalent diameter')
        ax3.set_ylabel('Betweenness centrality')
        ax3.set_ylim([0,1.05*np.max(bc_all)])
        #ax3.set_xlim([-0.0005,0.005])
    plt.tight_layout()
    
    fig1 = plt.figure(num=88, figsize=(11,4))
    plt.clf()
    ax1 = fig1.add_subplot(1,3,1)
    ax1.plot(pn['pore.coords'][pn['pore.internal'],0], bc_all, 'b.')
    ax2 = fig1.add_subplot(1,3,2)
    ax2.plot(pn['pore.coords'][pn['pore.internal'],1], bc_all, 'b.')
    ax3 = fig1.add_subplot(1,3,3)
    ax3.plot(pn['pore.coords'][pn['pore.internal'],2], bc_all, 'b.')    
    plt.suptitle('Betweenness centrality vs. directions')