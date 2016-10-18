# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 13:40:03 2016

@author: vlchaplin@gmail.com
"""

#import h5py
#import sys
import numpy as np
#from scipy import io
#from scipy.signal import hilbert as hilbert
from math import *;
#
#import re
#import glob
#import os
import matplotlib.pyplot as plt
from matplotlib import image
import matplotlib.gridspec as gridspec
#from matplotlib.backends.backend_pdf import PdfPages

from scipy.interpolate import interp1d

import pandas


plt_attr={}
plt_attr["single"]=dict(color=(0,0,1))
plt_attr["multi_1"]=dict(color=(0.9,0.1,0))
plt_attr["multi_2"]=dict(color=(0.5,0.5,0))

def grnumaxis(nslc):
    sqr=np.sqrt(nslc)
    if np.abs(np.mod( sqr, 1.0))<1e-16:
        #is perfect square
        nnc=int(sqr)
    else:
        nnc=floor(sqr)+1
    
    nnr = ceil(nslc/nnc)

    return (nnr,nnc)
def grax_ij(slidx, nnr,nnc):
    j_ax = np.mod(slidx,nnc)
    i_ax = floor(slidx/nnc)
    return (i_ax,j_ax)


#%%


mrisets = [ 
    dict(use=1, path='/Users/Vandiver/Data/Verasonics/sonalleve_20160709/',
            pars='HifuScanParameters.csv', temp='scans_batch_DT_unwr_20160709.csv'),
    dict(use=1, path='/Users/Vandiver/Data/sonalleve/HifuCav20160810/',
            pars='HifuScanParams.csv', temp='scans_batch_DT_unwr_20160810.csv')]
            
tempfile='/Users/Vandiver/Data/Verasonics/MulitfocusCav/cav_scans_batch_DT_unwr2.csv'

merged_frames=[]
for mi in range(len(mrisets)):
    
    if not mrisets[mi]['use']:
        continue
    
    parsfile=mrisets[mi]['path']+mrisets[mi]['pars']
    scanHifuParams = pandas.read_csv(parsfile)
            
    #tempfile=mrisets[mi]['path']+mrisets[mi]['temp']
    scanAnalysisDT = pandas.read_csv(tempfile)
    
    
    merged = scanAnalysisDT.merge(scanHifuParams)
    merged_frames.append(merged)
    
merged = pandas.concat(merged_frames)
#%%

grouped = merged.query("material=='phantom'").groupby(['tag','material'])
#computed = grouped['finalT'].agg([np.min, np.max])

fig=plt.figure(figsize=(6,6))

for gidx,rowids in grouped:
    
    case=gidx[0]
    g=grouped.get_group(gidx).groupby('power')
    
    powers = list(g.groups.keys())  
    
    avgTr = g['finalT'].agg([np.min, np.max, np.mean])
    #numr = g['n5'].agg([np.min, np.max])
    #n10r = g['n10'].agg([np.min, np.max])
    #n15r = g['n15'].agg([np.min, np.max])

    plt.plot(avgTr.index, avgTr['mean'], 'o-',**plt_attr[case])
    plt.xlabel('Power (W)')
    plt.ylabel('ROI $<T>$ ($^o$C)')
    
fig=plt.figure(figsize=(6,6))

for gidx,rowids in grouped:
    
    case=gidx[0]
    g=grouped.get_group(gidx).groupby('power')
    
    TwSize=g['momR'].agg([np.mean])

    plt.plot(TwSize.index, TwSize['mean'], 'o-',**plt_attr[case])
    plt.xlabel('Power (W)')
    plt.ylabel('T-weighted moment (mm)')
    
#%%
    
#Tvalues = np.array(list( map ( lambda crv: max(map(float,crv[1:-1].split())), subset.maxT)  ))
        
def vectorstr2arr(panda_vec_string,dtype=np.float):
    return np.array(list(map(float,panda_vec_string[1:-1].split())), dtype=dtype)
    

powers=[5,10,20,40,60,80]
powAxLookup={'val':powers,'axi':list(range(len(powers)))}
powAxLookup = pandas.DataFrame(powAxLookup,index=powAxLookup['val'])

(nnr,nnc) = grnumaxis(len(powers))
gsPows=gridspec.GridSpec( nnr, nnc, wspace=0.5,hspace=0.25)

cases=["single","multi_1", "multi_2"]
gsCases=gridspec.GridSpec( 1, 3,wspace=0.5,hspace=0.3)


#fig=plt.figure(figsize=(9,6))
fig2=plt.figure(figsize=(14,10))

ax2s=[]
for pi in range(len(powers)):
    (ii,jj) = grax_ij(pi,nnr,nnc)
    ax2s.append( fig2.add_subplot(gsPows[ii, jj])  )



#dataset=merged.query("material=='phantom'").query("tag=='single'")

#cases=['single']
for ci in range(0,len(cases)):
    
    case=cases[ci]

    dataset=merged.query("material=='phantom'").query("tag=='%s'"%(case)).sort_values(by='power')

    if dataset.size==0:
        continue
    #row=dataset.iloc[0]
    
    pgrouped=dataset.groupby('power')
    powlist=sorted( pgrouped.groups.keys() )
    
    numTrialsPerPower=pgrouped.size()
    #concat the groups if more than trial exists
    tempvecsVsPower=[]
    for pi in powlist:
        subframe= pgrouped.get_group(pi).maxTdata
        
        temperatures = vectorstr2arr( subframe.iloc[0] )
        
        for gi in range(1,numTrialsPerPower[pi]):
            temperatures = np.concatenate( [temperatures, vectorstr2arr( subframe.iloc[gi] )] )
        
        tempvecsVsPower.append(temperatures)

#    ax=fig.add_subplot(gsCases[0,ci])
#
#    vplot=ax.violinplot(tempvecsVsPower,vert=False)
#    
#    for vpc in vplot['bodies']:
#        vpc.set_color(plt_attr[case]['color'])
#        
#    ax.set_xlim([0,25])
    
    for pi in range(len(powlist)):
        p=powlist[pi]
        try:
            pax=powAxLookup.loc[p].axi
        except KeyError:
            continue
        
        
        (hy,hx) = np.histogram(tempvecsVsPower[pi],bins=30,range=[-5,35],normed=True)        
        
        hxc = (hx[0:-1] + hx[1:] )/2
        hinterp = interp1d( hxc, hy, kind='cubic', fill_value=0,bounds_error=False)
        
        xip = np.arange(hx[0], hx[-1], .1)
        kwarg=dict(linewidth=2.0)
        
        (xx,yy) = (xip, hinterp(xip))
        yy[yy<0]=0
        
        ax2s[pax].fill_between(xx,yy,alpha=0.2, color=plt_attr[case]['color'], step=None,**kwarg) 
        ax2s[pax].plot(hxc,hy,'.',color=plt_attr[case]['color'])
        ax2s[pax].tick_params(labelsize=14)
        ax2s[pax].set_xlabel('$\Delta T$ ($^o$C)',fontsize=16)
        ax2s[pax].set_ylabel('Fraction of ROI',fontsize=16)
        ax2s[pax].text(0.9,0.9,'%dW'%p, transform=ax2s[pax].transAxes, horizontalalignment='right', fontsize=18,color='k')

        
#        vplot=ax2s[pax].violinplot(tempvecsVsPower[pi] ,positions=[ci],vert=False,widths=0.8)
#        for vpc in vplot['bodies']:
#            vpc.set_color(plt_attr[case]['color']) 
#        #ax2s[pax].set_xlim([0,25])
#        ax2s[pax].yaxis.set_ticks([])
    

#%%
#
#powers=[5,10,20,40,60,80]
#
#gsPow=gridspec.GridSpec( 1, 3,wspace=0.05,hspace=0.25)
#fig=plt.figure(figsize=(9,6))
#
##vs power
#for pi in range(len(powers)):
#    power=powers[pi]
#
#    dataset=merged.query("material=='phantom'").query("tag=='%s'"%(case)).sort(columns='power')
