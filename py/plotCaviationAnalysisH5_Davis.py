# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 21:45:49 2016

@author: vlchaplin@gmail.com
"""

import h5py
import sys
import numpy as np
from scipy import io
from scipy.signal import hilbert as hilbert
from math import *;

import re
import glob
import os
import matplotlib.pyplot as plt
from matplotlib import image

#hfile='C:/Users/Vandiver/Data/DavisData/May12_ImageFormed.h5'
#
#
#f = h5py.File(hfile,'r')
#
#folders = list(f.keys())
#
#datasets = list(f['May_12'].keys())
#rmsVsFrame = f['May_12']['bubble_center_6V_2']['rms'].value
#
#f.close()
#
#
#plt.imshow(rmsVsFrame,interpolation='None')


#%%
hfile='C:/Users/Vandiver/Data/DavisData/May13_BgSubImageFormed.h5'
#hfile='C:/Users/Vandiver/Data/DavisData/May13_ImageFormed.h5'

f = h5py.File(hfile,'r')

folders = list(f.keys())

xscandatasets = list(f['x_axis_scans'].keys())
yscandatasets = list(f['y_axis_scans'].keys())
#rmsVsFrame = f['May_13']['bubble_center_6V_2']['rms'].value

#f.close()

xsetA = [fn for fn in xscandatasets if re.match('bu.*linescan_\dV', fn)]
xsetB = [fn for fn in xscandatasets if re.match('bu.*linescan_point_source', fn)]

xsetA=sorted(xsetA)
xsetB=sorted(xsetB)


ysetA = [fn for fn in yscandatasets if re.match('bu.*linescan_\dV', fn)]
ysetB = [fn for fn in yscandatasets if re.match('bu.*linescan_point_source', fn)]

ysetA=sorted(ysetA)
ysetB=sorted(ysetB)

xbkg = [fn for fn in xscandatasets if re.match('back.*linescan_', fn)]
ybkg = [fn for fn in yscandatasets if re.match('back.*linescan_', fn)]

frame_pos = np.arange(-5.0, 5.001, 0.5)
axlabels = ["%0.1f" % x for x in frame_pos]

for i in range(1,len(axlabels),2):
    axlabels[i]=''

#%%
(l,b,w,h)=(0.1, 0.1, 0.5, 0.8)
def plot_page(extension, label, pdf=None,positions=frame_pos,labels=axlabels,xlabel='mm'):
    plt.figure(figsize=(10,6))
    plt.axes([l,b,w,h])
    plt.boxplot(extension['rms'].value.transpose(), positions=positions, labels=labels, widths= 0.25)
    plt.text(0.5,0.9, label, horizontalalignment='center',transform=plt.gca().transAxes,fontsize=16)
    plt.xlabel(xlabel,fontsize=14)
    plt.ylabel("Source strength",fontsize=14)
    
    ax2=plt.axes([w+1.5*l,b,1-w-2.5*l,h])        
    
    dz = extension['map'].attrs.get('dz')*1e2
    dx = extension['map'].attrs.get('dx')*1e2
    (nz,nx)=extension['map'].shape
    extent=[ -dx*nx/2,dx*nx/2,nz*dz,0]
    im=ax2.imshow(  extension['map'].value, extent=extent,cmap=image.cm.gray )
    ax2.xaxis.set_ticks([-2.0, 0.0, 2.0])
    plt.xlabel("cm",fontsize=14)
    cb=plt.colorbar(mappable=im,orientation='vertical')
    
    if pdf is not None:
        pdf.savefig()
        plt.close()
    
#%%
from matplotlib.backends.backend_pdf import PdfPages

with PdfPages('/Users/Vandiver/Data/DavisData/may13_xscans_scans.pdf') as pdf:
    ext=f['x_axis_scans']
    for dset in xsetA:        
        plot_page(ext[dset], dset + " (A)", pdf=pdf)
        
    
    for dset in xsetB:
        plot_page(ext[dset], dset + " (B)", pdf=pdf)

#%%
with PdfPages('/Users/Vandiver/Data/DavisData/may13_yscans_scans.pdf') as pdf:
    ext=f['y_axis_scans']
    for dset in ysetA:        
        plot_page(ext[dset], dset + " (A)", pdf=pdf)
        
    for dset in ysetB:
        plot_page(ext[dset], dset + " (B)", pdf=pdf) 


#%%
if len(xbkg)>0 or len(ybkg)>0:
    with PdfPages('/Users/Vandiver/Data/DavisData/may13_BKGscans.pdf') as pdf:
        ext=f['x_axis_scans']
        for dset in xbkg:
            plot_page(ext[dset], dset, pdf=pdf)
            
        ext=f['y_axis_scans']
        for dset in ybkg:
            plot_page(ext[dset], dset, pdf=pdf)
        
f.close()



#%%

hfile='C:/Users/Vandiver/Data/DavisData/May13_BgSubOffAxImageFormed.h5'
#hfile='C:/Users/Vandiver/Data/DavisData/May13_ImageFormed.h5'

f = h5py.File(hfile,'r')

dsets=list(f['first_off_axis'].keys())
with PdfPages('/Users/Vandiver/Data/DavisData/may13_OffAxis.pdf') as pdf:
    for dset in dsets:
        plot_page(f['first_off_axis'][dset], dset, pdf=pdf,positions=None,labels=None,xlabel='Frame')
        
f.close()


#%%
hfile='C:/Users/Vandiver/Data/DavisData/May13_Ytube.h5'
#hfile='C:/Users/Vandiver/Data/DavisData/May13_ImageFormed.h5'

f = h5py.File(hfile,'r')

dsets=list(f['tube_in_plane_correctly_aligned'].keys())
with PdfPages('/Users/Vandiver/Data/DavisData/may13_YscanTubeInPlane.pdf') as pdf:
    for dset in dsets:
        plot_page(f['tube_in_plane_correctly_aligned'][dset], dset, pdf=pdf,positions=None,labels=None,xlabel='Frame')
        
f.close()