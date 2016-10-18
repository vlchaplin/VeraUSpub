# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 09:15:47 2016

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
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

import pandas

hfile='/Users\Vandiver\Data\Verasonics\sonalleve_20160709\pcdmaps20160709_BBfilt.h5'
f = h5py.File(hfile,'r')

extensions = ('single','multi_1','multi_2')
probes = ("probe1","probe2")
powers = [5,10,20,40,60,80]
#powers=[10,60]
#%% same image acquistion was used for all the images
dsetAttr1=f['multi_2/multi_2_10W/probe1/map'].attrs
dsetAttr2=f['multi_2/multi_2_10W/probe2/map'].attrs
extent1=dsetAttr1['extent']
extent2=dsetAttr2['extent']

(dz,dx,dt,c)=list(map( lambda x:dsetAttr1[x], ['dz','dx','dt','c']))

xaxis1 = np.arange(extent1[0],extent1[1]+0e-7,dx*1e2)*1e-2
zaxis1 = np.arange(extent1[3],extent1[2]+1e-7,dz*1e2)*1e-2 
xaxis2 = np.arange(extent2[0],extent2[1]+0e-7,dx*1e2)*1e-2 
zaxis2 = np.arange(extent2[3],extent2[2]+1e-7,dz*1e2)*1e-2 

g1 = np.meshgrid(zaxis1,xaxis1,indexing='ij')
g2 = np.meshgrid(zaxis2,xaxis2,indexing='ij')

dist1 = np.sqrt( g1[0]**2 + g1[1]**2)
dist2 = np.sqrt( g2[0]**2 + g2[1]**2)

maxMom2_1 = np.max( f['multi_2/multi_2_80W/probe1/mom2'] )
maxMom2_2 = np.max( f['multi_2/multi_2_80W/probe2/mom2'] )


(l,b,w,h)=(0.1, 0.1, 0.5, 0.8)

    
def plot_page(extension, label, pdf=None,xlabel='frame'):
    plt.figure(figsize=(10,8))
    plt.axes([l,b,w,h])
    
    (numf,numa) = extension['rms'].shape    
    noise = extension['BB'].value.reshape([numf,numa])
    plt.boxplot(noise.transpose() )
    plt.text(0.5,0.9, label, horizontalalignment='center',transform=plt.gca().transAxes,fontsize=16)
    plt.xlabel(xlabel,fontsize=14)
    plt.ylabel("Source strength",fontsize=14)
    
    ax2=plt.axes([w+1.5*l,b,1-w-2.5*l,h])        
    
    dz = extension['map'].attrs.get('dz')*1e2
    dx = extension['map'].attrs.get('dx')*1e2
    (nz,nx)=extension['map'].shape
    extent=[ -dx*nx/2,dx*nx/2,nz*dz,0]
    im=ax2.imshow(  extension['mom2'][20], extent=extent,cmap=image.cm.gray )
    ax2.xaxis.set_ticks([-2.0, 0.0, 2.0])
    plt.xlabel("cm",fontsize=14)
    cb=plt.colorbar(mappable=im,orientation='vertical')
    
    if pdf is not None:
        pdf.savefig()
        plt.close()    

def plot_probenoise(base_extension, label=None,acq_list=None, axs=None,acq_x=None,xlabel='acq',ylabel="Source strength"):

    if label==None:
        default_label=True
    else:
        default_label=False
    if axs==None:
        fig=plt.figure(figsize=(10,8))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        axs=(ax1,ax2)
    else:
        (ax1,ax2)=axs
    for pi in range(0,len(probes)):
        
        extension = base_extension[probes[pi]]
        (numf,numa) = extension['rms'].shape    
        noise = extension['BB'].value.reshape([numf,numa])
        if default_label:
            strlabel = extension.name
        else:
            strlabel=label + " "+probes[pi]
        
        axlabels = ["" for x in range(1,numf+1)]
        for i in range(5,len(axlabels)+1,5):
            axlabels[i-1]="%d" % (i)
        
            
        #axs[pi].boxplot(noise.transpose(), labels=axlabels )
        #noise[0,0]=noise[0,1] #the first acquisition usually has a large signal
        
        
        if acq_list is not None:
            if acq_x is None:
                acq_x = acq_list
            axs[pi].plot( acq_x, noise.flatten()[acq_list] )
        else:
            noiseflat=noise.flatten()
            if acq_x is None:
                axs[pi].plot( acq_x, noiseflat )
            else:
                axs[pi].plot( noiseflat )
            
            
        
        axs[pi].text(0.5,0.9, strlabel, horizontalalignment='center',transform=axs[pi].transAxes,fontsize=16)
        axs[pi].set_xlabel(xlabel,fontsize=16)
        axs[pi].set_ylabel(ylabel,fontsize=16)
        axs[pi].tick_params(labelsize=16)
    
    y1 = ax1.get_ylim()[1]
    y2 = ax2.get_ylim()[1]

    if y1 > y2:
        ax2.set_ylim(ax1.get_ylim())
    else:
        ax1.set_ylim(ax2.get_ylim())
    

plot_probenoise(f['single/single_80W'],acq_list=range(0,200))
#%%
#with PdfPages('/Users/Vandiver/Data/Verasonics/sonalleve_20160709/plots_20160709.pdf') as pdf:
#    plot_page(f['single/single_80W/probe1'],'single',pdf=pdf)
#    plot_page(f['multi_1/multi_1_80W/probe1'],'multi_1',pdf=pdf)
#    plot_page(f['multi_2/multi_2_80W/probe1'],'multi_2',pdf=pdf)

#%%

def make_extname (hifutype, power, probe): 
    return "%s/%s_%dW/%s"%(hifutype, hifutype,power, probe)


row_it = probes
col_it = extensions

plt_attr={}
plt_attr["single"]=dict(color=(0,0,1))
plt_attr["multi_1"]=dict(color=(0.9,0.1,0))
plt_attr["multi_2"]=dict(color=(0.5,0.5,0))

gs=gridspec.GridSpec(3,3,wspace=0.05,hspace=0.35)
gsp=gridspec.GridSpec(3,3,wspace=0.05,hspace=0.1)
vsPowerCurves = {}
vsPowerCurves['cav']= np.zeros([len(row_it),len(col_it),len(powers)])
vsPowerCurves['cav_img']= np.zeros([len(col_it),len(powers)])

distmats = (dist1,dist2)
#powi=3

num_img_frames=1

acq_list = np.arange(1,580,1,dtype=int)
acq_x = acq_list*50e-3

class NoneContext:
    def __init__(self):
        pass
    def __enter__(self):
        return None
    def __exit__(self, type, value, traceback):
        return
#with PdfPages('/Users/Vandiver/Data/Verasonics/sonalleve_20160709/maps_20160709_4frames.pdf') as pdf:    
#with NoneContext() as pdf:

#pdf = PdfPages('/Users/Vandiver/Data/Verasonics/sonalleve_20160709/maps_20160709_filtHarm.pdf')
pdf=None   
for powi in range(0,len(powers)):
    
    fig=plt.figure(figsize=(11,11))
    axList1 = []
    powval=powers[powi]        
    
    for n in range(len(col_it)):
        (ax1,ax2)=(fig.add_subplot(gsp[0,n]), fig.add_subplot(gsp[1,n]) )
    
        ext=make_extname( col_it[n], powval, "" )
        plot_probenoise(f[ext], label=col_it[n], xlabel="sec", acq_list=acq_list,acq_x=acq_x, axs=(ax1,ax2))
        axList1.append(ax1)
        axList1.append(ax2)
        if n==0:
            ymax=ax1.get_ylim()[1]
            
        ax1.axes.xaxis.set_ticks([])
        ax1.set_xlabel('')
        if n>0:
            if ax1.get_ylim()[1]>ymax:
                ymax=ax1.get_ylim()[1]
                
            
            ax1.axes.yaxis.set_ticks([])
            ax2.axes.yaxis.set_ticks([])
            
            ax1.set_ylabel('')
            ax2.set_ylabel('')
            
        
    fig.text(0.5,0.95, "%dW"%powval, horizontalalignment='center',fontsize=24,color='r',fontweight='bold')

    for ax in axList1:
        ax.set_ylim((0,ymax))

    if pdf is not None:        
        pdf.savefig(fig)
        plt.close(fig)

for powi in range(0,len(powers),1):
    fig=plt.figure(figsize=(11,11))
    axList = []
    powval=powers[powi]
    maxIm=-9999
    
    for n in range(len(col_it)):
        
        probeSumM1 =0
        probeSumM2 =0
            
        for m in range(len(row_it)):
            
            if m==0:
                extent=extent1
            else:
                extent=extent1                
            
            axList.append(fig.add_subplot(gs[m,n]))
            ext=make_extname( col_it[n], powval, row_it[m] )
            
            label=col_it[n]            
            (numf,numa) = f[ext]['rms'].shape              
            
                        
            
            imdata2=( np.mean( f[ext+'/mom2'][0:num_img_frames], axis=0)  )
            imdata1=( np.mean( f[ext+'/mom1'][0:num_img_frames], axis=0)  )
            imdata = (imdata2 - imdata1**2)
            
 
            #imdata = imdata1
            #imdata/=np.sum(imdata)
            
            #imdata = imdata1
            if m==0:
                imprev=imdata.copy()
            
                probeSumM1=imdata1[31:,0:-5]    
                probeSumM2=imdata2[31:,0:-5]
                probeCross=f[ext+'/mom1'][0:num_img_frames,31:,0:-5]
                (znmax,xnmax)=probeSumM1.shape
                
            if m==1:
                imdata1 =  np.flipud(np.fliplr(imdata1))
                imdata2 =  np.flipud(np.fliplr(imdata2))
                #imdata*=np.sum(imprev)/np.sum(imdata)
                mom1frames=f[ext+'/mom1'][0:num_img_frames,0:znmax,0:xnmax]                
                
                probeCross*=mom1frames[:,::-1,::-1]
                probeSumM1+=imdata1[0:znmax,0:xnmax]        
                probeSumM2+=imdata2[0:znmax,0:xnmax]
                
                sumimg = probeSumM2 - probeSumM1**2             
                #sumimg = probeSumM2 -  2*np.mean(probeCross,axis=0)                
                
                #sumimg = np.zeros_like(imdata)
                #subsetim1 = imprev[31:,0:-5]
                #subsetim1=imprev
                #sumimg[0:subsetim1.shape[0],0:subsetim1.shape[1]] = subsetim1
                                    
                #sumimg += (imdata+imprev)**2/4
                #imdata.
            
            imax=0.8*np.max(imdata) 
            if maxIm<imax:
                maxIm=imax
                
            
            #vsPowerCurves['cav'][m,n,powi]=np.mean(f[ext+'/rms'][0:6,:])
            vsPowerCurves['cav'][m,n,powi]=np.mean(f[ext+'/BB'][0:(num_img_frames*numa)])
            axList[-1].imshow((imdata))
            plt.text(0.5,0.9, label, horizontalalignment='center',transform=plt.gca().transAxes,fontsize=16, color='r',fontweight='bold')
            #axList[-1].imshow( np.log10( np.abs(hilbert( f[ext+'/mom1'][0]) ) ) )

            #endfor    

        #maxIm=0.7*np.max(sumimg) 
        #if maxIm<imax:
        #    maxIm=imax
            
        vsPowerCurves['cav_img'][n,powi]=np.sum(sumimg)
        axList.append(fig.add_subplot(gs[2,n]))
        axList[-1].imshow( sumimg)


    fig.text(0.5,0.95, "%dW"%powval, horizontalalignment='center',fontsize=24,color='r',fontweight='bold')
    
    for ax in axList:
        im=ax.get_images()[0]
        im.set_clim(vmin=0,vmax=(maxIm) )    
        im.set_cmap(image.cm.jet)
    
    if pdf is not None:        
        pdf.savefig(fig)
        plt.close(fig)

#%%
pdir='/Users/Vandiver/Data/Verasonics/sonalleve_20160709/'
pfile='/Users/Vandiver/Data/Verasonics/sonalleve_20160709/HifuScanParameters.csv'
scanHifuParams=pandas.read_csv(pfile)
    
scanAnalysisDT=pandas.read_csv(pdir+'scans_batch_DT_unwr.csv')

merged=scanAnalysisDT.merge(scanHifuParams)

#%%

gsPow=gridspec.GridSpec(3,2,wspace=0.25,hspace=0.35)

#with PdfPages('/Users/Vandiver/Data/Verasonics/sonalleve_20160709/plots_20160709_4frames.pdf') as pdf:
#with NoneContext() as pdf:
fig2=plt.figure(figsize=(10,10))

cases=('single','multi_1','multi_2')

for n in range(2):
    ax=fig2.add_subplot(gsPow[0,n])
    
    ax.plot(powers, vsPowerCurves['cav'][n,0,:],'-o', **plt_attr[extensions[0]])
    ax.plot(powers, vsPowerCurves['cav'][n,1,:],'-o', **plt_attr[extensions[1]])
    ax.plot(powers, vsPowerCurves['cav'][n,2,:],'-o', **plt_attr[extensions[2]])
    ax.set_xlabel('Input Power (W)',fontsize=18)
    ax.set_ylabel('Cav. (BB sum)',fontsize=18)
    ax.tick_params(labelsize=16)
        #plt.text(0.5,0.9, powval, horizontalalignment='center',transform=plt.gca().transAxes,fontsize=16)
    
(axCavim, axT1)=(fig2.add_subplot(gsPow[1,0]) , fig2.add_subplot(gsPow[1,1]))
(axP1,axP2)=(fig2.add_subplot(gsPow[2,0]), fig2.add_subplot(gsPow[2,1]))

for n in range(len(cases)):
    axCavim.plot( powers, vsPowerCurves['cav_img'][n,:], 'o-',**plt_attr[cases[n]])

axCavim.set_xlabel('Input Power (W)',fontsize=18)
axCavim.set_ylabel('Cav. (Img+probe sum)',fontsize=18)
axCavim.tick_params(labelsize=16)


for n in range(len(cases)):
    subset=merged.query('tag=="'+cases[n]+'"').sort_values(by='power')


    #Tvalues = np.array(list( map ( lambda crv: max(map(float,crv[1:-1].split())), subset.maxT)  ))
    Tvalues = subset.finalT
    #Tvalue = np.max(subset.maxT)
    axT1.plot( subset.power, Tvalues, 'o-',**plt_attr[cases[n]])

    #axT2.plot( subset.power, subset.finalT, 'o-',**plt_attr[cases[n]])
    axP1.plot( vsPowerCurves['cav'][0,n,:] + vsPowerCurves['cav'][1,n,:], subset.finalT, 'o-',**plt_attr[cases[n]])
    axP2.plot( vsPowerCurves['cav_img'][n,:], Tvalues, 'o-',**plt_attr[cases[n]])
    #lt.plot( vsPowerCurves['cav'][0,n,:], subset.finalT,'o-')
    
axP1.set_xlabel('Cav. (BB sum)',fontsize=18)
axP2.set_xlabel('Cav. (Img+probe sum)',fontsize=18)
for ax in (axP1,axP2):    
        ax.set_ylabel('$\Delta$T ($^o$C)',fontsize=18)
        ax.tick_params(labelsize=16)

axT1.tick_params(labelsize=16)
axT1.set_xlabel('Input Power (W)',fontsize=18)
axT1.set_ylabel('$\Delta$T ($^o$C)',fontsize=18)

if pdf is not None:
    pdf.savefig(fig2)
    pdf.close()