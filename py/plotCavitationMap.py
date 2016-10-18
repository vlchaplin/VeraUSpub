# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 15:21:20 2016

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

import argparse

sys.path.append('C:\\Users\\vchaplin\\Documents\\HiFU\\code\\myPy')  
import stringmisc

#sys.path.append('C:\\Users\\vchaplin\\Documents\\HiFU\\code\\AblationSims')

#check if in interactive mode

try:
    sys.ps1
    interactive=True
except AttributeError:
    interactive=False


def getParamFile(file):
    """
    """
    datapath = os.path.dirname(file)
    (dataname,ext) = os.path.basename(file).split(".")

    return (datapath +"/" + dataname +"_params.mat",ext)

def load_cavitation_params(params_file):
    """
    """
    pardict=io.loadmat(params_file,squeeze_me=True,struct_as_record=False)
    params = pardict['params']
    return params
    
def load_veradata_bin(datafile,params=None):
    """
    """
    if params is None:
        params = load_cavitation_params(params_file)
    
    ns=params.numRcvSamples
    Nchan = 128
    
    rf_data = np.zeros([ns, Nchan, params.numframes*params.numacq])     
    
    with open(datafile, 'rb') as fid:
        #chandat = np.fromfile(fid, np.int16).reshape([ns, params.numacq, Nchan, params.numframes])
        #dimensions are reversed from matlab since python is row major
        chandat = np.fromfile(fid, np.int16).reshape([params.numframes,Nchan,params.numacq,ns])
        
    for f in range(0, params.numframes):
        for a in range(0, params.numacq):
            #rf_data[:,:,f*params.numacq + a] = chandat[:,a,:,f]
            rf_data[:,:,f*params.numacq + a] = chandat[f,:,a,:].transpose()
            
    return (rf_data)
    
def load_veradata_mat(datafile,params=None,removeNsamples=0):
    """
    """
    if params is None:
        params = load_cavitation_params(params_file)
    
    ns=params.numRcvSamples
    Nchan = 128
    
    rf_data = np.zeros([ns, Nchan, params.numframes*params.numacq])  
        
    matdict=io.loadmat(datafile,squeeze_me=False,struct_as_record=False)
    
    for f in range(0, params.numframes):
        for a in range(0, params.numacq):
            rf_data[:,:,f*params.numacq + a] = matdict['tmpDat'][(a*ns):(a*ns+ns),:,f]

    rf_data = rf_data[removeNsamples: ,: ,:]
    ns=rf_data.shape[0]
    
    return (rf_data)
    
def getBkgFile(datafile):
    path=os.path.dirname(datafile)
    datafile_basename=os.path.basename(datafile)
    candidates = [os.path.basename(fn) for fn in glob.glob(path+"/"+"back*.mat" ) if not re.match(".*params.*",fn)]
    bkgdfile=''
    if (len(candidates) > 0) :
        lcsmatches=[stringmisc.longest_common_substring(fn,datafile_basename) for fn in candidates]
        matchlengths = list(map(len, lcsmatches))
        idx=  np.argmax(matchlengths) 
    
        bkgdfile = candidates[idx]
    
    return path+"/"+bkgdfile
    
def rescaleIm(inArray,newmin=0.0,newmax=1.0,vmin=None,vmax=None,trunc=False):
    """
    Rescale input array values to the range [newmin, newmax]
    Keywords & defaults: newmin=0.0, newmin=1.0
        vmin=  The value in input array that will be mapped to newmin. Default=min(inArray)
        vmax=  The value in input array that will be mapped to newmax. Default=max(inArray)
        trunc=  If true, truncate the output range. Values less than newmin will be set to newmin. Likewise for newmax.
    """
    oldmin = np.min(inArray)
    oldmax = np.max(inArray)
    
    if vmin is None:
        vmin = oldmin
    if vmax is None:
        vmax = oldmax
        
    newim = (inArray - vmin)*(newmax - newmin) / (vmax - vmin) + newmin
    
    if trunc:
        newim[newim < newmin] = newmin
        newim[newim > newmax] = newmax
        
    return newim
    
#%%
    
delayVsFrame = None
#h5outfile="/Users/vchaplin/Data/Davis Visit Cavitation Mapping/May 12/May12_ImageFormed.h5"

#fileList = glob.glob("/Users/vchaplin/Data/Davis Visit Cavitation Mapping/May 12/*.mat")

#datafile="/Users/vchaplin/Data/Davis Visit Cavitation Mapping/May 12/bubble_center_6V_1.mat"

#fileList = [fn for fn in fileList if not re.match(".*params.*",fn)]
path="/Users/vchaplin/Data/DavisData/May_13/first_off_axis/"
path="/Users/vchaplin/Data/Davisdata/May_13/x_axis_scans/"
#files=[path+"bubbles_center_6V.mat", path+"bubbles_left5mm_6V.mat", path+"bubbles_right5mm_6V.mat"]

files = ["/Users/vchaplin/Data/Verasonics/sonalleve_0607/3foc-agar-1000av.bin"]
files = [path+"bubble_linescan_6V_5c.mat"]

datafile=files[0]
datafile_basename=os.path.basename(datafile)

 

path=os.path.dirname(datafile)
params_file,ext = getParamFile(datafile)

params = load_cavitation_params(params_file)

bkg_data=None

        
numf = params.numframes
numa = params.numacq
    
Nchan=params.numRcvChannels
ns = params.numRcvSamples

#%% construct delay array on first file only
    
if 1:
    c=1480; 
    Trans={}
    Trans['spacing']=params.pitch;
    Trans['frequency']=params.fs/4*1e-6;
    #wavelength = c/(Trans['frequency']*1e6);
    
    Fs = 4*Trans['frequency']*1e6;
    dt = 1/Fs;
    
    dx = 2e-4
    dz = 2e-4
    
    ducer_width = params.pitch*Nchan
    
    Nx = round(ducer_width/dx)
    #Nz = round(ns*c*dt/dz);
    #Nz = round(6e-2/dz)
    
    xpnts = np.linspace(-0.5,0.5,Nx)*Nx*dx
    #zpnts = np.linspace(0,1,Nz)*Nz*dz
    
    zpnts = np.arange(0.03,0.06,dz)
    Nz=len(zpnts)    
    
    #sensor positions
    ux = (np.linspace(-0.5,0.5,Nchan)*Nchan + 0.5)*params.pitch     
            
    #image and resampling grids
    ndZ,ndX,ndux = np.meshgrid(zpnts,xpnts,ux, indexing='ij')
    
    distances = np.sqrt( (ndX-ndux)**2 + ndZ**2 )
    delayinds = np.round( distances / (c*dt)).astype(int)
    inbounds3 = (delayinds < ns)
    ii=np.nonzero(inbounds3)

#distflat = distances[inbounds3]

#%%image formation
delayed=np.zeros([Nz,Nx,Nchan])

superFrames = np.arange(0,numf)
#○superFrames=[19]
acqs = np.arange(1,numa)
numSF = len(superFrames)

acqs=[1,2,3]

dataOut=[]
L2sqOut=[]
for datafile in files:
    print('')
    print('Processing: %s'% datafile)   
    #removeNsamples = 829
    #rf_data = load_veradata_mat(datafile, params=params, removeNsamples=removeNsamples)
    #(ns,Nchan,totframes)=rf_data.shape    
    
    
#    if ext=="mat":
#        #if mat file
#        removeNsamples = 829
#        
#        rf_data = load_veradata_mat(datafile, params=params, removeNsamples=removeNsamples)
#        (ns,Nchan,totframes)=rf_data.shape
#        
#        # look for background file in dirrectory
#        bkgdfile = getBkgFile(datafile)
#        if len(bkgdfile)>0:
#            print ("BACKGROUND file: ", bkgdfile, flush=True)
#            bkg_data = load_veradata_mat(bkgdfile,removeNsamples=removeNsamples)        
#    else:
#        #if binary file
#        rf_data = load_veradata_bin(datafile, params=params)
#        (ns,Nchan,totframes)=rf_data.shape    
    
    
    
    if re.match(".*linescan.*", os.path.basename(datafile)):
        delayVsFrame = np.abs(np.round(72.0/5*np.arange(-5.0, 5.001,0.5)).astype(int))
    elif re.match(".*_center_.*", os.path.basename(datafile)):
        delayVsFrame = np.zeros(numf,dtype=int)
    elif re.match(".*left5mm.*", os.path.basename(datafile)):
        delayVsFrame = np.zeros(numf,dtype=int) + 72
    elif re.match(".*right5mm.*", os.path.basename(datafile)):
        delayVsFrame = np.zeros(numf,dtype=int) + 72
    elif delayVsFrame is not None:
        delayVsFrame = np.zeros(numf,dtype=int)    
    
    
    imgf = np.zeros([numSF,Nz,Nx])
    chsumf = np.zeros_like(imgf)
    #imga = np.zeros([numSF,Nz,Nx])
    imgBkg = np.zeros_like(imgf)
    rmsPerImage = np.zeros([numSF,params.numacq])
    
    for sfi in range(0,numSF):
        sf = superFrames[sfi]
        
        if delayVsFrame is not None:
            
            delayinds = np.round( distances / (c*dt)).astype(int) +delayVsFrame[sf]
            inbounds3 = (delayinds < ns)
            ii=np.nonzero(inbounds3)
        
        
        for a in acqs:
            #delayed[:]=0
            delayed[ii[0],ii[1],ii[2]] = bkg_data[delayinds[inbounds3],ii[2],sf*params.numacq + a]
    
            summed_passivemap = np.sum(distances*delayed,axis=2)
            L2sq=np.sum( (distances*delayed)**2,axis=2)
            rms=summed_passivemap**2
            
            
            rmsPerImage[sfi,a] = np.sum(rms,axis=(0,1))
            imgf[sfi]+=rms
            chsumf[sfi]+=L2sq
            
            print("\r%d/%d"%(sf*params.numacq + a, superFrames[-1]*params.numacq + acqs[-1]),end='', flush=True)
        
#        if bkg_data is not None:
#            for a in acqs:        
#                delayed[ii[0],ii[1],ii[2]] = bkg_data[delayinds[inbounds3],ii[2],sf*params.numacq + a]
#                imgBkg[sfi]+=np.sum(distances*delayed,axis=2)**2
#            imgBkg[sfi]/=len(acqs)
            
        imgf[sfi]/=len(acqs)
        chsumf[sfi]/=len(acqs)
    dataOut.append(imgf)
    L2sqOut.append(chsumf)
    #img=np.sum(imgf,axis=0)

#%% bmode ish

delayedBinds = np.round( distances / (2*1480*dt)).astype(int) -817
inbb = (delayedBinds < ns)
iibb=np.nonzero(inbb)
delayed[iibb[0],iibb[1],iibb[2]] = rf_data[delayedBinds[iibb],iibb[2],0*params.numacq + 0]

bmode = np.sum(delayed,axis=2)
plt.imshow(np.abs(hilbert(bmode)),extent=extent)

#%%
#img = np.log10( np.abs(hilbert(img)) )

#sfi=0
#pcImg = imgf[sfi]
#bakImg=imgBkg[sfi]



#pcmin=np.min(pcImg)
#pcmax=np.max(pcImg)

extent = [ 100*xpnts[0], 100*xpnts[-1], 100*zpnts[-1], 100*zpnts[0] ]

fig=plt.figure(figsize=(8,10))
ax=plt.gca()
imB=plt.imshow(np.abs(hilbert(bmode)), cmap=image.cm.gray,extent=extent, interpolation='none')
ax.set_xlabel('cm',fontsize=14)
ax.set_ylabel('cm',fontsize=14)
ax.tick_params(labelsize=14)  

#vmin=1e3
#vmax=1e6

(vmin,vmax)=(None,None)

(vmin,vmax)=(10.97,2183237.72030)

for di in range(0,len(dataOut)):
    pcImg = np.mean(dataOut[di] ,axis=0)
    #pcImg = dataOut[di][0] - 3*L2sqOut[di][0]
    #pcImg = np.mean(dataOut[di] - L2sqOut[di] ,axis=0)
    if vmin is None:
        vmin=np.min(pcImg)
    if vmax is None:
        vmax=np.max(pcImg)    
    
    newmin=0
    newmax=800
    scaledMatrix=rescaleIm(pcImg,vmin=vmin,vmax=vmax, newmin=newmin,newmax=newmax, trunc=True)
    vmin=1
    vmax=newmax
    scaledMatrix=rescaleIm(scaledMatrix, vmin=vmin, vmax=vmax, trunc=False)

    rgbaIm1 = image.cm.jet(scaledMatrix )
    mask1 = scaledMatrix < 0.0
    
    rgbaIm1[mask1,3]=0
    
    im=plt.imshow(rgbaIm1,extent=extent,interpolation='none')
    
cb = plt.colorbar(mappable=im)

tickp = np.arange(0,1.01,.25)
tickstr = list(map(lambda s:"%0.e" % s,tickp*(vmax-vmin) + vmin))
cb.set_ticks(tickp)
cb.set_ticklabels(tickstr)