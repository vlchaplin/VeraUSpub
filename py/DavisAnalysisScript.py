# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 11:26:01 2016

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
            
    return params
    
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
    
    
    
    
#%%
    
parser = argparse.ArgumentParser()
    
parser.add_argument('-f', metavar='matfiles', help="Name .mat files to process. Filenames containing string 'params' are ignored", type=str, nargs='+')
parser.add_argument("-o", metavar='output_HDF5',help="Name of output hdf5 file to create. It is opened in append mode", type=str )
parser.add_argument("-maxf", metavar='num',help="end number of super frames (default will use all)", type=int )
parser.add_argument("-maxa", metavar='num',help="end number of acquisitions to use for each frame (default will use all)", type=int )
parser.add_argument("--plot",help="make plots", action="store_true" )
parser.add_argument("--each",help="Recompute image space and delays for each file (default uses the delays from the first file)", action="store_true" )
parser.add_argument("--nobg",help="Do not look for a matching background_* file for background subtraction", action="store_true" )

args = parser.parse_args()    
    
if args.f :
    fileList = args.f
else:
    print("No files passed",flush=True)
    sys.exit()
if args.o :
    h5outfile = args.o
else:
    h5outfile= None
    
    
delayVsFrame = None
#h5outfile="/Users/vchaplin/Data/Davis Visit Cavitation Mapping/May 12/May12_ImageFormed.h5"

#fileList = glob.glob("/Users/vchaplin/Data/Davis Visit Cavitation Mapping/May 12/*.mat")

#datafile="/Users/vchaplin/Data/Davis Visit Cavitation Mapping/May 12/bubble_center_6V_1.mat"

fileList = [fn for fn in fileList if not re.match(".*params.*",fn)]

fnum=0
for datafile in fileList:
#bkgdfile="/Users/vchaplin/Data/Davis Visit Cavitation Mapping/May 12/background_center_6V.mat"
#datafile="/Users/vchaplin/Data/Davis Visit Cavitation Mapping/May 12/background_center_6V.mat"
#datafile="/Users/Vandiver/Data/Verasonics/sonalleve_0607/3foc-agar-1000av.bin"
    
    datafile_basename=os.path.basename(datafile)
    print('')
    print('Processing: %s'% datafile_basename)    
    
    path=os.path.dirname(datafile)
    params_file,ext = getParamFile(datafile)
    
    params = load_cavitation_params(params_file)
    
    
    if ext=="mat":
        #if mat file
        removeNsamples = 829
        
        rf_data = load_veradata_mat(datafile, params=params, removeNsamples=removeNsamples)
        (ns,Nchan,totframes)=rf_data.shape
        
        
        # look for background file in dirrectory
        candidates = [os.path.basename(fn) for fn in glob.glob(path+"/"+"back*.mat" ) if not re.match(".*params.*",fn)]
        bkgdfile=''
        if (len(candidates) > 0) and (not args.nobg):
            lcsmatches=[stringmisc.longest_common_substring(fn,datafile_basename) for fn in candidates]
            matchlengths = list(map(len, lcsmatches))
            idx=  np.argmax(matchlengths) 
        
            bkgdfile = candidates[idx]
            
            print ("BACKGROUND file: ", bkgdfile, flush=True)
        
            bkg_data = load_veradata_mat(path+"/"+bkgdfile,removeNsamples=removeNsamples)        
                
            rf_data = rf_data - bkg_data
            
            bkgdfile=path+"/"+bkgdfile
            
            
    else:
        #if binary file
        rf_data = load_veradata_bin(datafile, params=params)
        (ns,Nchan,totframes)=rf_data.shape  
            
    if args.maxf:
        numf=args.maxf
    else:            
        numf = params.numframes
        
    if args.maxa:
        numa=args.maxa
    else:            
        numa = params.numacq
        
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
    
    #%% construct delay array on first file only
        
    if fnum==0 or args.each:
        c=1490; 
        Trans={}
        Trans['spacing']=params.pitch;
        Trans['frequency']=params.fs/4*1e-6;
        #wavelength = c/(Trans['frequency']*1e6);
        
        Fs = 4*Trans['frequency']*1e6;
        dt = 1/Fs;
        
        dx = 2e-4
        dz = 1e-4
        
        ducer_width = params.pitch*Nchan
        
        Nx = round(ducer_width/dx)
        Nz = round(ns*c*dt/dz);
        #Nz = round(6e-2/dz)
        
        xpnts = np.linspace(-0.5,0.5,Nx)*Nx*dx
        zpnts = np.linspace(0,1,Nz)*Nz*dz
        
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
    acqs = np.arange(1,numa)
    numSF = len(superFrames)
    
    imgf = np.zeros([numSF,Nz,Nx])
    rmsPerImage = np.zeros([numSF,params.numacq])
    
    for sfi in range(0,numSF):
        sf = superFrames[sfi]
        
        if delayVsFrame is not None:
            
            delayinds = np.round( distances / (c*dt)).astype(int) +delayVsFrame[sf]
            inbounds3 = (delayinds < ns)
            ii=np.nonzero(inbounds3)
        
        
        for a in acqs:
            #delayed[:]=0
            delayed[ii[0],ii[1],ii[2]] = rf_data[delayinds[inbounds3],ii[2],sf*params.numacq + a]
    
            summed_passivemap = np.sum(distances*delayed,axis=2)**2
            
            rmsPerImage[sfi,a] = np.sum(summed_passivemap,axis=(0,1))
            imgf[sfi]+=summed_passivemap
            
            print("\r%d/%d"%(sf*params.numacq + a, superFrames[-1]*params.numacq + acqs[-1]),end='', flush=True)
        
        imgf[sfi]/=len(acqs)
        
    img=np.sum(imgf,axis=0)
    
    #%%
    #img = np.log10( np.abs(hilbert(img)) )
    
    pcmin=np.min(img)
    pcmax=np.max(img)
    
    extent = [ 100*xpnts[0], 100*xpnts[-1], 100*zpnts[-1], 100*zpnts[0] ]
    #try:
    #    ax=fig.gca()
    #    im.set_data(img)
    #    im.set_clim(vmin=pcmin,vmax=pcmax)
    #    im.set_extent(extent)
    #    cbar.set_array(img)
    #    plt.draw_all()
    #    fig.show()
    #except NameError:
    
#    fig=plt.figure(figsize=(8,10))
#    ax=plt.gca()
#    im=ax.imshow(img,cmap=image.cm.gray,extent=extent,interpolation='none')
#    cbar=plt.colorbar(mappable=im)
#    ax.set_xlabel('cm',fontsize=14)
#    ax.set_ylabel('cm',fontsize=14)
#    ax.tick_params(labelsize=14)
    
    path=os.path.dirname(datafile)
    folder = path.split('/')[-1].replace(' ','_')
    (dataname,ext) = os.path.basename(datafile).split(".")
    
    if args.plot:
        if not interactive:
            figName = path+"/"+dataname+".png"
            
            print("Writing: ", figName,flush=True)
            plt.ioff()
        fig=plt.figure(figsize=(8,10))
        ax=plt.gca()
        im=ax.imshow(img,cmap=image.cm.gray,extent=extent,interpolation='none')
        cbar=plt.colorbar(mappable=im)
        ax.set_xlabel('cm',fontsize=14)
        ax.set_ylabel('cm',fontsize=14)
        ax.tick_params(labelsize=14)        
        
        if not interactive:
            plt.savefig(figName)    
            plt.close()
        
    if h5outfile is not None:
        f = h5py.File(h5outfile, mode='a')
        
        dsetExtension = folder+"/"+dataname
        #check if extension is already present. If so, remove it before writing the new data.
        try:
            dset = f[dsetExtension]
            del f[dsetExtension]
            f.flush()
            
            print("Writing extension (existing deleted): ", dsetExtension,flush=True)
        except KeyError:
            print("Writing extension: ", dsetExtension,flush=True)
            #pass
        
        dset = f.create_dataset(dsetExtension+"/map",data=img)
        dset.attrs["dz"]=dz
        dset.attrs["dx"]=dx
        dset.attrs["dt"]=dt
        dset.attrs["c"]=c
        dset.attrs['file']=datafile
        dset.attrs['bkg']=bkgdfile
        f.flush()
        
        dset = f.create_dataset(dsetExtension+"/rms",data=rmsPerImage)
        f.flush()
        f.close()

    fnum+=1