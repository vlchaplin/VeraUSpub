# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 14:45:07 2016

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

try:
    sys.ps1
    interactive=True
except AttributeError:
    interactive=False
    
def oneLogisticNotchFiltFunc(fx,f_lo,f_hi,scale):
    return 1.0 / (1 + np.exp(scale*(fx-f_lo)) ) + 1.0 / (1 + np.exp(-scale*(fx-f_hi)) )
    
def LogisticNotchFilter(fx,f0,window,nharm,harmstep=1.0, scale=200.0):
    """
    Recursive construction of filter.
    Actual filter is [fn + window[0], fn + window[1]] about each step
    """

    fn = f0*nharm
    f_lo = fn+window[0]
    f_hi = fn+window[1]
    F = oneLogisticNotchFiltFunc(fx,f_lo, f_hi, scale)        
        
    if (nharm-harmstep)<=1e-9:
        #lowest recursion depth (fundamental or lowest subharmonic):
        return F
    else:    
        Fb = LogisticNotchFilter(fx, f0, window, nharm-harmstep, harmstep=harmstep, scale=scale)
    
    return Fb*F
    
def getParamFile(file):
    """
    """
    datapath = os.path.dirname(file)
    (dataname,ext) = os.path.basename(file).split(".")

    return (os.path.join(datapath , dataname +"_params.mat"),ext)

def load_cavitation_params(params_file):
    """
    """
    pardict=io.loadmat(params_file,squeeze_me=True,struct_as_record=False)
    params = pardict['params']
    
    return params    
    
def load_veradata_2x128_bin(datafile,params=None):
    """
    Load a 256-recieve channel data set, formed from two 128-channel probes.
    Output rf_data is split into two 128-channel arrays of shape (numsamples, 128, numframes)
    
    (rf1, rf2) = load_veradata_2x128_bin(file)    
    
    """
    if params is None:
        params = load_cavitation_params(params_file)
    
    ns=params.numRcvSamples
    Nchan = params.numRcvChannels
    
    if Nchan!=256:
        raise ValueError("Expected 256-channel data (params file "+params_file+")")
    
    bmodes = params.num_final_bmode_acs
    
    rf_data = np.zeros([ns, Nchan, params.numframes*params.numacq])     
    
    bmodeStartidx = ns*Nchan*params.numframes*params.numacq
    
    with open(datafile, 'rb') as fid:
        chandat = np.fromfile(fid, np.int16)
        if bmodes>0:
            bmodedat = chandat[bmodeStartidx:].reshape([1,Nchan,bmodes,ns])
        else:
            bmodedat=None
            
        chandat = chandat[0:bmodeStartidx].reshape([params.numframes,Nchan,params.numacq,ns])
        
    for f in range(0, params.numframes):
        for a in range(0, params.numacq):
            #rf_data[:,:,f*params.numacq + a] = chandat[:,a,:,f]
            rf_data[:,:,f*params.numacq + a] = chandat[f,:,a,:].transpose()
   
    (probe1, probe2) = np.array_split(rf_data,2,axis=1)
    
    return (probe1, probe2, bmodedat)
    
#%%
parser = argparse.ArgumentParser()
    
parser.add_argument('-f', metavar='data files', help="Name .bin files to process. Filenames containing string 'params' are ignored", type=str, nargs='+')
parser.add_argument("-o", metavar='output_HDF5',help="Name of output hdf5 file to create. It is opened in append mode", type=str )
parser.add_argument("-maxf", metavar='num',help="end number of super frames (default will use all)", type=int )
parser.add_argument("-maxa", metavar='num',help="end number of acquisitions to use for each frame (default will use all)", type=int )
parser.add_argument("--plot",help="make plots", action="store_true" )
parser.add_argument("--each",help="Recompute image space and delays for each file (default uses the delays from the first file)", action="store_true" )
parser.add_argument("--nobg",help="Do not look for a matching background_* file for background subtraction", action="store_true" )

args = parser.parse_args()    

#args.plot=True    
#args.f =["/Users/Vandiver/Data/Verasonics/sonalleve_20160709/multi_2/multi_2_80W.bin"]    
#h5outfile = "/Users/Vandiver/Data/Verasonics/sonalleve_20160709/pcdmaps20160709.h5"
if args.f :
    fileList = args.f
elif ~interactive:
    print("No files passed",flush=True)
    sys.exit()
if args.o :
    h5outfile = args.o
else:
    h5outfile= None
    
    
delayVsFrame = None


#fileList = glob.glob("/Users/vchaplin/Data/Davis Visit Cavitation Mapping/May 12/*.mat")

#datafile="/Users/vchaplin/Data/Davis Visit Cavitation Mapping/May 12/bubble_center_6V_1.mat"

fileList = [fn for fn in fileList if not re.match(".*params.*",fn)]



NFFT=0    

do_fft_filt=True
bkgdfile=""
fnum=0

for datafile in fileList:
    datafile_basename=os.path.basename(datafile)
    print('')
    print('Processing: %s'% datafile_basename,flush=True)    
    
    path=os.path.dirname(os.path.abspath(datafile))
    params_file,ext = getParamFile(datafile)
    
    params = load_cavitation_params(params_file)
    
    
    #if binary file
    (rf_probe1,rf_probe2, BmodeData) = load_veradata_2x128_bin(datafile, params=params)
    (ns,Nchan,totframes)=rf_probe1.shape  
            
    if args.maxf:
        numf=args.maxf
    else:            
        numf = params.numframes
        
    if args.maxa:
        numa=args.maxa
    else:            
        numa = params.numacq
        


    c=1540; 
    Trans={}
    Trans['spacing']=params.pitch;
    Trans['frequency']=params.fs/4*1e-6;
    #wavelength = c/(Trans['frequency']*1e6);
    
    Fs = 4*Trans['frequency']*1e6;
    dt = 1/Fs;
    
    if NFFT!=ns:
        NFFT = 2**ceil(log2(ns)) #next power of 2
        endidx=floor(NFFT/2)+1
        f0MHz=1.2
        maxFiltFreq=20.0
        nharm=ceil(maxFiltFreq/f0MHz)
        fx = 1e-6*np.linspace(0,1,NFFT)*Fs
        win=[-0.1, 0.1]
        harmonicFilter = LogisticNotchFilter(fx,f0MHz,[-f0MHz, 0.1],1,harmstep=1.0,scale=200)
        harmonicFilter *= LogisticNotchFilter(fx,f0MHz,[-0.05, 0.05],nharm,harmstep=1.0,scale=200)        
        harmonicFilter *= LogisticNotchFilter(fx,f0MHz*0.5,[-0.03, 0.03],nharm,harmstep=1.0,scale=200)
        #harmonicFilter *= LogisticNotchFilter(fx,f0MHz*0.25,[-0.03, 0.03],nharm,harmstep=1.0,scale=200)
        #print("HF max",np.max(harmonicFilter),flush=True)        
        
    probenum=0
    for rf_data in (rf_probe1, rf_probe2):
        
        
        
        probenum+=1
        summed_raw_ch = np.sum(rf_data,axis=1)
        rf_fft = np.fft.fft( summed_raw_ch,NFFT,axis=0)
        
        #efficiently multiply the filter along the first dimension
        filtFFT = rf_fft*harmonicFilter[:,np.newaxis]
        
        bbNoiseLevel_vs_acq = np.sum( np.abs(filtFFT[0:endidx]),axis=0)
        
        if h5outfile is not None:
            
            folder = os.path.basename(path).replace(' ','_')
            (dataname,ext) = os.path.basename(datafile).split(".")            
            
            f = h5py.File(h5outfile, mode='a')
            
            dsetExtension = folder+"/"+dataname+"/probe%d"%probenum+"/BB"
            
            try:
                dset = f[dsetExtension]
                del f[dsetExtension]
                f.flush()
                
                print("Writing extension (existing deleted): ", dsetExtension,flush=True)
                
            except KeyError:
                print("Writing extension: ", dsetExtension,flush=True)
                
                #pass
            
            dset = f.create_dataset(dsetExtension,data=bbNoiseLevel_vs_acq)
            dset.attrs["fMHz"]=fx
            dset.attrs["filter"]=harmonicFilter
            dset.attrs["NFFT"]=NFFT
            dset.attrs["endidx"]=endidx
            f.flush()
            f.close()
            
        if args.plot:
            if not interactive:
                (dataname,ext) = os.path.basename(datafile).split(".") 
                figName = path+"/"+dataname+("_probe%d"%probenum)+"_filtspec100acq.png"
                
                print("Writing: ", figName,flush=True)
                plt.ioff()
            fig=plt.figure(figsize=(9,7))
            ax=plt.gca()
            
            filtSpec = np.sum(np.abs(filtFFT[0:endidx,0:100]),axis=1)        
            removedSpec= (1-harmonicFilter[0:endidx])*np.sum(np.abs(rf_fft[0:endidx,0:100]),axis=1)
            ax.plot(fx[0:endidx],filtSpec,'b')
            ax.plot(fx[0:endidx],removedSpec,color=(0.8,0.8,0.8))
                        
            ax.set_xlabel('MHZ',fontsize=14)
            ax.set_ylabel('amplitude',fontsize=14)
            ax.tick_params(labelsize=14)        
            
            if not interactive:
                plt.savefig(figName)    
                plt.close(fig)