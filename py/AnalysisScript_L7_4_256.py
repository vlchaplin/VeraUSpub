# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 11:13:53 2016

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
from matplotlib.backends.backend_pdf import PdfPages

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
    Actual filter is [fn + window[0], fn + window[1]] around each step
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
    
args.f =["/Users/Vandiver/Data/Verasonics/sonalleve_20160709/multi_test/multi_test_30W_rot=45.bin"]    
args.plot=True

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

do_fft_filt=False
#args.plot=True
args.maxf=2
args.maxa=4
bkgdfile=""
fnum=0

#%%
for datafile in fileList:
    datafile_basename=os.path.basename(datafile)
    print('')
    print('Processing: %s'% datafile_basename)    
    
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
        
        
        
    
    #construct delay array on first file only
        
    if fnum==0 or args.each:
        c=1540; 
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
        #Nz = round(ns/2*c*dt/dz);
        #Nz = round(9e-2/dz)
        
        xpnts = np.linspace(-0.5,0.5,Nx)*Nx*dx
        #zpnts = np.linspace(0,1,Nz)*Nz*dz
        zpnts = np.arange(4e-2,9e-2,dz)
        Nz=len(zpnts)
        
        #sensor positions
        ux = (np.linspace(-0.5,0.5,Nchan)*Nchan + 0.5)*params.pitch     
                
        #image and resampling grids
        ndZ,ndX,ndux = np.meshgrid(zpnts,xpnts,ux, indexing='ij')
        
        distances = np.sqrt( (ndX-ndux)**2 + ndZ**2 )
        delayinds = np.round( distances / (c*dt)).astype(int)
        inbounds3 = (delayinds < ns)
        ii=np.nonzero(inbounds3)
        
        if do_fft_filt and NFFT!=ns:
            NFFT = 2**ceil(log2(ns)) #next power of 2
            endidx=floor(NFFT/2)+1
            f0MHz=1.2
            maxFiltFreq=20.0
            nharm=ceil(maxFiltFreq/f0MHz)
            fx = 1e-6*np.linspace(0,1,NFFT)*Fs
            win=[-0.1, 0.1]
            harmonicFilter = LogisticNotchFilter(fx,f0MHz,[-f0MHz, 0.1],1,harmstep=1.0,scale=200)
            harmonicFilter *= LogisticNotchFilter(fx,f0MHz,[-0.1, 0.1],nharm,harmstep=1.0,scale=200)        
            harmonicFilter *= LogisticNotchFilter(fx,f0MHz*0.5,[-0.1, 0.1],nharm,harmstep=1.0,scale=200)
            harmonicFilter=1-harmonicFilter
    
    #distflat = distances[inbounds3]
    
    #%%image formation.  Repeat for each probe data
    probenum=0
    for rf_data in (rf_probe1, rf_probe2):
        probenum+=1        
        
        delayed=np.zeros([Nz,Nx,Nchan])
        
        superFrames = np.arange(0,numf)
        #superFrames = np.arange(27,27+numf)
        
        acqs = np.arange(0,numa)
        numSF = len(superFrames)
        
        Moment1Imgf = np.zeros([numSF,Nz,Nx])
        Moment2Imgf = np.zeros([numSF,Nz,Nx])
        rmsPerImage = np.zeros([numSF,params.numacq])
        bbNoisePerAcq = np.zeros([numSF*params.numacq])
        
        if do_fft_filt:
            spectrumPerFrame = np.zeros([numSF,endidx])        
        
        for sfi in range(0,numSF):
            sf = superFrames[sfi]
            
            if delayVsFrame is not None:
                
                delayinds = np.round( distances / (c*dt)).astype(int) +delayVsFrame[sf]
                inbounds3 = (delayinds < ns)
                ii=np.nonzero(inbounds3)
            
            
            for a in acqs:
                #delayed[:]=0
            
                if do_fft_filt:
                    rf_fft = np.fft.fft( rf_data[:,:,sf*params.numacq + a],NFFT,axis=0)
        
                    #efficiently multiply the filter along the first dimension
                    filtFFT = rf_fft*harmonicFilter[:,np.newaxis]
                    chSumSpec = np.sum(np.abs(2*rf_fft[0:endidx]),axis=1)
                    bbNoiseLevel = np.sum( 2*np.abs(filtFFT[0:endidx])) 
                    
                    delayed[ii[0],ii[1],ii[2]] = np.real ( np.fft.ifft(filtFFT,n=NFFT,axis=0) )[delayinds[inbounds3],ii[2]]
                     
                    bbNoisePerAcq[sfi*params.numacq + a]=bbNoiseLevel
                    spectrumPerFrame[sfi] += chSumSpec               
                else:
                    
                    delayed[ii[0],ii[1],ii[2]] = rf_data[delayinds[inbounds3],ii[2],sf*params.numacq + a]
        
                chandelaysum = np.sum(distances*delayed,axis=2)
                summed_passivemap = chandelaysum**2
                
                rmsPerImage[sfi,a] = np.sum(summed_passivemap,axis=(0,1))
                
                Moment1Imgf[sfi]+=chandelaysum
                Moment2Imgf[sfi]+=summed_passivemap
                
                print("\r%d/%d"%(sf*params.numacq + a, superFrames[-1]*params.numacq + acqs[-1]),end='', flush=True)
            
            Moment1Imgf[sfi]/=len(acqs)        
            Moment2Imgf[sfi]/=len(acqs)
            
        img1=np.mean(Moment1Imgf,axis=0)
        img2=np.mean(Moment2Imgf,axis=0)
        
        
        #img = np.log10( np.abs(hilbert(img)) )
        
        pcmin=np.min(img2)
        pcmax=np.max(img2)
        
        extent = [ 100*xpnts[0], 100*xpnts[-1], 100*zpnts[-1], 100*zpnts[0] ]
    
        
        
        #folder = path.split('/')[-1].replace(' ','_')
        folder = os.path.basename(path).replace(' ','_')
        (dataname,ext) = os.path.basename(datafile).split(".")
        
        if args.plot:
            #make image plots
            if not interactive:
                figName = path+"/"+dataname+("_probe%d"%probenum)+".png"
                
                print("Writing: ", figName,flush=True)
                plt.ioff()
            fig=plt.figure(figsize=(8,10))
            ax=plt.gca()
            im=ax.imshow(img2,cmap=image.cm.gray,extent=extent,interpolation='none')
            cbar=plt.colorbar(mappable=im)
            ax.set_xlabel('cm',fontsize=14)
            ax.set_ylabel('cm',fontsize=14)
            ax.tick_params(labelsize=14)        
            
            if not interactive:
                plt.savefig(figName)    
                plt.close()
                
                
                
            #spectrum plots
            if do_fft_filt:
                if not interactive:
                    
                    pdfName = path+"/"+dataname+("_probe%d"%probenum)+"_filtspec.pdf"
                    pdf = PdfPages(pdfName)
                    print("Writing: ", pdfName,flush=True)
                    plt.ioff()
                
                for sfi in range(numSF):
                    sf=superFrames[sfi]
                    fig=plt.figure(figsize=(9,7))
                    ax=plt.gca()
                    rawSpec = spectrumPerFrame[sfi]
                    usedSpec = (harmonicFilter[0:endidx])*rawSpec       
                    removedSpec= (1-harmonicFilter[0:endidx])*rawSpec
                    ax.plot(fx[0:endidx],usedSpec,'b')
                    ax.plot(fx[0:endidx],removedSpec,color=(0.8,0.8,0.8))
                                
                    ax.set_xlabel('MHZ',fontsize=18)
                    ax.set_ylabel('amplitude',fontsize=18)
                    ax.tick_params(labelsize=18)
                    label = 'Acqs. %d - %d (probe %d)'%(sf*params.numacq+acqs[0]+1, sf*params.numacq+acqs[-1]+1, probenum)
                    ax.text(0.5,0.95, label, horizontalalignment='center',transform=ax.transAxes,fontsize=16, color='r',fontweight='bold')
                    ax.set_ylim([0.0, 2e9])
                
                    if not interactive:
                        pdf.savefig(fig)
                        plt.close(fig)
                if not interactive:
                    pdf.close()
            
        if h5outfile is not None:
            f = h5py.File(h5outfile, mode='a')
            
            dsetExtension = folder+"/"+dataname+"/probe%d"%probenum
            #check if extension is already present. If so, remove it before writing the new data.
            try:
                dset = f[dsetExtension]
                del f[dsetExtension]
                f.flush()
                
                print("Writing extension (existing deleted): ", dsetExtension,flush=True)
            except KeyError:
                print("Writing extension: ", dsetExtension,flush=True)
                #pass
            
            dset = f.create_dataset(dsetExtension+"/map",data=img2)
            dset.attrs["extent"]=extent
            dset.attrs["dz"]=dz
            dset.attrs["dx"]=dx
            dset.attrs["dt"]=dt
            dset.attrs["c"]=c
            dset.attrs['file']=datafile
            dset.attrs['bkg']=bkgdfile
            dset.attrs['framelist']=superFrames
            dset.attrs['acqlist']=acqs
            f.flush()
            
            dset = f.create_dataset(dsetExtension+"/mom1",data=Moment1Imgf)
            f.flush()
            dset = f.create_dataset(dsetExtension+"/mom2",data=Moment2Imgf)
            f.flush()
            
            dset = f.create_dataset(dsetExtension+"/rms",data=rmsPerImage)
            f.flush()
            
            
            if do_fft_filt:
                
                dset = f.create_dataset(dsetExtension+"/BB",data=bbNoisePerAcq)
                dset.attrs["fMHz"]=fx
                dset.attrs["filter"]=harmonicFilter
                dset.attrs["NFFT"]=NFFT
                dset.attrs["endidx"]=endidx
                f.flush()
                
                dset = f.create_dataset(dsetExtension+"/specra",data=spectrumPerFrame)
                f.flush()
            
            
            f.close()
            
    #end of probe for-loop

    fnum+=1