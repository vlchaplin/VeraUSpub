# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 11:06:18 2016

@author: caskeylab
"""


import numpy as np
from scipy import io


import os

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
    
def load_veradata_128_bin(datafile,params=None):
    """
    Load a 128-recieve channel data set
    Output rf_data is split into a 128-channel arrays of shape (numsamples, 128, numframes)
    
    
    """
    if params is None:
        params = load_cavitation_params(params_file)
    
    ns=params.numRcvSamples
    Nchan = params.numRcvChannels
    
    if Nchan!=128:
        raise ValueError("Expected 128-channel data (params file "+params_file+")")
    
    bmodes = params.num_final_bmode_acs
    
    rf_data = np.zeros([ns, Nchan, params.numframes*params.numacq],dtype=np.int16)     
    
    bmodeStartidx = ns*Nchan*params.numframes*params.numacq
    
    with open(datafile, 'rb') as fid:
        chandat = np.fromfile(fid, np.int16)
        if bmodes>0:
            bmodedat = np.zeros([ns,128,bmodes],dtype=np.int16)
            bchandat = chandat[bmodeStartidx:].reshape([1,128,bmodes,ns])
            for bi in range(bmodes):
                print(bi)
                bmodedat[:,:,bi]=bchandat[0,:,bi,:].transpose()
                            
            bmodedat=np.array_split(bmodedat,2,axis=1)
        else:
            bmodedat=None
            
        chandat = chandat[0:bmodeStartidx].reshape([params.numframes,Nchan,params.numacq,ns])
        
    for f in range(0, params.numframes):
        for a in range(0, params.numacq):
            #rf_data[:,:,f*params.numacq + a] = chandat[:,a,:,f]
            rf_data[:,:,f*params.numacq + a] = chandat[f,:,a,:].transpose()
   
    
    return (rf_data, bmodedat)
    
    
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
    
    rf_data = np.zeros([ns, Nchan, params.numframes*params.numacq],dtype=np.int16)     
    
    bmodeStartidx = ns*Nchan*params.numframes*params.numacq
    
    with open(datafile, 'rb') as fid:
        chandat = np.fromfile(fid, np.int16)
        if bmodes>0:
            bmodedat = np.zeros([ns,256,bmodes],dtype=np.int16)
            bchandat = chandat[bmodeStartidx:].reshape([1,256,bmodes,ns])
            for bi in range(bmodes):
                print(bi)
                bmodedat[:,:,bi]=bchandat[0,:,bi,:].transpose()
                            
            bmodedat=np.array_split(bmodedat,2,axis=1)
        else:
            bmodedat=None
            
        chandat = chandat[0:bmodeStartidx].reshape([params.numframes,Nchan,params.numacq,ns])
        
    for f in range(0, params.numframes):
        for a in range(0, params.numacq):
            #rf_data[:,:,f*params.numacq + a] = chandat[:,a,:,f]
            rf_data[:,:,f*params.numacq + a] = chandat[f,:,a,:].transpose()
   
    (probe1, probe2) = np.array_split(rf_data,2,axis=1)
    
    return (probe1, probe2, bmodedat)
    
    
    
    
def getProbe2DelaySet(uz2,ux2, zpnts, xpnts, c,dt,ns,yplane=0.0,trips=1,fnum=0):
    ndZ2,ndX2,ndux2 = np.meshgrid(zpnts,xpnts,ux2, indexing='ij')
    distances2 = trips*np.sqrt( (ndX2-ndux2)**2 + yplane**2 + (ndZ2-uz2)**2 )
    delayinds2 = np.round( distances2 / (c*dt)).astype(int)
    inbounds2 = (delayinds2 < ns)
    ii2=np.nonzero(inbounds2)
    ape2 = np.abs((ndZ2-uz2) / (2*(ndX2-ndux2)))>=fnum
        
    return  (distances2,delayinds2, inbounds2, ii2,ape2)

def probe_ccf(delaySets,frames=[0],acqs=[1],do_fft_filt=False,rfSets=[],parSets=[]):
    probenum=0
    numSF=len(frames)
    
    Nz,Nx,Nchan=delaySets[0][0].shape    
    
    ProbeCrossCorr = np.ones([numSF,len(acqs),Nz,Nx])
    
    mom2 = np.zeros([Nz,Nx])
    mom1 = np.zeros([Nz,Nx])
    for rf_data in rfSets:

        #delayed*=0
        (distances, delayinds, inbounds3, ii,ape) = delaySets[probenum]
        
        params = parSets[probenum]
        
        delayed = np.zeros_like(distances)

        for sfi in range(0,numSF):
            sf = frames[sfi]

            for ai in range(len(acqs)):
                a=acqs[ai]
                delayed[:]=0

                if do_fft_filt:
                    #need to pass NFFT, harmonicFilt, etc.)
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

                beamformed=ape*distances*delayed
                chandelaysum = np.sum(beamformed,axis=2)
                chandelaysum2 = np.sum( (beamformed)**2,axis=2)

                ProbeCrossCorr[sfi,ai]*=chandelaysum

                mom1+=chandelaysum**2
                mom2+=chandelaysum2


                print("\r%d/%d"%(sf*params.numacq + a, frames[-1]*params.numacq + acqs[-1]),end='', flush=True)

        probenum+=1
    
    avgCorr = np.mean(ProbeCrossCorr,axis=(0,1))

    return (avgCorr, mom1, mom2)