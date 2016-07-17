# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 10:40:39 2016

@author: vlchaplin@gmail.com
"""

import h5py
import sys
import numpy as np
from scipy import io
from math import *;

import glob
from scipy.signal import hilbert as hilbert
import matplotlib.pyplot as plt
from matplotlib import image

sys.path.append('C:\\Users\\vchaplin\\Documents\\HiFU\\code\\myPy')  
sys.path.append('C:\\Users\\vchaplin\\Documents\\HiFU\\code\\AblationSims')
sys.path.append('C:\\Users\\Vandiver\\Documents\\HiFU\\code\\BioHeatCpp\\PBHEswig\\x64')


simDir="/Users/vchaplin/Data/Verasonics/PSF_simulation/runs0706cw/"
simFiles=glob.glob(simDir + "*.mat")

outH5file="/Users/vchaplin/Data/Verasonics/PSF_simulation/PSF_point_0609.h5"


#%% load template file
file="/Users/vchaplin/Data/Verasonics/PSF_simulation/runs0706cw/simulation_x335_y335.mat"

matdict=io.loadmat(file,squeeze_me=True,struct_as_record=False)

#load data
pr_data=matdict['sensor_data'].p
p_final=matdict['sensor_data'].p_final

Nchan=pr_data.shape[0]
Nt=pr_data.shape[1]
dt=matdict['dt']
times = np.arange(0, Nt*dt,dt)
simDims = (matdict['Nx'], matdict['Ny'] )

d = matdict['d']

pitch=matdict['pitch']*matdict['d']

c0 = np.mean(matdict['medium'].sound_speed)

#target image space
ducer_width = pitch*Nchan

dx = 2e-4
dz = 2e-4
Nx = round(ducer_width/dx)
Nz = round(6.0e-2 /dz )

xpnts = np.linspace(-0.5,0.5,Nx)*Nx*dx
zpnts = np.linspace(0,1,Nz)*Nz*dz
#sensor positions
ux = (np.linspace(-0.5,0.5,Nchan)*Nchan + 0.5)*pitch


#image and resampling grids
ndZ,ndX,ndux = np.meshgrid(zpnts,xpnts,ux, indexing='ij')

distances = np.sqrt( (ndX-ndux)**2 + ndZ**2 )
delayinds = np.round( distances / (c0*dt)).astype(int)
inbounds3 = (delayinds <= Nt)
ii=np.nonzero(inbounds3)

#%% determine input sampling
PSFsourcePlacementMask = np.zeros(simDims,dtype=bool)
uniqZ = {}
uniqX = {}

# "x" in the sim file name maps to depth (z) dimension
for fi in range(0,len(simFiles)):
    matdict=io.loadmat(simFiles[fi],squeeze_me=True,struct_as_record=False,variable_names='sourceloc')
    zk,xk=matdict['sourceloc'][[0,1]]
    PSFsourcePlacementMask[zk,xk]=True
    
    uniqZ["%d"%zk]=1
    uniqX["%d"%xk]=1
#%% PSF

npointsZ = len(uniqZ)
npointsX = len(uniqX)

#large dimensionality, store as 32 bit
PSF = np.zeros([npointsZ,npointsX,Nz,Nx],dtype=np.int32)

#%%image formation
delayed=np.zeros([Nz,Nx,Nchan])
delayed[ii[0],ii[1],ii[2]] = pr_data[ii[2], delayinds[inbounds3]]
img=np.sum(delayed*distances,axis=2)

#%%plot
plt.figure(figsize=(12,12))
ext=[xpnts[0],xpnts[-1],zpnts[-1],zpnts[0]]
plt.imshow( img, extent=ext, cmap=image.cm.gray,interpolation='None')


