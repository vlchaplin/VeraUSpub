# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 09:54:58 2016

@author: vlchaplin@gmail.com
"""

import h5py
import sys
import numpy as np
from scipy import io

from math import *;

import re
import glob
import os
import matplotlib.pyplot as plt
from matplotlib import image


import pandas as pd

geom_file = "/Users/Vandiver/Data/DavisData/pos_elmts_Imasonic.mat"

geomat=io.loadmat(geom_file,squeeze_me=True,struct_as_record=False)


plt.plot(geomat['pos_elmts'][:,0], geomat['pos_elmts'][:,1],'o')
plt.plot(geomat['pos_elmts'][:,0], geomat['pos_elmts'][:,2],'o')


zz=np.arange(45,60, 1.0)
stdrr = list(map(lambda z: np.std( np.sqrt(np.sum( (geomat['pos_elmts'] - [0., 0., z])**2, axis=1))), zz ))
plt.plot(zz,stdrr)

xscans = glob.glob("/Users/Vandiver/Data/DavisData/May_13/x_axis_scans/*.mat")
xscansA = [fn for fn in xscans if not ( re.match(".*params.*",fn) or re.match(".*point.*",fn) )]
xscansB = [fn for fn in xscans if  ( not re.match(".*params.*",fn) and re.match(".*point.*",fn) )]



#%%

#check if pulse is shorter than travel length disparity
xp=5
yp=55
x1=50
y1=34
xN=-50
yN=34

d1 = sqrt( (x1-xp)**2 + (y1-yp)**2)
dN = sqrt( (xN-xp)**2 + (yN-yp)**2)

cycleLength =1e3*(1490/1e6)*10

dN-d1 > cycleLength
