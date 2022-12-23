#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Generate Obspy Stream per Station with Seismic Data
#
# John Diaz December 2022

import os
import warnings

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import obspy
from obspy.core import UTCDateTime
from obspy.core.stream import Stream

import pickle

from IPython import get_ipython
get_ipython().magic('reset -sf')
warnings.filterwarnings("ignore")
plt.close('all')

os.system('clear')

# DG folder
# DGFolder = "DGrun/Amatrice_ID1_1.0Hz/"
# DGFolder = "DGrun/Amatrice_ID2_1.0Hz/"
DGFolder = "DGrun/Amatrice_ID3_1.0Hz/"
# Event time
teven = "2016-08-24T01:36:32"

# Stations name in the same order as GEODG3D.acquisition DG file.
Stats = [ "AMT","NRC","CSC","RM33","MSC","SPD","LSS","PZI1","TERO","FEMA","ANT"
         ,"TRL","ASP","SNO","SPM","MNF","TRE","CLF","FOS"]

nstats = int(len(Stats))

#nt = 822
#dt = 0.0486502

nt = 813
dt = 0.0491877

print("  ")
print(" START PROGRAM ")
print("  ")

fileFSx = DGFolder + "VX_1"
fileFSy = DGFolder + "VY_1"
fileFSz = DGFolder + "VZ_1"

print(" Loading DG velocity files:\n ")
print(f" {fileFSx}")
print(f" {fileFSy}")
print(f" {fileFSz}")

Vsyn_FSx = np.reshape(np.fromfile(fileFSx, dtype=np.float32),(nstats,nt), order ='F')
Vsyn_FSy = np.reshape(np.fromfile(fileFSy, dtype=np.float32),(nstats,nt), order ='F')
Vsyn_FSz = np.reshape(np.fromfile(fileFSz, dtype=np.float32),(nstats,nt), order ='F')

print()
print(' Loading Syntetic Velocity files:')
print()

for istat in range (0,nstats):
    # DG Data Proccesing
    velox_DG = Vsyn_FSx[istat,:]
    veloy_DG = Vsyn_FSy[istat,:]
    veloz_DG = Vsyn_FSz[istat,:]
    
    trx = obspy.Trace()
    trx.stats.station = Stats[istat]
    trx.stats.starttime = UTCDateTime(teven)
    trx.stats.delta     = dt
    trx.stats.npts      = nt
    trx.stats.channel   = "DGx"
    trx.data            = velox_DG*100
    trx.plot()
    
    tryy = obspy.Trace()
    tryy.stats.station = Stats[istat]
    tryy.stats.starttime = UTCDateTime(teven)
    tryy.stats.delta     = dt
    tryy.stats.npts      = nt
    tryy.stats.channel   = "DGy"
    tryy.data            = veloy_DG*100
    tryy.plot()
    
    trz = obspy.Trace()
    trz.stats.station = Stats[istat]
    trz.stats.starttime = UTCDateTime(teven)
    trz.stats.delta     = dt
    trz.stats.npts      = nt
    trz.stats.channel   = "DGz"
    trz.data            = veloz_DG*100
    trz.plot()
    
    st = Stream([trx,tryy,trz])
    stfile = DGFolder+Stats[istat]+"_DGVEL.pickle" 
    print(f" Writing Stream in file: {stfile} ")
    st.write(stfile, format = "PICKLE")

