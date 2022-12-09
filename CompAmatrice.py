#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Compare velocity for Amatrice 2017
#
# John Diaz October 2022

# To Do List:
# Add Data
    
import os 
import warnings

import obspy

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from scipy import signal
from scipy.io import FortranFile

from IPython import get_ipython
get_ipython().magic('reset -sf')
warnings.filterwarnings("ignore")
plt.close('all') 
os.system('clear')

# DG folder
DGFolder = Path('DGrun/Amatrice_ID1_1.0Hz/')

# Data folder
DataFolder = Path('Data/ASCII')

Stats = [ "AMT","NRC","CSC","RM33","MSC","SPD","LSS","PZI1","TERO","FEMA","ANT"
         ,"TRL","ASP","SNO","SPM","MNF","TRE","CLF","FOS"]

offset = np.array([ 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
it0 = np.array([100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
itf = it0 + 8000
nstats = int(len(Stats))

# Filtering parameters
lowcut  = 0.5                              # low pass cut frecuency
highcut = 0.06                              # cut high pass frecuency
forder  = 4                                # Butterworth order

nt_DG  = 822
dt_DG  = 0.0486502

# initial time for plotting components (s)
tiy = 5;  tix = 50;  tiz = 95;


print("  ")
print(" START PROGRAM ")
print("  ")

fileFSx = fnameDGx = DGFolder.joinpath('VX_1')
fileFSy = fnameDGx = DGFolder.joinpath('VY_1')
fileFSz = fnameDGx = DGFolder.joinpath('VZ_1')

print(" Loading DG velocity files:\n ")
print(f" {fileFSx}")
print(f" {fileFSy}")
print(f" {fileFSz}")

Vsyn_FSx = np.reshape(np.fromfile(fileFSx, dtype=np.float32),(nstats,nt_DG), order ='F')
Vsyn_FSy = np.reshape(np.fromfile(fileFSy, dtype=np.float32),(nstats,nt_DG), order ='F')
Vsyn_FSz = np.reshape(np.fromfile(fileFSz, dtype=np.float32),(nstats,nt_DG), order ='F')

# DG Time Vector
time_DG = np.linspace(0, dt_DG*(nt_DG-1), nt_DG)

# Coefs for DG filtering
fs_DG   = 1/dt_DG
wlow_DG = lowcut/(fs_DG/2)                         # Normalize the frequency
whp_DG  = highcut/(fs_DG/2)
blow_DG, alow_DG = signal.butter(forder, wlow_DG, 'low')
bhp_DG, ahp_DG   = signal.butter(forder, whp_DG, 'hp')

print()
print(' Loading Syntetic Velocity files:')
print()


    
for istat in range (0,nstats):
   # DG Data Proccesing
    velox_DG = Vsyn_FSx[istat,:]
    veloy_DG = Vsyn_FSy[istat,:]
    veloz_DG = Vsyn_FSz[istat,:]
    # low pass filtering
    DGvelox_low = signal.filtfilt(blow_DG, alow_DG, velox_DG)
    DGveloy_low = signal.filtfilt(blow_DG, alow_DG, veloy_DG)
    DGveloz_low = signal.filtfilt(blow_DG, alow_DG, veloz_DG)
    # high pass filtering
    DGvelox = signal.filtfilt(bhp_DG, ahp_DG, DGvelox_low)
    DGveloy = signal.filtfilt(bhp_DG, ahp_DG, DGveloy_low)
    DGveloz = signal.filtfilt(bhp_DG, ahp_DG, DGveloz_low)

    di = (istat-1)*0.1

    fig = plt.figure(1)
    plt.plot(time_DG+tiy,DGveloy-di,color='k')
    plt.plot(time_DG+tix,DGvelox-di,color='k')
    plt.plot(time_DG+tiz,DGveloz-di,color='k')


for istat in range (0,nstats):
    # Synthetic Data Proccesing
    preStat = '/IT.'
    posStat = '..HGE.D.EMSC-20160824_0000006.VEL.MP.ASC'
    if Stats[istat] == 'RM33' or Stats[istat] == 'TERO' or Stats[istat] == 'FEMA':
        preStat = '/IV.'
        posStat = '..HNE.D.EMSC-20160824_0000006.VEL.MP.ASC'
        
    DataFileVE = DataFolder.joinpath( Stats[istat] + preStat + Stats[istat] + posStat)
    print(DataFileVE)
    with open(DataFileVE, 'r') as f: lines = f.readlines()
    dt_SYN = float(lines[28].split(':')[1])
    nt_SYN = int(lines[29].split(':')[1])-1
    velox_SYN = np.array(lines[64:-1]).astype(np.float)/100
    time_SYN = np.linspace(0, dt_SYN*(nt_SYN-1), nt_SYN)
    # Coefs for SYN filtering
    fs_SYN = 1/dt_SYN
    wlow_SYN = lowcut/(fs_SYN/2)                         # Normalize the frequency
    whp_SYN  = highcut/(fs_SYN/2)
    blow_SYN, alow_SYN = signal.butter(forder, wlow_SYN, 'low')
    bhp_SYN, ahp_SYN = signal.butter(forder, whp_SYN, 'hp')
    # low pass filtering
    SYNvelox_low = signal.filtfilt(blow_SYN, alow_SYN, velox_SYN)
    # high pass filtering
    SYNvelox = signal.filtfilt(bhp_SYN, ahp_SYN, SYNvelox_low)
    
    # DG Data Proccesing
    velox_DG = Vsyn_FSx[istat,:]
    veloy_DG = Vsyn_FSy[istat,:]
    veloz_DG = Vsyn_FSz[istat,:]
    # low pass filtering
    DGvelox_low = signal.filtfilt(blow_DG, alow_DG, velox_DG)
    DGveloy_low = signal.filtfilt(blow_DG, alow_DG, veloy_DG)
    DGveloz_low = signal.filtfilt(blow_DG, alow_DG, veloz_DG)
    # high pass filtering
    DGvelox = signal.filtfilt(bhp_DG, ahp_DG, DGvelox_low)
    DGveloy = signal.filtfilt(bhp_DG, ahp_DG, DGveloy_low)
    DGveloz = signal.filtfilt(bhp_DG, ahp_DG, DGveloz_low)

    
    di = (istat-1)*0.1
    print(f" dt: {dt_SYN} nt: {nt_SYN} "  )
    i0 = it0[istat]
    timesyn = time_SYN[i0:]-offset[istat]
    synvx = SYNvelox[i0:]
    fig = plt.figure(2)
    plt.plot(timesyn+tix,synvx-di,color='r')
    plt.plot(time_DG+tix,DGvelox-di,color='k')

















