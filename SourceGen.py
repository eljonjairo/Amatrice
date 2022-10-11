#!/usr/bin/env python3

#
# Generate sliprates and fault coordinates files for DG
#
# John Diaz July 2022

# To Do List:
#    Filter sliprates 
#    Sliprate animation
#   
    
import os 
import warnings
    
import numpy as np
import pickle
import matplotlib.pyplot as plt

from scipy import signal
import matplotlib.cm as cm
from scipy.io import FortranFile
from pathlib import Path
from scipy import interpolate


from IPython import get_ipython
get_ipython().magic('reset -sf')
warnings.filterwarnings("ignore")
plt.close('all') 
os.system('clear')

inDir  = Path('Input/')
outDir = Path('Outputs/3DFaults/')

infile = 's2016AMATRI01PIZZ.fsp'  # Input file with fault data
name = "AmatricePizzi2017_dhF"  # Output files name

# Name of the Fault file
ObjName = 'AmatricePizzi2017_dhF500m.pickle'

# Source Inputs
dt = 0.01;
lowcut  = 2.0                              # srate low pass cut frecuency
highcut = 0.01                             # srate cut high pass frecuency
fs = 1/dt                                  # Sample rate

rake = -85

# ID of simulation
ID = 1

# Dimension of input fault
nstkin = 30
ndipin = 12
# Number of time windows
ntwin = 20
# dt between time windows  
dtwin = 0.4

# m to Km
m = 1000

print("  ")
print(" START PROGRAM ")
print("  ")

ObjName = outDir.joinpath(ObjName)

# Load the Fault object
with open(ObjName, 'rb') as handle:
    Fault = pickle.load(handle)

#Load Fault input file
infile = inDir.joinpath(infile)

inData = np.loadtxt(infile, skiprows = 50)

SRateinMat = np.zeros([ndipin,nstkin,ntwin])
SlipFinMat = np.zeros([ndipin,nstkin])

twin = np.arange(6,46,2)
itw = 0
for it in twin:
    SRatetmp = np.flipud(np.reshape(inData[:,it], (ndipin, nstkin)))
    SlipFinMat = SlipFinMat + SRatetmp
    SRateinMat[:,:,itw] = SRatetmp
    itw += 1

# Load Fault parameters
stkMat = Fault['stkMat']
dipMat = Fault['dipMat']
stkVec = Fault['stkVec']
dipVec = Fault['dipVec']
stkinMat = Fault['stkinMat']
dipinMat = Fault['dipinMat']
stkinVec = Fault['stkinVec']
dipinVec = Fault['dipinVec']
SlipMat = Fault['SlipMat']
nstk = Fault['nstk']
ndip = Fault['ndip']
hypoistk = Fault['hypoistk']
hypoidip = Fault['hypoidip']

# Comparison Interpolated Slip and Sum of Srate time windows inputs.
fig = plt.figure()
ax = fig.subplots(1,2)
sin = ax[0].pcolormesh(stkMat,dipMat,SlipMat)
plt.colorbar(sin,location='bottom',label="Slip (m)",shrink=.9)
sinf = ax[1].pcolormesh(stkinMat,dipinMat,SlipFinMat)
plt.colorbar(sinf,location='bottom',label="Slip (m)",shrink=.9)

# Spatial interpolation of srates
SRateMatEsp = np.zeros([ndip,nstk,ntwin])
SlipFMatEsp = np.zeros([ndip,nstk])
for itw in range (0,ntwin):
    SRatetmp = SRateinMat[:,:,itw] 
    SRateF = interpolate.interp2d(stkinVec,dipinVec,SRatetmp, kind = "linear")
    SRateMatEsp[:,:,itw] = SRateF(stkVec, dipVec)
    SlipFMatEsp = SlipFMatEsp + SRateMatEsp[:,:,itw]

# Temporal interpolation of srates
inTime = np.arange(0,8,dtwin)
tmax = np.max(inTime);
Time = np.arange(0,tmax,dt)
nt = int(Time.size)
SRateMat = np.zeros([ndip,nstk,nt])
SlipFMat = np.zeros([ndip,nstk])
for istk in range (0,nstk):
   for idip in range (0,ndip ):
       SRateTime  = SRateMatEsp[idip,istk,:]
       SRateTimeF = interpolate.interp1d(inTime,SRateTime)
       SRateMat[idip,istk,:] = SRateTimeF(Time)

# Check temporal interpolation in the hypocenter
fig = plt.figure()
plt.scatter(inTime,SRateMatEsp[hypoidip,hypoistk,:])

plt.plot(Time,SRateMat[hypoidip,hypoistk,:].transpose())
       
for it in range(0,nt):       
       SlipFMat = SlipFMat + SRateMat[:,:,it]
       
# Comparison Interpolated Slip and Sum of Srate time windows inputs.
fig = plt.figure()
ax = fig.subplots(1,3)
sin = ax[0].pcolormesh(stkMat,dipMat,SlipMat)
plt.colorbar(sin,location='bottom',label="Slip (m)",shrink=.9)
sinE = ax[1].pcolormesh(stkMat,dipMat,SlipFMatEsp)
plt.colorbar(sinE,location='bottom',label="Slip (m)",shrink=.9)
sinT = ax[2].pcolormesh(stkMat,dipMat,SlipFMat)
plt.colorbar(sinT,location='bottom',label="Slip (m)",shrink=.9)
       
print("  ")
print(" END PROGRAM ")
print("  ")




