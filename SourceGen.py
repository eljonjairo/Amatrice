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

SlipinMat = np.zeros([ndipin,nstkin,ntwin])
SlipFinMat = np.zeros([ndipin,nstkin])

twin = np.arange(6,46,2)
itw = 0
for it in twin:
    Sliptmp = np.flipud(np.reshape(inData[:,it], (ndipin, nstkin)))
    SlipFinMat = SlipFinMat + Sliptmp
    SlipinMat[:,:,itw] = Sliptmp
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
Slip = Fault['SlipMat']
nstk = Fault['nstk']
ndip = Fault['ndip']
hypoistk = Fault['hypoistk']
hypoidip = Fault['hypoidip']

# Comparison Interpolated Slip and Sum of Srate time windows inputs.
fig = plt.figure()
ax = fig.subplots(1,2)
sin = ax[0].pcolormesh(stkMat,dipMat,Slip,cmap='hsv')
plt.colorbar(sin,location='bottom',label="Slip (m)",shrink=.9)
sinf = ax[1].pcolormesh(stkinMat,dipinMat,SlipFinMat,cmap='hsv')
plt.colorbar(sinf,location='bottom',label="Slip (m)",shrink=.9)

# Spatial interpolation of srates
SlipMatEsp = np.zeros([ndip,nstk,ntwin])
SlipFMatEsp = np.zeros([ndip,nstk])
for itw in range (0,ntwin):
    Sliptmp = SlipinMat[:,:,itw] 
    SlipF = interpolate.interp2d(stkinVec,dipinVec,Sliptmp, kind = "linear")
    SlipMatEsp[:,:,itw] = SlipF(stkVec, dipVec)
    SlipFMatEsp = SlipFMatEsp + SlipMatEsp[:,:,itw]

# Temporal interpolation of slips
inTime = np.arange(0,8,dtwin)
tmax = np.max(inTime);
Time = np.arange(0,tmax,dt)
nt = int(Time.size)
SlipMat = np.zeros([ndip,nstk,nt])

for istk in range (0,nstk):
    for idip in range (0,ndip ):
        SlipTime  = SlipMatEsp[idip,istk,:]
        SlipTimeF = interpolate.interp1d(inTime,SlipTime,kind="cubic")
        SlipMat[idip,istk,:] = SlipTimeF(Time)

# Check temporal interpolation in the hypocenter
fig = plt.figure()
plt.scatter(inTime,SlipMatEsp[hypoidip,hypoistk,:])
plt.plot(Time,SlipMat[hypoidip,hypoistk,:].transpose())

# Calculate SlipRates
SRate = np.zeros([ndip,nstk,nt])
for it in range(1,nt):
    for istk in range (0,nstk):
        for idip in range (0,ndip ):
            SRate[idip,istk,it] = ( SlipMat[idip,istk,it]-SlipMat[idip,istk,it-1])/dt
            
itws = [0,1,2,3,4,5,6,7]
fig = plt.figure()
ax = fig.subplots(8,1)
for itw in itws:
    it = int(itw/dtwin)
    itime = ("%3.1f" %(it*dtwin) )
    SRatefig = SRate[:,:,it]
    sin = ax[itw].pcolormesh(stkMat,dipMat,SRatefig,cmap='hsv',vmin=0,vmax=0.3)
    ax[itw].text(32,7,itime,fontsize=8,fontweight='bold',color='white')

# plt.colorbar(sin,location='bottom',label="Slip (m)",shrink=.9)     

print("  ")
print(" END PROGRAM ")
print("  ")




