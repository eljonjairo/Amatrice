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
import matplotlib.colors
from matplotlib.colors import Normalize

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
dt = 0.1;
lowcut  = 1.0                              # srate low pass cut frecuency
highcut = 0.01                             # srate cut high pass frecuency
fn = 1/dt/2                                # Nyquist frecuency (Sample rate/2)
 

rake = -85

# ID of simulation
ID = 3

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

# Load Slip in each time window
DeltaSlipinMat = np.zeros([ndipin,nstkin,ntwin])
twin = np.arange(6,46,2)
itw = 0
for it in twin:
    #DSliptmp = np.flipud(np.reshape(inData[:,it], (ndipin, nstkin)))
    DSliptmp = np.reshape(inData[:,it], (ndipin, nstkin))
    DeltaSlipinMat[:,:,itw] = DSliptmp
    itw += 1

# Calculate the accumulated slip in each time window
# Acumulated Slip
SlipinMat = np.zeros([ndipin,nstkin,ntwin])
for it in range (1,itw):
    SlipinMat[:,:,it] = SlipinMat[:,:,it-1] + DeltaSlipinMat[:,:,it] 

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

# Spatial interpolation of accumulated  slip
SlipMatEsp = np.zeros([ndip,nstk,ntwin])

for itw in range (0,ntwin):
    Sliptmp = SlipinMat[:,:,itw] 
    SlipF = interpolate.interp2d(stkinVec,dipinVec,Sliptmp, kind = "linear")
    SlipMatEsp[:,:,itw] = SlipF(stkVec, dipVec)
   

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

# Comparison of loaded Interpolated Final Slip and last accumulated Slip 
# spatial and temporal interpolated window.
fig = plt.figure()
ax = fig.subplots(1,2)
sin = ax[0].pcolormesh(stkMat,dipMat,Slip,cmap='hsv')
plt.colorbar(sin,location='bottom',label="Slip (m)",shrink=.9)
sinf = ax[1].pcolormesh(stkMat,dipMat,SlipMat[:,:,-1],cmap='hsv')
plt.colorbar(sinf,location='bottom',label="Slip (m)",shrink=.9)

# Check temporal interpolation in the hypocenter
fig = plt.figure()
plt.scatter(inTime,SlipMatEsp[hypoidip,hypoistk,:])
plt.plot(Time,SlipMat[hypoidip,hypoistk,:].transpose())

# Sliprate calculation
SRate = np.zeros([ndip,nstk,nt])
for it in range(1,nt):
    for istk in range (0,nstk):
        for idip in range (0,ndip ):
            SRate[idip,istk,it] = ( SlipMat[idip,istk,it]-SlipMat[idip,istk,it-1])/dt
            
# Filtered Sliprate
b, a = signal.butter(4, lowcut/fn)
SRateFilt = SRate
for istk in range (0,nstk):
    for idip in range (0,ndip ):
        SRatetmp = SRate[idip,istk,:]                
        SRateFilt[idip,istk,:] = signal.filtfilt(b, a, SRatetmp)

tfig = np.array([0,1,2,3,4,5,7])
ntfig = int(tfig.size)
itfig = np.zeros([ntfig,])
for it in range (0,ntfig):
    itfig[it] = int(tfig[it]/dt)
 

fig, axes = plt.subplots(nrows=8, ncols=1,sharex='col', sharey='row')
for i,ax in enumerate(axes.flat):
    it = int(i/dt)
    itime = ("%3.1f" %(it*dt) )
    SRatefig = SRate[:,:,it]
    im=ax.pcolormesh(stkMat,dipMat,SRatefig,cmap="hot_r",vmin=0,vmax=0.3)
    ax.text(32,7,itime,fontsize=8,fontweight='bold',color='black')
fig.colorbar(im, ax=axes.ravel().tolist())


fig, axes = plt.subplots(nrows=8, ncols=1,sharex='col', sharey='row')
for i,ax in enumerate(axes.flat):
    it = int(i/dt)
    itime = ("%3.1f" %(it*dt) )
    SRatefig = SRateFilt[:,:,it]
    im=ax.pcolormesh(stkMat,dipMat,SRatefig,cmap="hot_r",vmin=0,vmax=0.3)
    ax.text(32,7,itime,fontsize=8,fontweight='bold',color='black')
fig.colorbar(im, ax=axes.ravel().tolist())


# Slip positive in the direction of dip and in the direction of strike
SRdip = (-SRate.flatten(order='F'))*np.sin(np.deg2rad(rake))
SRstk = (SRate.flatten(order='F'))*np.cos(np.deg2rad(rake))
SRdip = np.float32(SRdip)
SRstk = np.float32(SRstk)
Nsr = np.arange(0,SRdip.size,1)

outname = Fault['outname']

SRdipName = 'Outputs/3DFaults/srate_dip_'+outname+'_ID_'+str(ID)
SRstkName = 'Outputs/3DFaults/srate_str_'+outname+'_ID_'+str(ID)

fsrd = FortranFile(SRdipName, 'w')
fsrd.write_record(SRdip)
fsrs = FortranFile(SRstkName, 'w')
fsrs.write_record(SRstk)

# Write fcoor file
nflt = nstk*ndip
fcoorHeader = "%d  %d %4.2f " %(nflt, nt, dt)
fcoor = np.array(Fault['fcoor'])*m
fcoorName = 'Outputs/3DFaults/fcoor_'+outname+'_ID_'+str(ID)+'.in'

with open(fcoorName,'wb') as f:
    np.savetxt(f, fcoor, header=fcoorHeader, comments=' ',fmt = '%9.4f')

print(f" Coordinates saved in file: {fcoorName}" )
print(f" SlipRate dip saved in file: {SRdipName}" )
print(f" SlipRate stk saved in file: {SRstkName}" )


print("  ")
print(" END PROGRAM ")
print("  ")





 



