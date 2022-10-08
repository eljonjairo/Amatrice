#!/home/jon/anaconda3/bin/python

#
# Generate Fault and Topography output files
#
# John Diaz July 2022

# To Do List:
#    Add index from hypocenter coordinates
#    Gempy.
#

import os
import warnings

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.spatial import Delaunay
import pickle
from pathlib import Path
import pygmt
import utm
import pandas as pd
import matplotlib.colors
from matplotlib import cm


warnings.filterwarnings("ignore")
os.system('clear')


dhF = 0.5  # Output subfaults size in Km
inDir  = Path('Input/')
outDir = Path('Outputs/3DFaults/')

infile = 's2016AMATRI01PIZZ.fsp'  # Input file with fault data
name = "AmatricePizzi2017_dhF"  # Output files name

# Dimension of input fault
nstkin = 30
ndipin = 12

# Lon-Lat limits for Topo Download
Latmin = 42.1
Latmax = 43.3
Lonmin = 12.4
Lonmax = 13.9

# X and Y UTM limits (Km) for Model
xmin = 300.0
xmax = 400.0
ymin = 4680.0
ymax = 4780.0
zmin = -60.0

# Hypocenter coordinates (Km) 
hypolon = 13.25
hypolat = 42.70
hypoz =-5.0

# m to Km
m = 1000

print("  ")
print(" START PROGRAM ")
print("  ")

# Load Earth relief data for the entire globe and a subset region
region = [Lonmin,Lonmax,Latmin,Latmax]
grid = pygmt.datasets.load_earth_relief(resolution="03s", region=region)

TopoZMat = grid.data
TopoLat = grid.lat.data
TopoLon = grid.lon.data

zmax = np.max(TopoZMat)/m

TopoLonMat, TopoLatMat = np.meshgrid(TopoLon,TopoLat)
TopoXMat,TopoYMat, tmp1, tmp2 = utm.from_latlon(TopoLatMat, TopoLonMat,33,'N')
TopoXMat = TopoXMat/m
TopoYMat = TopoYMat/m

#Load Stations input file
inStats =  pd.read_csv(inDir.joinpath("stations.dat"), delimiter= '\s+')

inStatsX, inStatsY, tmp1, tmp2 = utm.from_latlon(inStats.Lat, inStats.Lon,33,'N')
inStats['X'] = inStatsX/m
inStats['Y'] = inStatsY/m

fig = pygmt.Figure()
pygmt.makecpt(cmap="geo", series=[-10000, 10000])
fig.grdimage(grid=grid, region=region,
             frame=['WSrt+t" Topografia Italia Central"', "xa0.4", "ya0.4"])
fig.colorbar(frame=["xa2000f500+lElevación ", "y+lm"])
#stations
fig.plot(x=inStats.Lon, y=inStats.Lat,style="t0.5c",color='red', pen="black") 
fig.text(textfiles=None,x=inStats.Lon-0.05, y=inStats.Lat+0.05, position=None,text=inStats.Name,  
 angle=0, font='9p,Helvetica-Bold,black', justify='LM')

fig.show()

cvals  = [ 0, 100, 1000, 2000, 3000, 4000, ]
colors = ["darkgreen","green","forestgreen","yellowgreen",
          "orange","maroon","sienna","brown","white"]

norm=plt.Normalize(min(cvals),max(cvals))
tuples = list(zip(map(norm,cvals), colors))
cTopo = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)

#Load Fault input file
infile = inDir.joinpath(infile)

inData = np.loadtxt(infile, skiprows = 50)
LatFinMat = np.reshape(inData[:,0], (ndipin, nstkin))
LonFinMat = np.reshape(inData[:,1], (ndipin, nstkin))
ZFinMat = np.reshape(inData[:,4], (ndipin, nstkin))
SlipinMat = np.reshape(inData[:,5], (ndipin, nstkin))

XFinMat,YFinMat, tmp1, tmp2 = utm.from_latlon(LatFinMat, LonFinMat,33,'N')
XFinMat = XFinMat/m
YFinMat = YFinMat/m



fig = plt.figure(figsize = (10,10))
ax = fig.subplots(1,1)
mp=ax.pcolormesh(TopoXMat,TopoYMat,TopoZMat*m, cmap=cTopo)
plt.colorbar(mp,location='bottom',label="Elevación (m)",shrink=.6)
fp=ax.pcolormesh(XFinMat,YFinMat,SlipinMat, cmap='viridis')
plt.colorbar(fp,location='right',label="Slip (m)",shrink=.6)
ax.scatter(inStats.X,inStats.Y,c ="red",linewidths = 0.5,marker ="^",edgecolor ="black",s = 80)
for i in range(0,inStats.X.size):
    ax.text(inStats.X[i]+0.7, inStats.Y[i]+0.5, inStats.Name[i],fontsize=10,fontweight='bold')

ax.set_xlabel(" X (Km)")
ax.set_ylabel(" Y (Km)")
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_aspect('equal',adjustable='box')
ax.set_title('Topografía Italia Central Input Slip')

# Interpolation of fault plane
dstk = np.array([ XFinMat[0,1]-XFinMat[0,0], YFinMat[0,1]-YFinMat[0,0],
                  ZFinMat[0,1]-ZFinMat[0,0] ])
dstkin = np.linalg.norm(dstk)

ddip = np.array([ XFinMat[1,0]-XFinMat[0,0], YFinMat[1,0]-YFinMat[0,0],
                  ZFinMat[1,0]-ZFinMat[0,0] ])
ddipin = np.linalg.norm(ddip)

dstk = dstk*dhF
ddip = ddip*dhF

# Calculate the strike and dip unitary vetors
univec_stk = np.linalg.norm(dstk)
univec_dip = np.linalg.norm(ddip)

stk = round((nstkin-1)*dstkin)
dip = round((ndipin-1)*ddipin)

hypox, hypoy, tmp1, tmp2 = utm.from_latlon(hypolat,hypolon,33,'N')
hypox = hypox/m
hypoy = hypoy/m
print()
print(" Original Fault Dimensions:")
print(" Strike (Km): %6.2f nstk: %d dstk (Km): %6.2f " %(stk,nstkin,dstkin) )
print(" Dip (Km): %6.2f ndip: %d ddip (Km): %6.2f" %(dip,ndipin,ddipin) )
print(" Hypocenter Coordinates x, y and z (Km): %6.2f %6.2f %6.2f " %(hypox,hypoy,hypoz) )

dipinVec = np.linspace(0, dip, ndipin)
stkinVec = np.linspace(0, stk, nstkin)
stkinMat, dipinMat = np.meshgrid(stkinVec,dipinVec)

#interpolation
nstk = int(stk/dhF)+1
ndip = int(dip/dhF)+1
stkVec = np.linspace(0, stk, nstk)
dipVec = np.linspace(0, dip, ndip)
stkMat, dipMat = np.meshgrid(stkVec,dipVec)

# Slip Interpolation
SlipF = interpolate.interp2d(stkinVec,dipinVec,SlipinMat, kind = "linear")
SlipMat = SlipF(stkVec, dipVec)

#Coordinates Interpolation
inivec = np.array([ XFinMat[0,0], YFinMat[0,0], ZFinMat[0,0] ])

XFMat = np.zeros((ndip,nstk))
YFMat = np.zeros((ndip,nstk))
ZFMat = np.zeros((ndip,nstk))

for istk in range (0,nstk):
    delta_stk = istk*dstk
    for idip in range (0,ndip ):
        delta_dip = idip*ddip
        XFMat[idip,istk] = inivec[0] + delta_stk[0] + delta_dip[0]
        YFMat[idip,istk] = inivec[1] + delta_stk[1] + delta_dip[1]
        ZFMat[idip,istk] = inivec[2] + delta_stk[2] + delta_dip[2]

# From matrix to column vector following fortran 
XF3D = XFMat.flatten(order='F').transpose()
YF3D = YFMat.flatten(order='F').transpose()
ZF3D = ZFMat.flatten(order='F').transpose()


fig = plt.figure()
ax1 = fig.add_subplot(221, projection='3d')
surf = ax1.plot_surface( XFinMat, YFinMat, ZFinMat, facecolors=cm.hsv(SlipinMat), linewidth=0,
                        antialiased=False )
#plt.colorbar(surf,location='top',label="Slip (m)",shrink=.6)
ax1.set_title("Input Slip")
ax1.set_xlabel(" X (Km)")
ax1.set_ylabel(" Y (Km)")
ax1.set_zlabel(" Z (Km)")
ax1.set_xlim(340,380)
ax1.set_ylim(4720,4760)
ax1.set_zlim(-15.0,15.0)
ax1.set_aspect('equal',adjustable='box')
ax1.azim = -120
ax1.dist = 10
ax1.elev = 10

ax2 = fig.add_subplot(222, projection='3d')
surf = ax2.plot_surface( XFMat, YFMat, ZFMat, facecolors=cm.hsv(SlipMat), linewidth=0,
                        antialiased=False )
#plt.colorbar(surf,location='top',label="Slip (m)",shrink=.6)
ax2.set_title("Interpolated Slip")
ax2.set_xlabel(" X (Km)")
ax2.set_ylabel(" Y (Km)")
ax2.set_zlabel(" Z (Km)")
ax2.set_xlim(340,380)
ax2.set_ylim(4720,4760)
ax2.set_zlim(-15.0,15.0)
ax2.set_aspect('equal',adjustable='box')
ax2.azim = -120
ax2.dist = 10
ax2.elev = 10

ax3 = fig.add_subplot(223, projection='3d')
surf = ax3.plot_surface( XFinMat, YFinMat, ZFinMat, facecolors=cm.hsv(SlipinMat), linewidth=0,
                        antialiased=False )
#plt.colorbar(surf,location='top',label="Slip (m)",shrink=.6)
ax3.set_xlabel(" X (Km)")
ax3.set_ylabel(" Y (Km)")
ax3.set_xlim(340,380)
ax3.set_ylim(4700,4750)
ax3.set_zlim(-15.0,15.0)
ax3.set_aspect('equal',adjustable='box')
ax3.azim = -90
ax3.dist = 10
ax3.elev = 90

ax4 = fig.add_subplot(224, projection='3d')
surf = ax4.plot_surface( XFMat, YFMat, ZFMat, facecolors=cm.hsv(SlipMat), linewidth=0,
                        antialiased=False )
#plt.colorbar(surf,location='top',label="Slip (m)",shrink=.6)
ax4.set_xlabel(" X (Km)")
ax4.set_ylabel(" Y (Km)")
ax4.set_xlim(340,380)
ax4.set_ylim(4700,4750)
ax4.set_zlim(-15.0,15.0)
ax4.set_aspect('equal',adjustable='box')
ax4.azim = -90
ax4.dist = 10
ax4.elev = 90
plt.show()

print("  ")
print(" END PROGRAM ")
print("  ")

