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
ax.set_title('Topografía Italia Central')

fig = plt.figure(figsize = (10,10))
ax = fig.subplots(1,3)
ax[0].pcolormesh(XFinMat,YFinMat,SlipinMat, cmap="viridis")
ax[0].set_aspect('equal',adjustable='box')
ax[1].pcolormesh(XFinMat,ZFinMat,SlipinMat, cmap="viridis")
ax[1].set_aspect('equal',adjustable='box')
ax[2].pcolormesh(ZFinMat,YFinMat,SlipinMat, cmap="viridis")
ax[2].set_aspect('equal',adjustable='box')

fig = plt.figure()
ax = plt.axes(projection ='3d')

surf = ax.plot_surface( XFinMat, YFinMat, ZFinMat, facecolors=cm.hsv(SlipinMat), linewidth=0,
                        antialiased=False )
#plt.colorbar(surf,location='top',label="Slip (m)",shrink=.6)
ax.set_title("Input Slip")
ax.set_xlabel(" X (Km)")
ax.set_ylabel(" Y (Km)")
ax.set_zlabel(" Z (Km)")
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_aspect('equal',adjustable='box')
ax.azim = 0
ax.dist = 5
ax.elev = 0
plt.show()

fig = plt.figure()
ax = plt.axes(projection ='3d')

surf = ax.plot_surface( XFinMat, YFinMat, ZFinMat, facecolors=cm.hsv(SlipinMat), linewidth=0,
                        antialiased=False )
#plt.colorbar(surf,location='top',label="Slip (m)",shrink=.6)
ax.set_title("Input Slip")
ax.set_xlabel(" X (Km)")
ax.set_ylabel(" Y (Km)")
ax.set_zlabel(" Z (Km)")
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_aspect('equal',adjustable='box')
ax.azim = -90
ax.dist = 3
ax.elev = 0
plt.show()


fig = plt.figure()
ax = plt.axes(projection ='3d')
surf = ax.plot_surface( XFinMat, YFinMat, ZFinMat, facecolors=cm.hsv(SlipinMat), linewidth=0,
                        antialiased=False )
#plt.colorbar(surf,location='top',label="Slip (m)",shrink=.6)
ax.set_title("Input Slip")
ax.set_xlabel(" X (Km)")
ax.set_ylabel(" Y (Km)")
ax.set_zlabel(" Z (Km)")
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_aspect('equal',adjustable='box')
ax.azim = 0
ax.dist = 5
ax.elev = -90
plt.show()


fig = plt.figure()
ax1 = fig.add_subplot(121, projection='3d')
surf = ax1.plot_surface( XFinMat, YFinMat, ZFinMat, facecolors=cm.hsv(SlipinMat), linewidth=0,
                        antialiased=False )
#plt.colorbar(surf,location='top',label="Slip (m)",shrink=.6)
ax1.set_title("Input Slip")
ax1.set_xlabel(" X (Km)")
ax1.set_ylabel(" Y (Km)")
ax1.set_zlabel(" Z (Km)")
ax1.set_xlim(xmin,xmax)
ax1.set_ylim(ymin,ymax)
ax1.set_zlim(-15.0,15.0)
ax1.set_aspect('equal',adjustable='box')
ax1.azim = -120
ax1.dist = 5
ax1.elev = 10

ax2 = fig.add_subplot(122, projection='3d')
surf = ax2.plot_surface( XFinMat, YFinMat, -ZFinMat, facecolors=cm.hsv(SlipinMat), linewidth=0,
                        antialiased=False )
#plt.colorbar(surf,location='top',label="Slip (m)",shrink=.6)
ax2.set_title("Input Slip")
ax2.set_xlabel(" X (Km)")
ax2.set_ylabel(" Y (Km)")
ax2.set_zlabel(" Z (Km)")
ax2.set_xlim(xmin,xmax)
ax2.set_ylim(ymin,ymax)
ax2.set_zlim(-15.0,15.0)
ax2.set_aspect('equal',adjustable='box')
ax2.azim = -120
ax2.dist = 5
ax2.elev = 10
plt.show()


plt.show()

print("  ")
print(" END PROGRAM ")
print("  ")

