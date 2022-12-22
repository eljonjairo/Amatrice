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
DGFolder = "DGrun/Amatrice_ID1_1.0Hz/"

# Data folder
DataFolder = "Data/Streams/"

# Stats = [ "AMT","NRC","CSC","RM33","MSC","SPD","LSS","PZI1","TERO","FEMA","ANT"
#         ,"TRL","ASP","SNO","SPM","MNF","TRE","CLF","FOS"]

Stats  = [ "AMT","NRC","MNF","TRL"]


# Filtering parameters
lowcut  = 0.5                              # low pass cut frecuency
highcut = 0.06                             # cut high pass frecuency

# initial time for plotting components (s)
tiy = 30;  tix = 80;  tiz = 130;
xmin = 0; xmax = 205
ymin = -10; ymax = 0

print("  ")
print(" START PROGRAM ")
print("  ")

di=-1

for stat in Stats:
    # DG Synthetic Data Reading
    DGfile = DGFolder + stat + "_DGVEL.pickle"
    print(f" Reading DG file {DGfile}")
    DGst = obspy.read(DGfile)    
    DGvx = DGst[0]
    DGvy = DGst[1]
    DGvz = DGst[2]
    # Data Reading
    Datafile = DataFolder + stat + "_VEL.pickle"
    print(f" Reading Data file {Datafile}")
    Datast = obspy.read(Datafile)
    t = DGvx.stats.starttime
    Datavx = Datast[0].slice(t,t+40)
    Datavy = Datast[1].slice(t,t+40)
    Datavz = Datast[2].slice(t,t+40)
    di -= 1.1
    
    OBStimex = Datavx.times(reftime=DGvx.stats.starttime)
    SYNtimex = DGvx.times()
    maxv =  max([abs(ele) for ele in Datast.max()])
    Datavx.normalize(norm=maxv)
    DGvx.normalize(norm=maxv)
    OBSvx = Datavx.filter("lowpass",freq=lowcut).data
    SYNvx = DGvx.filter("lowpass",freq=lowcut).data
    OBStimey = Datavy.times(reftime=DGvy.stats.starttime)
    SYNtimey = DGvy.times()
    Datavy.normalize(norm=maxv)
    DGvy.normalize(norm=maxv)
    OBSvy = Datavy.filter("lowpass",freq=lowcut).data
    SYNvy = DGvy.filter("lowpass",freq=lowcut).data
    OBStimez = Datavz.times(reftime=DGvz.stats.starttime)
    SYNtimez = DGvz.times()
    Datavz.normalize(norm=maxv)
    DGvz.normalize(norm=maxv)
    OBSvz = Datavz.filter("lowpass",freq=lowcut).data
    SYNvz = DGvz.filter("lowpass",freq=lowcut).data

    
    # Velocity Comparison
    fig = plt.figure(1)
    plt.title(" Amatrice ")    
    plt.plot(OBStimey+tiy,OBSvy+di,color='r')
    plt.plot(SYNtimey+tiy,SYNvy+di,color='k')
    plt.plot(OBStimex+tix,OBSvx+di,color='r')
    plt.plot(SYNtimex+tix,SYNvx+di,color='k')
    plt.plot(OBStimez+tiz,OBSvz+di,color='r')
    plt.plot(SYNtimez+tiz,SYNvz+di,color='k')
    plt.text(4, di, stat,fontsize=10,fontweight='bold')
    plt.text(2, -1, 'Station',fontsize=10,fontweight='bold')
    plt.text(40, -1, ' Vy',fontsize=10,fontweight='bold')
    plt.text(90, -1, ' Vx',fontsize=10,fontweight='bold')
    plt.text(145, -1, ' Vz',fontsize=10,fontweight='bold')
    plt.text(175, -0.5, ' max vel ',fontsize=10,fontweight='bold')
    plt.text(175, -1, ' (cm/s)',fontsize=10,fontweight='bold')
    plt.text(177, di, '{:5.2f}'.format(maxv),fontsize=10,fontweight='bold')
    plt.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
    plt.legend(['Data','DG'],loc='lower center')
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    













