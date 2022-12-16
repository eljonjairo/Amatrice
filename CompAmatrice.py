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

Stats = [ "AMT","NRC","CSC","RM33","MSC"]

# Filtering parameters
lowcut  = 0.5                              # low pass cut frecuency
highcut = 0.06                             # cut high pass frecuency

# initial time for plotting components (s)
tiy = 20;  tix = 65;  tiz = 110;

print("  ")
print(" START PROGRAM ")
print("  ")

di=0

for stat in Stats:
    # DG Synthetic Data Reading
    DGfile = DGFolder + stat + "_DGVEL.pickle"
    print(f" Reading DG file {DGfile}")
    DGst = obspy.read(DGfile)    
    DGx = DGst[0]
    DGy = DGst[1]
    DGz = DGst[2]
    # Data Reading
    Datafile = DataFolder + stat + "_VEL.pickle"
    print(f" Reading Data file {Datafile}")
    Datast = obspy.read(Datafile)
    Datax = Datast[0]
    Datay = Datast[1]
    Dataz = Datast[2]
    di -= max(abs(np.array(Datast.max())))*0.8
    
    OBStimex = Datax.times(reftime=DGx.stats.starttime)
    SYNtimex = DGx.times()
    OBSvx = di + Datax.filter("lowpass",freq=lowcut).data
    SYNvx = di + DGx.filter("lowpass",freq=lowcut).data
    OBStimey = Datay.times(reftime=DGy.stats.starttime)
    SYNtimey = DGy.times()
    OBSvy = di + Datay.filter("lowpass",freq=lowcut).data
    SYNvy = di + DGy.filter("lowpass",freq=lowcut).data
    OBStimez = Dataz.times(reftime=DGz.stats.starttime)
    SYNtimez = DGz.times()
    OBSvz = di + Dataz.filter("lowpass",freq=lowcut).data
    SYNvz = di + DGz.filter("lowpass",freq=lowcut).data
        
    plt.plot(OBStimey+tiy,OBSvy,color='r')
    plt.plot(SYNtimey+tiy,SYNvy,color='k')
    plt.plot(OBStimex+tix,OBSvx,color='r')
    plt.plot(SYNtimex+tix,SYNvx,color='k')
    plt.plot(OBStimez+tiz,OBSvz,color='r')
    plt.plot(SYNtimez+tiz,SYNvz,color='k')
    plt.text(-2, di, stat,fontsize=10,fontweight='bold')
    plt.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)















