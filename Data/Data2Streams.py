
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

# Data folder
DataFolder = Path('ASCII')

# Object Folder
StFolder = "Streams/"

Stats = [ "AMT","NRC","CSC","RM33","MSC","SPD","LSS","PZI1","TERO","FEMA","ANT"
         ,"TRL","ASP","SNO","SPM","MNF","TRE","CLF","FOS"]

nstats = int(len(Stats))

# Field DISP: displacement, VEL: velocity, ACC: acceleration
field = "ACC"


def File2Trace(file,station,comp):
    tr = obspy.Trace()
    tr.stats.station = station
    with open(file, 'r') as f: lines = f.readlines()
    stime = lines[26].split(':')[1].split('_')
    time = stime[0]+'T'+stime[1]
    tr.stats.starttime = UTCDateTime(time)
    tr.stats.delta     = float(lines[28].split(':')[1])
    tr.stats.npts      = int(lines[29].split(':')[1])
    tr.stats.channel   = comp
    tr.data            = np.array(lines[64:-1],dtype=float)
 
    return tr

  
for istat in range (0,nstats):
    # Synthetic Data Proccesing
    print(f" Station {Stats[istat]} ")
    st = Stream()
    stfile = StFolder+Stats[istat]+"_"+field+".pickle" 
    for comp in ["E","N","Z"]:
 
        print(f"    Component: { comp } Field: { field } ")
        preStat = '/IT.'
        posStat = "..HG"+ comp +".D.EMSC-20160824_0000006."+field+".MP.ASC"
        if Stats[istat] == 'RM33' or Stats[istat] == 'TERO' or Stats[istat] == 'FEMA':
            preStat = '/IV.'
            posStat = "..HN"+ comp +".D.EMSC-20160824_0000006."+field+".MP.ASC"
    
        DataFileV = DataFolder.joinpath( Stats[istat] + preStat + Stats[istat] + posStat)
        print(f"    Loading file : {DataFileV} ")
        tr = File2Trace(DataFileV, Stats[istat], comp)
        tr.plot()
    
        st.append(tr)
        
    print(f" Writing Stream in file: {stfile} ")
    st.write(stfile, format = "PICKLE")
