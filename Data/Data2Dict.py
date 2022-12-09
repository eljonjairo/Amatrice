
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Generate Python Dict with Seismic Data
#
# John Diaz December 2022

import os
import warnings

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

from IPython import get_ipython
get_ipython().magic('reset -sf')
warnings.filterwarnings("ignore")
plt.close('all')

os.system('clear')

# Data folder
DataFolder = Path('ASCII')

Stats = [ "AMT","NRC","CSC","RM33","MSC","SPD","LSS","PZI1","TERO","FEMA","ANT"
         ,"TRL","ASP","SNO","SPM","MNF","TRE","CLF","FOS"]

nstats = int(len(Stats))
nt = 8000
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
    dataini = (0 for element in range(nt))
    dataList = list(dataini)
    data = np.zeros([nt])
    teven = float(lines[3].split(':')[1])
    tini = float(lines[26].split(':')[1])-20160824000000
    dtSYN = float(lines[28].split(':')[1])
    ntSYN = int(lines[29].split(':')[1])-1
    ntini = abs(int((tini-teven)//dtSYN))
    ntfin = nt-ntini
    synthe = list(lines[64:-1])
    dataList.insert(ntini,synthe)
    del dataList[nt:]
    data = np.asarray(dataList[:nt], dtype=np.float64) #.astype(float)/100 
    time = np.linspace(0, dtSYN*(nt-1), nt)
    print(teven)
    print(tini)
    di = (istat-1)*0.1
    #fig = plt.figure(1)
    #plt.plot(time,data,color='r')
    