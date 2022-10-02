
#
# Generate Fault and Topography output files
#
# John Diaz July 2022

# To Do List:
#    Add index from hypocenter coordinates
#    Test python dealunay triangulation.
#

import os
import warnings

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.spatial import Delaunay
import pickle

import pygmt
import utm

warnings.filterwarnings("ignore")
os.system('clear')


dhF = 0.5  # Output subfaults size in Km
infile = 'Input/s2016AMATRI01PIZZ.fsp'  # Input file with fault data
name = "AmatricePizzi2017_dhF"  # Output files name

# Dimension of input fault
nstkin = 30
ndipin = 12

# X and Y UTM limits (Km) for Model
xmin = 290.0
xmax = 410.0
ymin = 4670.0
ymax = 4790.0
zmin = -60.0

# Hypocenter coordinates (Km) 
hypolon = 13.25
hypolat = 42.70
hypoz =-5.0

print("  ")
print(" START PROGRAM ")
print("  ")








