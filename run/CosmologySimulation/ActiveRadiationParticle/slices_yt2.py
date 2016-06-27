from matplotlib import use; use('Agg')
from yt.mods import *
import matplotlib.pyplot as plt
import numpy as np

first = 0
last = 3
center = [0.5]*3
Width = 500.0 #in kpc

########################################################################
# Radial profiles for some outputs
########################################################################
Fields = ["Density", "Temperature", "HI_Fraction", "HII_Fraction", "HeI_Fraction", 
          "HeII_Fraction", "HI_kph", "HeI_kph", "HeII_kph", "Electron_Fraction"]

outputs = [1,2,3]

for outp in outputs:
    amrfile = "RD%4.4d/RD%4.4d" % (outp, outp)
    print "Load file %s" % (amrfile)
    ds = load(amrfile)
    
    for field in Fields:
        slc = SlicePlot(ds, 1, field, center=center, width=(Width, 'kpc'))
        slc.save("%s_slice_%d.png" % (field, outp))
