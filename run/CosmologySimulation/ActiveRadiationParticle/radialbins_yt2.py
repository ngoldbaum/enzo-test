from matplotlib import use; use('Agg')
from yt.mods import *
import matplotlib.pyplot as plt
import numpy as np


first = 0
last = 4
center = [0.5]*3


########################################################################
# Radial profiles for some outputs
########################################################################
YFields = []
Fields = ["Density", "Temperature", "HI_Fraction", "HII_Fraction", "HeI_Fraction",
          "HeII_Fraction", "HI_kph", "HeI_kph", "HeII_kph", "Electron_Fraction"]
FieldsFigure = {}
FieldsAxes = {}
for elem in Fields:
    #Setup Figures
    tmpfigure = plt.figure()
    FieldsFigure[elem] = tmpfigure
    FieldsAxes[elem] = tmpfigure.add_subplot(111)
    YFields.append(elem)

print "Yfields = ", YFields
outputs = [1,2,3]

for outp in outputs:
    amrfile = "RD%4.4d/RD%4.4d" % (outp, outp)
    print "Load file %s" % (amrfile)
    ds = load(amrfile)
    sphere = ds.h.sphere(center, (250, 'kpc'))
    for elem in YFields:
        rp = ProfilePlot(sphere, 'Radiuskpc',  elem, n_bins=32,
                         weight_field='Density')
        rp.save()
