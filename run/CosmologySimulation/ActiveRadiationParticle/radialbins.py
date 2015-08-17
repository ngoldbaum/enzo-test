from matplotlib import use; use('Agg')
import yt
import matplotlib.pyplot as plt
import numpy as np


first = 0
last = 4
center = [0.5]*3


########################################################################
# Radial profiles for some outputs
########################################################################
YFields = []
Fields = ["Density", "Temperature", "H_fraction", "H_p1_fraction", "He_p1_fraction", 
          "He_p2_fraction", "HI_kph", "HeI_kph", "HeII_kph"]
FieldsFigure = {}
FieldsAxes = {}
for elem in Fields:
    #Setup Figures
    tmpfigure = plt.figure()
    FieldsFigure[elem] = tmpfigure
    FieldsAxes[elem] = tmpfigure.add_subplot(111)
    YFields.append(elem)

print "Yfields = ", YFields
outputs = [1,2,3,4]

for outp in outputs:
    amrfile = "RD%4.4d/RD%4.4d" % (outp, outp)
    print "Load file %s" % (amrfile)
    ds = yt.load(amrfile)
    sphere = ds.h.sphere(center, (250, 'kpc'))
    rp = yt.create_profile(sphere, 'radius',  YFields, n_bins=32,
                           units = {'radius': 'kpc', 'HI_kph' : '1/s',
                                    'HeI_kph' : '1/s','HeII_kph' : '1/s'},
                           weight_field='density')
    
    for elem in YFields:
            FieldsAxes[elem].semilogy(rp.x.value, rp[elem], 
                                      label="z = %.2f" % (ds.current_redshift))
            FieldsAxes[elem].set_xlabel("Radius (kpc)")
            FieldsAxes[elem].set_ylabel(elem)
            FieldsAxes[elem].legend(loc='best')
            
for elem in YFields:
    FieldsFigure[elem].savefig("%s_RadialProfile.png" % (elem))
    
