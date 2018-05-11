import matplotlib
matplotlib.use('Agg')
import yt
import numpy as np
import sys
import glob
import matplotlib.pyplot as plt
from yt.units import G, kboltz, mh
from yt.units.yt_array import YTQuantity, YTArray
yt.enable_parallelism()
Mu = 3.0
CENTRE = [0.5, 0.5, 0.5]

AccretionRadius = 2
BASE = "./"
fns = glob.glob(BASE + "DD001*/DD001*.hierarchy")
fns.sort()

WIDTH = 8000


ActiveParticle = "AccretingParticle"


for f in fns:
    ff = f
    ds = yt.load(ff)
    print("Stats = ", ds.index.get_smallest_dx().in_units("au"))
    dx = ds.index.get_smallest_dx().in_units("au")
    print("dx = ", dx.d)
    dd = ds.all_data()
    centre = CENTRE
    print("Time = ", ds.current_time.in_units("yr"))
    #Now create a slice of the density field through the sink particle
    NumAPs = 0
    flist = dir(ds.fields)
    ap = flist[0]
    try:
        NumAPs = len(dd[ap, "level"])
        print("Number of %s active particles = %d" % (flist[0], NumAPs))
    except:
        print("No active particles found")
    
    slc = yt.SlicePlot(ds, "z", "density", center=centre, width=(WIDTH, 'au'))
    if(NumAPs > 0):
        for i in range(NumAPs):
            pos = dd[ActiveParticle, "particle_position"][i]
            slc.annotate_sphere(pos, radius=(dx.d*AccretionRadius, 'au'),
                                circle_args={'color':'black', 'fill':True})
    slc.annotate_timestamp(corner='upper_left', redshift=False, draw_inset_box=True)
    slc.save("Slice_z_%s.png" % (ds))
    


