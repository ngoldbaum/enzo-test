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

DO_MASSPLOT = True

Mu = 3.0
CENTRE = [0.5, 0.5, 0.5]

fns = glob.glob(BASE + "DD*/*.hierarchy")
fns.sort()
#print fns
Files = fns
Start = 0
End = 63
step = 1
Files = []
for i in range(Start, End, step):
    Files.append("DD%04d/DD%04d" % (i, i))

CT   = {}

TimeArray = []
MassDict = {}
for f in Files:
    NumAPs = 0
    ff = BASE + f
    ds = yt.load(ff)
    dd = ds.all_data()
    flist = dir(ds.fields)
    ap = flist[0]
    
    try:
        NumAPs = len(dd[ap, "level"])
        print "Number of %s active particles = %d" % (flist[0], NumAPs)
    except:
        print "No active particles found"
    if(NumAPs > 0):
        print(dir(ds.fields.SmartStar))
        for i in range(NumAPs):
            print "\n%s #%d" % (ap, i)
            print "Particle Position = ", dd["SmartStar", "particle_position"][i]
            print "ID = ", dd["SmartStar", "identifier"][i]
            print "Mass = ", dd["SmartStar", "particle_mass"][i].in_units("Msun")
            print "Accretion Rate = ", dd["SmartStar", "AccretionRate"][i]
            print "Creation Time = ", dd["SmartStar", "creation_time"][i].in_units("kyr")
            print "Age = ", dd["SmartStar", "age"][i].in_units("kyr")
            print "Refinement Level = ", dd["SmartStar", "level"][i]
            print "Metallicity = ", dd["SmartStar", "metallicity"][i]
            Time = dd["SmartStar", "creation_time"][i].in_units("kyr").d
            Mass = dd["SmartStar", "particle_mass"][i].in_units("Msun").d
            if Time not in CT.keys():
                #print dd["SmartStar", "creation_time"][i]
                CT[float(Time)] = Mass
            #if(Time == 42.2660115218954):
            #    print "Mass = ", Mass

if DO_MASSPLOT:
    for time in CT.keys():
        print "\n\n\n\nLooking for AP with CT = ", time
        MassArray = []
        AgeArray = []
        for f in Files:
            NumAPs = 0
            ff = BASE + f
            ds = yt.load(ff)
            dd = ds.all_data()
            flist = dir(ds.fields)
            ap = flist[0] 
        
            try:
                NumAPs = len(dd[ap, "level"])
                print "Number of %s active particles = %d" % (flist[0], NumAPs)
            except:
                print "No active particles found"
            if(NumAPs > 0):
                for i in range(NumAPs):
                    if(time == dd["SmartStar", "creation_time"][i].in_units("kyr").d):
                        print "Creation Time = ", dd["SmartStar", "creation_time"][i].in_units("kyr")
                        print "Mass = ", dd["SmartStar", "particle_mass"][i].in_units("Msun")
                        print "Age = ", dd["SmartStar", "age"][i].in_units("kyr")
                        Mass = dd["SmartStar", "particle_mass"][i].in_units("Msun").d
                        MassArray.append(float(Mass))
                        Age = time + dd["SmartStar", "age"][i].in_units("kyr").d
                        AgeArray.append(float(Age))

                        print "Mass = ", MassArray
                        print AgeArray
        if(len(MassArray) > 5):
            plt.plot(AgeArray, MassArray, label="T = %1.3f" % time, marker='x')
            plt.legend(loc = "upper left")
            plt.xlabel("Time [kyrs]")
            plt.ylabel("Mass [Msun]")
            plt.savefig("Mass.png")
