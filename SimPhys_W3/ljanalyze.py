from numpy import *
from matplotlib.pyplot import *
import sys, os
import pickle

# get filename on command line
if len(sys.argv) != 2:
    print "Usage: python %s FILE" % sys.argv[0]
    sys.exit(2)
datafilename = sys.argv[1]

# check whether data file exists
if not os.path.exists(datafilename):
    print "ERROR: %s doesn't exist."
    sys.exit(1)

print "Reading data from %s." % datafilename
datafile = open(datafilename, 'r')
ts, Es = pickle.load(datafile)
datafile.close()

ts = array(ts)
Es = array(Es)

plot(ts, Es)
show()