import matplotlib.pyplot as plt
import numpy
import h5py
import sys

if len(sys.argv) != 3:
  print "Usage:"
  print "  python %s <data file> <linetype>" % sys.argv[0]
  sys.exit()

data_file  = sys.argv[1]
linetype   = sys.argv[2]

t   = []
TKE = []

f = open(data_file)

for line in f:
    line = line.strip().split()
    if line and line[0].isdigit():
        t.append(float(line[2]))
        TKE.append(float(line[12]))

f.close()

t = numpy.array(t)
TKE = numpy.array(TKE)

plt.figure(1)
plt.plot(t,TKE,linetype,fillstyle='none')
plt.xlabel(r'$t$', fontsize=20)
plt.ylabel(r'$TKE$', fontsize=20)