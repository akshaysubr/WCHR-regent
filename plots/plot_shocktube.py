import matplotlib.pyplot as plt
import numpy
import h5py
import sys

if len(sys.argv) != 3:
  print "Usage:"
  print "  python %s <coordinate file> <data file>" % sys.argv[0]
  sys.exit()

coord_file = sys.argv[1]
data_file  = sys.argv[2]

f = h5py.File(coord_file)
x = numpy.squeeze(f['y_c'].value)
f.close()

f = h5py.File(data_file)
rho = numpy.squeeze(f['rho'].value)
u   = numpy.squeeze(f['v'].value)
p   = numpy.squeeze(f['p'].value)
f.close()

plt.figure()
plt.plot(x,rho,'ok',fillstyle='none')
plt.xlim((-0.5,0.5))
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$\rho$', fontsize=20)
plt.show(block=False)

plt.figure()
plt.plot(x,u,'ok',fillstyle='none')
plt.xlim((-0.5,0.5))
plt.ylim((-0.1,1.0))
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$u$', fontsize=20)
plt.show(block=False)

plt.figure()
plt.plot(x,p,'ok',fillstyle='none')
plt.xlim((-0.5,0.5))
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$p$', fontsize=20)
plt.show()

