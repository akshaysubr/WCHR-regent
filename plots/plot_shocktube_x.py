import matplotlib.pyplot as plt
import numpy
import h5py
import sys

if len(sys.argv) != 4:
  print "Usage:"
  print "  python %s <coordinate file> <data file> <linetype>" % sys.argv[0]
  sys.exit()

coord_file = sys.argv[1]
data_file  = sys.argv[2]
linetype   = sys.argv[3]

f = h5py.File(coord_file)
x = (f['x_c'].value)[0,0,:]
f.close()

f = h5py.File(data_file)
rho = (f['rho'].value)[0,0,:]
u   = (f['u'].value)[0,0,:]
p   = (f['p'].value)[0,0,:]
f.close()

plt.figure(1)
plt.plot(x,rho,linetype,fillstyle='none')
plt.xlim((-0.5,0.5))
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$\rho$', fontsize=20)
# plt.show(block=False)

plt.figure(2)
plt.plot(x,u,linetype,fillstyle='none')
plt.xlim((-0.5,0.5))
plt.ylim((-0.1,1.0))
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$v$', fontsize=20)
# plt.show(block=False)

plt.figure(3)
plt.plot(x,p,linetype,fillstyle='none')
plt.xlim((-0.5,0.5))
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$p$', fontsize=20)
# plt.show()

