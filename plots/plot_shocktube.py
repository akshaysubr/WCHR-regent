import matplotlib.pyplot as plt
import numpy
import h5py

f = h5py.File('../src/cell_coords.h5')
x = f['x_c'].value
f.close()

f = h5py.File('../src/cell_primitive0400.h5')
rho = f['rho'].value
u = f['u'].value
p = f['p'].value
f.close()

plt.figure()
plt.plot(x[0,0,:],rho[0,0,:],'ok',fillstyle='none')
plt.xlim((-0.5,0.5))
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$\rho$', fontsize=20)
plt.show(block=False)

plt.figure()
plt.plot(x[0,0,:],u[0,0,:],'ok',fillstyle='none')
plt.xlim((-0.5,0.5))
plt.ylim((-0.1,1.0))
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$u$', fontsize=20)
plt.show(block=False)

plt.figure()
plt.plot(x[0,0,:],p[0,0,:],'ok',fillstyle='none')
plt.xlim((-0.5,0.5))
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$p$', fontsize=20)
plt.show()

