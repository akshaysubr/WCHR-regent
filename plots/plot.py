import matplotlib.pyplot as plt
import numpy
import h5py

f = h5py.File('../tests/cell_coords.h5')

x = f['x_c'].value
y = f['y_c'].value
z = f['z_c'].value

dx = x[0,0,1] - x[0,0,0]

f = h5py.File('../tests/cell_primitive0000.h5')
rho = f['rho'].value
u   = f['u'].value
v   = f['v'].value
w   = f['w'].value
p   = f['p'].value
f.close()

f = h5py.File('../tests/edge_primitive_l_x0000.h5')
rho_l = f['rho'].value
u_l   = f['u'].value
v_l   = f['v'].value
w_l   = f['w'].value
p_l   = f['p'].value
f.close()

f = h5py.File('../tests/rhs_l_x0000.h5')
rho_rhs = f['rho'].value
f.close()

plt.plot(x[0,0,:], rho[0,0,:], 'ok-')

x_e = numpy.hstack((x[0,0,:]-0.5*dx, x[0,0,-1]+0.5*dx))
plt.plot(x_e, rho_l[0,0,:], 'or-')
plt.plot(x_e, rho_rhs[0,0,:], 'ob-')
plt.show()
