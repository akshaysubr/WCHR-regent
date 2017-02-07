import matplotlib.pyplot as plt
import numpy
import h5py

f = h5py.File('../src/cell_coords.h5')

x = f['x_c'].value
y = f['y_c'].value
z = f['z_c'].value

dx = x[0,0,1] - x[0,0,0]

f = h5py.File('../src/cell_primitive.h5')
rho = f['rho'].value
u   = f['u'].value
v   = f['v'].value
w   = f['w'].value
p   = f['p'].value

f = h5py.File('../src/edge_primitive_l_x.h5')
rho_l = f['rho'].value
u_l   = f['u'].value
v_l   = f['v'].value
w_l   = f['w'].value
p_l   = f['p'].value

plt.plot(x[0,0,:], rho[0,0,:], 'ok-')
plt.plot(x[0,0,:]-dx/2., rho_l[0,0,:-1], 'or')
plt.show()
