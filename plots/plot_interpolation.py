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

x_e = numpy.hstack((x[0,0,:]-0.5*dx, x[0,0,-1]+0.5*dx))

f = h5py.File('../tests/edge_primitive_l_x0000.h5')
rho_l_x = f['rho'].value
u_l_x   = f['u'].value
v_l_x   = f['v'].value
w_l_x   = f['w'].value
p_l_x   = f['p'].value
f.close()

f = h5py.File('../tests/edge_primitive_r_x0000.h5')
rho_r_x = f['rho'].value
u_r_x   = f['u'].value
v_r_x   = f['v'].value
w_r_x   = f['w'].value
p_r_x   = f['p'].value
f.close()

f = h5py.File('../tests/edge_primitive_l_y0000.h5')
rho_l_y = f['rho'].value
u_l_y   = f['u'].value
v_l_y   = f['v'].value
w_l_y   = f['w'].value
p_l_y   = f['p'].value
f.close()

f = h5py.File('../tests/edge_primitive_r_y0000.h5')
rho_r_y = f['rho'].value
u_r_y   = f['u'].value
v_r_y   = f['v'].value
w_r_y   = f['w'].value
p_r_y   = f['p'].value
f.close()

plt.figure()
plt.plot(x[0,0,:], rho[0,0,:], 'ok-')
plt.plot(x_e, rho_l_x[0,0,:], 'ob-', fillstyle='none')
plt.plot(x_e, rho_r_x[0,0,:], 'or-', fillstyle='none')
plt.plot(x_e, rho_l_y[0,:,0], 'sb', fillstyle='none')
plt.plot(x_e, rho_r_y[0,:,0], 'sr', fillstyle='none')
plt.show()

# plt.figure()
# plt.plot(x[0,0,:], u[0,0,:], 'ok-')
# plt.plot(x_e, u_l[0,0,:], 'ob-')
# plt.show()
# 
# plt.figure()
# plt.plot(x[0,0,:], v[0,0,:], 'ok-')
# plt.plot(x_e, v_l[0,0,:], 'ob-')
# plt.show()
# 
# plt.figure()
# plt.plot(x[0,0,:], w[0,0,:], 'ok-')
# plt.plot(x_e, w_l[0,0,:], 'ob-')
# plt.show()
# 
# plt.figure()
# plt.plot(x[0,0,:], p[0,0,:], 'ok-')
# plt.plot(x_e, p_l[0,0,:], 'ob-')
# plt.show()
