import matplotlib.pyplot as plt
import numpy

coordsfile = '../src/interpolation_x_c_coords_px0000_pz0000.dat'
edgecoordsfile = '../src/interpolation_x_l_coords_px0000_pz0000.dat'
f = open(coordsfile)
line0 = f.readline()
line1 = f.readline()
f.close()

nx = int(line1.split()[0]) - int(line0.split()[0]) + 1
ny = int(line1.split()[1]) - int(line0.split()[1]) + 1
nz = int(line1.split()[2]) - int(line0.split()[2]) + 1

x, y, z = numpy.loadtxt(coordsfile, skiprows=2, unpack=True)
x = x.reshape((nx,ny,nz), order='F')
y = y.reshape((nx,ny,nz), order='F')
z = z.reshape((nx,ny,nz), order='F')

dx = x[1,0,0] - x[0,0,0]
dy = y[0,1,0] - y[0,0,0]
dz = z[0,0,1] - z[0,0,0]

cellfile = '../src/interpolation_x_c_0000_px0000_pz0000.dat'
rho, u, v, w, p = numpy.loadtxt(cellfile, skiprows=2, unpack=True)
rho = rho.reshape((nx,ny,nz), order='F')
u = u.reshape((nx,ny,nz), order='F')
v = v.reshape((nx,ny,nz), order='F')
w = w.reshape((nx,ny,nz), order='F')
p = p.reshape((nx,ny,nz), order='F')

edgefile_l_x = '../src/interpolation_x_l_0000_px0000_pz0000.dat'
rho_l_x, u_l_x, v_l_x, w_l_x, p_l_x = numpy.loadtxt(edgefile_l_x, skiprows=2, unpack=True)
rho_l_x = rho_l_x.reshape((nx-6+1,ny,nz), order='F')
u_l_x = u_l_x.reshape((nx-6+1,ny,nz), order='F')
v_l_x = v_l_x.reshape((nx-6+1,ny,nz), order='F')
w_l_x = w_l_x.reshape((nx-6+1,ny,nz), order='F')
p_l_x = p_l_x.reshape((nx-6+1,ny,nz), order='F')

edgefile_r_x = '../src/interpolation_x_r_0000_px0000_pz0000.dat'
rho_r_x, u_r_x, v_r_x, w_r_x, p_r_x = numpy.loadtxt(edgefile_r_x, skiprows=2, unpack=True)
rho_r_x = rho_r_x.reshape((nx-6+1,ny,nz), order='F')
u_r_x = u_r_x.reshape((nx-6+1,ny,nz), order='F')
v_r_x = v_r_x.reshape((nx-6+1,ny,nz), order='F')
w_r_x = w_r_x.reshape((nx-6+1,ny,nz), order='F')
p_r_x = p_r_x.reshape((nx-6+1,ny,nz), order='F')

# edgefile_l_y = '../src/interpolation_y_l_0000_px0000_pz0000.dat'
# rho_l_y, u_l_y, v_l_y, w_l_y, p_l_y = numpy.loadtxt(edgefile_l_y, skiprows=2, unpack=True)
# rho_l_y = rho_l_y.reshape((nx,ny+1,nz), order='F')
# u_l_y = u_l_y.reshape((nx,ny+1,nz), order='F')
# v_l_y = v_l_y.reshape((nx,ny+1,nz), order='F')
# w_l_y = w_l_y.reshape((nx,ny+1,nz), order='F')
# p_l_y = p_l_y.reshape((nx,ny+1,nz), order='F')
# 
# edgefile_r_y = '../src/interpolation_y_r_0000_px0000_pz0000.dat'
# rho_r_y, u_r_y, v_r_y, w_r_y, p_r_y = numpy.loadtxt(edgefile_r_y, skiprows=2, unpack=True)
# rho_r_y = rho_r_y.reshape((nx,ny+1,nz), order='F')
# u_r_y = u_r_y.reshape((nx,ny+1,nz), order='F')
# v_r_y = v_r_y.reshape((nx,ny+1,nz), order='F')
# w_r_y = w_r_y.reshape((nx,ny+1,nz), order='F')
# p_r_y = p_r_y.reshape((nx,ny+1,nz), order='F')
# 
# edgefile_l_z = '../src/interpolation_z_l_0000_px0000_pz0000.dat'
# rho_l_z, u_l_z, v_l_z, w_l_z, p_l_z = numpy.loadtxt(edgefile_l_z, skiprows=2, unpack=True)
# rho_l_z = rho_l_z.reshape((nx,ny,nz+1), order='F')
# u_l_z = u_l_z.reshape((nx,ny,nz+1), order='F')
# v_l_z = v_l_z.reshape((nx,ny,nz+1), order='F')
# w_l_z = w_l_z.reshape((nx,ny,nz+1), order='F')
# p_l_z = p_l_z.reshape((nx,ny,nz+1), order='F')
# 
# edgefile_r_z = '../src/interpolation_z_r_0000_px0000_pz0000.dat'
# rho_r_z, u_r_z, v_r_z, w_r_z, p_r_z = numpy.loadtxt(edgefile_r_z, skiprows=2, unpack=True)
# rho_r_z = rho_r_z.reshape((nx,ny,nz+1), order='F')
# u_r_z = u_r_z.reshape((nx,ny,nz+1), order='F')
# v_r_z = v_r_z.reshape((nx,ny,nz+1), order='F')
# w_r_z = w_r_z.reshape((nx,ny,nz+1), order='F')
# p_r_z = p_r_z.reshape((nx,ny,nz+1), order='F')

# y_e = numpy.hstack((y[0,:,0]-0.5*dy, y[0,-1,0]+0.5*dy))
x_e = numpy.loadtxt(edgecoordsfile, skiprows=2, unpack=True, usecols=(1,))
x_e = x_e.reshape((nx-6+1,ny,nz), order='F')

plt.figure()
plt.plot(x[:,0,0], rho[:,0,0], 'k-')

plt.plot(x_e[:,0,0], rho_l_x[:,0,0], 'b--',  fillstyle='none', label=r'$X-left $')
plt.plot(x_e[:,0,0], rho_r_x[:,0,0], 'r--',  fillstyle='none', label=r'$X-right$')

# plt.plot(x_e, rho_l_y[0,:,0], '^b',  fillstyle='none', label=r'$Y-left $')
# plt.plot(x_e, rho_r_y[0,:,0], '^r',  fillstyle='none', label=r'$Y-right$')
# 
# plt.plot(x_e, rho_l_z[0,0,:], 'ob',  fillstyle='none', label=r'$Z-left $')
# plt.plot(x_e, rho_r_z[0,0,:], 'or',  fillstyle='none', label=r'$Z-right$')

plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$\rho$', fontsize=20)
plt.title(r'Characteristic Interpolation', fontsize=20)
plt.legend(loc='upper right')
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
