import "regent"

local c       = regentlib.c
local cmath   = terralib.includec("math.h")
local PI      = cmath.M_PI

require("fields")
require("IO")
require("partition")
local interpolation = require("interpolation")

-- Grid dimensions
local NX = 64
local NY = 64
local NZ = 64

-- Domain size
local LX = 2.0*math.pi
local LY = 2.0*math.pi
local LZ = 2.0*math.pi

local X1 = 0.0
local Y1 = 0.0
local Z1 = 0.0

-- Grid spacing
local DX = LX / NX
local DY = LY / NY
local DZ = LZ / NZ

local ONEBYDX = 1.0 / DX
local ONEBYDY = 1.0 / DY
local ONEBYDZ = 1.0 / DZ

local periodic_x = true
local periodic_y = true
local periodic_z = true

-- Make the node and midpoint-node differencing tasks (Using pentadiagonal solver for this instead of tridiagonal solver)
require("derivatives")
local alpha10d1 = 1.0/3.0
local beta10d1  = 0.0
local a06d1 = ( 14.0/ 9.0)/2.0
local b06d1 = (  1.0/ 9.0)/4.0
local c06d1 = (  0.0/100.0)/6.0

-- Compact MND finite difference
-- local alpha06MND = -1.0/12.0
-- local beta06MND  = 0.0
-- local a06MND = 16.0/9.0
-- local b06MND = (-17.0/18.0)/2.0
-- local c06MND = (0.0)/3.0

-- Compact staggered finite difference
local alpha06MND = 9.0/62.0
local beta06MND  = 0.0
local a06MND = 63.0/62.0
local b06MND = (0.0/18.0)/2.0
local c06MND = (17.0/62.0)/3.0

local r_prim       = regentlib.newsymbol(region(ispace(int3d), primitive), "r_prim")
local r_prim_e     = regentlib.newsymbol(region(ispace(int3d), primitive), "r_prim_e")
local r_derivative = regentlib.newsymbol(region(ispace(int3d), primitive), "r_derivative")

-- local ddx     = make_ddx(r_prim, "rho", r_derivative, "rho", NX, NY, NZ, ONEBYDX, a06d1, b06d1, c06d1)
-- local ddy     = make_ddy(r_prim, "rho", r_derivative, "rho", NX, NY, NZ, ONEBYDX, a06d1, b06d1, c06d1)
-- local ddz     = make_ddz(r_prim, "rho", r_derivative, "rho", NX, NY, NZ, ONEBYDX, a06d1, b06d1, c06d1)

local ddx_MND = make_ddx_MND(r_prim, r_prim_e, "rho", r_derivative, "rho", NX, NY, NZ, ONEBYDX, periodic_x)
local ddy_MND = make_ddy_MND(r_prim, r_prim_e, "rho", r_derivative, "rho", NX, NY, NZ, ONEBYDY, periodic_y)
local ddz_MND = make_ddz_MND(r_prim, r_prim_e, "rho", r_derivative, "rho", NX, NY, NZ, ONEBYDZ, periodic_z)
-----------------------------------------------------

task initialize( coords     : region(ispace(int3d), coordinates),
                 coords_x   : region(ispace(int3d), coordinates),
                 coords_y   : region(ispace(int3d), coordinates),
                 coords_z   : region(ispace(int3d), coordinates),
                 r_prim_c   : region(ispace(int3d), primitive),
                 r_prim_l_x : region(ispace(int3d), primitive),
                 r_prim_l_y : region(ispace(int3d), primitive),
                 r_prim_l_z : region(ispace(int3d), primitive),
                 dx         : double,
                 dy         : double,
                 dz         : double,
                 n_ghosts   : int64 )
where
  reads writes(coords, coords_x, coords_y, coords_z, r_prim_c, r_prim_l_x, r_prim_l_y, r_prim_l_z)
do
  for i in coords.ispace do
    coords[i].x_c = X1 + (i.x - n_ghosts + 0.5) * dx
    coords[i].y_c = Y1 + (i.y - n_ghosts + 0.5) * dy
    coords[i].z_c = Z1 + (i.z - n_ghosts + 0.5) * dz

    r_prim_c[i].rho = cmath.sin(coords[i].x_c) * cmath.cos(coords[i].y_c) * cmath.cos(coords[i].z_c)
  end

  for i in r_prim_l_x do
    var x_c : double = X1 + (i.x + 0.0) * dx
    var y_c : double = Y1 + (i.y - n_ghosts + 0.5) * dy
    var z_c : double = Z1 + (i.z - n_ghosts + 0.5) * dz

    coords_x[i].x_c = x_c
    coords_x[i].y_c = y_c
    coords_x[i].z_c = z_c

    r_prim_l_x[i].rho = cmath.sin(x_c) * cmath.cos(y_c) * cmath.cos(z_c)
  end

  for i in r_prim_l_y do
    var x_c : double = X1 + (i.x - n_ghosts + 0.5) * dx
    var y_c : double = Y1 + (i.y + 0.0) * dy
    var z_c : double = Z1 + (i.z - n_ghosts + 0.5) * dz

    coords_y[i].x_c = x_c
    coords_y[i].y_c = y_c
    coords_y[i].z_c = z_c

    r_prim_l_y[i].rho = cmath.sin(x_c) * cmath.cos(y_c) * cmath.cos(z_c)
  end

  for i in r_prim_l_z do
    var x_c : double = X1 + (i.x - n_ghosts + 0.5) * dx
    var y_c : double = Y1 + (i.y - n_ghosts + 0.5) * dy
    var z_c : double = Z1 + (i.z + 0.0) * dz

    coords_z[i].x_c = x_c
    coords_z[i].y_c = y_c
    coords_z[i].z_c = z_c

    r_prim_l_z[i].rho = cmath.sin(x_c) * cmath.cos(y_c) * cmath.cos(z_c)
  end

  return 1
end

task check_ddx( coords : region(ispace(int3d), coordinates),
                r_der  : region(ispace(int3d), primitive) )
where
  reads writes(coords, r_der)
do
  var err : double = 0.0
  for i in coords.ispace do
    var err_t : double = cmath.fabs( r_der[i].rho - cmath.cos(coords[i].x_c) * cmath.cos(coords[i].y_c) * cmath.cos(coords[i].z_c) )
    if err_t > err then
      err = err_t
    end
  end
  return err
end

task check_ddy( coords : region(ispace(int3d), coordinates),
                r_der  : region(ispace(int3d), primitive) )
where
  reads writes(coords, r_der)
do
  var err : double = 0.0
  for i in coords.ispace do
    var err_t : double = cmath.fabs( r_der[i].rho + cmath.sin(coords[i].x_c) * cmath.sin(coords[i].y_c) * cmath.cos(coords[i].z_c) )
    if err_t > err then
      err = err_t
    end
  end
  return err
end

task check_ddz( coords : region(ispace(int3d), coordinates),
                r_der  : region(ispace(int3d), primitive) )
where
  reads writes(coords, r_der)
do
  var err : double = 0.0
  for i in coords.ispace do
    var err_t : double = cmath.fabs( r_der[i].rho + cmath.sin(coords[i].x_c) * cmath.cos(coords[i].y_c) * cmath.sin(coords[i].z_c) )
    if err_t > err then
      err = err_t
    end
  end
  return err
end

terra wait_for(x : int)
  return x
end

task main()
  var Nx : int64 = NX
  var Ny : int64 = NY
  var Nz : int64 = NZ

  var n_ghosts = interpolation.n_ghosts
  var Nx_g : int64 = NX + 2*n_ghosts
  var Ny_g : int64 = NY + 2*n_ghosts
  var Nz_g : int64 = NZ + 2*n_ghosts

  var Lx : double = LX
  var Ly : double = LY
  var Lz : double = LZ

  var dx : double = DX
  var dy : double = DY
  var dz : double = DZ

  c.printf("================ Problem parameters ================\n")
  c.printf("           Nx, Ny, Nz = %d, %d, %d\n", Nx, Ny, Nz)
  c.printf("           Lx, Ly, Lz = %f, %f, %f\n", Lx, Ly, Lz)
  c.printf("           dx, dy, dz = %f, %f, %f\n", dx, dy, dz)
  c.printf("====================================================\n")

  --------------------------------------------------------------------------------------------
  --                       DATA STRUCTURES
  --------------------------------------------------------------------------------------------
  var grid_c     = ispace(int3d, {x = Nx_g, y = Ny_g, z = Nz_g})  -- Cell center index space

  var grid_e_x   = ispace(int3d, {x = Nx+1, y = Ny_g, z = Nz_g})  -- x cell edge index space
  var grid_e_y   = ispace(int3d, {x = Nx_g, y = Ny+1, z = Nz_g})  -- y cell edge index space
  var grid_e_z   = ispace(int3d, {x = Nx_g, y = Ny_g, z = Nz+1})  -- z cell edge index space

  var coords     = region(grid_c,   coordinates)  -- Coordinates of cell center
  var coords_x   = region(grid_e_x, coordinates)  -- Coordinates of cell center
  var coords_y   = region(grid_e_y, coordinates)  -- Coordinates of cell center
  var coords_z   = region(grid_e_z, coordinates)  -- Coordinates of cell center

  var r_prim_c   = region(grid_c,   primitive)  -- Primitive variables at cell center
  var r_prim_l_x = region(grid_e_x, primitive)  -- Primitive variables at left x cell edge
  var r_prim_l_y = region(grid_e_y, primitive)  -- Primitive variables at left y cell edge
  var r_prim_l_z = region(grid_e_z, primitive)  -- Primitive variables at left z cell edge
  
  var r_der      = region(grid_c,   primitive)  -- RHS for time stepping at cell center

  var mat_x      = region(ispace(int3d, {x = Nx, y = 1+2, z = 1+2}), LU_coeffs) -- Data structure to hold x derivative LHS matrix
  var mat_y      = region(ispace(int3d, {x = Ny, y = 1+2, z = 1+2}), LU_coeffs) -- Data structure to hold y derivative LHS matrix
  var mat_z      = region(ispace(int3d, {x = Nz, y = 1+2, z = 1+2}), LU_coeffs) -- Data structure to hold z derivative LHS matrix

  var LU_x       = region(ispace(int3d, {x = Nx, y = 1+2, z = 1+2}), LU_struct) -- Data structure to hold x derivative LU decomposition
  var LU_y       = region(ispace(int3d, {x = Ny, y = 1+2, z = 1+2}), LU_struct) -- Data structure to hold y derivative LU decomposition
  var LU_z       = region(ispace(int3d, {x = Nz, y = 1+2, z = 1+2}), LU_struct) -- Data structure to hold z derivative LU decomposition

  var pencil          = ispace(int2d, int2d {1+2, 1+2}) -- All pencil partitions including the ghost pencils
  var pencil_interior = ispace(int2d, int2d {1,   1}, int2d {1, 1}) -- Only the interior pencil partitions

  var p_coords_x = partition_xpencil_coords(coords, n_ghosts, false, pencil)
  var p_coords_y = partition_ypencil_coords(coords, n_ghosts, false, pencil)
  var p_coords_z = partition_zpencil_coords(coords, n_ghosts, false, pencil)

  var p_prim_c_x = partition_xpencil_prim(r_prim_c, n_ghosts, true, pencil)
  var p_prim_c_y = partition_ypencil_prim(r_prim_c, n_ghosts, true, pencil)
  var p_prim_c_z = partition_zpencil_prim(r_prim_c, n_ghosts, true, pencil)

  var p_prim_l_x = partition_xpencil_prim(r_prim_l_x, n_ghosts, true, pencil)
  var p_prim_l_y = partition_ypencil_prim(r_prim_l_y, n_ghosts, true, pencil)
  var p_prim_l_z = partition_zpencil_prim(r_prim_l_z, n_ghosts, true, pencil)

  var p_der_x    = partition_xpencil_prim(r_der, n_ghosts, false, pencil)
  var p_der_y    = partition_ypencil_prim(r_der, n_ghosts, false, pencil)
  var p_der_z    = partition_zpencil_prim(r_der, n_ghosts, false, pencil)

  var p_LU_x     = partition_LU(LU_x, pencil)
  var p_LU_y     = partition_LU(LU_y, pencil)
  var p_LU_z     = partition_LU(LU_z, pencil)

  var p_mat_x    = partition_mat(mat_x, pencil)
  var p_mat_y    = partition_mat(mat_y, pencil)
  var p_mat_z    = partition_mat(mat_z, pencil)

  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  var token = initialize(coords, coords_x, coords_y, coords_z, r_prim_c, r_prim_l_x, r_prim_l_y, r_prim_l_z, dx, dy, dz, n_ghosts)
  wait_for(token)

  var err = 0.0

  -- fill(r_der.rho, 0.0)
  -- get_derivative_matrix(mat_x, beta10d1, alpha10d1, 1.0, alpha10d1, beta10d1)
  -- get_LU_decomposition(LU_x, mat_x)
  -- token += ddx(r_prim_c, r_der, LU_x)
  -- wait_for(token)
  -- err = check_ddx(coords, r_der)
  -- c.printf("Error in ddx     = %g\n", err)
  -- regentlib.assert( err <= 1.e-9, "Derivative test failed for task ddx")

  -- fill(r_der.rho, 0.0)
  -- get_derivative_matrix(mat_y, beta10d1, alpha10d1, 1.0, alpha10d1, beta10d1)
  -- get_LU_decomposition(LU_y, mat_y)
  -- token += ddy(r_prim_c, r_der, LU_y)
  -- wait_for(token)
  -- err = check_ddy(coords, r_der)
  -- c.printf("Error in ddy     = %g\n", err)
  -- regentlib.assert( err <= 1.e-9, "Derivative test failed for task ddy")

  -- fill(r_der.rho, 0.0)
  -- get_derivative_matrix(mat_z, beta10d1, alpha10d1, 1.0, alpha10d1, beta10d1)
  -- get_LU_decomposition(LU_z, mat_z)
  -- token += ddz(r_prim_c, r_der, LU_z)
  -- wait_for(token)
  -- err = check_ddz(coords, r_der)
  -- c.printf("Error in ddz     = %g\n", err)
  -- regentlib.assert( err <= 1.e-9, "Derivative test failed for task ddz")

  -- c.printf("\n")

  fill(r_der.rho, 0.0)
  for i in pencil_interior do
    get_MND_matrix(p_mat_x[i], periodic_x)
    get_LU_decomposition(p_LU_x[i], p_mat_x[i])
  end
  for i in pencil_interior do
    token += ddx_MND(p_prim_c_x[i], p_prim_l_x[i], p_der_x[i], p_LU_x[i])
  end
  wait_for(token)
  for i in pencil_interior do
    err = check_ddx(p_coords_x[i], p_der_x[i])
  end
  c.printf("Error in ddx_MND = %g\n", err)
  regentlib.assert( err <= 1.e-09, "Derivative test failed for task ddx_MND")

  fill(r_der.rho, 0.0)
  for i in pencil_interior do
    get_MND_matrix(p_mat_y[i], periodic_y)
    get_LU_decomposition(p_LU_y[i], p_mat_y[i])
  end
  for i in pencil_interior do
    token += ddy_MND(p_prim_c_y[i], p_prim_l_y[i], p_der_y[i], p_LU_y[i])
  end
  wait_for(token)
  for i in pencil_interior do
    err = check_ddy(p_coords_y[i], p_der_y[i])
  end
  c.printf("Error in ddy_MND = %g\n", err)
  regentlib.assert( err <= 1.e-09, "Derivative test failed for task ddy_MND")

  fill(r_der.rho, 0.0)
  for i in pencil_interior do
    get_MND_matrix(p_mat_z[i], periodic_z)
    get_LU_decomposition(p_LU_z[i], p_mat_z[i])
  end
  for i in pencil_interior do
    token += ddz_MND(p_prim_c_z[i], p_prim_l_z[i], p_der_z[i], p_LU_z[i])
  end
  wait_for(token)
  for i in pencil_interior do
    err = check_ddz(p_coords_z[i], p_der_z[i])
  end
  c.printf("Error in ddz_MND = %g\n", err)
  regentlib.assert( err <= 1.e-09, "Derivative test failed for task ddz_MND")

  c.printf("\nAll tests passed! :)\n")

  -- var IOtoken = 0
  -- IOtoken += write_coords(coords, "derivatives_x_c_", {0,0})
  -- wait_for(IOtoken)
  -- IOtoken += write_primitive(r_prim_c, "derivatives_x_c_", 0, {0,0})
  -- wait_for(IOtoken)
  -- IOtoken += write_coords(coords_x, "derivatives_x_e_", {0,0})
  -- wait_for(IOtoken)
  -- IOtoken += write_primitive(r_prim_l_x, "derivatives_x_e_", 0, {0,0})
  -- wait_for(IOtoken)
  -- IOtoken += write_coords(coords, "derivatives_x_der_", {0,0})
  -- wait_for(IOtoken)
  -- IOtoken += write_primitive(r_der, "derivatives_x_der_", 0, {0,0})
  -- wait_for(IOtoken)

end

regentlib.start(main)
