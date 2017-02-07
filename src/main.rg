import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")
require("IO")
local superlu = require("superlu_util")

-- Grid dimensions
local NX = 200
local NY = 8
local NZ = 8

-- Domain size
-- local LX = 2.0*math.pi
-- local LY = 2.0*math.pi
-- local LZ = 2.0*math.pi
local LX = 10.0
local LY = 1.0
local LZ = 1.0

local X1 = -5.0
local Y1 = 0.0
local Z1 = 0.0

-- Grid spacing
local DX = LX / NX
local DY = LY / NY
local DZ = LZ / NZ

local ONEBYDX = 1.0 / DX
local ONEBYDY = 1.0 / DY
local ONEBYDZ = 1.0 / DZ

-- Make the node and midpoint-node differencing tasks (Using pentadiagonal solver for this instead of tridiagonal solver)
-- require("derivatives")
-- local alpha10d1 = 1.0/3.0
-- local beta10d1  = 0.0
-- local a06d1 = ( 14.0/ 9.0)/2.0
-- local b06d1 = (  1.0/ 9.0)/4.0
-- local c06d1 = (  0.0/100.0)/6.0
-- 
-- local alpha06MND = -1.0/12.0
-- local beta06MND  = 0.0
-- local a06MND = 16.0/9.0
-- local b06MND = (-17.0/18.0)/2.0
-- local c06MND = (0.0)/3.0
-- 
-- local r_flux   = regentlib.newsymbol(region(ispace(int3d), conserved), "r_flux")
-- local r_flux_e = regentlib.newsymbol(region(ispace(int3d), conserved), "r_flux_e")
-- local r_cnsr   = regentlib.newsymbol(region(ispace(int3d), conserved), "r_cnsr")
-- 
-- local ddx     = make_ddx(r_flux, "rho", r_cnsr, "rho", NX, NY, NZ, ONEBYDX, a06d1, b06d1, c06d1)
-- local ddx_MND = make_ddx_MND(r_flux, r_flux_e, "rho", r_cnsr, "rho", NX, NY, NZ, ONEBYDX, a06MND, b06MND, c06MND)
-----------------------------------------------------

task initialize( coords     : region(ispace(int3d), coordinates),
                 r_prim_c   : region(ispace(int3d), primitive),
                 r_prim_l_x : region(ispace(int3d), primitive),
                 r_prim_l_y : region(ispace(int3d), primitive),
                 r_prim_l_z : region(ispace(int3d), primitive),
                 dx         : double,
                 dy         : double,
                 dz         : double )
where
  reads writes(coords, r_prim_c, r_prim_l_x, r_prim_l_y, r_prim_l_z)
do
  for i in coords.ispace do
    coords[i].x_c = X1 + (i.x + 0.5) * dx
    coords[i].y_c = Y1 + (i.y + 0.5) * dy
    coords[i].z_c = Z1 + (i.z + 0.5) * dz

    if (coords[i].x_c < -4.0) then
      r_prim_c[i].rho = 3.857143
      r_prim_c[i].u   = 2.62936
      r_prim_c[i].v   = 0.0 
      r_prim_c[i].w   = 0.0
      r_prim_c[i].p   = 10.33339
    else
      r_prim_c[i].rho = 1.0 + 0.2*cmath.sin(5.0*coords[i].x_c)
      r_prim_c[i].u   = 0.0
      r_prim_c[i].v   = 0.0 
      r_prim_c[i].w   = 0.0
      r_prim_c[i].p   = 1.0
    end
  end

  for i in r_prim_l_x do
    var x_c : double = X1 + (i.x) * dx
    var y_c : double = Y1 + (i.y) * dy
    var z_c : double = Z1 + (i.z) * dz
    if (x_c < -4.0) then
      r_prim_l_x[i].rho = 3.857143
      r_prim_l_x[i].u   = 2.62936
      r_prim_l_x[i].v   = 0.0 
      r_prim_l_x[i].w   = 0.0
      r_prim_l_x[i].p   = 10.33339
    else
      r_prim_l_x[i].rho = 1.0 + 0.2*cmath.sin(5.0*x_c)
      r_prim_l_x[i].u   = 0.0
      r_prim_l_x[i].v   = 0.0 
      r_prim_l_x[i].w   = 0.0
      r_prim_l_x[i].p   = 1.0
    end
  end

  return 1
end

task check_ddx( coords     : region(ispace(int3d), coordinates),
                r_cnsr     : region(ispace(int3d), conserved) )
where
  reads writes(coords, r_cnsr)
do
  var err : double = 0.0
  for i in coords.ispace do
    var err_t : double = cmath.fabs( r_cnsr[i].rho - cmath.cos(coords[i].x_c) * cmath.cos(coords[i].y_c) * cmath.cos(coords[i].z_c) )
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

  generate_hdf5_file_coords("cell_coords.h5", Nx, Ny, Nz)
  generate_hdf5_file("cell_primitive.h5", Nx, Ny, Nz)
  generate_hdf5_file("edge_primitive_l_x.h5", Nx+1, Ny, Nz)

  --------------------------------------------------------------------------------------------
  --                       DATA STUCTURES
  --------------------------------------------------------------------------------------------
  var grid_c     = ispace(int3d, {x = Nx,   y = Ny,   z = Nz  })  -- Cell center index space

  var grid_e_x   = ispace(int3d, {x = Nx+1, y = Ny,   z = Nz  })  -- x cell edge index space
  var grid_e_y   = ispace(int3d, {x = Nx,   y = Ny+1, z = Nz  })  -- y cell edge index space
  var grid_e_z   = ispace(int3d, {x = Nx,   y = Ny,   z = Nz+1})  -- z cell edge index space

  var coords     = region(grid_c, coordinates)  -- Coordinates of cell center

  var r_cnsr     = region(grid_c,   conserved)  -- Conserved variables at cell center
  var r_aux      = region(grid_c,   auxiliary)  -- Conserved variables at cell center

  var r_prim_c   = region(grid_c,   primitive)  -- Primitive variables at cell center
  var r_prim_l_x = region(grid_e_x, primitive)  -- Primitive variables at left x cell edge
  var r_prim_l_y = region(grid_e_y, primitive)  -- Primitive variables at left y cell edge
  var r_prim_l_z = region(grid_e_z, primitive)  -- Primitive variables at left z cell edge
  var r_prim_r_x = region(grid_e_x, primitive)  -- Primitive variables at right x cell edge
  var r_prim_r_y = region(grid_e_y, primitive)  -- Primitive variables at right y cell edge
  var r_prim_r_z = region(grid_e_z, primitive)  -- Primitive variables at right z cell edge
  
  var r_vect_x   = region(grid_e_x, soe_vector) -- RHS variables at cell edge for characteristic interpolation
  var r_vect_y   = region(grid_e_y, soe_vector) -- RHS variables at cell edge for characteristic interpolation
  var r_vect_z   = region(grid_e_z, soe_vector) -- RHS variables at cell edge for characteristic interpolation

  var r_flux_c   = region(grid_c,   conserved)  -- Flux at cell center
  var r_flux_e_x = region(grid_e_x, conserved)  -- Flux at x cell edge
  var r_flux_e_y = region(grid_e_y, conserved)  -- Flux at y cell edge
  var r_flux_e_z = region(grid_e_z, conserved)  -- Flux at z cell edge
  
  var r_frhs_c_x = region(grid_c,   double)     -- RHS while getting x flux derivative
  var r_frhs_c_y = region(grid_c,   double)     -- RHS while getting y flux derivative
  var r_frhs_c_z = region(grid_c,   double)     -- RHS while getting z flux derivative
  
  var r_fder_c_x = region(grid_c,   double)     -- x flux derivative
  var r_fder_c_y = region(grid_c,   double)     -- y flux derivative
  var r_fder_c_z = region(grid_c,   double)     -- z flux derivative
  
  var r_rhs      = region(grid_c,   conserved)  -- RHS for time stepping at cell center

  var LU_x       = region(ispace(int3d, {x = Nx, y = 1, z = 1}), LU_struct) -- Data structure to hold x derivative LU decomposition
  var LU_y       = region(ispace(int3d, {x = Ny, y = 1, z = 1}), LU_struct) -- Data structure to hold y derivative LU decomposition
  var LU_z       = region(ispace(int3d, {x = Nz, y = 1, z = 1}), LU_struct) -- Data structure to hold z derivative LU decomposition

  var pgrid_x    = ispace(int2d, {x = 1, y = 1}) -- Processor grid in x
  var pgrid_y    = ispace(int2d, {x = 1, y = 1}) -- Processor grid in y
  var pgrid_z    = ispace(int2d, {x = 1, y = 1}) -- Processor grid in z

  var slu_x      = region(pgrid_x, superlu.c.superlu_vars_t) -- Super LU data structure for x interpolation
  var slu_y      = region(pgrid_y, superlu.c.superlu_vars_t) -- Super LU data structure for y interpolation
  var slu_z      = region(pgrid_z, superlu.c.superlu_vars_t) -- Super LU data structure for z interpolation

  var matrix_x_l = region(pgrid_x, superlu.CSR_matrix) -- matrix data structure for x left interpolation
  var matrix_x_r = region(pgrid_x, superlu.CSR_matrix) -- matrix data structure for x right interpolation
  var matrix_y_l = region(pgrid_y, superlu.CSR_matrix) -- matrix data structure for y left interpolation
  var matrix_y_r = region(pgrid_y, superlu.CSR_matrix) -- matrix data structure for y right interpolation
  var matrix_z_l = region(pgrid_z, superlu.CSR_matrix) -- matrix data structure for z left interpolation
  var matrix_z_r = region(pgrid_z, superlu.CSR_matrix) -- matrix data structure for z right interpolation
  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  attach(hdf5, coords.{x_c, y_c, z_c}, "cell_coords.h5", regentlib.file_read_write)
  attach(hdf5, r_prim_c.{rho, u, v, w, p}, "cell_primitive.h5", regentlib.file_read_write)
  attach(hdf5, r_prim_l_x.{rho, u, v, w, p}, "edge_primitive_l_x.h5", regentlib.file_read_write)
  -- acquire(r_prim_c)
  -- acquire(r_prim_l_x)

  var token = initialize(coords, r_prim_c, r_prim_l_x, r_prim_l_y, r_prim_l_z, dx, dy, dz)
  wait_for(token)

  -- release(r_prim_c)
  -- release(r_prim_l_x)
  detach(hdf5, coords.{x_c, y_c, z_c})
  detach(hdf5, r_prim_c.{rho, u, v, w, p})
  detach(hdf5, r_prim_l_x.{rho, u, v, w, p})

end

regentlib.start(main)
