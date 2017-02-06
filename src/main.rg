import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")
require("derivatives")
local superlu = require("superlu_util")

-- Grid dimensions
local NX = 64
local NY = 64
local NZ = 64

-- Domain size
local LX = 2.0*math.pi
local LY = 2.0*math.pi
local LZ = 2.0*math.pi

-- Grid spacing
local DX = LX / NX
local DY = LY / NY
local DZ = LZ / NZ

local ONEBYDX = 1.0 / DX
local ONEBYDY = 1.0 / DY
local ONEBYDZ = 1.0 / DZ





-- local r_flux = regentlib.newsymbol(region(ispace(int3d), conserved), "r_flux")
-- local r_cnsr = regentlib.newsymbol(region(ispace(int3d), conserved), "r_cnsr")
-- 
-- local privileges_r_flux = regentlib.privilege(regentlib.reads,  r_flux, "rho")
-- local privileges_r_cnsr = regentlib.privilege(regentlib.writes, r_cnsr, "rho")
-- 
-- local ComputeXRHS  = make_stencil_x(r_flux, privileges_r_flux, "rho", r_cnsr, privileges_r_cnsr, "rho", NX, NY, NZ, ONEBYDX, a10d1, b10d1, c10d1, 1)
-- local ComputeYRHS  = make_stencil_y(r_flux, privileges_r_flux, "rho", r_cnsr, privileges_r_cnsr, "rho", NX, NY, NZ, ONEBYDY, a10d1, b10d1, c10d1, 1)
-- local ComputeZRHS  = make_stencil_z(r_flux, privileges_r_flux, "rho", r_cnsr, privileges_r_cnsr, "rho", NX, NY, NZ, ONEBYDZ, a10d1, b10d1, c10d1, 1)





task initialize( coords     : region(ispace(int3d), coordinates),
                 r_flux_c   : region(ispace(int3d), conserved),
                 r_flux_e_x : region(ispace(int3d), conserved),
                 r_flux_e_y : region(ispace(int3d), conserved),
                 r_flux_e_z : region(ispace(int3d), conserved),
                 dx         : double,
                 dy         : double,
                 dz         : double )
where
  reads writes(coords, r_flux_c, r_flux_e_x, r_flux_e_y, r_flux_e_z)
do
  for i in coords.ispace do
    coords[i].x_c = (i.x + 0.5) * dx
    coords[i].y_c = (i.y + 0.5) * dy
    coords[i].z_c = (i.z + 0.5) * dz
    r_flux_c[i].rho = cmath.sin(coords[i].x_c) * cmath.cos(coords[i].y_c) * cmath.cos(coords[i].z_c)
  end

  for i in r_flux_e_x do
    var x_c : double = (i.x + 0.0) * dx
    var y_c : double = (i.y + 0.5) * dy
    var z_c : double = (i.z + 0.5) * dz
    r_flux_e_x[i].rho = cmath.sin(x_c) * cmath.cos(y_c) * cmath.cos(z_c)
  end

  for i in r_flux_e_y do
    var x_c : double = (i.x + 0.5) * dx
    var y_c : double = (i.y + 0.0) * dy
    var z_c : double = (i.z + 0.5) * dz
    r_flux_e_y[i].rho = cmath.sin(x_c) * cmath.cos(y_c) * cmath.cos(z_c)
  end

  for i in r_flux_e_z do
    var x_c : double = (i.x + 0.5) * dx
    var y_c : double = (i.y + 0.5) * dy
    var z_c : double = (i.z + 0.0) * dz
    r_flux_e_z[i].rho = cmath.sin(x_c) * cmath.cos(y_c) * cmath.cos(z_c)
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
  c.printf("           dx, dy, dz = %f, %f, %f\n", dx, dy, dz)
  c.printf("====================================================\n")

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

  var token = initialize(coords, r_flux_c, r_flux_e_x, r_flux_e_y, r_flux_e_z, dx, dy, dz)
  wait_for(token)

  -- var alpha10d1 : double = 1.0/2.0
  -- var beta10d1  : double = 1.0/20.0
  var alpha10d1 : double = 1.0/3.0
  var beta10d1  : double = 0.0
  get_LU_decomposition(LU_x, beta10d1, alpha10d1, 1.0, alpha10d1, beta10d1)
  token += ddx(r_flux_c, r_cnsr, LU_x)
  wait_for(token)
  c.printf("Error in ddx = %g\n", check_ddx(coords, r_cnsr))

  var alpha06MND : double = -1.0/12.0
  var beta06MND  : double = 0.0
  get_LU_decomposition(LU_x, beta06MND, alpha06MND, 1.0, alpha06MND, beta06MND)
  token += ddx_MND(r_flux_c, r_flux_e_x, r_cnsr, LU_x)
  wait_for(token)
  c.printf("Error in ddx_MND = %g\n", check_ddx(coords, r_cnsr))
end

regentlib.start(main)
