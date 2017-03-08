import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")
require("IO")
require("interpolation")
local superlu = require("superlu_util")

local csuperlu_mapper = require("superlu_mapper")

-- Grid dimensions
local NX = 5
local NY = 80
local NZ = 5

-- Domain size
-- local LX = 2.0*math.pi
-- local LY = 2.0*math.pi
-- local LZ = 2.0*math.pi
local LX = 1.0
local LY = 10.0
local LZ = 1.0

local X1 = 0.0
local Y1 = -5.0
local Z1 = 0.0

-- Grid spacing
local DX = LX / NX
local DY = LY / NY
local DZ = LZ / NZ

local ONEBYDX = 1.0 / DX
local ONEBYDY = 1.0 / DY
local ONEBYDZ = 1.0 / DZ

local alpha06CI = 3.0/16.0
local beta06CI  = 5.0/8.0
local gamma06CI = 3.0/16.0

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

    if (coords[i].y_c < -4.0) then
      r_prim_c[i].rho = 27./7.
      r_prim_c[i].u   = 0.0 
      r_prim_c[i].v   = 4.0*cmath.sqrt(35.0)/9.0
      r_prim_c[i].w   = 0.0
      r_prim_c[i].p   = 31./3.
    else
      r_prim_c[i].rho = 1.0 + 0.2*cmath.sin(5.0*coords[i].y_c)
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
    if (y_c < -4.0) then
      r_prim_l_x[i].rho = 27./7.
      r_prim_l_x[i].u   = 0.0 
      r_prim_l_x[i].v   = 4.0*cmath.sqrt(35.0)/9.0
      r_prim_l_x[i].w   = 0.0
      r_prim_l_x[i].p   = 31./3.
    else
      r_prim_l_x[i].rho = 1.0 + 0.2*cmath.sin(5.0*y_c)
      r_prim_l_x[i].u   = 0.0
      r_prim_l_x[i].v   = 0.0 
      r_prim_l_x[i].w   = 0.0
      r_prim_l_x[i].p   = 1.0
    end
  end

  return 1
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

  --------------------------------------------------------------------------------------------
  --                       DATA STUCTURES
  --------------------------------------------------------------------------------------------
  var grid_c     = ispace(int3d, {x = Nx,   y = Ny,   z = Nz  })  -- Cell center index space

  var grid_e_x   = ispace(int3d, {x = Nx+1, y = Ny,   z = Nz  })  -- x cell edge index space
  var grid_e_y   = ispace(int3d, {x = Nx,   y = Ny+1, z = Nz  })  -- y cell edge index space
  var grid_e_z   = ispace(int3d, {x = Nx,   y = Ny,   z = Nz+1})  -- z cell edge index space

  var coords     = region(grid_c, coordinates)  -- Coordinates of cell center

  var r_prim_c   = region(grid_c,   primitive)  -- Primitive variables at cell center
  var r_prim_l_x = region(grid_e_x, primitive)  -- Primitive variables at left x cell edge
  var r_prim_r_x = region(grid_e_x, primitive)  -- Primitive variables at right x cell edge
  var r_prim_l_y = region(grid_e_y, primitive)  -- Primitive variables at left y cell edge
  var r_prim_r_y = region(grid_e_y, primitive)  -- Primitive variables at right y cell edge
  var r_prim_l_z = region(grid_e_z, primitive)  -- Primitive variables at left z cell edge
  var r_prim_r_z = region(grid_e_z, primitive)  -- Primitive variables at right z cell edge
  
  var r_rhs_l_y  = region(grid_e_y, primitive)  -- RHS for left biased interpolation y cell edge
  var r_rhs_r_y  = region(grid_e_y, primitive)  -- RHS for right biased interpolation y cell edge

  var pgrid_y    = ispace(int2d, {x = 1, y = 1}) -- Processor grid in y

  var slu_y      = region(pgrid_y, superlu.c.superlu_vars_t) -- Super LU data structure for y interpolation

  var matrix_l_y = region(pgrid_y, superlu.CSR_matrix) -- matrix data structure for y left interpolation
  var matrix_r_y = region(pgrid_y, superlu.CSR_matrix) -- matrix data structure for y right interpolation
  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  -- Initialize SuperLU stuff
  -- matrix_l_y[{0,0}] = superlu.initialize_matrix_char_y(alpha06CI, beta06CI, gamma06CI, Nx, Ny, Nz)
  -- matrix_r_y[{0,0}] = superlu.initialize_matrix_char_y(alpha06CI, beta06CI, gamma06CI, Nx, Ny, Nz)

  -- superlu.initialize_superlu_vars( matrix_l_y[{0,0}], 5*Nx*(Ny+1)*Nz, r_rhs_l_y, r_prim_l_y, slu_y ) 

  superlu.initialize_matrix_char_y(matrix_l_y, alpha06CI, beta06CI, gamma06CI, Nx, Ny, Nz)
  superlu.initialize_matrix_char_y(matrix_r_y, alpha06CI, beta06CI, gamma06CI, Nx, Ny, Nz)

  fill( r_rhs_l_y.{rho,u,v,w,p}, 0.0 )
  superlu.init_superlu_vars( matrix_l_y, 5*Nx*(Ny+1)*Nz, r_rhs_l_y, r_prim_l_y, slu_y )

  var token = initialize(coords, r_prim_c, r_prim_l_x, r_prim_l_y, r_prim_l_z, dx, dy, dz)
  wait_for(token)

  var t_start = c.legion_get_current_time_in_micros()
  token += WCHR_interpolation_y( r_prim_c, r_prim_l_y, r_prim_r_y, r_rhs_l_y, r_rhs_r_y, matrix_l_y, matrix_r_y, slu_y )
  wait_for(token)
  var t_WCHR = c.legion_get_current_time_in_micros() - t_start
  c.printf("Time to get the WCHR interpolation: %12.5e\n", (t_WCHR)*1e-6)

  -- write_coords(coords)
  -- write_primitive(r_prim_c, "cell_primitive", 0)
  -- write_primitive(r_prim_l_y, "edge_primitive_l_y", 0)
  write_primitive(r_prim_r_y, "edge_primitive_r_y", 0)
end

regentlib.start(main, csuperlu_mapper.register_mappers)
