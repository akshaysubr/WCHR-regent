import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")
require("IO")
local interpolation = require("interpolation")
local problem = require("problem")

-- Grid dimensions
local NX = 64
local NY = 64
local NZ = 64

-- Domain size
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

    if (coords[i].y_c < -4.0) then
      r_prim_c[i].rho = 27./7.
      r_prim_c[i].u   = 1.0 
      r_prim_c[i].v   = 4.0*cmath.sqrt(35.0)/9.0
      r_prim_c[i].w   = 1.0
      r_prim_c[i].p   = 31./3.
    else
      r_prim_c[i].rho = 1.0 + 0.2*cmath.sin(5.0*coords[i].y_c)
      r_prim_c[i].u   = 0.0
      r_prim_c[i].v   = 0.0 
      r_prim_c[i].w   = 0.0
      r_prim_c[i].p   = 1.0
    end
  end

  for i in coords_x.ispace do
    coords_x[i].x_c = X1 + (i.x                 ) * dx
    coords_x[i].y_c = Y1 + (i.y - n_ghosts + 0.5) * dy
    coords_x[i].z_c = Z1 + (i.z - n_ghosts + 0.5) * dz
  end

  for i in coords_y.ispace do
    coords_y[i].x_c = X1 + (i.x - n_ghosts + 0.5) * dx
    coords_y[i].y_c = Y1 + (i.y                 ) * dy
    coords_y[i].z_c = Z1 + (i.z - n_ghosts + 0.5) * dz
  end

  for i in coords_z.ispace do
    coords_z[i].x_c = X1 + (i.x - n_ghosts + 0.5) * dx
    coords_z[i].y_c = Y1 + (i.y - n_ghosts + 0.5) * dy
    coords_z[i].z_c = Z1 + (i.z                 ) * dz
  end

  for i in r_prim_l_y do
    var x_c : double = X1 + (i.x) * dx
    var y_c : double = Y1 + (i.y) * dy
    var z_c : double = Z1 + (i.z) * dz
    if (y_c < -4.0) then
      r_prim_l_y[i].rho = 27./7.
      r_prim_l_y[i].u   = 1.0 
      r_prim_l_y[i].v   = 4.0*cmath.sqrt(35.0)/9.0
      r_prim_l_y[i].w   = 0.0
      r_prim_l_y[i].p   = 31./3.
    else
      r_prim_l_y[i].rho = 1.0 + 0.2*cmath.sin(5.0*y_c)
      r_prim_l_y[i].u   = 0.0
      r_prim_l_y[i].v   = 0.0 
      r_prim_l_y[i].w   = 0.0
      r_prim_l_y[i].p   = 1.0
    end
  end

  return 1
end

task periodic_ghost_cells_y( r_prim_c   : region(ispace(int3d), primitive),
                             n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Ny_g = bounds_c.hi.y + 1
  var Ny   = Ny_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z+1 do
    for i = bounds_c.lo.x, bounds_c.hi.x+1 do
      for j = 0,n_ghosts do
        var ghost_l = int3d {i,j,k}
        var int_r   = int3d {i,Ny+j,k}
        r_prim_c[ghost_l].rho = r_prim_c[int_r].rho
        r_prim_c[ghost_l].u   = r_prim_c[int_r].u
        r_prim_c[ghost_l].v   = r_prim_c[int_r].v
        r_prim_c[ghost_l].w   = r_prim_c[int_r].w
        r_prim_c[ghost_l].p   = r_prim_c[int_r].p

        var ghost_r = int3d {i,Ny+n_ghosts+j,k}
        var int_l   = int3d {i,n_ghosts+j,k}
        r_prim_c[ghost_r].rho = r_prim_c[int_l].rho
        r_prim_c[ghost_r].u   = r_prim_c[int_l].u
        r_prim_c[ghost_r].v   = r_prim_c[int_l].v
        r_prim_c[ghost_r].w   = r_prim_c[int_l].w
        r_prim_c[ghost_r].p   = r_prim_c[int_l].p
      end
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

  var n_ghosts : int64 = interpolation.n_ghosts
  var Nx_g : int64 = Nx + 2*n_ghosts
  var Ny_g : int64 = Ny + 2*n_ghosts
  var Nz_g : int64 = Nz + 2*n_ghosts

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

  var coords     = region(grid_c, coordinates)  -- Coordinates of cell center

  var coords_x   = region(grid_e_x, coordinates) -- Coordinates of x cell edge
  var coords_y   = region(grid_e_y, coordinates) -- Coordinates of y cell edge
  var coords_z   = region(grid_e_z, coordinates) -- Coordinates of z cell edge

  var r_prim_c   = region(grid_c,   primitive)  -- Primitive variables at cell center
  var r_prim_l_x = region(grid_e_x, primitive)  -- Primitive variables at left x cell edge
  var r_prim_r_x = region(grid_e_x, primitive)  -- Primitive variables at right x cell edge
  var r_prim_l_y = region(grid_e_y, primitive)  -- Primitive variables at left y cell edge
  var r_prim_r_y = region(grid_e_y, primitive)  -- Primitive variables at right y cell edge
  var r_prim_l_z = region(grid_e_z, primitive)  -- Primitive variables at left z cell edge
  var r_prim_r_z = region(grid_e_z, primitive)  -- Primitive variables at right z cell edge
  
  var alpha_l = region( grid_e_y, coeffs )
  var beta_l  = region( grid_e_y, coeffs )
  var gamma_l = region( grid_e_y, coeffs )

  var alpha_r = region( grid_e_y, coeffs )
  var beta_r  = region( grid_e_y, coeffs )
  var gamma_r = region( grid_e_y, coeffs )

  var rho_avg = region( grid_e_y, double )
  var sos_avg = region( grid_e_y, double )

  var block_d    = region( grid_e_y, double[9] )
  var block_Uinv = region( grid_e_y, double[9] )
  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  fill( r_prim_l_y.{rho,u,v,w,p}, 0.0 )
  fill( r_prim_r_y.{rho,u,v,w,p}, 0.0 )

  var token = initialize(coords, coords_x, coords_y, coords_z, r_prim_c, r_prim_l_x, r_prim_l_y, r_prim_l_z, dx, dy, dz, n_ghosts)
  if (problem.periodic_y) then
    token += periodic_ghost_cells_y(r_prim_c, n_ghosts)
  end
  wait_for(token)

  var t_start = c.legion_get_current_time_in_micros()
  token += WCHR_interpolation_y( r_prim_c, r_prim_l_y, r_prim_r_y, alpha_l, beta_l, gamma_l, 
                                 alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv )
  wait_for(token)
  var t_WCHR = c.legion_get_current_time_in_micros() - t_start
  c.printf("Time to get the WCHR interpolation: %12.5e\n", (t_WCHR)*1e-6)

  var IOtoken = 0
  IOtoken += write_coords(coords, "interpolation_y_c_", 0, {1,1})
  wait_for(IOtoken)
  IOtoken += write_coords(coords_y, "interpolation_y_l_", 0, {1,1})
  wait_for(IOtoken)
  IOtoken += write_coords(coords_y, "interpolation_y_r_", 0, {1,1})
  wait_for(IOtoken)
  IOtoken += write_primitive(r_prim_c, "interpolation_y_c_", 0, 0, {1,1})
  wait_for(IOtoken)
  IOtoken += write_primitive(r_prim_l_y, "interpolation_y_l_", 0, 0, {1,1})
  wait_for(IOtoken)
  IOtoken += write_primitive(r_prim_r_y, "interpolation_y_r_", 0, 0, {1,1})
  wait_for(IOtoken)
end

regentlib.start(main)
