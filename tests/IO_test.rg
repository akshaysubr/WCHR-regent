import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")
require("IO")

-- Grid dimensions
local NX = 200
local NY = 8
local NZ = 8

-- Domain size
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
  var r_prim_l_y = region(grid_e_y, primitive)  -- Primitive variables at left y cell edge
  var r_prim_l_z = region(grid_e_z, primitive)  -- Primitive variables at left z cell edge
  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  var token = initialize(coords, r_prim_c, r_prim_l_x, r_prim_l_y, r_prim_l_z, dx, dy, dz)
  wait_for(token)

  write_coords(coords, "cell_")
  write_primitive(r_prim_c, "cell_", 0)
  -- write_primitive(r_prim_l_x, "cell_primitive_l_x", 0)

end

regentlib.start(main)
