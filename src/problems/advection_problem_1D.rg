import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")

local problem = {}

-- Grid dimensions
problem.NX = 200
problem.NY = 1
problem.NZ = 1

-- Domain size
problem.LX = 2.0
problem.LY = 1.0
problem.LZ = 1.0

problem.X1 = -1.0
problem.Y1 = -0.5
problem.Z1 = -0.5

-- Grid spacing
problem.DX = problem.LX / problem.NX
problem.DY = problem.LY / problem.NY
problem.DZ = problem.LZ / problem.NZ

problem.ONEBYDX = 1.0 / problem.DX
problem.ONEBYDY = 1.0 / problem.DY
problem.ONEBYDZ = 1.0 / problem.DZ

problem.dt    = 0.2 * problem.DX
problem.tstop = 2.0
problem.tviz  = 1.0

task problem.initialize( coords     : region(ispace(int3d), coordinates),
                         r_prim_c   : region(ispace(int3d), primitive),
                         dx         : double,
                         dy         : double,
                         dz         : double )
where
  reads writes(coords, r_prim_c)
do
  for i in coords.ispace do
    coords[i].x_c = problem.X1 + (i.x + 0.5) * dx
    coords[i].y_c = problem.Y1 + (i.y + 0.5) * dy
    coords[i].z_c = problem.Z1 + (i.z + 0.5) * dz

    r_prim_c[i].rho = 1.0 + 0.5*cmath.sin(PI*coords[i].x_c)
    r_prim_c[i].u   = 1.0
    r_prim_c[i].v   = 0.0 
    r_prim_c[i].w   = 0.0
    r_prim_c[i].p   = 1.0
  end

  return 1
end

task problem.get_errors( coords     : region(ispace(int3d), coordinates),
                         r_prim_c   : region(ispace(int3d), primitive),
                         tsim       : double )
where
  reads(coords, r_prim_c)
do

  var errors : double[5] = array(0.0, 0.0, 0.0, 0.0, 0.0)

  for i in r_prim_c do
    var err : double

    err = cmath.fabs( r_prim_c[i].rho - (1.0 + 0.5*cmath.sin(PI*(coords[i].x_c - tsim))) )
    if err > errors[0] then
      errors[0] = err
    end

    err = cmath.fabs( r_prim_c[i].u   - 1.0 )
    if err > errors[1] then
      errors[1] = err
    end

    err = cmath.fabs( r_prim_c[i].v   - 0.0 )
    if err > errors[2] then
      errors[2] = err
    end

    err = cmath.fabs( r_prim_c[i].w   - 0.0 )
    if err > errors[3] then
      errors[3] = err
    end

    err = cmath.fabs( r_prim_c[i].p   - 1.0 )
    if err > errors[4] then
      errors[4] = err
    end

  end

  return errors
end

return problem
