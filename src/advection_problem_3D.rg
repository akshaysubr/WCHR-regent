import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")

local problem = {}

-- Grid dimensions
problem.NX = 64
problem.NY = 64
problem.NZ = 64

-- Domain size
problem.LX = 2.0
problem.LY = 2.0
problem.LZ = 2.0

problem.X1 = -1.0
problem.Y1 = -1.0
problem.Z1 = -1.0

-- Grid spacing
problem.DX = problem.LX / problem.NX
problem.DY = problem.LY / problem.NY
problem.DZ = problem.LZ / problem.NZ

problem.ONEBYDX = 1.0 / problem.DX
problem.ONEBYDY = 1.0 / problem.DY
problem.ONEBYDZ = 1.0 / problem.DZ

problem.dt    = 0.2 * cmath.fmin(problem.DX, problem.DY)
problem.tstop = 0.05

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

    r_prim_c[i].rho = 1.0 + 0.5*cmath.exp(-cmath.pow((coords[i].x_c/0.2), 2) - cmath.pow((coords[i].y_c/0.2), 2) - cmath.pow((coords[i].z_c/0.2), 2))
    r_prim_c[i].u   = 1.0
    r_prim_c[i].v   = 1.0 
    r_prim_c[i].w   = 1.0
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

    var x0 : double = coords[i].x_c - tsim
    x0 = x0 - cmath.nearbyint(x0/problem.LX)*problem.LX

    var y0 : double = coords[i].y_c - tsim
    y0 = y0 - cmath.nearbyint(y0/problem.LY)*problem.LY

    var z0 : double = coords[i].z_c - tsim
    z0 = z0 - cmath.nearbyint(z0/problem.LZ)*problem.LZ

    err = cmath.fabs( r_prim_c[i].rho - (1.0 + 0.5*cmath.exp( -cmath.pow(( x0/0.2), 2) - cmath.pow(( y0/0.2), 2) - cmath.pow(( z0/0.2), 2) )) )
    if err > errors[0] then
      errors[0] = err
    end

    err = cmath.fabs( r_prim_c[i].u   - 1.0 )
    if err > errors[1] then
      errors[1] = err
    end

    err = cmath.fabs( r_prim_c[i].v   - 1.0 )
    if err > errors[2] then
      errors[2] = err
    end

    err = cmath.fabs( r_prim_c[i].w   - 1.0 )
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
