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

problem.X1 = -0.5
problem.Y1 = -0.5
problem.Z1 = -0.5

-- Grid spacing
problem.DX = problem.LX / problem.NX
problem.DY = problem.LY / problem.NY
problem.DZ = problem.LZ / problem.NZ

problem.ONEBYDX = 1.0 / problem.DX
problem.ONEBYDY = 1.0 / problem.DY
problem.ONEBYDZ = 1.0 / problem.DZ

problem.dt    = 0.001
problem.tstop = 0.2
problem.tviz  = 0.1

task problem.initialize( coords     : region(ispace(int3d), coordinates),
                         r_prim_c   : region(ispace(int3d), primitive),
                         dx         : double,
                         dy         : double,
                         dz         : double )
where
  reads writes(coords, r_prim_c)
do

  var rhoL : double = 1.0
  var rhoR : double = 0.125
  var pL   : double = 1.0
  var pR   : double = 0.1

  var thick : double = 1.e-10

  for i in coords.ispace do
    coords[i].x_c = problem.X1 + (i.x + 0.5) * dx
    coords[i].y_c = problem.Y1 + (i.y + 0.5) * dy
    coords[i].z_c = problem.Z1 + (i.z + 0.5) * dz

    var tmp : double = cmath.tanh( (coords[i].x_c - 1.e-6)/(thick*dx) ) - cmath.tanh( (coords[i].x_c - (1.-1.e-6))/(thick*dx) ) - 1.

    r_prim_c[i].rho = 0.5*(rhoL + rhoR) - 0.5*(rhoL - rhoR)*tmp
    r_prim_c[i].u   = 0.0
    r_prim_c[i].v   = 0.0 
    r_prim_c[i].w   = 0.0
    r_prim_c[i].p   = 0.5*(pL + pR) - 0.5*(pL - pR)*tmp
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

  return errors
end

return problem
