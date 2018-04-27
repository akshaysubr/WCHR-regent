import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")

local problem = {}

-- Problem specific parameters
problem.gamma = 1.4  -- Ratio of specific heats
problem.Rgas  = 1.0

-- Grid dimensions
problem.NX = 200
problem.NY = 1
problem.NZ = 1

-- Periodicity
problem.periodic_x = true
problem.periodic_y = true
problem.periodic_z = true

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

problem.timestepping_setting = "CONSTANT_CFL_NUM" -- "CONSTANT_TIME_STEP" / "CONSTANT_CFL_NUM"
problem.dt_or_CFL_num        = 0.6
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

task problem.get_transport_coeffs( r_prim : region(ispace(int3d), primitive),
                                   r_aux  : region(ispace(int3d), auxiliary),
                                   r_visc : region(ispace(int3d), transport_coeffs) )
where
  reads(r_prim.{}, r_aux.T), writes(r_visc)
do
  for i in r_visc do
    r_visc[i].mu_s  = 0.
    r_visc[i].mu_b  = 0.
    r_visc[i].kappa = 0.
  end
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
