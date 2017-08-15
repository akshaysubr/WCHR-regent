import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")

local problem = {}

-- Problem specific parameters
problem.gamma = 1.4
problem.Rgas  = 1.0
problem.Re    = 100.
problem.Pr    = 1.

-- Grid dimensions
problem.NX = 64 
problem.NY = 64 
problem.NZ = 64 

-- Periodicity
problem.periodic_x = true
problem.periodic_y = true
problem.periodic_z = true

-- Domain size
problem.LX = 2.0*PI
problem.LY = 2.0*PI
problem.LZ = 2.0*PI

problem.X1 = 0.0
problem.Y1 = 0.0
problem.Z1 = 0.0

-- Grid spacing
problem.DX = problem.LX / problem.NX
problem.DY = problem.LY / problem.NY
problem.DZ = problem.LZ / problem.NZ

problem.ONEBYDX = 1.0 / problem.DX
problem.ONEBYDY = 1.0 / problem.DY
problem.ONEBYDZ = 1.0 / problem.DZ

problem.timestepping_setting = "CONSTANT_CFL_NUM" -- "CONSTANT_TIME_STEP" / "CONSTANT_CFL_NUM"
problem.dt_or_CFL_num        = 0.6
problem.tstop                = 10.0
problem.tviz                 = 0.1

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

    r_prim_c[i].rho = 1.0
    r_prim_c[i].u   = cmath.sin(coords[i].x_c) * cmath.cos(coords[i].y_c) * cmath.cos(coords[i].z_c)
    r_prim_c[i].v   =-cmath.cos(coords[i].x_c) * cmath.sin(coords[i].y_c) * cmath.cos(coords[i].z_c) 
    r_prim_c[i].w   = 0.0
    r_prim_c[i].p   = 100.0/problem.gamma + (1.0/16.0)*( (cmath.cos(2.0*coords[i].z_c) +2.0)*(cmath.cos(2.0*coords[i].x_c) + cmath.cos(2.0*coords[i].y_c)) - 2.0)
  end

  return 1
end

task problem.get_transport_coeffs( r_prim : region(ispace(int3d), primitive),
                                   r_aux  : region(ispace(int3d), auxiliary),
                                   r_visc : region(ispace(int3d), transport_coeffs) )
where
  reads(r_prim.{}, r_aux.T), writes(r_visc)
do
  var mu_s  : double = 1. / problem.Re
  var kappa : double = problem.Pr * (problem.gamma / (problem.gamma - 1.)) * problem.Rgas * mu_s
  for i in r_visc do
    r_visc[i].mu_s  = mu_s
    r_visc[i].mu_b  = 0.
    r_visc[i].kappa = kappa
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

task problem.TKE( r_prim_c : region(ispace(int3d), primitive) )
where
  reads(r_prim_c)
do
  var TKE : double = 0.0
  for i in r_prim_c do
    TKE += 0.5 * r_prim_c[i].rho * (r_prim_c[i].u*r_prim_c[i].u + r_prim_c[i].v*r_prim_c[i].v + r_prim_c[i].w*r_prim_c[i].w)
  end
  return TKE
end

task problem.enstrophy( r_duidxj : region(ispace(int3d), tensor2) )
where
  reads(r_duidxj)
do
  var enstrophy : double = 0.0
  for i in r_duidxj do
    var omega_x : double = r_duidxj[i]._32 - r_duidxj[i]._23
    var omega_y : double = r_duidxj[i]._13 - r_duidxj[i]._31
    var omega_z : double = r_duidxj[i]._21 - r_duidxj[i]._12
    enstrophy += omega_x*omega_x + omega_y*omega_y + omega_z*omega_z
  end
  return enstrophy
end

return problem
