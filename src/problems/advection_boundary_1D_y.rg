import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")

local problem = {}

-- Problem specific parameters
problem.gamma = 1.4  -- Ratio of specific heats
problem.Rgas  = 1.0
problem.Re    = 100.
problem.Pr    = 1.
problem.viscous = false

problem.epsilon  = 0.1
problem.width    = 0.05
problem.velocity = 0.5

-- Grid dimensions
problem.NX = 1
problem.NY = 128
problem.NZ = 1

-- Periodicity
problem.periodic_x = true
problem.periodic_y = false
problem.periodic_z = true

-- Boundary (if not periodic)
-- condition: DIRICHLET, EXTRAPOLATION, SUBSONIC_INFLOW, SUBSONIC_OUTFLOW
problem.boundary_l_x = { condition="EXTRAPOLATION", rho=1., u=0., v=0., w=0., p=1. }
problem.boundary_r_x = { condition="EXTRAPOLATION", rho=1., u=0., v=0., w=0., p=1. }
problem.boundary_l_y = { condition="EXTRAPOLATION", rho=1., u=0., v=0., w=0., p=1. }
problem.boundary_r_y = { condition="EXTRAPOLATION", rho=1., u=0., v=0., w=0., p=1. }
problem.boundary_l_z = { condition="EXTRAPOLATION", rho=1., u=0., v=0., w=0., p=1. }
problem.boundary_r_z = { condition="EXTRAPOLATION", rho=1., u=0., v=0., w=0., p=1. }

-- Domain size
problem.LX = 1.0
problem.LY = 1.0
problem.LZ = 1.0

problem.X1 = -0.5
problem.Y1 =  0.0
problem.Z1 = -0.5

-- Grid spacing
problem.DX = problem.LX / problem.NX
problem.DY = problem.LY / problem.NY
problem.DZ = problem.LZ / problem.NZ

problem.ONEBYDX = 1.0 / problem.DX
problem.ONEBYDY = 1.0 / problem.DY
problem.ONEBYDZ = 1.0 / problem.DZ

problem.timestepping_setting = "CONSTANT_TIME_STEP" -- "CONSTANT_TIME_STEP" / "CONSTANT_CFL_NUM"
problem.dt_or_CFL_num        = 2.0e-3
problem.tstop                = 2.0e0
problem.tviz                 = 0.1

task problem.initialize( coords     : region(ispace(int3d), coordinates),
                         r_prim_c   : region(ispace(int3d), primitive),
                         dx         : double,
                         dy         : double,
                         dz         : double,
                         n_ghosts   : int64 )
where
  reads writes(coords, r_prim_c)
do
  for i in coords.ispace do
    var idx = int3d {x = i.x - n_ghosts, y = i.y - n_ghosts, z = i.z - n_ghosts}
    coords[i].x_c = problem.X1 + (idx.x + 0.5) * dx
    coords[i].y_c = problem.Y1 + (idx.y + 0.5) * dy
    coords[i].z_c = problem.Z1 + (idx.z + 0.5) * dz

    var exponent = (coords[i].y_c - problem.Y1 - 0.5*problem.LY)/problem.width
    r_prim_c[i].rho = 1.0 + problem.epsilon * cmath.exp(-exponent*exponent)
    r_prim_c[i].u   = 0.0
    r_prim_c[i].v   = problem.velocity 
    r_prim_c[i].w   = 0.0
    r_prim_c[i].p   = 1.0
  end

  return 1
end

task problem.get_transport_coeffs( r_prim : region(ispace(int3d), primitive),
                                   r_aux  : region(ispace(int3d), auxiliary),
                                   r_visc : region(ispace(int3d), transport_coeffs) )
where
  writes(r_visc)
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

  for i in r_prim_c do
    var err : double

    var exponent = (coords[i].y_c - problem.velocity*tsim - problem.Y1 - 0.5*problem.LY)/problem.width
    var rho_exact = 1.0 + problem.epsilon * cmath.exp(-exponent*exponent)
    err = cmath.fabs( r_prim_c[i].rho - rho_exact )
    if err > errors[0] then
      errors[0] = err
    end

    err = cmath.fabs( r_prim_c[i].u   - 0.0 )
    if err > errors[1] then
      errors[1] = err
    end

    err = cmath.fabs( r_prim_c[i].v   - problem.velocity )
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

return problem
