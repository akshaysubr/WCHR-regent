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
problem.boundary_l_y = { condition="EXTRAPOLATION", rho=1., u=0., v=0., w=0., p=1. }
problem.boundary_r_y = { condition="EXTRAPOLATION", rho=1., u=0., v=0., w=0., p=1. }

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


-- DEFAULT BOUNDARY CONDITIONS --
if problem.boundary_l_x           == nil then problem.boundary_l_x           = {}              end
if problem.boundary_l_x.condition == nil then problem.boundary_l_x.condition = "EXTRAPOLATION" end
if problem.boundary_l_x.rho       == nil then problem.boundary_l_x.rho       = 1.              end
if problem.boundary_l_x.u         == nil then problem.boundary_l_x.u         = 0.              end
if problem.boundary_l_x.v         == nil then problem.boundary_l_x.v         = 0.              end
if problem.boundary_l_x.w         == nil then problem.boundary_l_x.w         = 0.              end
if problem.boundary_l_x.p         == nil then problem.boundary_l_x.p         = 1.              end
if problem.boundary_l_x.L_x       == nil then problem.boundary_l_x.L_x       = 0.1             end
if problem.boundary_l_x.sigma     == nil then problem.boundary_l_x.sigma     = 0.005           end
if problem.boundary_l_x.eta_1     == nil then problem.boundary_l_x.eta_1     = 2.0             end
if problem.boundary_l_x.eta_2     == nil then problem.boundary_l_x.eta_2     = 2.0             end
if problem.boundary_l_x.eta_3     == nil then problem.boundary_l_x.eta_3     = 2.0             end
if problem.boundary_l_x.eta_4     == nil then problem.boundary_l_x.eta_4     = 2.0             end
if problem.boundary_l_x.eta_5     == nil then problem.boundary_l_x.eta_5     = 2.0             end

if problem.boundary_r_x           == nil then problem.boundary_r_x           = {}              end
if problem.boundary_r_x.condition == nil then problem.boundary_r_x.condition = "EXTRAPOLATION" end
if problem.boundary_r_x.rho       == nil then problem.boundary_r_x.rho       = 1.              end
if problem.boundary_r_x.u         == nil then problem.boundary_r_x.u         = 0.              end
if problem.boundary_r_x.v         == nil then problem.boundary_r_x.v         = 0.              end
if problem.boundary_r_x.w         == nil then problem.boundary_r_x.w         = 0.              end
if problem.boundary_r_x.p         == nil then problem.boundary_r_x.p         = 1.              end
if problem.boundary_r_x.L_x       == nil then problem.boundary_r_x.L_x       = 0.1             end
if problem.boundary_r_x.sigma     == nil then problem.boundary_r_x.sigma     = 0.005           end
if problem.boundary_r_x.eta_1     == nil then problem.boundary_r_x.eta_1     = 2.0             end
if problem.boundary_r_x.eta_2     == nil then problem.boundary_r_x.eta_2     = 2.0             end
if problem.boundary_r_x.eta_3     == nil then problem.boundary_r_x.eta_3     = 2.0             end
if problem.boundary_r_x.eta_4     == nil then problem.boundary_r_x.eta_4     = 2.0             end
if problem.boundary_r_x.eta_5     == nil then problem.boundary_r_x.eta_5     = 2.0             end

if problem.boundary_l_y           == nil then problem.boundary_l_y           = {}              end
if problem.boundary_l_y.condition == nil then problem.boundary_l_y.condition = "EXTRAPOLATION" end
if problem.boundary_l_y.rho       == nil then problem.boundary_l_y.rho       = 1.              end
if problem.boundary_l_y.u         == nil then problem.boundary_l_y.u         = 0.              end
if problem.boundary_l_y.v         == nil then problem.boundary_l_y.v         = 0.              end
if problem.boundary_l_y.w         == nil then problem.boundary_l_y.w         = 0.              end
if problem.boundary_l_y.p         == nil then problem.boundary_l_y.p         = 1.              end
if problem.boundary_l_y.L_x       == nil then problem.boundary_l_y.L_x       = 0.1             end
if problem.boundary_l_y.sigma     == nil then problem.boundary_l_y.sigma     = 0.005           end
if problem.boundary_l_y.eta_1     == nil then problem.boundary_l_y.eta_1     = 2.0             end
if problem.boundary_l_y.eta_2     == nil then problem.boundary_l_y.eta_2     = 2.0             end
if problem.boundary_l_y.eta_3     == nil then problem.boundary_l_y.eta_3     = 2.0             end
if problem.boundary_l_y.eta_4     == nil then problem.boundary_l_y.eta_4     = 2.0             end
if problem.boundary_l_y.eta_5     == nil then problem.boundary_l_y.eta_5     = 2.0             end

if problem.boundary_r_y           == nil then problem.boundary_r_y           = {}              end
if problem.boundary_r_y.condition == nil then problem.boundary_r_y.condition = "EXTRAPOLATION" end
if problem.boundary_r_y.rho       == nil then problem.boundary_r_y.rho       = 1.              end
if problem.boundary_r_y.u         == nil then problem.boundary_r_y.u         = 0.              end
if problem.boundary_r_y.v         == nil then problem.boundary_r_y.v         = 0.              end
if problem.boundary_r_y.w         == nil then problem.boundary_r_y.w         = 0.              end
if problem.boundary_r_y.p         == nil then problem.boundary_r_y.p         = 1.              end
if problem.boundary_r_y.L_x       == nil then problem.boundary_r_y.L_x       = 0.1             end
if problem.boundary_r_y.sigma     == nil then problem.boundary_r_y.sigma     = 0.005           end
if problem.boundary_r_y.eta_1     == nil then problem.boundary_r_y.eta_1     = 2.0             end
if problem.boundary_r_y.eta_2     == nil then problem.boundary_r_y.eta_2     = 2.0             end
if problem.boundary_r_y.eta_3     == nil then problem.boundary_r_y.eta_3     = 2.0             end
if problem.boundary_r_y.eta_4     == nil then problem.boundary_r_y.eta_4     = 2.0             end
if problem.boundary_r_y.eta_5     == nil then problem.boundary_r_y.eta_5     = 2.0             end

if problem.boundary_l_z           == nil then problem.boundary_l_z           = {}              end
if problem.boundary_l_z.condition == nil then problem.boundary_l_z.condition = "EXTRAPOLATION" end
if problem.boundary_l_z.rho       == nil then problem.boundary_l_z.rho       = 1.              end
if problem.boundary_l_z.u         == nil then problem.boundary_l_z.u         = 0.              end
if problem.boundary_l_z.v         == nil then problem.boundary_l_z.v         = 0.              end
if problem.boundary_l_z.w         == nil then problem.boundary_l_z.w         = 0.              end
if problem.boundary_l_z.p         == nil then problem.boundary_l_z.p         = 1.              end
if problem.boundary_l_z.L_x       == nil then problem.boundary_l_z.L_x       = 0.1             end
if problem.boundary_l_z.sigma     == nil then problem.boundary_l_z.sigma     = 0.005           end
if problem.boundary_l_z.eta_1     == nil then problem.boundary_l_z.eta_1     = 2.0             end
if problem.boundary_l_z.eta_2     == nil then problem.boundary_l_z.eta_2     = 2.0             end
if problem.boundary_l_z.eta_3     == nil then problem.boundary_l_z.eta_3     = 2.0             end
if problem.boundary_l_z.eta_4     == nil then problem.boundary_l_z.eta_4     = 2.0             end
if problem.boundary_l_z.eta_5     == nil then problem.boundary_l_z.eta_5     = 2.0             end

if problem.boundary_r_z           == nil then problem.boundary_r_z           = {}              end
if problem.boundary_r_z.condition == nil then problem.boundary_r_z.condition = "EXTRAPOLATION" end
if problem.boundary_r_z.rho       == nil then problem.boundary_r_z.rho       = 1.              end
if problem.boundary_r_z.u         == nil then problem.boundary_r_z.u         = 0.              end
if problem.boundary_r_z.v         == nil then problem.boundary_r_z.v         = 0.              end
if problem.boundary_r_z.w         == nil then problem.boundary_r_z.w         = 0.              end
if problem.boundary_r_z.p         == nil then problem.boundary_r_z.p         = 1.              end
if problem.boundary_r_z.L_x       == nil then problem.boundary_r_z.L_x       = 0.1             end
if problem.boundary_r_z.sigma     == nil then problem.boundary_r_z.sigma     = 0.005           end
if problem.boundary_r_z.eta_1     == nil then problem.boundary_r_z.eta_1     = 2.0             end
if problem.boundary_r_z.eta_2     == nil then problem.boundary_r_z.eta_2     = 2.0             end
if problem.boundary_r_z.eta_3     == nil then problem.boundary_r_z.eta_3     = 2.0             end
if problem.boundary_r_z.eta_4     == nil then problem.boundary_r_z.eta_4     = 2.0             end
if problem.boundary_r_z.eta_5     == nil then problem.boundary_r_z.eta_5     = 2.0             end

return problem
