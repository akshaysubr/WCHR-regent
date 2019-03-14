import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")

local problem = {}

-- Problem specific parameters
problem.gamma   = 1.4  -- Ratio of specific heats
problem.Rgas    = 1.0
problem.viscous = false

-- Grid dimensions
problem.NX = 1
problem.NY = 1
problem.NZ = 201

-- Periodicity
problem.periodic_x = true
problem.periodic_y = true
problem.periodic_z = false

-- Domain size
problem.LX = 1.0
problem.LY = 1.0
problem.LZ = 4.0

problem.X1 = -0.5
problem.Y1 = -0.5
problem.Z1 = 0.0

-- Grid spacing
problem.DX = problem.LX / problem.NX
problem.DY = problem.LY / problem.NY
problem.DZ = problem.LZ / problem.NZ

problem.ONEBYDX = 1.0 / problem.DX
problem.ONEBYDY = 1.0 / problem.DY
problem.ONEBYDZ = 1.0 / problem.DZ

problem.interpolation_scheme     = "WCHR"
problem.Riemann_solver           = "HLLC"
problem.use_flux_difference_form = true
problem.use_positivity_limiter   = true
problem.timestepping_setting     = "CONSTANT_TIME_STEP" -- "CONSTANT_TIME_STEP" / "CONSTANT_CFL_NUM"
problem.dt_or_CFL_num            = 1.0e-6
problem.tstop                    = 1.0e-3
problem.tviz                     = 1.0e-4

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
    
    if (coords[i].z_c > (problem.LZ/2.0 - problem.DZ)) and (coords[i].z_c < (problem.LZ/2.0 + problem.DZ)) then
      r_prim_c[i].rho = 1.0
      r_prim_c[i].u   = 0.0
      r_prim_c[i].v   = 0.0 
      r_prim_c[i].w   = 0.0
      r_prim_c[i].p   = 3200000/problem.DZ*(problem.gamma - 1.0)
    else
      r_prim_c[i].rho = 1.0
      r_prim_c[i].u   = 0.0
      r_prim_c[i].v   = 0.0 
      r_prim_c[i].w   = 0.0
      r_prim_c[i].p   = 4.0e-13
    end
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

-- DEFAULT SCHEME TO USE --
if problem.interpolation_scheme == nil then problem.interpolation_scheme = "WCHR" end

-- DEFAULT POSITIVITY LIMITER LIMITS --
if problem.positivity == nil then
  problem.positivity = { epsilon_rho = 1.e-13, epsilon_p = 1.e-13 }
end

-- DEFAULT VISCOUS TERM FORMULATION TO USE --
if problem.conservative_viscous_terms == nil then problem.conservative_viscous_terms = false end

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
if problem.boundary_l_x.beta      == nil then problem.boundary_l_x.beta      = 0.5             end
if problem.boundary_l_x.eta_1     == nil then problem.boundary_l_x.eta_1     = 0.01            end
if problem.boundary_l_x.eta_2     == nil then problem.boundary_l_x.eta_2     = 0.01            end
if problem.boundary_l_x.eta_3     == nil then problem.boundary_l_x.eta_3     = 0.01            end
if problem.boundary_l_x.eta_4     == nil then problem.boundary_l_x.eta_4     = 0.01            end
if problem.boundary_l_x.eta_5     == nil then problem.boundary_l_x.eta_5     = 0.01            end
if problem.boundary_l_x.condition ~= "CUSTOM" then
  task problem.fill_ghost_cells_l_x( coords   : region(ispace(int3d), coordinates),
                                     r_prim_c : region(ispace(int3d), primitive),
                                     tsim     : double,
                                     n_ghosts : int64 )
  end
end

if problem.boundary_r_x           == nil then problem.boundary_r_x           = {}              end
if problem.boundary_r_x.condition == nil then problem.boundary_r_x.condition = "EXTRAPOLATION" end
if problem.boundary_r_x.rho       == nil then problem.boundary_r_x.rho       = 1.              end
if problem.boundary_r_x.u         == nil then problem.boundary_r_x.u         = 0.              end
if problem.boundary_r_x.v         == nil then problem.boundary_r_x.v         = 0.              end
if problem.boundary_r_x.w         == nil then problem.boundary_r_x.w         = 0.              end
if problem.boundary_r_x.p         == nil then problem.boundary_r_x.p         = 1.              end
if problem.boundary_r_x.L_x       == nil then problem.boundary_r_x.L_x       = 0.1             end
if problem.boundary_r_x.sigma     == nil then problem.boundary_r_x.sigma     = 0.005           end
if problem.boundary_r_x.beta      == nil then problem.boundary_r_x.beta      = 0.5             end
if problem.boundary_r_x.eta_1     == nil then problem.boundary_r_x.eta_1     = 0.01            end
if problem.boundary_r_x.eta_2     == nil then problem.boundary_r_x.eta_2     = 0.01            end
if problem.boundary_r_x.eta_3     == nil then problem.boundary_r_x.eta_3     = 0.01            end
if problem.boundary_r_x.eta_4     == nil then problem.boundary_r_x.eta_4     = 0.01            end
if problem.boundary_r_x.eta_5     == nil then problem.boundary_r_x.eta_5     = 0.01            end
if problem.boundary_r_x.condition ~= "CUSTOM" then
  task problem.fill_ghost_cells_r_x( coords   : region(ispace(int3d), coordinates),
                                     r_prim_c : region(ispace(int3d), primitive),
                                     tsim     : double,
                                     n_ghosts : int64 )
  end
end

if problem.boundary_l_y           == nil then problem.boundary_l_y           = {}              end
if problem.boundary_l_y.condition == nil then problem.boundary_l_y.condition = "EXTRAPOLATION" end
if problem.boundary_l_y.rho       == nil then problem.boundary_l_y.rho       = 1.              end
if problem.boundary_l_y.u         == nil then problem.boundary_l_y.u         = 0.              end
if problem.boundary_l_y.v         == nil then problem.boundary_l_y.v         = 0.              end
if problem.boundary_l_y.w         == nil then problem.boundary_l_y.w         = 0.              end
if problem.boundary_l_y.p         == nil then problem.boundary_l_y.p         = 1.              end
if problem.boundary_l_y.L_x       == nil then problem.boundary_l_y.L_x       = 0.1             end
if problem.boundary_l_y.sigma     == nil then problem.boundary_l_y.sigma     = 0.005           end
if problem.boundary_l_y.beta      == nil then problem.boundary_l_y.beta      = 0.5             end
if problem.boundary_l_y.eta_1     == nil then problem.boundary_l_y.eta_1     = 0.01            end
if problem.boundary_l_y.eta_2     == nil then problem.boundary_l_y.eta_2     = 0.01            end
if problem.boundary_l_y.eta_3     == nil then problem.boundary_l_y.eta_3     = 0.01            end
if problem.boundary_l_y.eta_4     == nil then problem.boundary_l_y.eta_4     = 0.01            end
if problem.boundary_l_y.eta_5     == nil then problem.boundary_l_y.eta_5     = 0.01            end
if problem.boundary_l_y.condition ~= "CUSTOM" then
  task problem.fill_ghost_cells_l_y( coords   : region(ispace(int3d), coordinates),
                                     r_prim_c : region(ispace(int3d), primitive),
                                     tsim     : double,
                                     n_ghosts : int64 )
  end
end

if problem.boundary_r_y           == nil then problem.boundary_r_y           = {}              end
if problem.boundary_r_y.condition == nil then problem.boundary_r_y.condition = "EXTRAPOLATION" end
if problem.boundary_r_y.rho       == nil then problem.boundary_r_y.rho       = 1.              end
if problem.boundary_r_y.u         == nil then problem.boundary_r_y.u         = 0.              end
if problem.boundary_r_y.v         == nil then problem.boundary_r_y.v         = 0.              end
if problem.boundary_r_y.w         == nil then problem.boundary_r_y.w         = 0.              end
if problem.boundary_r_y.p         == nil then problem.boundary_r_y.p         = 1.              end
if problem.boundary_r_y.L_x       == nil then problem.boundary_r_y.L_x       = 0.1             end
if problem.boundary_r_y.sigma     == nil then problem.boundary_r_y.sigma     = 0.005           end
if problem.boundary_r_y.beta      == nil then problem.boundary_r_y.beta      = 0.5             end
if problem.boundary_r_y.eta_1     == nil then problem.boundary_r_y.eta_1     = 0.01            end
if problem.boundary_r_y.eta_2     == nil then problem.boundary_r_y.eta_2     = 0.01            end
if problem.boundary_r_y.eta_3     == nil then problem.boundary_r_y.eta_3     = 0.01            end
if problem.boundary_r_y.eta_4     == nil then problem.boundary_r_y.eta_4     = 0.01            end
if problem.boundary_r_y.eta_5     == nil then problem.boundary_r_y.eta_5     = 0.01            end
if problem.boundary_r_y.condition ~= "CUSTOM" then
  task problem.fill_ghost_cells_r_y( coords   : region(ispace(int3d), coordinates),
                                     r_prim_c : region(ispace(int3d), primitive),
                                     tsim     : double,
                                     n_ghosts : int64 )
  end
end

if problem.boundary_l_z           == nil then problem.boundary_l_z           = {}              end
if problem.boundary_l_z.condition == nil then problem.boundary_l_z.condition = "EXTRAPOLATION" end
if problem.boundary_l_z.rho       == nil then problem.boundary_l_z.rho       = 1.              end
if problem.boundary_l_z.u         == nil then problem.boundary_l_z.u         = 0.              end
if problem.boundary_l_z.v         == nil then problem.boundary_l_z.v         = 0.              end
if problem.boundary_l_z.w         == nil then problem.boundary_l_z.w         = 0.              end
if problem.boundary_l_z.p         == nil then problem.boundary_l_z.p         = 1.              end
if problem.boundary_l_z.L_x       == nil then problem.boundary_l_z.L_x       = 0.1             end
if problem.boundary_l_z.sigma     == nil then problem.boundary_l_z.sigma     = 0.005           end
if problem.boundary_l_z.beta      == nil then problem.boundary_l_z.beta      = 0.5             end
if problem.boundary_l_z.eta_1     == nil then problem.boundary_l_z.eta_1     = 0.01            end
if problem.boundary_l_z.eta_2     == nil then problem.boundary_l_z.eta_2     = 0.01            end
if problem.boundary_l_z.eta_3     == nil then problem.boundary_l_z.eta_3     = 0.01            end
if problem.boundary_l_z.eta_4     == nil then problem.boundary_l_z.eta_4     = 0.01            end
if problem.boundary_l_z.eta_5     == nil then problem.boundary_l_z.eta_5     = 0.01            end
if problem.boundary_l_z.condition ~= "CUSTOM" then
  task problem.fill_ghost_cells_l_z( coords   : region(ispace(int3d), coordinates),
                                     r_prim_c : region(ispace(int3d), primitive),
                                     tsim     : double,
                                     n_ghosts : int64 )
  end
end

if problem.boundary_r_z           == nil then problem.boundary_r_z           = {}              end
if problem.boundary_r_z.condition == nil then problem.boundary_r_z.condition = "EXTRAPOLATION" end
if problem.boundary_r_z.rho       == nil then problem.boundary_r_z.rho       = 1.              end
if problem.boundary_r_z.u         == nil then problem.boundary_r_z.u         = 0.              end
if problem.boundary_r_z.v         == nil then problem.boundary_r_z.v         = 0.              end
if problem.boundary_r_z.w         == nil then problem.boundary_r_z.w         = 0.              end
if problem.boundary_r_z.p         == nil then problem.boundary_r_z.p         = 1.              end
if problem.boundary_r_z.L_x       == nil then problem.boundary_r_z.L_x       = 0.1             end
if problem.boundary_r_z.sigma     == nil then problem.boundary_r_z.sigma     = 0.005           end
if problem.boundary_r_z.beta      == nil then problem.boundary_r_z.beta      = 0.5             end
if problem.boundary_r_z.eta_1     == nil then problem.boundary_r_z.eta_1     = 0.01            end
if problem.boundary_r_z.eta_2     == nil then problem.boundary_r_z.eta_2     = 0.01            end
if problem.boundary_r_z.eta_3     == nil then problem.boundary_r_z.eta_3     = 0.01            end
if problem.boundary_r_z.eta_4     == nil then problem.boundary_r_z.eta_4     = 0.01            end
if problem.boundary_r_z.eta_5     == nil then problem.boundary_r_z.eta_5     = 0.01            end
if problem.boundary_r_z.condition ~= "CUSTOM" then
  task problem.fill_ghost_cells_r_z( coords   : region(ispace(int3d), coordinates),
                                     r_prim_c : region(ispace(int3d), primitive),
                                     tsim     : double,
                                     n_ghosts : int64 )
  end
end

return problem
