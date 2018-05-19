import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")

local problem = {}

-- Problem specific parameters
problem.gamma    = 1.4
problem.Rgas  = 1.0
problem.viscous = false

problem.x_shock  = 0.
problem.M_v      = 1.
problem.x_vortex = 4.
problem.y_vortex = 0.
problem.R_vortex = 1.

local Mach    = 1.2

local rho_pre = 1.0
local p_pre   = 1.0 / problem.gamma
local u_pre   = - Mach * math.sqrt( problem.gamma * p_pre / rho_pre )

local rho_post = rho_pre * (problem.gamma + 1.)*Mach*Mach / ( 2. + (problem.gamma - 1)*Mach*Mach )
local u_post   = rho_pre * u_pre / rho_post
local p_post   = p_pre * ( 1. + 2*problem.gamma/(problem.gamma+1) * (Mach*Mach - 1.) )

-- Grid dimensions
problem.NX = 512
problem.NY = 512
problem.NZ = 1

-- Periodicity
problem.periodic_x = false
problem.periodic_y = true
problem.periodic_z = true

-- Boundary (if not periodic)
-- condition: DIRICHLET, EXTRAPOLATION, SUBSONIC_INFLOW, SUBSONIC_OUTFLOW
problem.boundary_l_x = { condition="DIRICHLET", rho=rho_post, u=u_post, v=0., w=0., p=p_post }
problem.boundary_r_x = { condition="DIRICHLET", rho=rho_pre,  u=u_pre,  v=0., w=0., p=p_pre  }
problem.boundary_l_y = { condition="EXTRAPOLATION", rho=1., u=0., v=0., w=0., p=1. }
problem.boundary_r_y = { condition="EXTRAPOLATION", rho=1., u=0., v=0., w=0., p=1. }
problem.boundary_l_z = { condition="EXTRAPOLATION", rho=1., u=0., v=0., w=0., p=1. }
problem.boundary_r_z = { condition="EXTRAPOLATION", rho=1., u=0., v=0., w=0., p=1. }

-- Domain size
problem.LX = 40.
problem.LY = 40.
problem.LZ = 1.0

problem.X1 = -30.
problem.Y1 = -20.
problem.Z1 = -0.5

-- Grid spacing
problem.DX = problem.LX / problem.NX
problem.DY = problem.LY / problem.NY
problem.DZ = problem.LZ / problem.NZ

problem.ONEBYDX = 1.0 / problem.DX
problem.ONEBYDY = 1.0 / problem.DY
problem.ONEBYDZ = 1.0 / problem.DZ

problem.timestepping_setting = "CONSTANT_CFL_NUM" -- "CONSTANT_TIME_STEP" / "CONSTANT_CFL_NUM"
problem.dt_or_CFL_num        = 0.5
problem.tstop                = 16.0
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
  c.printf("======================\n")
  c.printf("Pre-shock state:\n")
  c.printf("    rho = %g\n", rho_pre)
  c.printf("    u   = %g\n", u_pre)
  c.printf("    p   = %g\n", p_pre)
  c.printf("Post-shock state:\n")
  c.printf("    rho = %g\n", rho_post)
  c.printf("    u   = %g\n", u_post)
  c.printf("    p   = %g\n", p_post)
  c.printf("======================\n\n")

  for i in coords.ispace do
    var idx = int3d {x = i.x - n_ghosts, y = i.y - n_ghosts, z = i.z - n_ghosts}
    coords[i].x_c = problem.X1 + (idx.x + 0.5) * dx
    coords[i].y_c = problem.Y1 + (idx.y + 0.5) * dy
    coords[i].z_c = problem.Z1 + (idx.z + 0.5) * dz

    if coords[i].x_c <= problem.x_shock then
      r_prim_c[i].rho = rho_post
      r_prim_c[i].u   = u_post
      r_prim_c[i].v   = 0.
      r_prim_c[i].w   = 0.
      r_prim_c[i].p   = p_post
    else
      var x_v : double = coords[i].x_c - problem.x_vortex
      var y_v : double = coords[i].y_c - problem.y_vortex
      var rad : double = cmath.sqrt(x_v*x_v + y_v*y_v)
      var expfactor : double = cmath.exp( 1 - (rad / problem.R_vortex)*(rad / problem.R_vortex) )
      var expfactor_half : double = cmath.sqrt( expfactor )
      var rho_vortex = rho_pre * cmath.pow( 1. - 0.5*( problem.gamma - 1. ) * problem.M_v*problem.M_v * expfactor, 1./(problem.gamma - 1.) )

      r_prim_c[i].rho = rho_vortex
      r_prim_c[i].u   = u_pre - problem.M_v * expfactor_half * y_v
      r_prim_c[i].v   = 0.    + problem.M_v * expfactor_half * x_v
      r_prim_c[i].w   = 0.
      r_prim_c[i].p   = p_pre * cmath.pow(rho_vortex, problem.gamma)
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

return problem
