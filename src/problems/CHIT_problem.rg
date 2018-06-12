import "regent"

local c      = regentlib.c
local cmath  = terralib.includec("math.h")
local PI     = cmath.M_PI

require("fields")

local problem = {}

-- Problem specific parameters
problem.gamma = 1.4  -- Ratio of specific heats
problem.Rgas  = 1.0  -- Gas constant
problem.Mt = 0.6     -- Initial turbulent Mach number
problem.k0 = 4.0     -- Inital peak energy wavenumber
problem.Re = 100.    -- Initial Taylor scale Reynolds number
problem.Pr = 0.7     -- Prandtl number
problem.datafile = "/home/akshays/Data/WCHR/CHIT/setup/CHIT-velocity-grid-independent-k04-N0064.dat"  -- Data file containing initial velocities
problem.viscous = true
problem.conservative_viscous_terms = false

problem.u_rms0  = problem.Mt / cmath.sqrt(3)        -- Initial RMS velocity
problem.lambda0 = 2.0 / problem.k0                  -- Initial Taylor microscale
problem.tau     = problem.lambda0 / problem.u_rms0  -- Eddy turnover time scale
problem.T_ref   = (1.0 / problem.gamma) / problem.Rgas
problem.mu_ref  =  problem.u_rms0 * problem.lambda0 / problem.Re -- Reference viscosity

-- Grid dimensions
problem.NX = 64
problem.NY = 64
problem.NZ = 64

-- Periodicity
problem.periodic_x = true
problem.periodic_y = true
problem.periodic_z = true

-- Boundary (if not periodic)
-- condition: DIRICHLET, EXTRAPOLATION, SUBSONIC_INFLOW, SUBSONIC_OUTFLOW

-- Interpolation scheme to use
-- WCHR, WCNS-LD, WCNS-Z, WCNS-JS
problem.interpolation_scheme = "WCHR"

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
problem.dt_or_CFL_num        = 0.5
problem.tstop                = 4.0 * problem.tau
problem.tviz                 = 0.02 * problem.tau

terra read_grid( nx : &int64, ny : &int64, nz : &int64, f : &c.FILE )
  c.fscanf(f, "%d", nx)
  c.fscanf(f, "%d", ny)
  c.fscanf(f, "%d", nz)
end

terra read_data( ix : &int64, iy : &int64, iz : &int64,
                 u  : &double, v : &double, w : &double,
                 f : &c.FILE )
  c.fscanf(f, "%d", ix)
  c.fscanf(f, "%d", iy)
  c.fscanf(f, "%d", iz)
  c.fscanf(f, "%lf", u)
  c.fscanf(f, "%lf", v)
  c.fscanf(f, "%lf", w)
end

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
    coords[i].x_c = problem.X1 + (i.x - n_ghosts) * dx
    coords[i].y_c = problem.Y1 + (i.y - n_ghosts) * dy
    coords[i].z_c = problem.Z1 + (i.z - n_ghosts) * dz

    r_prim_c[i].rho = 1.0
    r_prim_c[i].p   = 1.0 / problem.gamma
  end

  var nx : int64[1]
  var ny : int64[1]
  var nz : int64[1]
  nx[0], ny[0], nz[0] = 0, 0, 0

  var ix : int64[1]
  var iy : int64[1]
  var iz : int64[1]
  ix[0], iy[0], iz[0] = 0, 0, 0

  var u_dat : double[1]
  var v_dat : double[1]
  var w_dat : double[1]
  u_dat[0], v_dat[0], w_dat[0] = 0., 0., 0.

  var file_handle = c.fopen(problem.datafile, "r")
  read_grid(nx, ny, nz, file_handle)

  regentlib.assert( (nx[0] == problem.NX), "Grid size in datafile does not match problem's grid size in x" )
  regentlib.assert( (ny[0] == problem.NY), "Grid size in datafile does not match problem's grid size in y" )
  regentlib.assert( (nz[0] == problem.NZ), "Grid size in datafile does not match problem's grid size in z" )

  var bound = r_prim_c.bounds
  var counter = 0

  while counter < (bound.hi.x-bound.lo.x+1)*(bound.hi.y-bound.lo.y+1)*(bound.hi.z-bound.lo.z+1) do
    read_data(ix, iy, iz, u_dat, v_dat, w_dat, file_handle)
    if ( (iz[0]+n_ghosts >= bound.lo.z) and (iz[0]+n_ghosts <= bound.hi.z) ) then
      if ( (iy[0]+n_ghosts >= bound.lo.y) and (iy[0]+n_ghosts <= bound.hi.y) ) then
        if ( (ix[0]+n_ghosts >= bound.lo.x) and (ix[0]+n_ghosts <= bound.hi.x) ) then
          r_prim_c[{ix[0]+n_ghosts,iy[0]+n_ghosts,iz[0]+n_ghosts}].u = problem.u_rms0 * u_dat[0]
          r_prim_c[{ix[0]+n_ghosts,iy[0]+n_ghosts,iz[0]+n_ghosts}].v = problem.u_rms0 * v_dat[0]
          r_prim_c[{ix[0]+n_ghosts,iy[0]+n_ghosts,iz[0]+n_ghosts}].w = problem.u_rms0 * w_dat[0]
          counter += 1
        end
      end
    end
  end

  c.fclose(file_handle)

  return 1
end

task problem.get_transport_coeffs( r_prim : region(ispace(int3d), primitive),
                                   r_aux  : region(ispace(int3d), auxiliary),
                                   r_visc : region(ispace(int3d), transport_coeffs) )
where
  reads(r_prim.{}, r_aux.T), writes(r_visc)
do
  for i in r_visc do
    var mu_s : double  = problem.mu_ref * cmath.pow( r_aux[i].T / problem.T_ref, 3./4. ) -- Power law for mu
    r_visc[i].mu_s  = mu_s -- Power law for mu
    r_visc[i].mu_b  = 0.
    r_visc[i].kappa = (1.0/problem.Pr) * (problem.gamma / (problem.gamma - 1.)) * problem.Rgas * mu_s
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
if problem.boundary_l_x.eta_1     == nil then problem.boundary_l_x.eta_1     = 0.5             end
if problem.boundary_l_x.eta_2     == nil then problem.boundary_l_x.eta_2     = 0.5             end
if problem.boundary_l_x.eta_3     == nil then problem.boundary_l_x.eta_3     = 0.5             end
if problem.boundary_l_x.eta_4     == nil then problem.boundary_l_x.eta_4     = 0.5             end
if problem.boundary_l_x.eta_5     == nil then problem.boundary_l_x.eta_5     = 0.5             end

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
if problem.boundary_r_x.eta_1     == nil then problem.boundary_r_x.eta_1     = 0.5             end
if problem.boundary_r_x.eta_2     == nil then problem.boundary_r_x.eta_2     = 0.5             end
if problem.boundary_r_x.eta_3     == nil then problem.boundary_r_x.eta_3     = 0.5             end
if problem.boundary_r_x.eta_4     == nil then problem.boundary_r_x.eta_4     = 0.5             end
if problem.boundary_r_x.eta_5     == nil then problem.boundary_r_x.eta_5     = 0.5             end

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
if problem.boundary_l_y.eta_1     == nil then problem.boundary_l_y.eta_1     = 0.5             end
if problem.boundary_l_y.eta_2     == nil then problem.boundary_l_y.eta_2     = 0.5             end
if problem.boundary_l_y.eta_3     == nil then problem.boundary_l_y.eta_3     = 0.5             end
if problem.boundary_l_y.eta_4     == nil then problem.boundary_l_y.eta_4     = 0.5             end
if problem.boundary_l_y.eta_5     == nil then problem.boundary_l_y.eta_5     = 0.5             end

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
if problem.boundary_r_y.eta_1     == nil then problem.boundary_r_y.eta_1     = 0.5             end
if problem.boundary_r_y.eta_2     == nil then problem.boundary_r_y.eta_2     = 0.5             end
if problem.boundary_r_y.eta_3     == nil then problem.boundary_r_y.eta_3     = 0.5             end
if problem.boundary_r_y.eta_4     == nil then problem.boundary_r_y.eta_4     = 0.5             end
if problem.boundary_r_y.eta_5     == nil then problem.boundary_r_y.eta_5     = 0.5             end

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
if problem.boundary_l_z.eta_1     == nil then problem.boundary_l_z.eta_1     = 0.5             end
if problem.boundary_l_z.eta_2     == nil then problem.boundary_l_z.eta_2     = 0.5             end
if problem.boundary_l_z.eta_3     == nil then problem.boundary_l_z.eta_3     = 0.5             end
if problem.boundary_l_z.eta_4     == nil then problem.boundary_l_z.eta_4     = 0.5             end
if problem.boundary_l_z.eta_5     == nil then problem.boundary_l_z.eta_5     = 0.5             end

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
if problem.boundary_r_z.eta_1     == nil then problem.boundary_r_z.eta_1     = 0.5             end
if problem.boundary_r_z.eta_2     == nil then problem.boundary_r_z.eta_2     = 0.5             end
if problem.boundary_r_z.eta_3     == nil then problem.boundary_r_z.eta_3     = 0.5             end
if problem.boundary_r_z.eta_4     == nil then problem.boundary_r_z.eta_4     = 0.5             end
if problem.boundary_r_z.eta_5     == nil then problem.boundary_r_z.eta_5     = 0.5             end

return problem
