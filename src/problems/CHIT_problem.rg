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
problem.datafile = "/home/akshays/Data/WCHR/CHIT/setup/CHIT-velocity-k04-N0064.dat"  -- Data file containing initial velocities

problem.u_rms0  = problem.Mt / cmath.sqrt(3)        -- Initial RMS velocity
problem.lambda0 = 2.0 / problem.k0                  -- Initial Taylor microscale
problem.tau     = problem.lambda0 / problem.u_rms0  -- Eddy turnover time scale
problem.T_ref   = (1.0 / problem.gamma) / problem.Rgas
problem.mu_ref  =  problem.u_rms0 * problem.lambda0 / problem.Re -- Reference viscosity

-- Grid dimensions
problem.NX = 64 
problem.NY = 64 
problem.NZ = 64 

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
problem.tviz                 = 0.1 * problem.tau

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
                         dz         : double )
where
  reads writes(coords, r_prim_c)
do
  for i in coords.ispace do
    coords[i].x_c = problem.X1 + (i.x) * dx
    coords[i].y_c = problem.Y1 + (i.y) * dy
    coords[i].z_c = problem.Z1 + (i.z) * dz

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
    if ( (iz[0] >= bound.lo.z) and (iz[0] <= bound.hi.z) ) then
      if ( (iy[0] >= bound.lo.y) and (iy[0] <= bound.hi.y) ) then
        if ( (ix[0] >= bound.lo.x) and (ix[0] <= bound.hi.x) ) then
          r_prim_c[{ix[0],iy[0],iz[0]}].u = problem.u_rms0 * u_dat[0]
          r_prim_c[{ix[0],iy[0],iz[0]}].v = problem.u_rms0 * v_dat[0]
          r_prim_c[{ix[0],iy[0],iz[0]}].w = problem.u_rms0 * w_dat[0]
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
    r_visc[i].kappa = problem.Pr * (problem.gamma / (problem.gamma - 1.)) * problem.Rgas * mu_s
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
