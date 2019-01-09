import "regent"

local c       = regentlib.c
local cmath   = terralib.includec("math.h")
local PI      = cmath.M_PI
local cstring = terralib.includec("string.h")
local min     = regentlib.fmin

local mapper  = require("load_mapper")

require("fields")
require("derivatives")
require("SOE")
require("RHS")
require("partition")
require("boundary")
local interpolation = require("interpolation")
local use_io = require("IO")

local problem = require("problem")
local Config  = require("config")


terra wait_for(x : int)
  return x
end

terra wait_for_double(x : double)
  return x
end

task main()
  var Nx : int64 = problem.NX
  var Ny : int64 = problem.NY
  var Nz : int64 = problem.NZ

  var n_ghosts : int64 = interpolation.n_ghosts
  var Nx_g : int64 = Nx + 2*n_ghosts
  var Ny_g : int64 = Ny + 2*n_ghosts
  var Nz_g : int64 = Nz + 2*n_ghosts

  var Lx : double = problem.LX
  var Ly : double = problem.LY
  var Lz : double = problem.LZ

  var dx : double = problem.DX
  var dy : double = problem.DY
  var dz : double = problem.DZ
  
  -- Get time settings.
  var timestepping_setting = problem.timestepping_setting
  
  var dt      : double = 0.0
  var dt_fix  : double = 0.0
  var CFL_num : double = 0.0

  var useCFL : bool = true
  if cstring.strcmp(timestepping_setting, "CONSTANT_TIME_STEP") == 0 then
    useCFL = false
    dt_fix = problem.dt_or_CFL_num
  elseif cstring.strcmp(timestepping_setting, "CONSTANT_CFL_NUM") == 0 then
    CFL_num = problem.dt_or_CFL_num
  else
    regentlib.assert( false, "Unknown time stepping setting! Choose between \"CONSTANT_TIME_STEP\" or \"CONSTANT_CFL_NUM\"!" )
  end


  var tstop    : double = problem.tstop
  var tsim     : double = 0.0
  var step     : int64  = 0
  var tviz     : double = problem.tviz
  var vizcount : int    = 0
  var vizcond  : bool   = false

  var config : Config
  config:initialize_from_command(Nx, Ny, Nz)

  c.printf("======================= Problem parameters ======================\n")
  c.printf("           Nx, Ny, Nz           = %d, %d, %d\n", Nx, Ny, Nz)
  c.printf("           Lx, Ly, Lz           = %f, %f, %f\n", Lx, Ly, Lz)
  c.printf("           dx, dy, dz           = %f, %f, %f\n", dx, dy, dz)
  c.printf("           Time stepping method = ")
  c.printf(timestepping_setting)
  c.printf("\n")
  if useCFL then
    c.printf("           CFL_num              = %f\n", CFL_num)
  else
    c.printf("           dt                   = %f\n", dt_fix)
  end
  
  c.printf("           prow, pcol           = %d, %d\n", config.prow, config.pcol)
  if use_io then
    c.printf("           fileIO               = true\n")
    c.printf("           prefix               = %s\n", config.filename_prefix)
  else
    c.printf("           fileIO               = false\n")
  end
  if config.restart then
    c.printf("          restart               = true\n")
  else
    c.printf("          restart               = false\n")
  end

  c.printf("================================================================\n")

  --------------------------------------------------------------------------------------------
  --                       DATA STRUCTURES
  --------------------------------------------------------------------------------------------
  
  var grid_c     = ispace(int3d, {x = Nx_g, y = Ny_g, z = Nz_g})  -- cell center index space (including ghost cells now)

  var grid_e_x   = ispace(int3d, {x = Nx+1, y = Ny_g, z = Nz_g})  -- x cell edge index space (assuming no ghost edges)
  var grid_e_y   = ispace(int3d, {x = Nx_g, y = Ny+1, z = Nz_g})  -- y cell edge index space (assuming no ghost edges)
  var grid_e_z   = ispace(int3d, {x = Nx_g, y = Ny_g, z = Nz+1})  -- z cell edge index space (assuming no ghost edges)

  var coords     = region(grid_c, coordinates)  -- coordinates of cell center

  var r_cnsr     = region(grid_c,   conserved)  -- conserved variables at cell center

  var r_prim_c   = region(grid_c,   primitive)  -- primitive variables at cell center

  var r_aux_c    = region(grid_c,   auxiliary)         -- auxiliary variables at cell center
  var r_visc     = region(grid_c,   transport_coeffs)  -- Transport coefficients at cell center

  var r_gradu    = region(grid_c,   tensor2)      -- Velocity gradient tensor at the cell center
  var r_grad2u   = region(grid_c,   tensor2)      -- Velocity second-derivative tensor at the cell center (d^2 u_i / d x_j^2)
  var r_tauij    = region(grid_c,   tensor2symm)  -- Viscous stress tensor at the cell center

  var r_gradrho  = region(grid_c,   vect)       -- Density gradient
  var r_gradp    = region(grid_c,   vect)       -- Pressure gradient

  var r_rhs      = region(grid_c,   conserved)  -- RHS for time stepping at cell center
  var r_qrhs     = region(grid_c,   conserved)  -- buffer for RK45 time stepping

  -- data structure to hold x derivative LU decomposition
  var mat_x      = region(ispace(int3d, {x = Nx, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For staggered finite difference
  var mat_N_x    = region(ispace(int3d, {x = Nx, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For nodal 1st order finite difference
  var mat2_N_x   = region(ispace(int3d, {x = Nx, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For nodal 2nd order finite difference
  
  var LU_x       = region(ispace(int3d, {x = Nx, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For staggered finite difference
  var LU_N_x     = region(ispace(int3d, {x = Nx, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For nodal 1st order finite difference
  var LU2_N_x    = region(ispace(int3d, {x = Nx, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For nodal 2nd order finite difference
  
  -- data structure to hold y derivative LU decomposition
  var mat_y      = region(ispace(int3d, {x = Ny, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For staggered finite difference
  var mat_N_y    = region(ispace(int3d, {x = Ny, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For nodal 1st order finite difference
  var mat2_N_y   = region(ispace(int3d, {x = Ny, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For nodal 2nd order finite difference
  
  var LU_y       = region(ispace(int3d, {x = Ny, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For staggered finite difference
  var LU_N_y     = region(ispace(int3d, {x = Ny, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For nodal 1st order finite difference
  var LU2_N_y    = region(ispace(int3d, {x = Ny, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For nodal 2nd order finite difference

  -- data structure to hold z derivative LU decomposition
  var mat_z      = region(ispace(int3d, {x = Nz, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For staggered finite difference
  var mat_N_z    = region(ispace(int3d, {x = Nz, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For nodal 1st order finite difference
  var mat2_N_z   = region(ispace(int3d, {x = Nz, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For nodal 2nd order finite difference
  
  var LU_z       = region(ispace(int3d, {x = Nz, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For staggered finite difference
  var LU_N_z     = region(ispace(int3d, {x = Nz, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For nodal 1st order finite difference
  var LU2_N_z    = region(ispace(int3d, {x = Nz, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For nodal 2nd order finite difference

  var pencil = ispace(int2d, int2d {config.prow+2, config.pcol+2}) -- All pencil partitions including the ghost pencils
  var pencil_interior = ispace(int2d, int2d {config.prow, config.pcol}, int2d {1, 1}) -- Only the interior pencil partitions

  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  --------------------------------------------------------------------------------------------
  --                       PARTITIONING
  --------------------------------------------------------------------------------------------
  
  var p_coords_x    = partition_xpencil_coords(coords,    n_ghosts, false, pencil)
  var p_coords_y    = partition_ypencil_coords(coords,    n_ghosts, false, pencil)
  var p_coords_z    = partition_zpencil_coords(coords,    n_ghosts, false, pencil)

  var p_cnsr_x      = partition_xpencil_cnsr(r_cnsr,      n_ghosts, false, pencil)
  var p_cnsr_y      = partition_ypencil_cnsr(r_cnsr,      n_ghosts, false, pencil)
  var p_cnsr_z      = partition_zpencil_cnsr(r_cnsr,      n_ghosts, false, pencil)

  var p_prim_c_x    = partition_xpencil_prim(r_prim_c,    n_ghosts, false, pencil)
  var p_prim_c_y    = partition_ypencil_prim(r_prim_c,    n_ghosts, false, pencil)
  var p_prim_c_z    = partition_zpencil_prim(r_prim_c,    n_ghosts, false, pencil)

  var p_prim_c_x_wg = partition_xpencil_prim(r_prim_c,    n_ghosts,  true, pencil)
  var p_prim_c_y_wg = partition_ypencil_prim(r_prim_c,    n_ghosts,  true, pencil)
  var p_prim_c_z_wg = partition_zpencil_prim(r_prim_c,    n_ghosts,  true, pencil)

  var p_aux_c_x     = partition_xpencil_aux (r_aux_c,     n_ghosts, false, pencil)
  var p_aux_c_y     = partition_ypencil_aux (r_aux_c,     n_ghosts, false, pencil)
  var p_aux_c_z     = partition_zpencil_aux (r_aux_c,     n_ghosts, false, pencil)

  var p_visc_x      = partition_xpencil_visc(r_visc,      n_ghosts, false, pencil)
  var p_visc_y      = partition_ypencil_visc(r_visc,      n_ghosts, false, pencil)
  var p_visc_z      = partition_zpencil_visc(r_visc,      n_ghosts, false, pencil)

  var p_gradrho_x   = partition_xpencil_vect(r_gradrho,   n_ghosts, false, pencil)
  var p_gradrho_y   = partition_ypencil_vect(r_gradrho,   n_ghosts, false, pencil)
  var p_gradrho_z   = partition_zpencil_vect(r_gradrho,   n_ghosts, false, pencil)

  var p_gradp_x     = partition_xpencil_vect(r_gradp,     n_ghosts, false, pencil)
  var p_gradp_y     = partition_ypencil_vect(r_gradp,     n_ghosts, false, pencil)
  var p_gradp_z     = partition_zpencil_vect(r_gradp,     n_ghosts, false, pencil)

  var p_gradu_x     = partition_xpencil_tnsr2(r_gradu,    n_ghosts, false, pencil)
  var p_gradu_y     = partition_ypencil_tnsr2(r_gradu,    n_ghosts, false, pencil)
  var p_gradu_z     = partition_zpencil_tnsr2(r_gradu,    n_ghosts, false, pencil)

  var p_grad2u_x    = partition_xpencil_tnsr2(r_grad2u,   n_ghosts, false, pencil)
  var p_grad2u_y    = partition_ypencil_tnsr2(r_grad2u,   n_ghosts, false, pencil)
  var p_grad2u_z    = partition_zpencil_tnsr2(r_grad2u,   n_ghosts, false, pencil)

  var p_tauij_x     = partition_xpencil_tnsr2symm(r_tauij, n_ghosts, false, pencil)
  var p_tauij_y     = partition_ypencil_tnsr2symm(r_tauij, n_ghosts, false, pencil)
  var p_tauij_z     = partition_zpencil_tnsr2symm(r_tauij, n_ghosts, false, pencil)

  var p_rhs_x       = partition_xpencil_cnsr(r_rhs,      n_ghosts, false, pencil)
  var p_rhs_y       = partition_ypencil_cnsr(r_rhs,      n_ghosts, false, pencil)
  var p_rhs_z       = partition_zpencil_cnsr(r_rhs,      n_ghosts, false, pencil)

  -- var p_qrhs_x      = partition_xpencil_cnsr(r_qrhs,     n_ghosts, false, pencil)
  var p_qrhs_y      = partition_ypencil_cnsr(r_qrhs,     n_ghosts, false, pencil)
  -- var p_qrhs_z      = partition_zpencil_cnsr(r_qrhs,     n_ghosts, false, pencil)

  var p_mat_x       = partition_mat(mat_x, pencil)
  var p_mat_y       = partition_mat(mat_y, pencil)
  var p_mat_z       = partition_mat(mat_z, pencil)

  var p_LU_x        = partition_LU(LU_x, pencil)
  var p_LU_y        = partition_LU(LU_y, pencil)
  var p_LU_z        = partition_LU(LU_z, pencil)

  var p_mat_N_x     = partition_mat(mat_N_x, pencil)
  var p_mat_N_y     = partition_mat(mat_N_y, pencil)
  var p_mat_N_z     = partition_mat(mat_N_z, pencil)

  var p_LU_N_x      = partition_LU(LU_N_x, pencil)
  var p_LU_N_y      = partition_LU(LU_N_y, pencil)
  var p_LU_N_z      = partition_LU(LU_N_z, pencil)

  var p_mat2_N_x    = partition_mat(mat2_N_x, pencil)
  var p_mat2_N_y    = partition_mat(mat2_N_y, pencil)
  var p_mat2_N_z    = partition_mat(mat2_N_z, pencil)

  var p_LU2_N_x     = partition_LU(LU2_N_x, pencil)
  var p_LU2_N_y     = partition_LU(LU2_N_y, pencil)
  var p_LU2_N_z     = partition_LU(LU2_N_z, pencil)

  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  var token : int = 0

  -- Initialize MND derivatives stuff.
  __demand(__parallel)
  for i in pencil_interior do
    get_MND_matrix(p_mat_x[i], problem.periodic_x)
  end
  __demand(__parallel)
  for i in pencil_interior do
    get_MND_matrix(p_mat_y[i], problem.periodic_y)
  end
  __demand(__parallel)
  for i in pencil_interior do
    get_MND_matrix(p_mat_z[i], problem.periodic_z)
  end

  __demand(__parallel)
  for i in pencil_interior do
    token += get_LU_decomposition(p_LU_x[i], p_mat_x[i])
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_LU_decomposition(p_LU_y[i], p_mat_y[i])
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_LU_decomposition(p_LU_z[i], p_mat_z[i])
  end

  -- Initialize nodal first derivative stuff
  __demand(__parallel)
  for i in pencil_interior do
    get_compact_matrix(p_mat_N_x[i], 1, problem.periodic_x)
  end
  __demand(__parallel)
  for i in pencil_interior do
    get_compact_matrix(p_mat_N_y[i], 1, problem.periodic_y)
  end
  __demand(__parallel)
  for i in pencil_interior do
    get_compact_matrix(p_mat_N_z[i], 1, problem.periodic_z)
  end

  __demand(__parallel)
  for i in pencil_interior do
    token += get_LU_decomposition(p_LU_N_x[i], p_mat_N_x[i])
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_LU_decomposition(p_LU_N_y[i], p_mat_N_y[i])
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_LU_decomposition(p_LU_N_z[i], p_mat_N_z[i])
  end

  -- Initialize nodal second derivative stuff
  __demand(__parallel)
  for i in pencil_interior do
    get_compact_matrix(p_mat2_N_x[i], 2, problem.periodic_x)
  end
  __demand(__parallel)
  for i in pencil_interior do
    get_compact_matrix(p_mat2_N_y[i], 2, problem.periodic_y)
  end
  __demand(__parallel)
  for i in pencil_interior do
    get_compact_matrix(p_mat2_N_z[i], 2, problem.periodic_z)
  end

  __demand(__parallel)
  for i in pencil_interior do
    token += get_LU_decomposition(p_LU2_N_x[i], p_mat2_N_x[i])
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_LU_decomposition(p_LU2_N_y[i], p_mat2_N_y[i])
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_LU_decomposition(p_LU2_N_z[i], p_mat2_N_z[i])
  end

  wait_for(token)
  c.printf("Finished LU initialization\n")

  if config.restart then
    -- Restart the simulation from the latest viz dump
    c.printf("Restarting using viz dump\n")
    __demand(__parallel)
    for i in pencil_interior do
      read_coords(p_coords_y[i], config.filename_prefix, n_ghosts, i)
    end
    c.printf("Finished reading in coords\n")

    vizcount = config.restart_count
    __demand(__parallel)
    for i in pencil_interior do
      read_primitive(p_prim_c_y[i], config.filename_prefix, vizcount, n_ghosts, i)
    end
    c.printf("Finished restarting using viz dump %d\n", vizcount)
    tsim = tviz*vizcount
    vizcount += 1
  else
    -- If not restarting, then initialize from the problem initialization task
    __demand(__parallel)
    for i in pencil_interior do
      -- Initialize everything in y decomposition.
      token += problem.initialize(p_coords_y[i], p_prim_c_y[i], dx, dy, dz, n_ghosts)
    end
    c.printf("Finished initialization\n")
  end

  __demand(__parallel)
  for i in pencil_interior do
    -- Get temperature from the initial primitive variables
    token += get_temperature_r( p_prim_c_y[i], p_aux_c_y[i] )
  end

  -- Get the velocity derivatives at initial condition 
  __demand(__parallel)
  for i in pencil_interior do
    token += get_velocity_x_derivatives( p_prim_c_x[i], p_gradu_x[i], p_grad2u_x[i], p_LU_N_x[i], p_LU2_N_x[i] )
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_velocity_y_derivatives( p_prim_c_y[i], p_gradu_y[i], p_grad2u_y[i], p_LU_N_y[i], p_LU2_N_y[i] )
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_velocity_z_derivatives( p_prim_c_z[i], p_gradu_z[i], p_grad2u_z[i], p_LU_N_z[i], p_LU2_N_z[i] )
  end
 
  -- Get the density derivatives at initial condition 
  __demand(__parallel)
  for i in pencil_interior do
    token += get_density_x_derivatives( p_prim_c_x[i], p_gradrho_x[i], p_LU_N_x[i] )
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_density_y_derivatives( p_prim_c_y[i], p_gradrho_y[i], p_LU_N_y[i] )
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_density_z_derivatives( p_prim_c_z[i], p_gradrho_z[i], p_LU_N_z[i] )
  end
 
  -- Get the pressure derivatives at initial condition 
  __demand(__parallel)
  for i in pencil_interior do
    token += get_pressure_x_derivatives( p_prim_c_x[i], p_gradp_x[i], p_LU_N_x[i] )
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_pressure_y_derivatives( p_prim_c_y[i], p_gradp_y[i], p_LU_N_y[i] )
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_pressure_z_derivatives( p_prim_c_z[i], p_gradp_z[i], p_LU_N_z[i] )
  end
 
  -- Fill ghost cells in non-periodic directions first
  if not problem.periodic_x then
    __demand(__parallel)
    for i in pencil_interior do
      -- Fill in ghost cells
      nonperiodic_ghost_cells_x(p_coords_x[i], p_prim_c_x_wg[i], p_gradrho_x[i], p_gradu_x[i], p_gradp_x[i], tsim, n_ghosts)
    end
  end
  if not problem.periodic_y then
    __demand(__parallel)
    for i in pencil_interior do
      -- Fill in ghost cells
      nonperiodic_ghost_cells_y(p_coords_y[i], p_prim_c_y_wg[i], tsim, n_ghosts)
    end
  end
  if not problem.periodic_z then
    __demand(__parallel)
    for i in pencil_interior do
      -- Fill in ghost cells
      nonperiodic_ghost_cells_z(p_coords_z[i], p_prim_c_z_wg[i], tsim, n_ghosts)
    end
  end

  -- Fill ghost cells in periodic directions next
  if problem.periodic_x then
    __demand(__parallel)
    for i in pencil_interior do
      -- Fill in ghost cells
      periodic_ghost_cells_x(p_prim_c_x_wg[i], n_ghosts)
    end
  end
  if problem.periodic_y then
    __demand(__parallel)
    for i in pencil_interior do
      -- Fill in ghost cells
      periodic_ghost_cells_y(p_prim_c_y_wg[i], n_ghosts)
    end
  end
  if problem.periodic_z then
    __demand(__parallel)
    for i in pencil_interior do
      -- Fill in ghost cells
      periodic_ghost_cells_z(p_prim_c_z_wg[i], n_ghosts)
    end
  end

  var TKE0 : double = 1.0e-16
  do
    var t : double = 0.0
    __demand(__parallel)
    for i in pencil_interior do
      t += problem.TKE(p_prim_c_y[i])
    end
    TKE0 = wait_for_double(t)
  end

  var enstrophy0 : double = 1.0e-16
  -- do
  --   var e : double = 0.0
  --   __demand(__parallel)
  --   for i in pencil do
  --     e += problem.enstrophy( p_gradu_y[i] )
  --   end
  --   enstrophy0 = wait_for_double(e)
  -- end
  c.printf("TKE, Enstrophy = %g, %g", TKE0, enstrophy0)

  var IOtoken = 0
  if (use_io and (not config.restart)) then
    __demand(__parallel)
    for i in pencil_interior do
      IOtoken += write_coords(p_coords_y[i], config.filename_prefix, n_ghosts, i)
    end
  end
  
  var A_RK45 = array(0.0,
                     -6234157559845.0/12983515589748.0,
                     -6194124222391.0/4410992767914.0,
                     -31623096876824.0/15682348800105.0,
                     -12251185447671.0/11596622555746.0 )

  var B_RK45 = array( 494393426753.0/4806282396855.0,
                      4047970641027.0/5463924506627.0,
                      9795748752853.0/13190207949281.0,
                      4009051133189.0/8539092990294.0,
                      1348533437543.0/7166442652324.0 )
  
  -- Get conserved variables after initialization.
  __demand(__parallel)
  for i in pencil_interior do
    token += get_conserved_r(p_prim_c_y[i], p_cnsr_y[i])
  end

  if (use_io and (not config.restart)) then
    wait_for(IOtoken)
    __demand(__parallel)
    for i in pencil_interior do
      IOtoken += write_primitive(p_prim_c_y[i], config.filename_prefix, vizcount, n_ghosts, i)
    end
    vizcount = vizcount + 1
  end
  
  wait_for(token)
  wait_for(IOtoken)
  var t_start = c.legion_get_current_time_in_micros()

  __demand(__spmd)
  while tsim < tstop*(1.0 - 1.0e-16) do

    if useCFL then
      -- Get stable dt.
      var new_dt = 1.0e+100
      if Nz >= 8 then
        __demand(__parallel)
        for i in pencil_interior do
          new_dt min= get_max_stable_dt_3d(p_prim_c_y[i], dx, dy, dz)
        end
      elseif Ny >= 8 then
        __demand(__parallel)
        for i in pencil_interior do
          new_dt min= get_max_stable_dt_2d(p_prim_c_y[i], dx, dy)
        end
      else
        __demand(__parallel)
        for i in pencil_interior do
          new_dt min= get_max_stable_dt_1d(p_prim_c_y[i], dx)
        end
      end
      dt = wait_for_double(new_dt) * CFL_num
    else
      dt = dt_fix
    end

    if use_io then
      -- Check if viz dump is imminent
      if tsim + dt >= tviz*vizcount then
        dt = tviz * vizcount - tsim
        vizcond = true
      end
    end

    var Q_t : double = 0.0

    __demand(__parallel)
    for i in pencil_interior do
      set_zero_cnsr( p_qrhs_y[i] )
    end

    --------------------------------------------------------------------------------------------
    -- Advance sub-steps.
    --------------------------------------------------------------------------------------------
    for isub = 0,5 do

      -- Set RHS to zero.
      __demand(__parallel)
      for i in pencil_interior do
        set_zero_cnsr( p_rhs_y[i] )
      end
      
      if problem.viscous then
        -- Get the transport coefficients.
        __demand(__parallel)
        for i in pencil_interior do
          problem.get_transport_coeffs( p_prim_c_y[i], p_aux_c_y[i], p_visc_y[i] )
        end

        -- Get the viscous stress tensor.
        __demand(__parallel)
        for i in pencil_interior do
          get_tauij( p_gradu_y[i], p_tauij_y[i], p_visc_y[i] )
        end
      end

      -- Add x-direction convective flux derivative to RHS.
      __demand(__parallel)
      for i in pencil_interior do
        add_xflux_der_to_rhs( p_prim_c_x_wg[i], p_rhs_x[i], p_LU_x[i] )
      end

      -- Add x-direction viscous flux derivative to RHS.
      if problem.conservative_viscous_terms then
        __demand(__parallel)
        for i in pencil_interior do
          add_viscous_xflux_der_to_rhs( p_prim_c_x[i], p_aux_c_x[i], p_visc_x[i], p_tauij_x[i], p_rhs_x[i], p_LU_N_x[i] )
        end
      else
        __demand(__parallel)
        for i in pencil_interior do
          add_nonconservative_viscous_xflux_der_to_rhs( p_prim_c_x[i], p_aux_c_x[i], p_visc_x[i], p_gradu_x[i], p_grad2u_x[i], 
                                                        p_tauij_x[i], p_rhs_x[i], p_LU_N_x[i], p_LU2_N_x[i] )
        end
      end

      -- Add y-direction convective flux derivative to RHS.
      __demand(__parallel)
      for i in pencil_interior do
        add_yflux_der_to_rhs( p_prim_c_y_wg[i], p_rhs_y[i], p_LU_y[i] )
      end

      -- Add y-direction viscous flux derivative to RHS.
      if problem.conservative_viscous_terms then
        __demand(__parallel)
        for i in pencil_interior do
          add_viscous_yflux_der_to_rhs( p_prim_c_y[i], p_aux_c_y[i], p_visc_y[i], p_tauij_y[i], p_rhs_y[i], p_LU_N_y[i] )
        end
      else
        __demand(__parallel)
        for i in pencil_interior do
          add_nonconservative_viscous_yflux_der_to_rhs( p_prim_c_y[i], p_aux_c_y[i], p_visc_y[i], p_gradu_y[i], p_grad2u_y[i], 
                                                        p_tauij_y[i], p_rhs_y[i], p_LU_N_y[i], p_LU2_N_y[i] )
        end
      end

      -- Add z-direction convective flux derivative to RHS.
      __demand(__parallel)
      for i in pencil_interior do
        add_zflux_der_to_rhs( p_prim_c_z_wg[i], p_rhs_z[i], p_LU_z[i] )
      end

      -- Add z-direction viscous flux derivative to RHS.
      if problem.conservative_viscous_terms then
        __demand(__parallel)
        for i in pencil_interior do
          add_viscous_zflux_der_to_rhs( p_prim_c_z[i], p_aux_c_z[i], p_visc_z[i], p_tauij_z[i], p_rhs_z[i], p_LU_N_z[i] )
        end
      else
        __demand(__parallel)
        for i in pencil_interior do
          add_nonconservative_viscous_zflux_der_to_rhs( p_prim_c_z[i], p_aux_c_z[i], p_visc_z[i], p_gradu_z[i], p_grad2u_z[i], 
                                                        p_tauij_z[i], p_rhs_z[i], p_LU_N_z[i], p_LU2_N_z[i] )
        end
      end

      -- Update solution in this substep.
      __demand(__parallel)
      for i in pencil_interior do
        update_substep( p_cnsr_y[i], p_rhs_y[i], p_qrhs_y[i], dt, A_RK45[isub], B_RK45[isub] )
      end

      -- Update simulation time as well.
      Q_t = dt + A_RK45[isub]*Q_t
      tsim += B_RK45[isub]*Q_t

      -- Update the primitive variables.
      __demand(__parallel)
      for i in pencil_interior do
        token += get_primitive_r(p_cnsr_y[i], p_prim_c_y[i])
      end

      -- Update temperature.
      __demand(__parallel)
      for i in pencil_interior do
        token += get_temperature_r( p_prim_c_y[i], p_aux_c_y[i] )
      end

      -- Update velocity gradient tensor.
      __demand(__parallel)
      for i in pencil_interior do
        token += get_velocity_x_derivatives( p_prim_c_x[i], p_gradu_x[i], p_grad2u_x[i], p_LU_N_x[i], p_LU2_N_x[i] )
      end
      __demand(__parallel)
      for i in pencil_interior do
        token += get_velocity_y_derivatives( p_prim_c_y[i], p_gradu_y[i], p_grad2u_y[i], p_LU_N_y[i], p_LU2_N_y[i] )
      end
      __demand(__parallel)
      for i in pencil_interior do
        token += get_velocity_z_derivatives( p_prim_c_z[i], p_gradu_z[i], p_grad2u_z[i], p_LU_N_z[i], p_LU2_N_z[i] )
      end

      -- Get the density derivatives
      __demand(__parallel)
      for i in pencil_interior do
        token += get_density_x_derivatives( p_prim_c_x[i], p_gradrho_x[i], p_LU_N_x[i] )
      end
      __demand(__parallel)
      for i in pencil_interior do
        token += get_density_y_derivatives( p_prim_c_y[i], p_gradrho_y[i], p_LU_N_y[i] )
      end
      __demand(__parallel)
      for i in pencil_interior do
        token += get_density_z_derivatives( p_prim_c_z[i], p_gradrho_z[i], p_LU_N_z[i] )
      end
 
      -- Get the pressure derivatives
      __demand(__parallel)
      for i in pencil_interior do
        token += get_pressure_x_derivatives( p_prim_c_x[i], p_gradp_x[i], p_LU_N_x[i] )
      end
      __demand(__parallel)
      for i in pencil_interior do
        token += get_pressure_y_derivatives( p_prim_c_y[i], p_gradp_y[i], p_LU_N_y[i] )
      end
      __demand(__parallel)
      for i in pencil_interior do
        token += get_pressure_z_derivatives( p_prim_c_z[i], p_gradp_z[i], p_LU_N_z[i] )
      end
 
      -- Fill ghost cells in non-periodic directions first
      if not problem.periodic_x then
        __demand(__parallel)
        for i in pencil_interior do
          -- Fill in ghost cells
          nonperiodic_ghost_cells_x(p_coords_x[i], p_prim_c_x_wg[i], p_gradrho_x[i], p_gradu_x[i], p_gradp_x[i], tsim, n_ghosts)
        end
      end
      if not problem.periodic_y then
        __demand(__parallel)
        for i in pencil_interior do
          -- Fill in ghost cells
          nonperiodic_ghost_cells_y(p_coords_y[i], p_prim_c_y_wg[i], tsim, n_ghosts)
        end
      end
      if not problem.periodic_z then
        __demand(__parallel)
        for i in pencil_interior do
          -- Fill in ghost cells
          nonperiodic_ghost_cells_z(p_coords_z[i], p_prim_c_z_wg[i], tsim, n_ghosts)
        end
      end

      -- Fill ghost cells in periodic directions next
      if problem.periodic_x then
        __demand(__parallel)
        for i in pencil_interior do
          -- Fill in ghost cells
          periodic_ghost_cells_x(p_prim_c_x_wg[i], n_ghosts)
        end
      end
      if problem.periodic_y then
        __demand(__parallel)
        for i in pencil_interior do
          -- Fill in ghost cells
          periodic_ghost_cells_y(p_prim_c_y_wg[i], n_ghosts)
        end
      end
      if problem.periodic_z then
        __demand(__parallel)
        for i in pencil_interior do
          -- Fill in ghost cells
          periodic_ghost_cells_z(p_prim_c_z_wg[i], n_ghosts)
        end
      end

    end

    -- Update time step.
    step = step + 1

    if (step-1)%(config.nstats*50) == 0 then
      c.printf("\n")
      c.printf("%6.6s |%16.16s |%16.16s |%16.16s |%16.16s\n", "Step","Time","Timestep","TKE","Enstrophy")
      c.printf("-------|-----------------|-----------------|-----------------|----------------\n")
    end

    if (step-1)%config.nstats == 0 then
      var TKE : double = 0.0
      __demand(__parallel)
      for i in pencil_interior do
        TKE += problem.TKE(p_prim_c_y[i])
      end

      var enstrophy : double = 0.0
      -- __demand(__parallel)
      -- for i in pencil do
      --   enstrophy += problem.enstrophy( p_gradu_y[i] )
      -- end

      do
        var TKE = wait_for_double(TKE)
        var enstrophy = wait_for_double(enstrophy)
        c.printf("%6d |%16.8e |%16.8e |%16.8e |%16.8e\n", step, tsim, dt, TKE/TKE0, enstrophy/enstrophy0)
      end
    end

    if use_io then
      if vizcond then
        wait_for(IOtoken)
        __demand(__parallel)
        for i in pencil_interior do
          IOtoken += write_primitive(p_prim_c_y[i], config.filename_prefix, vizcount, n_ghosts, i)
        end
        vizcount = vizcount + 1
        vizcond = false
      end
    end

  end
  
  wait_for(token)
  var t_simulation = c.legion_get_current_time_in_micros() - t_start
  
  var errors : double[5]
  for ierr = 0,5 do
    errors[ierr] = 0.0
  end
  for i in pencil_interior do
    var perrors = problem.get_errors(p_coords_y[i], p_prim_c_y[i], tsim)
    for ierr = 0,5 do
      if perrors[ierr] > errors[ierr] then
        errors[ierr] = perrors[ierr]
      end
    end
  end

  c.printf("\n")
  c.printf("Error in rho = %g\n", errors[0])
  c.printf("Error in u   = %g\n", errors[1])
  c.printf("Error in v   = %g\n", errors[2])
  c.printf("Error in w   = %g\n", errors[3])
  c.printf("Error in p   = %g\n", errors[4])

  c.printf("\n")
  c.printf("Average time per time step = %12.5e\n", (t_simulation)*1e-6/step)
  
end

if os.getenv('SAVEOBJ') == '1' then
  local root_dir = arg[0]:match(".*/") or "./"
  local link_flags = {}
  if use_io then
    local hdf_root = os.getenv('HDF_ROOT')
    local hdf_lib_dir = hdf_root .. "/lib"
    link_flags = {"-L" .. root_dir, "-L" .. hdf_lib_dir, "-lhdf5", "-lm", "-lblas", "-L./", "-lwchr"}
  else
    link_flags = {"-L" .. root_dir, "-lm", "-lblas", "-L./", "-lwchr"}
  end
  print("Saving executable to ./wchr")
  regentlib.saveobj(main, "wchr", "executable", mapper.register_mappers, link_flags)
else
  regentlib.start(main, mapper.register_mappers)
end
