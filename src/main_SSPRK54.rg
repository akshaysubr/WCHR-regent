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

  var r_cnsr     = region(grid_c, conserved)  -- conserved variables at cell center

  var r_prim_c   = region(grid_c, primitive)  -- primitive variables at cell center

  var r_aux_c    = region(grid_c, auxiliary)         -- auxiliary variables at cell center
  var r_visc     = region(grid_c, transport_coeffs)  -- Transport coefficients at cell center

  var r_gradu    = region(grid_c, tensor2)      -- Velocity gradient tensor at the cell center
  var r_grad2u   = region(grid_c, tensor2)      -- Velocity second-derivative tensor at the cell center (d^2 u_i / d x_j^2)
  var r_tauij    = region(grid_c, tensor2symm)  -- Viscous stress tensor at the cell center

  var r_gradrho  = region(grid_c, vect)       -- Density gradient
  var r_gradp    = region(grid_c, vect)       -- Pressure gradient

  var r_cnsr_1 = region(grid_c, conserved) -- conserved variables at cell center for SSP-RK(5,4)
  var r_cnsr_2 = region(grid_c, conserved) -- conserved variables at cell center for SSP-RK(5,4)
  var r_cnsr_3 = region(grid_c, conserved) -- conserved variables at cell center for SSP-RK(5,4)
  var r_cnsr_4 = region(grid_c, conserved) -- conserved variables at cell center for SSP-RK(5,4)

  var r_rhs_0 = region(grid_c, conserved)  -- RHS for time stepping at cell center for SSP-RK(5,4)
  var r_rhs_1 = region(grid_c, conserved)  -- RHS for time stepping at cell center for SSP-RK(5,4)
  var r_rhs_2 = region(grid_c, conserved)  -- RHS for time stepping at cell center for SSP-RK(5,4)
  var r_rhs_3 = region(grid_c, conserved)  -- RHS for time stepping at cell center for SSP-RK(5,4)
  var r_rhs_4 = region(grid_c, conserved)  -- RHS for time stepping at cell center for SSP-RK(5,4)

  -- data structure to hold x derivative LU decomposition
  var mat_x      = region(ispace(int3d, {x = Nx, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For staggered finite difference
  var mat_N_x    = region(ispace(int3d, {x = Nx, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For nodal 1st order finite difference
  var mat2_N_x   = region(ispace(int3d, {x = Nx, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For nodal 2nd order finite difference

  var mat_e_x    = region(ispace(int3d, {x = Nx+1, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For staggered finite difference in flux difference form

  var LU_x       = region(ispace(int3d, {x = Nx, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For staggered finite difference
  var LU_N_x     = region(ispace(int3d, {x = Nx, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For nodal 1st order finite difference
  var LU2_N_x    = region(ispace(int3d, {x = Nx, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For nodal 2nd order finite difference

  var LU_e_x     = region(ispace(int3d, {x = Nx+1, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For staggered finite difference in flux difference form

  -- data structure to hold y derivative LU decomposition
  var mat_y      = region(ispace(int3d, {x = Ny, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For staggered finite difference
  var mat_N_y    = region(ispace(int3d, {x = Ny, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For nodal 1st order finite difference
  var mat2_N_y   = region(ispace(int3d, {x = Ny, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For nodal 2nd order finite difference

  var mat_e_y    = region(ispace(int3d, {x = Ny+1, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For staggered finite difference in flux difference form

  var LU_y       = region(ispace(int3d, {x = Ny, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For staggered finite difference
  var LU_N_y     = region(ispace(int3d, {x = Ny, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For nodal 1st order finite difference
  var LU2_N_y    = region(ispace(int3d, {x = Ny, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For nodal 2nd order finite difference

  var LU_e_y     = region(ispace(int3d, {x = Ny+1, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For staggered finite difference in flux difference form

  -- data structure to hold z derivative LU decomposition
  var mat_z      = region(ispace(int3d, {x = Nz, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For staggered finite difference
  var mat_N_z    = region(ispace(int3d, {x = Nz, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For nodal 1st order finite difference
  var mat2_N_z   = region(ispace(int3d, {x = Nz, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For nodal 2nd order finite difference

  var mat_e_z    = region(ispace(int3d, {x = Nz+1, y = config.prow+2, z = config.pcol+2}), LU_coeffs) -- For staggered finite difference in flux difference form

  var LU_z       = region(ispace(int3d, {x = Nz, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For staggered finite difference
  var LU_N_z     = region(ispace(int3d, {x = Nz, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For nodal 1st order finite difference
  var LU2_N_z    = region(ispace(int3d, {x = Nz, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For nodal 2nd order finite difference

  var LU_e_z     = region(ispace(int3d, {x = Nz+1, y = config.prow+2, z = config.pcol+2}), LU_struct) -- For staggered finite difference in flux difference form

  var pencil          = ispace(int2d, int2d {config.prow+2, config.pcol+2}) -- All pencil partitions including the ghost pencils
  var pencil_interior = ispace(int2d, int2d {config.prow, config.pcol}, int2d {1, 1}) -- Only the interior pencil partitions

  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  --------------------------------------------------------------------------------------------
  --                       PARTITIONING
  --------------------------------------------------------------------------------------------
  
  var p_coords_x       = partition_xpencil_coords(coords, n_ghosts, false, pencil)
  var p_coords_y       = partition_ypencil_coords(coords, n_ghosts, false, pencil)
  var p_coords_z       = partition_zpencil_coords(coords, n_ghosts, false, pencil)

  var p_cnsr_x         = partition_xpencil_cnsr(r_cnsr, n_ghosts, false, pencil)
  var p_cnsr_y         = partition_ypencil_cnsr(r_cnsr, n_ghosts, false, pencil)
  var p_cnsr_z         = partition_zpencil_cnsr(r_cnsr, n_ghosts, false, pencil)

  var p_prim_c_x       = partition_xpencil_prim(r_prim_c, n_ghosts, false, pencil)
  var p_prim_c_y       = partition_ypencil_prim(r_prim_c, n_ghosts, false, pencil)
  var p_prim_c_z       = partition_zpencil_prim(r_prim_c, n_ghosts, false, pencil)

  var p_prim_c_x_wg    = partition_xpencil_prim(r_prim_c, n_ghosts,  true, pencil)
  var p_prim_c_y_wg    = partition_ypencil_prim(r_prim_c, n_ghosts,  true, pencil)
  var p_prim_c_z_wg    = partition_zpencil_prim(r_prim_c, n_ghosts,  true, pencil)

  var p_prim_c_x_wo_wg = partition_with_overlap_xpencil_prim(r_prim_c, 1, n_ghosts, true, pencil)
  var p_prim_c_y_wo_wg = partition_with_overlap_ypencil_prim(r_prim_c, 1, n_ghosts, true, pencil)
  var p_prim_c_z_wo_wg = partition_with_overlap_zpencil_prim(r_prim_c, 1, n_ghosts, true, pencil)

  var p_aux_c_x        = partition_xpencil_aux(r_aux_c, n_ghosts, false, pencil)
  var p_aux_c_y        = partition_ypencil_aux(r_aux_c, n_ghosts, false, pencil)
  var p_aux_c_z        = partition_zpencil_aux(r_aux_c, n_ghosts, false, pencil)

  var p_visc_x         = partition_xpencil_visc(r_visc, n_ghosts, false, pencil)
  var p_visc_y         = partition_ypencil_visc(r_visc, n_ghosts, false, pencil)
  var p_visc_z         = partition_zpencil_visc(r_visc, n_ghosts, false, pencil)

  var p_gradrho_x      = partition_xpencil_vect(r_gradrho, n_ghosts, false, pencil)
  var p_gradrho_y      = partition_ypencil_vect(r_gradrho, n_ghosts, false, pencil)
  var p_gradrho_z      = partition_zpencil_vect(r_gradrho, n_ghosts, false, pencil)

  var p_gradp_x        = partition_xpencil_vect(r_gradp, n_ghosts, false, pencil)
  var p_gradp_y        = partition_ypencil_vect(r_gradp, n_ghosts, false, pencil)
  var p_gradp_z        = partition_zpencil_vect(r_gradp, n_ghosts, false, pencil)

  var p_gradu_x        = partition_xpencil_tnsr2(r_gradu, n_ghosts, false, pencil)
  var p_gradu_y        = partition_ypencil_tnsr2(r_gradu, n_ghosts, false, pencil)
  var p_gradu_z        = partition_zpencil_tnsr2(r_gradu, n_ghosts, false, pencil)

  var p_grad2u_x       = partition_xpencil_tnsr2(r_grad2u, n_ghosts, false, pencil)
  var p_grad2u_y       = partition_ypencil_tnsr2(r_grad2u, n_ghosts, false, pencil)
  var p_grad2u_z       = partition_zpencil_tnsr2(r_grad2u, n_ghosts, false, pencil)

  var p_tauij_x        = partition_xpencil_tnsr2symm(r_tauij, n_ghosts, false, pencil)
  var p_tauij_y        = partition_ypencil_tnsr2symm(r_tauij, n_ghosts, false, pencil)
  var p_tauij_z        = partition_zpencil_tnsr2symm(r_tauij, n_ghosts, false, pencil)

  var p_cnsr_1_x       = partition_xpencil_cnsr(r_cnsr_1, n_ghosts, false, pencil)
  var p_cnsr_1_y       = partition_ypencil_cnsr(r_cnsr_1, n_ghosts, false, pencil)
  var p_cnsr_1_z       = partition_zpencil_cnsr(r_cnsr_1, n_ghosts, false, pencil)

  var p_cnsr_2_x       = partition_xpencil_cnsr(r_cnsr_2, n_ghosts, false, pencil)
  var p_cnsr_2_y       = partition_ypencil_cnsr(r_cnsr_2, n_ghosts, false, pencil)
  var p_cnsr_2_z       = partition_zpencil_cnsr(r_cnsr_2, n_ghosts, false, pencil)

  var p_cnsr_3_x       = partition_xpencil_cnsr(r_cnsr_3, n_ghosts, false, pencil)
  var p_cnsr_3_y       = partition_ypencil_cnsr(r_cnsr_3, n_ghosts, false, pencil)
  var p_cnsr_3_z       = partition_zpencil_cnsr(r_cnsr_3, n_ghosts, false, pencil)

  var p_cnsr_4_x       = partition_xpencil_cnsr(r_cnsr_4, n_ghosts, false, pencil)
  var p_cnsr_4_y       = partition_ypencil_cnsr(r_cnsr_4, n_ghosts, false, pencil)
  var p_cnsr_4_z       = partition_zpencil_cnsr(r_cnsr_4, n_ghosts, false, pencil)

  var p_rhs_0_x          = partition_xpencil_cnsr(r_rhs_0, n_ghosts, false, pencil)
  var p_rhs_0_y          = partition_ypencil_cnsr(r_rhs_0, n_ghosts, false, pencil)
  var p_rhs_0_z          = partition_zpencil_cnsr(r_rhs_0, n_ghosts, false, pencil)

  var p_rhs_1_x          = partition_xpencil_cnsr(r_rhs_1, n_ghosts, false, pencil)
  var p_rhs_1_y          = partition_ypencil_cnsr(r_rhs_1, n_ghosts, false, pencil)
  var p_rhs_1_z          = partition_zpencil_cnsr(r_rhs_1, n_ghosts, false, pencil)

  var p_rhs_2_x          = partition_xpencil_cnsr(r_rhs_2, n_ghosts, false, pencil)
  var p_rhs_2_y          = partition_ypencil_cnsr(r_rhs_2, n_ghosts, false, pencil)
  var p_rhs_2_z          = partition_zpencil_cnsr(r_rhs_2, n_ghosts, false, pencil)

  var p_rhs_3_x          = partition_xpencil_cnsr(r_rhs_3, n_ghosts, false, pencil)
  var p_rhs_3_y          = partition_ypencil_cnsr(r_rhs_3, n_ghosts, false, pencil)
  var p_rhs_3_z          = partition_zpencil_cnsr(r_rhs_3, n_ghosts, false, pencil)

  var p_rhs_4_x          = partition_xpencil_cnsr(r_rhs_4, n_ghosts, false, pencil)
  var p_rhs_4_y          = partition_ypencil_cnsr(r_rhs_4, n_ghosts, false, pencil)
  var p_rhs_4_z          = partition_zpencil_cnsr(r_rhs_4, n_ghosts, false, pencil)

  var p_mat_x          = partition_mat(mat_x, pencil)
  var p_mat_y          = partition_mat(mat_y, pencil)
  var p_mat_z          = partition_mat(mat_z, pencil)

  var p_mat_e_x        = partition_mat(mat_e_x, pencil)
  var p_mat_e_y        = partition_mat(mat_e_y, pencil)
  var p_mat_e_z        = partition_mat(mat_e_z, pencil)

  var p_LU_x           = partition_LU(LU_x, pencil)
  var p_LU_y           = partition_LU(LU_y, pencil)
  var p_LU_z           = partition_LU(LU_z, pencil)

  var p_mat_N_x        = partition_mat(mat_N_x, pencil)
  var p_mat_N_y        = partition_mat(mat_N_y, pencil)
  var p_mat_N_z        = partition_mat(mat_N_z, pencil)

  var p_LU_N_x         = partition_LU(LU_N_x, pencil)
  var p_LU_N_y         = partition_LU(LU_N_y, pencil)
  var p_LU_N_z         = partition_LU(LU_N_z, pencil)

  var p_mat2_N_x       = partition_mat(mat2_N_x, pencil)
  var p_mat2_N_y       = partition_mat(mat2_N_y, pencil)
  var p_mat2_N_z       = partition_mat(mat2_N_z, pencil)

  var p_LU2_N_x        = partition_LU(LU2_N_x, pencil)
  var p_LU2_N_y        = partition_LU(LU2_N_y, pencil)
  var p_LU2_N_z        = partition_LU(LU2_N_z, pencil)

  var p_LU_e_x         = partition_LU(LU_e_x, pencil)
  var p_LU_e_y         = partition_LU(LU_e_y, pencil)
  var p_LU_e_z         = partition_LU(LU_e_z, pencil)

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

  -- Initialize MND derivatives stuff.
  __demand(__parallel)
  for i in pencil_interior do
    get_MND_matrix_fd(p_mat_e_x[i], problem.periodic_x)
  end
  __demand(__parallel)
  for i in pencil_interior do
    get_MND_matrix_fd(p_mat_e_y[i], problem.periodic_y)
  end
  __demand(__parallel)
  for i in pencil_interior do
    get_MND_matrix_fd(p_mat_e_z[i], problem.periodic_z)
  end

  __demand(__parallel)
  for i in pencil_interior do
    token += get_LU_decomposition(p_LU_e_x[i], p_mat_e_x[i])
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_LU_decomposition(p_LU_e_y[i], p_mat_e_y[i])
  end
  __demand(__parallel)
  for i in pencil_interior do
    token += get_LU_decomposition(p_LU_e_z[i], p_mat_e_z[i])
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

    -- if useCFL then
    --   -- Get stable dt.
    --   var new_dt = 1.0e+100
    --   if Nz >= 8 then
    --     __demand(__parallel)
    --     for i in pencil_interior do
    --       new_dt min= get_max_stable_dt_3d(p_prim_c_y[i], dx, dy, dz)
    --     end
    --   elseif Ny >= 8 then
    --     __demand(__parallel)
    --     for i in pencil_interior do
    --       new_dt min= get_max_stable_dt_2d(p_prim_c_y[i], dx, dy)
    --     end
    --   else
    --     __demand(__parallel)
    --     for i in pencil_interior do
    --       new_dt min= get_max_stable_dt_1d(p_prim_c_y[i], dx)
    --     end
    --   end
    --   dt = wait_for_double(new_dt) * CFL_num
    -- else
    --   dt = dt_fix
    -- end

    var max_wave_speed_x : double = 0.0
    var max_wave_speed_y : double = 0.0
    var max_wave_speed_z : double = 0.0

    if Nx >= 8 then
      __demand(__parallel)
      for i in pencil_interior do
        max_wave_speed_x max= get_max_wave_speed_x(p_prim_c_y[i])
      end
    end

    if Ny >= 8 then
      __demand(__parallel)
      for i in pencil_interior do
        max_wave_speed_y max= get_max_wave_speed_y(p_prim_c_y[i])
      end
    end

    if Nz >= 8 then
      __demand(__parallel)
      for i in pencil_interior do
        max_wave_speed_z max= get_max_wave_speed_z(p_prim_c_y[i])
      end
    end

    if useCFL then
      var tau_x : double = max_wave_speed_x/dx
      var tau_y : double = max_wave_speed_y/dy
      var tau_z : double = max_wave_speed_z/dz

      dt = CFL_num/(tau_x + tau_y + tau_z)
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

    var lambda_x : double = 0.0
    var lambda_y : double = 0.0
    var lambda_z : double = 0.0

    if Nx >= 8 then
      lambda_x = dt/dx + dt/dy*max_wave_speed_y/max_wave_speed_x + dt/dz*max_wave_speed_z/max_wave_speed_x
    end

    if Ny >= 8 then
      lambda_y = dt/dx*max_wave_speed_x/max_wave_speed_y + dt/dy + dt/dz*max_wave_speed_z/max_wave_speed_y
    end

    if Nz >=8 then
      lambda_z = dt/dx*max_wave_speed_x/max_wave_speed_z + dt/dy*max_wave_speed_y/max_wave_speed_z + dt/dz
    end

    --------------------------------------------------------------------------------------------
    -- Advance sub-steps.
    --------------------------------------------------------------------------------------------

    -- Advance first sub-step.

    pre_substep( r_prim_c, p_prim_c_y, r_aux_c, p_aux_c_y, r_gradu, p_gradu_y, r_tauij, p_tauij_y, r_visc, p_visc_y, pencil_interior )
    
    get_RHS ( r_rhs_0, p_rhs_0_x, p_rhs_0_y, p_rhs_0_z,
              r_prim_c, p_prim_c_x, p_prim_c_y, p_prim_c_z, p_prim_c_x_wg, p_prim_c_y_wg, p_prim_c_z_wg, p_prim_c_x_wo_wg, p_prim_c_y_wo_wg, p_prim_c_z_wo_wg,
              r_aux_c, p_aux_c_x, p_aux_c_y, p_aux_c_z,
              r_visc, p_visc_x, p_visc_y, p_visc_z,
              r_gradu, p_gradu_x, p_gradu_y, p_gradu_z,
              r_grad2u, p_grad2u_x, p_grad2u_y, p_grad2u_z,
              r_tauij, p_tauij_x, p_tauij_y, p_tauij_z,
              LU_x, p_LU_x, LU_y, p_LU_y, LU_z, p_LU_z,
              LU_N_x, p_LU_N_x, LU_N_y, p_LU_N_y, LU_N_z, p_LU_N_z,
              LU2_N_x, p_LU2_N_x, LU2_N_y, p_LU2_N_y, LU2_N_z, p_LU2_N_z,
              LU_e_x, p_LU_e_x, LU_e_y, p_LU_e_y, LU_e_z, p_LU_e_z,
              pencil_interior, max_wave_speed_x, max_wave_speed_y, max_wave_speed_z, lambda_x, lambda_y, lambda_z )

    -- Update solution in this substep.
    __demand(__parallel)
    for i in pencil_interior do
      set_zero_cnsr( p_cnsr_1_y[i] )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_1_y[i], p_cnsr_y[i], 1.0 )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_1_y[i], p_rhs_0_y[i], 0.391752226571890*dt )
    end

    -- Update the primitive variables.
    __demand(__parallel)
    for i in pencil_interior do
      token += get_primitive_r( p_cnsr_1_y[i], p_prim_c_y[i] )
    end

    token += post_substep( coords, p_coords_x, p_coords_y, p_coords_z,
                           r_prim_c, p_prim_c_x, p_prim_c_y, p_prim_c_z, p_prim_c_x_wg, p_prim_c_y_wg, p_prim_c_z_wg,
                           r_aux_c, p_aux_c_y,
                           r_gradu, p_gradu_x, p_gradu_y, p_gradu_z,
                           r_grad2u, p_grad2u_x, p_grad2u_y, p_grad2u_z,
                           r_gradrho, p_gradrho_x, p_gradrho_y, p_gradrho_z,
                           r_gradp, p_gradp_x, p_gradp_y, p_gradp_z,
                           LU_N_x, p_LU_N_x, LU_N_y, p_LU_N_y, LU_N_z, p_LU_N_z,
                           LU2_N_x, p_LU2_N_x, LU2_N_y, p_LU2_N_y, LU2_N_z, p_LU2_N_z,
                           pencil_interior, tsim, n_ghosts )

    -- Advance second sub-step.

    pre_substep( r_prim_c, p_prim_c_y, r_aux_c, p_aux_c_y, r_gradu, p_gradu_y, r_tauij, p_tauij_y, r_visc, p_visc_y, pencil_interior )

    get_RHS ( r_rhs_1, p_rhs_1_x, p_rhs_1_y, p_rhs_1_z,
              r_prim_c, p_prim_c_x, p_prim_c_y, p_prim_c_z, p_prim_c_x_wg, p_prim_c_y_wg, p_prim_c_z_wg, p_prim_c_x_wo_wg, p_prim_c_y_wo_wg, p_prim_c_z_wo_wg,
              r_aux_c, p_aux_c_x, p_aux_c_y, p_aux_c_z,
              r_visc, p_visc_x, p_visc_y, p_visc_z,
              r_gradu, p_gradu_x, p_gradu_y, p_gradu_z,
              r_grad2u, p_grad2u_x, p_grad2u_y, p_grad2u_z,
              r_tauij, p_tauij_x, p_tauij_y, p_tauij_z,
              LU_x, p_LU_x, LU_y, p_LU_y, LU_z, p_LU_z,
              LU_N_x, p_LU_N_x, LU_N_y, p_LU_N_y, LU_N_z, p_LU_N_z,
              LU2_N_x, p_LU2_N_x, LU2_N_y, p_LU2_N_y, LU2_N_z, p_LU2_N_z,
              LU_e_x, p_LU_e_x, LU_e_y, p_LU_e_y, LU_e_z, p_LU_e_z,
              pencil_interior, max_wave_speed_x, max_wave_speed_y, max_wave_speed_z, lambda_x, lambda_y, lambda_z )

    -- Update solution in this substep.
    __demand(__parallel)
    for i in pencil_interior do
      set_zero_cnsr( p_cnsr_2_y[i] )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_2_y[i], p_cnsr_y[i], 0.444370493651235 )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_2_y[i], p_cnsr_1_y[i], 0.555629506348765 )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_2_y[i], p_rhs_1_y[i], 0.368410593050371*dt )
    end

    -- Update the primitive variables.
    __demand(__parallel)
    for i in pencil_interior do
      token += get_primitive_r( p_cnsr_2_y[i], p_prim_c_y[i] )
    end

    token += post_substep( coords, p_coords_x, p_coords_y, p_coords_z,
                           r_prim_c, p_prim_c_x, p_prim_c_y, p_prim_c_z, p_prim_c_x_wg, p_prim_c_y_wg, p_prim_c_z_wg,
                           r_aux_c, p_aux_c_y,
                           r_gradu, p_gradu_x, p_gradu_y, p_gradu_z,
                           r_grad2u, p_grad2u_x, p_grad2u_y, p_grad2u_z,
                           r_gradrho, p_gradrho_x, p_gradrho_y, p_gradrho_z,
                           r_gradp, p_gradp_x, p_gradp_y, p_gradp_z,
                           LU_N_x, p_LU_N_x, LU_N_y, p_LU_N_y, LU_N_z, p_LU_N_z,
                           LU2_N_x, p_LU2_N_x, LU2_N_y, p_LU2_N_y, LU2_N_z, p_LU2_N_z,
                           pencil_interior, tsim, n_ghosts )

    -- Advance third sub-step.

    pre_substep( r_prim_c, p_prim_c_y, r_aux_c, p_aux_c_y, r_gradu, p_gradu_y, r_tauij, p_tauij_y, r_visc, p_visc_y, pencil_interior )

    get_RHS ( r_rhs_2, p_rhs_2_x, p_rhs_2_y, p_rhs_2_z,
              r_prim_c, p_prim_c_x, p_prim_c_y, p_prim_c_z, p_prim_c_x_wg, p_prim_c_y_wg, p_prim_c_z_wg, p_prim_c_x_wo_wg, p_prim_c_y_wo_wg, p_prim_c_z_wo_wg,
              r_aux_c, p_aux_c_x, p_aux_c_y, p_aux_c_z,
              r_visc, p_visc_x, p_visc_y, p_visc_z,
              r_gradu, p_gradu_x, p_gradu_y, p_gradu_z,
              r_grad2u, p_grad2u_x, p_grad2u_y, p_grad2u_z,
              r_tauij, p_tauij_x, p_tauij_y, p_tauij_z,
              LU_x, p_LU_x, LU_y, p_LU_y, LU_z, p_LU_z,
              LU_N_x, p_LU_N_x, LU_N_y, p_LU_N_y, LU_N_z, p_LU_N_z,
              LU2_N_x, p_LU2_N_x, LU2_N_y, p_LU2_N_y, LU2_N_z, p_LU2_N_z,
              LU_e_x, p_LU_e_x, LU_e_y, p_LU_e_y, LU_e_z, p_LU_e_z,
              pencil_interior, max_wave_speed_x, max_wave_speed_y, max_wave_speed_z, lambda_x, lambda_y, lambda_z )

    -- Update solution in this substep.
    __demand(__parallel)
    for i in pencil_interior do
      set_zero_cnsr( p_cnsr_3_y[i] )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_3_y[i], p_cnsr_y[i], 0.620101851488403 )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_3_y[i], p_cnsr_2_y[i], 0.379898148511597 )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_3_y[i], p_rhs_2_y[i], 0.251891774271694*dt )
    end

    -- Update the primitive variables.
    __demand(__parallel)
    for i in pencil_interior do
      token += get_primitive_r( p_cnsr_3_y[i], p_prim_c_y[i] )
    end

    token += post_substep( coords, p_coords_x, p_coords_y, p_coords_z,
                           r_prim_c, p_prim_c_x, p_prim_c_y, p_prim_c_z, p_prim_c_x_wg, p_prim_c_y_wg, p_prim_c_z_wg,
                           r_aux_c, p_aux_c_y,
                           r_gradu, p_gradu_x, p_gradu_y, p_gradu_z,
                           r_grad2u, p_grad2u_x, p_grad2u_y, p_grad2u_z,
                           r_gradrho, p_gradrho_x, p_gradrho_y, p_gradrho_z,
                           r_gradp, p_gradp_x, p_gradp_y, p_gradp_z,
                           LU_N_x, p_LU_N_x, LU_N_y, p_LU_N_y, LU_N_z, p_LU_N_z,
                           LU2_N_x, p_LU2_N_x, LU2_N_y, p_LU2_N_y, LU2_N_z, p_LU2_N_z,
                           pencil_interior, tsim, n_ghosts )

    -- Advance fourth sub-step.

    pre_substep( r_prim_c, p_prim_c_y, r_aux_c, p_aux_c_y, r_gradu, p_gradu_y, r_tauij, p_tauij_y, r_visc, p_visc_y, pencil_interior )

    get_RHS ( r_rhs_3, p_rhs_3_x, p_rhs_3_y, p_rhs_3_z,
              r_prim_c, p_prim_c_x, p_prim_c_y, p_prim_c_z, p_prim_c_x_wg, p_prim_c_y_wg, p_prim_c_z_wg, p_prim_c_x_wo_wg, p_prim_c_y_wo_wg, p_prim_c_z_wo_wg,
              r_aux_c, p_aux_c_x, p_aux_c_y, p_aux_c_z,
              r_visc, p_visc_x, p_visc_y, p_visc_z,
              r_gradu, p_gradu_x, p_gradu_y, p_gradu_z,
              r_grad2u, p_grad2u_x, p_grad2u_y, p_grad2u_z,
              r_tauij, p_tauij_x, p_tauij_y, p_tauij_z,
              LU_x, p_LU_x, LU_y, p_LU_y, LU_z, p_LU_z,
              LU_N_x, p_LU_N_x, LU_N_y, p_LU_N_y, LU_N_z, p_LU_N_z,
              LU2_N_x, p_LU2_N_x, LU2_N_y, p_LU2_N_y, LU2_N_z, p_LU2_N_z,
              LU_e_x, p_LU_e_x, LU_e_y, p_LU_e_y, LU_e_z, p_LU_e_z,
              pencil_interior, max_wave_speed_x, max_wave_speed_y, max_wave_speed_z, lambda_x, lambda_y, lambda_z )

    -- Update solution in this substep.
    __demand(__parallel)
    for i in pencil_interior do
      set_zero_cnsr( p_cnsr_4_y[i] )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_4_y[i], p_cnsr_y[i], 0.178079954393132 )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_4_y[i], p_cnsr_3_y[i], 0.821920045606868 )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_4_y[i], p_rhs_3_y[i], 0.544974750228521*dt )
    end

    -- Update the primitive variables.
    __demand(__parallel)
    for i in pencil_interior do
      token += get_primitive_r( p_cnsr_4_y[i], p_prim_c_y[i] )
    end

    token += post_substep( coords, p_coords_x, p_coords_y, p_coords_z,
                           r_prim_c, p_prim_c_x, p_prim_c_y, p_prim_c_z, p_prim_c_x_wg, p_prim_c_y_wg, p_prim_c_z_wg,
                           r_aux_c, p_aux_c_y,
                           r_gradu, p_gradu_x, p_gradu_y, p_gradu_z,
                           r_grad2u, p_grad2u_x, p_grad2u_y, p_grad2u_z,
                           r_gradrho, p_gradrho_x, p_gradrho_y, p_gradrho_z,
                           r_gradp, p_gradp_x, p_gradp_y, p_gradp_z,
                           LU_N_x, p_LU_N_x, LU_N_y, p_LU_N_y, LU_N_z, p_LU_N_z,
                           LU2_N_x, p_LU2_N_x, LU2_N_y, p_LU2_N_y, LU2_N_z, p_LU2_N_z,
                           pencil_interior, tsim, n_ghosts )

    -- Advance fifth sub-step.

    pre_substep( r_prim_c, p_prim_c_y, r_aux_c, p_aux_c_y, r_gradu, p_gradu_y, r_tauij, p_tauij_y, r_visc, p_visc_y, pencil_interior )

    get_RHS ( r_rhs_4, p_rhs_4_x, p_rhs_4_y, p_rhs_4_z,
              r_prim_c, p_prim_c_x, p_prim_c_y, p_prim_c_z, p_prim_c_x_wg, p_prim_c_y_wg, p_prim_c_z_wg, p_prim_c_x_wo_wg, p_prim_c_y_wo_wg, p_prim_c_z_wo_wg,
              r_aux_c, p_aux_c_x, p_aux_c_y, p_aux_c_z,
              r_visc, p_visc_x, p_visc_y, p_visc_z,
              r_gradu, p_gradu_x, p_gradu_y, p_gradu_z,
              r_grad2u, p_grad2u_x, p_grad2u_y, p_grad2u_z,
              r_tauij, p_tauij_x, p_tauij_y, p_tauij_z,
              LU_x, p_LU_x, LU_y, p_LU_y, LU_z, p_LU_z,
              LU_N_x, p_LU_N_x, LU_N_y, p_LU_N_y, LU_N_z, p_LU_N_z,
              LU2_N_x, p_LU2_N_x, LU2_N_y, p_LU2_N_y, LU2_N_z, p_LU2_N_z,
              LU_e_x, p_LU_e_x, LU_e_y, p_LU_e_y, LU_e_z, p_LU_e_z,
              pencil_interior, max_wave_speed_x, max_wave_speed_y, max_wave_speed_z, lambda_x, lambda_y, lambda_z )

    -- Update solution in this substep.
    __demand(__parallel)
    for i in pencil_interior do
      set_zero_cnsr( p_cnsr_y[i] )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_y[i], p_cnsr_2_y[i], 0.517231671970585 )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_y[i], p_cnsr_3_y[i], 0.096059710526147 )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_y[i], p_rhs_3_y[i], 0.063692468666290*dt )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_y[i], p_cnsr_4_y[i], 0.386708617503269 )
    end
    __demand(__parallel)
    for i in pencil_interior do
      add_value_cnsr( p_cnsr_y[i], p_rhs_4_y[i], 0.226007483236906*dt )
    end

    -- Update the primitive variables.
    __demand(__parallel)
    for i in pencil_interior do
      token += get_primitive_r( p_cnsr_y[i], p_prim_c_y[i] )
    end

    token += post_substep( coords, p_coords_x, p_coords_y, p_coords_z,
                           r_prim_c, p_prim_c_x, p_prim_c_y, p_prim_c_z, p_prim_c_x_wg, p_prim_c_y_wg, p_prim_c_z_wg,
                           r_aux_c, p_aux_c_y,
                           r_gradu, p_gradu_x, p_gradu_y, p_gradu_z,
                           r_grad2u, p_grad2u_x, p_grad2u_y, p_grad2u_z,
                           r_gradrho, p_gradrho_x, p_gradrho_y, p_gradrho_z,
                           r_gradp, p_gradp_x, p_gradp_y, p_gradp_z,
                           LU_N_x, p_LU_N_x, LU_N_y, p_LU_N_y, LU_N_z, p_LU_N_z,
                           LU2_N_x, p_LU2_N_x, LU2_N_y, p_LU2_N_y, LU2_N_z, p_LU2_N_z,
                           pencil_interior, tsim, n_ghosts )



    -- Update time step.
    step = step + 1
    
    -- Update the simulation time.
    tsim += dt

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

