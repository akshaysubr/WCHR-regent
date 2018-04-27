import "regent"

local c       = regentlib.c
local cmath   = terralib.includec("math.h")
local PI      = cmath.M_PI
local cstring = terralib.includec("string.h")
local min     = regentlib.fmin

require("fields")
require("derivatives")
require("SOE")
require("RHS")
require("partition")
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
  
  var grid_c     = ispace(int3d, {x = Nx,   y = Ny,   z = Nz  })  -- cell center index space

  var ghost_x    = ispace(int3d, {x =  4,   y = Ny,   z = Nz  })  -- x ghost cells
  var ghost_y    = ispace(int3d, {x = Nx,   y =  4,   z = Nz  })  -- y ghost cells
  var ghost_z    = ispace(int3d, {x = Nx,   y = Ny,   z =  4  })  -- y ghost cells

  var grid_e_x   = ispace(int3d, {x = Nx+1, y = Ny,   z = Nz  })  -- x cell edge index space
  var grid_e_y   = ispace(int3d, {x = Nx,   y = Ny+1, z = Nz  })  -- y cell edge index space
  var grid_e_z   = ispace(int3d, {x = Nx,   y = Ny,   z = Nz+1})  -- z cell edge index space

  var coords     = region(grid_c, coordinates)  -- coordinates of cell center

  var r_cnsr     = region(grid_c,   conserved)  -- conserved variables at cell center

  var r_prim_c   = region(grid_c,   primitive)  -- primitive variables at cell center
  var r_prim_l_x = region(grid_e_x, primitive)  -- primitive variables at left x cell edge
  var r_prim_l_y = region(grid_e_y, primitive)  -- primitive variables at left y cell edge
  var r_prim_l_z = region(grid_e_z, primitive)  -- primitive variables at left z cell edge
  var r_prim_r_x = region(grid_e_x, primitive)  -- primitive variables at right x cell edge
  var r_prim_r_y = region(grid_e_y, primitive)  -- primitive variables at right y cell edge
  var r_prim_r_z = region(grid_e_z, primitive)  -- primitive variables at right z cell edge

  var gr_prim_l_x = region(ghost_x, primitive)  -- primitive variables at left x ghost cells
  var gr_prim_l_y = region(ghost_y, primitive)  -- primitive variables at left y ghost cells
  var gr_prim_l_z = region(ghost_z, primitive)  -- primitive variables at left z ghost cells
  var gr_prim_r_x = region(ghost_x, primitive)  -- primitive variables at right x ghost cells
  var gr_prim_r_y = region(ghost_y, primitive)  -- primitive variables at right y ghost cells
  var gr_prim_r_z = region(ghost_z, primitive)  -- primitive variables at right z ghost cells

  var r_aux_c    = region(grid_c,   auxiliary)         -- auxiliary variables at cell center
  var r_visc     = region(grid_c,   transport_coeffs)  -- Transport coefficients at cell center

  var r_gradu    = region(grid_c,   tensor2)      -- Velocity gradient tensor at the cell center
  var r_tauij    = region(grid_c,   tensor2symm)  -- Viscous stress tensor at the cell center
  var r_q        = region(grid_c,   vect)         -- Heat flux vector at the cell center

  var r_rhs_l_x  = region(grid_e_x, primitive)  -- store RHS for left interpolation in x
  var r_rhs_r_x  = region(grid_e_x, primitive)  -- store RHS for right interpolation in x
  var r_rhs_l_y  = region(grid_e_y, primitive)  -- store RHS for left interpolation in y
  var r_rhs_r_y  = region(grid_e_y, primitive)  -- store RHS for right interpolation in y
  var r_rhs_l_z  = region(grid_e_z, primitive)  -- store RHS for left interpolation in z
  var r_rhs_r_z  = region(grid_e_z, primitive)  -- store RHS for right interpolation in z

  var r_flux_c   = region(grid_c,   conserved)  -- flux at cell center
  var r_flux_e_x = region(grid_e_x, conserved)  -- flux at x cell edge
  var r_flux_e_y = region(grid_e_y, conserved)  -- flux at y cell edge
  var r_flux_e_z = region(grid_e_z, conserved)  -- flux at z cell edge
  
  var gr_flux_l_x = region(ghost_x, conserved)  -- flux at left x ghost cells
  var gr_flux_l_y = region(ghost_y, conserved)  -- flux at left y ghost cells
  var gr_flux_l_z = region(ghost_z, conserved)  -- flux at left z ghost cells
  var gr_flux_r_x = region(ghost_x, conserved)  -- flux at right x ghost cells
  var gr_flux_r_y = region(ghost_y, conserved)  -- flux at right y ghost cells
  var gr_flux_r_z = region(ghost_z, conserved)  -- flux at right z ghost cells

  var r_fder_c_x = region(grid_c,   conserved)  -- x flux derivative
  var r_fder_c_y = region(grid_c,   conserved)  -- y flux derivative
  var r_fder_c_z = region(grid_c,   conserved)  -- z flux derivative
  
  var r_rhs      = region(grid_c,   conserved)  -- RHS for time stepping at cell center
  var r_qrhs     = region(grid_c,   conserved)  -- buffer for RK45 time stepping

  -- data structure to hold x derivative LU decomposition
  var LU_x       = region(ispace(int3d, {x = Nx, y = config.prow, z = config.pcol}), LU_struct) -- For staggered finite difference
  var LU_N_x     = region(ispace(int3d, {x = Nx, y = config.prow, z = config.pcol}), LU_struct) -- For nodal finite difference
  
  -- data structure to hold y derivative LU decomposition
  var LU_y       = region(ispace(int3d, {x = Ny, y = config.prow, z = config.pcol}), LU_struct) -- For staggered finite difference
  var LU_N_y     = region(ispace(int3d, {x = Ny, y = config.prow, z = config.pcol}), LU_struct) -- For nodal finite difference

  -- data structure to hold z derivative LU decomposition
  var LU_z       = region(ispace(int3d, {x = Nz, y = config.prow, z = config.pcol}), LU_struct) -- For staggered finite difference
  var LU_N_z     = region(ispace(int3d, {x = Nz, y = config.prow, z = config.pcol}), LU_struct) -- For nodal finite difference

  var pencil = ispace(int2d, int2d {config.prow, config.pcol})

  var alpha_l_x = region( ispace(int3d, {Nx+1, Ny, Nz} ), coeffs )
  var beta_l_x  = region( ispace(int3d, {Nx+1, Ny, Nz} ), coeffs )
  var gamma_l_x = region( ispace(int3d, {Nx+1, Ny, Nz} ), coeffs )

  var alpha_r_x = region( ispace(int3d, {Nx+1, Ny, Nz} ), coeffs )
  var beta_r_x  = region( ispace(int3d, {Nx+1, Ny, Nz} ), coeffs )
  var gamma_r_x = region( ispace(int3d, {Nx+1, Ny, Nz} ), coeffs )

  var rho_avg_x = region( ispace(int3d, {Nx+1, Ny, Nz} ), double )
  var sos_avg_x = region( ispace(int3d, {Nx+1, Ny, Nz} ), double )

  var block_d_x    = region( ispace(int3d, {Nx+1, Ny, Nz} ), double[9] )
  var block_Uinv_x = region( ispace(int3d, {Nx+1, Ny, Nz} ), double[9] )

  var alpha_l_y = region( ispace(int3d, {Nx, Ny+1, Nz} ), coeffs )
  var beta_l_y  = region( ispace(int3d, {Nx, Ny+1, Nz} ), coeffs )
  var gamma_l_y = region( ispace(int3d, {Nx, Ny+1, Nz} ), coeffs )

  var alpha_r_y = region( ispace(int3d, {Nx, Ny+1, Nz} ), coeffs )
  var beta_r_y  = region( ispace(int3d, {Nx, Ny+1, Nz} ), coeffs )
  var gamma_r_y = region( ispace(int3d, {Nx, Ny+1, Nz} ), coeffs )

  var rho_avg_y = region( ispace(int3d, {Nx, Ny+1, Nz} ), double )
  var sos_avg_y = region( ispace(int3d, {Nx, Ny+1, Nz} ), double )

  var block_d_y    = region( ispace(int3d, {Nx, Ny+1, Nz} ), double[9] )
  var block_Uinv_y = region( ispace(int3d, {Nx, Ny+1, Nz} ), double[9] )

  var alpha_l_z = region( ispace(int3d, {Nx, Ny, Nz+1} ), coeffs )
  var beta_l_z  = region( ispace(int3d, {Nx, Ny, Nz+1} ), coeffs )
  var gamma_l_z = region( ispace(int3d, {Nx, Ny, Nz+1} ), coeffs )

  var alpha_r_z = region( ispace(int3d, {Nx, Ny, Nz+1} ), coeffs )
  var beta_r_z  = region( ispace(int3d, {Nx, Ny, Nz+1} ), coeffs )
  var gamma_r_z = region( ispace(int3d, {Nx, Ny, Nz+1} ), coeffs )

  var rho_avg_z = region( ispace(int3d, {Nx, Ny, Nz+1} ), double )
  var sos_avg_z = region( ispace(int3d, {Nx, Ny, Nz+1} ), double )

  var block_d_z    = region( ispace(int3d, {Nx, Ny, Nz+1} ), double[9] )
  var block_Uinv_z = region( ispace(int3d, {Nx, Ny, Nz+1} ), double[9] )

  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  --------------------------------------------------------------------------------------------
  --                       PARTITIONING
  --------------------------------------------------------------------------------------------
  
  var p_coords_x   = partition_xpencil_coords(coords,    pencil)
  var p_coords_y   = partition_ypencil_coords(coords,    pencil)
  var p_coords_z   = partition_zpencil_coords(coords,    pencil)

  var p_cnsr_x     = partition_xpencil_cnsr(r_cnsr,      pencil)
  var p_cnsr_y     = partition_ypencil_cnsr(r_cnsr,      pencil)
  var p_cnsr_z     = partition_zpencil_cnsr(r_cnsr,      pencil)

  var p_prim_c_x   = partition_xpencil_prim(r_prim_c,    pencil)
  var p_prim_c_y   = partition_ypencil_prim(r_prim_c,    pencil)
  var p_prim_c_z   = partition_zpencil_prim(r_prim_c,    pencil)

  var gp_prim_l_x  = partition_xpencil_prim(gr_prim_l_x, pencil)
  var gp_prim_l_y  = partition_ypencil_prim(gr_prim_l_y, pencil)
  var gp_prim_l_z  = partition_zpencil_prim(gr_prim_l_z, pencil)
  var gp_prim_r_x  = partition_xpencil_prim(gr_prim_r_x, pencil)
  var gp_prim_r_y  = partition_ypencil_prim(gr_prim_r_y, pencil)
  var gp_prim_r_z  = partition_zpencil_prim(gr_prim_r_z, pencil)

  var p_prim_l_x   = partition_xpencil_prim(r_prim_l_x,  pencil)
  var p_prim_l_y   = partition_ypencil_prim(r_prim_l_y,  pencil)
  var p_prim_l_z   = partition_zpencil_prim(r_prim_l_z,  pencil)

  var p_prim_r_x   = partition_xpencil_prim(r_prim_r_x,  pencil)
  var p_prim_r_y   = partition_ypencil_prim(r_prim_r_y,  pencil)
  var p_prim_r_z   = partition_zpencil_prim(r_prim_r_z,  pencil)

  var p_aux_c_x    = partition_xpencil_aux (r_aux_c,     pencil)
  var p_aux_c_y    = partition_ypencil_aux (r_aux_c,     pencil)
  var p_aux_c_z    = partition_zpencil_aux (r_aux_c,     pencil)

  var p_visc_x     = partition_xpencil_visc(r_visc,      pencil)
  var p_visc_y     = partition_ypencil_visc(r_visc,      pencil)
  var p_visc_z     = partition_zpencil_visc(r_visc,      pencil)

  var p_q_x        = partition_xpencil_vect(r_q,         pencil)
  var p_q_y        = partition_ypencil_vect(r_q,         pencil)
  var p_q_z        = partition_zpencil_vect(r_q,         pencil)

  var p_gradu_x    = partition_xpencil_tnsr2(r_gradu,    pencil)
  var p_gradu_y    = partition_ypencil_tnsr2(r_gradu,    pencil)
  var p_gradu_z    = partition_zpencil_tnsr2(r_gradu,    pencil)

  var p_tauij_x    = partition_xpencil_tnsr2symm(r_tauij, pencil)
  var p_tauij_y    = partition_ypencil_tnsr2symm(r_tauij, pencil)
  var p_tauij_z    = partition_zpencil_tnsr2symm(r_tauij, pencil)

  var p_rhs_l_x    = partition_xpencil_prim(r_rhs_l_x,  pencil)
  var p_rhs_l_y    = partition_ypencil_prim(r_rhs_l_y,  pencil)
  var p_rhs_l_z    = partition_zpencil_prim(r_rhs_l_z,  pencil)

  var p_rhs_r_x    = partition_xpencil_prim(r_rhs_r_x,  pencil)
  var p_rhs_r_y    = partition_ypencil_prim(r_rhs_r_y,  pencil)
  var p_rhs_r_z    = partition_zpencil_prim(r_rhs_r_z,  pencil)

  var p_flux_c_x   = partition_xpencil_cnsr(r_flux_c,   pencil)
  var p_flux_c_y   = partition_ypencil_cnsr(r_flux_c,   pencil)
  var p_flux_c_z   = partition_zpencil_cnsr(r_flux_c,   pencil)

  var gp_flux_l_x  = partition_xpencil_cnsr(gr_flux_l_x, pencil)
  var gp_flux_l_y  = partition_ypencil_cnsr(gr_flux_l_y, pencil)
  var gp_flux_l_z  = partition_zpencil_cnsr(gr_flux_l_z, pencil)
  var gp_flux_r_x  = partition_xpencil_cnsr(gr_flux_r_x, pencil)
  var gp_flux_r_y  = partition_ypencil_cnsr(gr_flux_r_y, pencil)
  var gp_flux_r_z  = partition_zpencil_cnsr(gr_flux_r_z, pencil)

  var p_flux_e_x   = partition_xpencil_cnsr(r_flux_e_x, pencil)
  var p_flux_e_y   = partition_ypencil_cnsr(r_flux_e_y, pencil)
  var p_flux_e_z   = partition_zpencil_cnsr(r_flux_e_z, pencil)

  var p_fder_c_x   = partition_xpencil_cnsr(r_fder_c_x, pencil)
  var p_fder_c_y   = partition_ypencil_cnsr(r_fder_c_y, pencil)
  var p_fder_c_z   = partition_zpencil_cnsr(r_fder_c_z, pencil)

  var p_rhs_x      = partition_xpencil_cnsr(r_rhs,      pencil)
  var p_rhs_y      = partition_ypencil_cnsr(r_rhs,      pencil)
  var p_rhs_z      = partition_zpencil_cnsr(r_rhs,      pencil)

  var p_qrhs_x     = partition_xpencil_cnsr(r_qrhs,     pencil)
  var p_qrhs_y     = partition_ypencil_cnsr(r_qrhs,     pencil)
  var p_qrhs_z     = partition_zpencil_cnsr(r_qrhs,     pencil)

  var p_LU_x       = partition_LU(LU_x, pencil)
  var p_LU_y       = partition_LU(LU_y, pencil)
  var p_LU_z       = partition_LU(LU_z, pencil)

  var p_LU_N_x     = partition_LU(LU_N_x, pencil)
  var p_LU_N_y     = partition_LU(LU_N_y, pencil)
  var p_LU_N_z     = partition_LU(LU_N_z, pencil)

  var p_alpha_l_x = partition_xpencil_coeffs(alpha_l_x, pencil)
  var p_beta_l_x  = partition_xpencil_coeffs(beta_l_x , pencil)
  var p_gamma_l_x = partition_xpencil_coeffs(gamma_l_x, pencil)

  var p_alpha_r_x = partition_xpencil_coeffs(alpha_r_x, pencil)
  var p_beta_r_x  = partition_xpencil_coeffs(beta_r_x , pencil)
  var p_gamma_r_x = partition_xpencil_coeffs(gamma_r_x, pencil)

  var p_rho_avg_x = partition_xpencil_double(rho_avg_x, pencil)
  var p_sos_avg_x = partition_xpencil_double(sos_avg_x, pencil)

  var p_block_d_x    = partition_xpencil_double9(block_d_x   , pencil)
  var p_block_Uinv_x = partition_xpencil_double9(block_Uinv_x, pencil)

  var p_alpha_l_y = partition_ypencil_coeffs(alpha_l_y, pencil)
  var p_beta_l_y  = partition_ypencil_coeffs(beta_l_y , pencil)
  var p_gamma_l_y = partition_ypencil_coeffs(gamma_l_y, pencil)

  var p_alpha_r_y = partition_ypencil_coeffs(alpha_r_y, pencil)
  var p_beta_r_y  = partition_ypencil_coeffs(beta_r_y , pencil)
  var p_gamma_r_y = partition_ypencil_coeffs(gamma_r_y, pencil)

  var p_rho_avg_y = partition_ypencil_double(rho_avg_y, pencil)
  var p_sos_avg_y = partition_ypencil_double(sos_avg_y, pencil)

  var p_block_d_y    = partition_ypencil_double9(block_d_y   , pencil)
  var p_block_Uinv_y = partition_ypencil_double9(block_Uinv_y, pencil)

  var p_alpha_l_z = partition_zpencil_coeffs(alpha_l_z, pencil)
  var p_beta_l_z  = partition_zpencil_coeffs(beta_l_z , pencil)
  var p_gamma_l_z = partition_zpencil_coeffs(gamma_l_z, pencil)

  var p_alpha_r_z = partition_zpencil_coeffs(alpha_r_z, pencil)
  var p_beta_r_z  = partition_zpencil_coeffs(beta_r_z , pencil)
  var p_gamma_r_z = partition_zpencil_coeffs(gamma_r_z, pencil)

  var p_rho_avg_z = partition_zpencil_double(rho_avg_z, pencil)
  var p_sos_avg_z = partition_zpencil_double(sos_avg_z, pencil)

  var p_block_d_z    = partition_zpencil_double9(block_d_z   , pencil)
  var p_block_Uinv_z = partition_zpencil_double9(block_Uinv_z, pencil)

  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  var token : int = 0

  -- Initialize characteristic interpolation matrices and SuperLU structs.
  if Nx >= 8 then
    __demand(__parallel)
    for i in pencil do
      token += set_zero_prim( p_rhs_l_x[i] )
    end

    __demand(__parallel)
    for i in pencil do
      token += set_zero_prim( p_rhs_r_x[i] )
    end
  end
  wait_for(token)

  if Ny >= 8 then
    __demand(__parallel)
    for i in pencil do
      token += set_zero_prim( p_rhs_l_y[i] )
    end

    __demand(__parallel)
    for i in pencil do
      token += set_zero_prim( p_rhs_r_y[i] )
    end
  end
  wait_for(token)

  if Nz >= 8 then
    __demand(__parallel)
    for i in pencil do
      token += set_zero_prim( p_rhs_l_z[i] )
    end

    __demand(__parallel)
    for i in pencil do
      token += set_zero_prim( p_rhs_r_z[i] )
    end
  end
  wait_for(token)
  
  -- Initialize derivatives stuff.
  __demand(__parallel)
  for i in pencil do
    token += get_LU_decomposition(p_LU_x[i], beta06MND, alpha06MND, 1.0, alpha06MND, beta06MND)
  end
  __demand(__parallel)
  for i in pencil do
    token += get_LU_decomposition(p_LU_y[i], beta06MND, alpha06MND, 1.0, alpha06MND, beta06MND)
  end
  __demand(__parallel)
  for i in pencil do
    token += get_LU_decomposition(p_LU_z[i], beta06MND, alpha06MND, 1.0, alpha06MND, beta06MND)
  end
  __demand(__parallel)
  for i in pencil do
    token += get_LU_decomposition(p_LU_N_x[i], beta06d1, alpha06d1, 1.0, alpha06d1, beta06d1)
  end
  __demand(__parallel)
  for i in pencil do
    token += get_LU_decomposition(p_LU_N_y[i], beta06d1, alpha06d1, 1.0, alpha06d1, beta06d1)
  end
  __demand(__parallel)
  for i in pencil do
    token += get_LU_decomposition(p_LU_N_z[i], beta06d1, alpha06d1, 1.0, alpha06d1, beta06d1)
  end
  wait_for(token)
  c.printf("Finished LU initialization\n")

  if config.restart then
    -- Restart the simulation from the latest viz dump
    c.printf("Restarting using viz dump\n")
    __demand(__parallel)
    for i in pencil do
      read_coords(p_coords_y[i], config.filename_prefix, i)
    end
    c.printf("Finished reading in coords\n")

    vizcount = config.restart_count
    __demand(__parallel)
    for i in pencil do
      read_primitive(p_prim_c_y[i], config.filename_prefix, vizcount, i)
    end
    c.printf("Finished restarting using viz dump %d\n", vizcount)
    tsim = tviz*vizcount
    vizcount += 1
  else
    -- If not restarting, then initialize from the problem initialization task
    __demand(__parallel)
    for i in pencil do
      -- Initialize everything in y decomposition.
      token += problem.initialize(p_coords_y[i], p_prim_c_y[i], dx, dy, dz)
    end
    c.printf("Finished initialization\n")
  end

  __demand(__parallel)
  for i in pencil do
    -- Get temperature from the initial primitive variables
    token += get_temperature_r( p_prim_c_y[i], p_aux_c_y[i] )
  end

  -- Get the velocity derivatives at initial condition 
  __demand(__parallel)
  for i in pencil do
    token += get_velocity_x_derivatives( p_prim_c_x[i], p_gradu_x[i], p_LU_N_x[i] )
  end
  __demand(__parallel)
  for i in pencil do
    token += get_velocity_y_derivatives( p_prim_c_y[i], p_gradu_y[i], p_LU_N_y[i] )
  end
  __demand(__parallel)
  for i in pencil do
    token += get_velocity_z_derivatives( p_prim_c_z[i], p_gradu_z[i], p_LU_N_z[i] )
  end
 
  var TKE0 : double = 1.0e-16
  -- do
  --   var t : double = 0.0
  --   __demand(__parallel)
  --   for i in pencil do
  --     t += problem.TKE(p_prim_c_y[i])
  --   end
  --   TKE0 = wait_for_double(t)
  -- end

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
    for i in pencil do
      IOtoken += write_coords(p_coords_y[i], config.filename_prefix, i)
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
  for i in pencil do
    token += get_conserved_r(p_prim_c_y[i], p_cnsr_y[i])
  end

  if (use_io and (not config.restart)) then
    wait_for(IOtoken)
    __demand(__parallel)
    for i in pencil do
      IOtoken += write_primitive(p_prim_c_y[i], config.filename_prefix, vizcount, i)
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
        for i in pencil do
          new_dt min= get_max_stable_dt_3d(p_prim_c_y[i], dx, dy, dz)
        end
      elseif Ny >= 8 then
        __demand(__parallel)
        for i in pencil do
          new_dt min= get_max_stable_dt_2d(p_prim_c_y[i], dx, dy)
        end
      else
        __demand(__parallel)
        for i in pencil do
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
    for i in pencil do
      set_zero_cnsr( p_qrhs_y[i] )
    end

    --------------------------------------------------------------------------------------------
    -- Advance sub-steps.
    --------------------------------------------------------------------------------------------
    for isub = 0,5 do

      -- Set RHS to zero.
      __demand(__parallel)
      for i in pencil do
        set_zero_cnsr( p_rhs_y[i] )
      end
      
      if problem.viscous then
        -- Get the transport coefficients.
        __demand(__parallel)
        for i in pencil do
          problem.get_transport_coeffs( p_prim_c_y[i], p_aux_c_y[i], p_visc_y[i] )
        end

        -- Get the viscous stress tensor.
        __demand(__parallel)
        for i in pencil do
          get_tauij( p_gradu_y[i], p_tauij_y[i], p_visc_y[i] )
        end
      end

      -- Add x-direction flux derivative to RHS.
      __demand(__parallel)
      for i in pencil do
        -- add_xflux_der_to_rhs( p_cnsr_x[i], p_prim_c_x[i], p_aux_c_x[i], p_visc_x[i], p_tauij_x[i], p_q_x[i],
        --                       p_prim_l_x[i], p_prim_r_x[i], p_rhs_l_x[i], p_rhs_r_x[i],
        --                       p_flux_c_x[i], p_flux_e_x[i], p_fder_c_x[i], p_rhs_x[i],
        --                       p_LU_x[i], p_slu_l_x[i], p_slu_r_x[i], p_matrix_l_x[i], p_matrix_r_x[i],
        --                       Nx, Ny, Nz )
        add_xflux_der_to_rhs( p_cnsr_x[i], p_prim_c_x[i], p_aux_c_x[i], p_visc_x[i], p_tauij_x[i], p_q_x[i],
                              p_prim_l_x[i], p_prim_r_x[i], p_flux_c_x[i], p_flux_e_x[i], p_fder_c_x[i], p_rhs_x[i],
                              p_alpha_l_x[i], p_beta_l_x[i], p_gamma_l_x[i],
                              p_alpha_r_x[i], p_beta_r_x[i], p_gamma_r_x[i], p_rho_avg_x[i], p_sos_avg_x[i], p_block_d_x[i],
                              p_block_Uinv_x[i], p_LU_x[i], Nx, Ny, Nz )
      end

      -- Add y-direction flux derivative to RHS.
      __demand(__parallel)
      for i in pencil do
        -- add_yflux_der_to_rhs( p_cnsr_y[i], p_prim_c_y[i], p_aux_c_y[i], p_visc_y[i], p_tauij_y[i], p_q_y[i], 
        --                       p_prim_l_y[i], p_prim_r_y[i], p_rhs_l_y[i], p_rhs_r_y[i],
        --                       p_flux_c_y[i], p_flux_e_y[i], p_fder_c_y[i], p_rhs_y[i],
        --                       p_LU_y[i], p_slu_l_y[i], p_slu_r_y[i], p_matrix_l_y[i], p_matrix_r_y[i],
        --                       Nx, Ny, Nz )
        add_yflux_der_to_rhs( p_cnsr_y[i], p_prim_c_y[i], p_aux_c_y[i], p_visc_y[i], p_tauij_y[i], p_q_y[i],
                              p_prim_l_y[i], p_prim_r_y[i], p_flux_c_y[i], p_flux_e_y[i], p_fder_c_y[i], p_rhs_y[i],
                              p_alpha_l_y[i], p_beta_l_y[i], p_gamma_l_y[i],
                              p_alpha_r_y[i], p_beta_r_y[i], p_gamma_r_y[i], p_rho_avg_y[i], p_sos_avg_y[i], p_block_d_y[i],
                              p_block_Uinv_y[i], p_LU_y[i], Nx, Ny, Nz )
      end

      -- Add z-direction flux derivative to RHS.
      __demand(__parallel)
      for i in pencil do
        -- add_zflux_der_to_rhs( p_cnsr_z[i], p_prim_c_z[i], p_aux_c_z[i], p_visc_z[i], p_tauij_z[i], p_q_z[i], 
        --                       p_prim_l_z[i], p_prim_r_z[i], p_rhs_l_z[i], p_rhs_r_z[i],
        --                       p_flux_c_z[i], p_flux_e_z[i], p_fder_c_z[i], p_rhs_z[i],
        --                       p_LU_z[i], p_slu_l_z[i], p_slu_r_z[i], p_matrix_l_z[i], p_matrix_r_z[i],
        --                       Nx, Ny, Nz )
        add_zflux_der_to_rhs( p_cnsr_z[i], p_prim_c_z[i], p_aux_c_z[i], p_visc_z[i], p_tauij_z[i], p_q_z[i],
                              p_prim_l_z[i], p_prim_r_z[i], p_flux_c_z[i], p_flux_e_z[i], p_fder_c_z[i], p_rhs_z[i],
                              p_alpha_l_z[i], p_beta_l_z[i], p_gamma_l_z[i],
                              p_alpha_r_z[i], p_beta_r_z[i], p_gamma_r_z[i], p_rho_avg_z[i], p_sos_avg_z[i], p_block_d_z[i],
                              p_block_Uinv_z[i], p_LU_z[i], Nx, Ny, Nz )
      end

       
      -- Update solution in this substep.
      __demand(__parallel)
      for i in pencil do
        update_substep( p_cnsr_y[i], p_rhs_y[i], p_qrhs_y[i], dt, A_RK45[isub], B_RK45[isub] )
      end

      -- Update simulation time as well.
      Q_t = dt + A_RK45[isub]*Q_t
      tsim += B_RK45[isub]*Q_t

      -- Update the primitive variables.
      __demand(__parallel)
      for i in pencil do
        token += get_primitive_r(p_cnsr_y[i], p_prim_c_y[i])
      end

      -- Update temperature.
      __demand(__parallel)
      for i in pencil do
        token += get_temperature_r( p_prim_c_y[i], p_aux_c_y[i] )
      end

      if problem.viscous then
        -- Update velocity gradient tensor.
        __demand(__parallel)
        for i in pencil do
          token += get_velocity_x_derivatives( p_prim_c_x[i], p_gradu_x[i], p_LU_N_x[i] )
        end
        __demand(__parallel)
        for i in pencil do
          token += get_velocity_y_derivatives( p_prim_c_y[i], p_gradu_y[i], p_LU_N_y[i] )
        end
        __demand(__parallel)
        for i in pencil do
          token += get_velocity_z_derivatives( p_prim_c_z[i], p_gradu_z[i], p_LU_N_z[i] )
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
      -- __demand(__parallel)
      -- for i in pencil do
      --   TKE += problem.TKE(p_prim_c_y[i])
      -- end

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
        for i in pencil do
          IOtoken += write_primitive(p_prim_c_y[i], config.filename_prefix, vizcount, i)
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
  for i in pencil do
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
    link_flags = {"-L" .. root_dir, "-L" .. hdf_lib_dir, "-lhdf5", "-lm", "-lblas"}
  else
    link_flags = {"-L" .. root_dir, "-lm", "-lblas"}
  end
  regentlib.saveobj(main, "wchr", "executable", link_flags)
else
  regentlib.start(main)
end
