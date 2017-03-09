import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")
require("derivatives")
require("IO")
require("RHS")
require("partition")

local superlu = require("superlu_util")
local problem = require("problem")
local Config  = require("config")

local csuperlu_mapper = require("superlu_mapper")

terra wait_for(x : int)
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

  var config : Config
  config:initialize_from_command(Nx, Ny, Nz)

  c.printf("================ Problem parameters ================\n")
  c.printf("           Nx, Ny, Nz = %d, %d, %d\n", Nx, Ny, Nz)
  c.printf("           Lx, Ly, Lz = %f, %f, %f\n", Lx, Ly, Lz)
  c.printf("           dx, dy, dz = %f, %f, %f\n", dx, dy, dz)
  c.printf("           prow, pcol = %d, %d\n", config.prow, config.pcol)

  if config.fileIO then
    c.printf("           fileIO     = true\n")
    c.printf("           prefix     = %s\n", config.filename_prefix)
  else
    c.printf("           fileIO     = false\n")
  end

  c.printf("====================================================\n")

  --------------------------------------------------------------------------------------------
  --                       DATA STUCTURES
  --------------------------------------------------------------------------------------------
  var grid_c     = ispace(int3d, {x = Nx,   y = Ny,   z = Nz  })  -- Cell center index space

  var grid_e_x   = ispace(int3d, {x = Nx+1, y = Ny,   z = Nz  })  -- x cell edge index space
  var grid_e_y   = ispace(int3d, {x = Nx,   y = Ny+1, z = Nz  })  -- y cell edge index space
  var grid_e_z   = ispace(int3d, {x = Nx,   y = Ny,   z = Nz+1})  -- z cell edge index space

  var coords     = region(grid_c, coordinates)  -- Coordinates of cell center

  var r_cnsr     = region(grid_c,   conserved)  -- Conserved variables at cell center

  var r_prim_c   = region(grid_c,   primitive)  -- Primitive variables at cell center
  var r_prim_l_x = region(grid_e_x, primitive)  -- Primitive variables at left x cell edge
  var r_prim_l_y = region(grid_e_y, primitive)  -- Primitive variables at left y cell edge
  var r_prim_l_z = region(grid_e_z, primitive)  -- Primitive variables at left z cell edge
  var r_prim_r_x = region(grid_e_x, primitive)  -- Primitive variables at right x cell edge
  var r_prim_r_y = region(grid_e_y, primitive)  -- Primitive variables at right y cell edge
  var r_prim_r_z = region(grid_e_z, primitive)  -- Primitive variables at right z cell edge

  var r_rhs_l_x  = region(grid_e_x, primitive)  -- Store RHS for left interpolation in x
  var r_rhs_r_x  = region(grid_e_x, primitive)  -- Store RHS for right interpolation in x
  var r_rhs_l_y  = region(grid_e_y, primitive)  -- Store RHS for left interpolation in y
  var r_rhs_r_y  = region(grid_e_y, primitive)  -- Store RHS for right interpolation in y
  var r_rhs_l_z  = region(grid_e_z, primitive)  -- Store RHS for left interpolation in z
  var r_rhs_r_z  = region(grid_e_z, primitive)  -- Store RHS for right interpolation in z

  var r_flux_c   = region(grid_c,   conserved)  -- Flux at cell center
  var r_flux_e_x = region(grid_e_x, conserved)  -- Flux at x cell edge
  var r_flux_e_y = region(grid_e_y, conserved)  -- Flux at y cell edge
  var r_flux_e_z = region(grid_e_z, conserved)  -- Flux at z cell edge
  
  var r_fder_c_x = region(grid_c,   conserved)     -- x flux derivative
  var r_fder_c_y = region(grid_c,   conserved)     -- y flux derivative
  var r_fder_c_z = region(grid_c,   conserved)     -- z flux derivative
  
  var r_rhs      = region(grid_c,   conserved)  -- RHS for time stepping at cell center
  var Q_rhs      = region(grid_c,   conserved)  -- Buffer for RK45 time stepping

  var LU_x       = region(ispace(int3d, {x = Nx, y = 1, z = 1}), LU_struct) -- Data structure to hold x derivative LU decomposition
  var LU_y       = region(ispace(int3d, {x = Ny, y = 1, z = 1}), LU_struct) -- Data structure to hold y derivative LU decomposition
  var LU_z       = region(ispace(int3d, {x = Nz, y = 1, z = 1}), LU_struct) -- Data structure to hold z derivative LU decomposition

  var pgrid_x    = ispace(int2d, {x = 1, y = 1}) -- Processor grid in x
  var pgrid_y    = ispace(int2d, {x = 1, y = 1}) -- Processor grid in y
  var pgrid_z    = ispace(int2d, {x = 1, y = 1}) -- Processor grid in z

  var slu_x      = region(pgrid_x, superlu.c.superlu_vars_t) -- Super LU data structure for x interpolation
  var slu_y      = region(pgrid_y, superlu.c.superlu_vars_t) -- Super LU data structure for y interpolation
  var slu_z      = region(pgrid_z, superlu.c.superlu_vars_t) -- Super LU data structure for z interpolation

  var matrix_l_x = region(pgrid_x, superlu.CSR_matrix) -- matrix data structure for x left interpolation
  var matrix_r_x = region(pgrid_x, superlu.CSR_matrix) -- matrix data structure for x right interpolation
  var matrix_l_y = region(pgrid_y, superlu.CSR_matrix) -- matrix data structure for y left interpolation
  var matrix_r_y = region(pgrid_y, superlu.CSR_matrix) -- matrix data structure for y right interpolation
  var matrix_l_z = region(pgrid_z, superlu.CSR_matrix) -- matrix data structure for z left interpolation
  var matrix_r_z = region(pgrid_z, superlu.CSR_matrix) -- matrix data structure for z right interpolation
  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  -- Initialize characteristic interpolation matrices and SuperLU structs
  if Nx >= 5 then
    superlu.initialize_matrix_char_x(matrix_l_x, alpha06CI, beta06CI, gamma06CI, Nx, Ny, Nz)
    superlu.initialize_matrix_char_x(matrix_r_x, alpha06CI, beta06CI, gamma06CI, Nx, Ny, Nz)

    fill( r_rhs_l_x.{rho,u,v,w,p}, 0.0 )
    superlu.init_superlu_vars( matrix_l_x, 5*(Nx+1)*Ny*Nz, r_rhs_l_x, r_prim_l_x, slu_x )
  end
  if Ny >= 5 then
    superlu.initialize_matrix_char_y(matrix_l_y, alpha06CI, beta06CI, gamma06CI, Nx, Ny, Nz)
    superlu.initialize_matrix_char_y(matrix_r_y, alpha06CI, beta06CI, gamma06CI, Nx, Ny, Nz)

    fill( r_rhs_l_y.{rho,u,v,w,p}, 0.0 )
    superlu.init_superlu_vars( matrix_l_y, 5*Nx*(Ny+1)*Nz, r_rhs_l_y, r_prim_l_y, slu_y )
  end
  if Nz >= 5 then
    superlu.initialize_matrix_char_z(matrix_l_z, alpha06CI, beta06CI, gamma06CI, Nx, Ny, Nz)
    superlu.initialize_matrix_char_z(matrix_r_z, alpha06CI, beta06CI, gamma06CI, Nx, Ny, Nz)

    fill( r_rhs_l_z.{rho,u,v,w,p}, 0.0 )
    superlu.init_superlu_vars( matrix_l_z, 5*Nx*Ny*(Nz+1), r_rhs_l_z, r_prim_l_z, slu_z )
  end
  
  -- Initialize derivatives stuff
  get_LU_decomposition(LU_x, beta06MND, alpha06MND, 1.0, alpha06MND, beta06MND)
  get_LU_decomposition(LU_y, beta06MND, alpha06MND, 1.0, alpha06MND, beta06MND)
  get_LU_decomposition(LU_z, beta06MND, alpha06MND, 1.0, alpha06MND, beta06MND)

  var token = problem.initialize(coords, r_prim_c, dx, dy, dz)
  wait_for(token)
  
  var IOtoken = 0
  if config.fileIO then
    IOtoken += write_coords(coords)
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

  var tstop : double = problem.tstop
  var dt    : double = problem.dt
  var tsim  : double = 0.0
  var step  : int64  = 0

  token += get_conserved_r(r_prim_c, r_cnsr) -- Get conserved variables after initialization

  wait_for(IOtoken)
  if config.fileIO then
    IOtoken += write_primitive(r_prim_c, "cell_primitive", step)
  end
  wait_for(IOtoken)
  
  wait_for(token)
  var t_start = c.legion_get_current_time_in_micros()

  __demand(__spmd)
  -- for step = 0,10 do -- Run for 10 time steps
  while tsim < tstop*(1.0 - 1.0e-16) do

    var Q_t : double = 0.0
    fill(Q_rhs.{rho, rhou, rhov, rhow, rhoE}, 0.0)
    for isub = 0,5 do
        -- Set RHS to zero
        set_rhs_zero( r_rhs )
        
        -- Add X direction flux derivatives to RHS
        add_xflux_der_to_rhs( r_cnsr, r_prim_c, r_prim_l_x, r_prim_r_x, r_rhs_l_x, r_rhs_r_x,
                              r_flux_c, r_flux_e_x, r_fder_c_x, r_rhs,
                              LU_x, slu_x, matrix_l_x, matrix_r_x )

        -- Add Y direction flux derivatives to RHS
        add_yflux_der_to_rhs( r_cnsr, r_prim_c, r_prim_l_y, r_prim_r_y, r_rhs_l_y, r_rhs_r_y,
                              r_flux_c, r_flux_e_y, r_fder_c_y, r_rhs,
                              LU_y, slu_y, matrix_l_y, matrix_r_y )

        -- Add Z direction flux derivatives to RHS
        add_zflux_der_to_rhs( r_cnsr, r_prim_c, r_prim_l_z, r_prim_r_z, r_rhs_l_z, r_rhs_r_z,
                              r_flux_c, r_flux_e_z, r_fder_c_z, r_rhs,
                              LU_z, slu_z, matrix_l_z, matrix_r_z )

        -- Update solution in this substep
        update_substep( r_cnsr, r_rhs, Q_rhs, dt, A_RK45[isub], B_RK45[isub] )

        -- Update simulation time as well
        Q_t = dt + A_RK45[isub]*Q_t
        tsim += B_RK45[isub]*Q_t

        token += get_primitive_r(r_cnsr, r_prim_c)
        c.printf("    Min rho r_cnsr   (main.rg) : %g\n", min_rho_c(r_cnsr) )
        c.printf("    Min rho r_prim_c (main.rg) : %g\n", min_rho_p(r_prim_c) )
    end
    step = step + 1

    c.printf("\n")

    c.printf("Simulation time: %g\n", tsim)
    c.printf("    Step: %d\n", step)
    c.printf("    Timestep: %g\n", dt)
  
    var errors = problem.get_errors(coords, r_prim_c, tsim)

    c.printf("    Error in rho = %g\n", errors[0])
    c.printf("    Error in u   = %g\n", errors[1])
    c.printf("    Error in v   = %g\n", errors[2])
    c.printf("    Error in w   = %g\n", errors[3])
    c.printf("    Error in p   = %g\n", errors[4])

    wait_for([int](errors[0]))
  end
  
  wait_for(token)
  var t_simulation = c.legion_get_current_time_in_micros() - t_start

  c.printf("Average time per time step = %12.5e\n", (t_simulation)*1e-6/step)
  
  if config.fileIO then
    write_primitive(r_prim_c, "cell_primitive", step)
  end

end

regentlib.start(main, csuperlu_mapper.register_mappers)
