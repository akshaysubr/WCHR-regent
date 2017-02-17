import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")
require("derivatives")
require("IO")
require("RHS")

local superlu = require("superlu_util")
local problem = require("problem")

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

  c.printf("================ Problem parameters ================\n")
  c.printf("           Nx, Ny, Nz = %d, %d, %d\n", Nx, Ny, Nz)
  c.printf("           Lx, Ly, Lz = %f, %f, %f\n", Lx, Ly, Lz)
  c.printf("           dx, dy, dz = %f, %f, %f\n", dx, dy, dz)
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

  -- Initialize SuperLU stuff
  matrix_l_x[{0,0}] = superlu.initialize_matrix_x(alpha06CI, beta06CI, gamma06CI, Nx, Ny, Nz)
  matrix_r_x[{0,0}] = superlu.initialize_matrix_x(alpha06CI, beta06CI, gamma06CI, Nx, Ny, Nz)
  superlu.initialize_superlu_vars( matrix_l_x[{0,0}], (Nx+1)*Ny*Nz, __physical(r_prim_r_x.rho), __fields(r_prim_r_x.rho),
                                   __physical(r_prim_l_x.rho), __fields(r_prim_l_x.rho), r_prim_l_x.bounds,
                                   __physical(slu_x)[0], __fields(slu_x)[0], slu_x.bounds)
  
  -- Initialize derivatives stuff
  get_LU_decomposition(LU_x, beta06MND, alpha06MND, 1.0, alpha06MND, beta06MND)

  var token = problem.initialize(coords, r_prim_c, dx, dy, dz)
  wait_for(token)
  
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
  wait_for(token)
  var t_start = c.legion_get_current_time_in_micros()

  __demand(__spmd)
  -- for step = 0,10 do -- Run for 10 time steps
  while tsim < tstop*(1.0 - 1.0e-16) do

    var Q_t : double = 0.0
    fill(Q_rhs.{rho, rhou, rhov, rhow, rhoE}, 0.0)
    for isub = 0,5 do
        fill(r_rhs.{rho, rhou, rhov, rhow, rhoE}, 0.0)
        add_xflux_der_to_rhs( r_cnsr, r_prim_c, r_prim_l_x, r_prim_r_x, r_rhs_l_x, r_rhs_r_x,
                              r_flux_c, r_flux_e_x, r_fder_c_x, r_rhs,
                              LU_x, slu_x, matrix_l_x, matrix_r_x )

        update_substep( r_cnsr, r_rhs, Q_rhs, dt, A_RK45[isub], B_RK45[isub] )

        -- Update simulation time as well
        Q_t = dt + A_RK45[isub]*Q_t
        tsim += B_RK45[isub]*Q_t

        token += get_primitive_r(r_cnsr, r_prim_c)
    end
    step = step + 1

    c.printf("Step: %d\n", step)
    c.printf("Simulation time: %g\n", tsim)
    c.printf("\n")
  end
  
  wait_for(token)
  var t_simulation = c.legion_get_current_time_in_micros() - t_start
  
  c.printf("Average time per time step = %12.5e\n", (t_simulation)*1e-6/step)

  var errors = problem.get_errors(coords, r_prim_c, tsim)

  c.printf("Error in rho = %g\n", errors[0])
  c.printf("Error in u   = %g\n", errors[1])
  c.printf("Error in v   = %g\n", errors[2])
  c.printf("Error in w   = %g\n", errors[3])
  c.printf("Error in p   = %g\n", errors[4])

end

regentlib.start(main)
