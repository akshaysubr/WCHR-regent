import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")

require("fields")
require("IO")
require("SOE")
require("block_tridiagonal")
local problem = require("problem")
local superlu = require("superlu_util")

local xi = 2.0/3.0
 
alpha06CI = - 45.*( xi - 1. ) / ( 16.*(xi + 5) )
beta06CI  = ( 53.*xi - 5. ) / ( 8.*(xi + 5) )
gamma06CI = - 45.*( xi - 1. ) / ( 16.*(xi + 5) )

local function v_index(i,is_left)
  if is_left then
    return i
  else
    return (6 - i - 1)
  end
end

local function nl_index(i,is_left)
  if is_left then
    return rexpr i end
  else
    return rexpr (4 - i - 1) end
  end
end

local function make_get_beta(is_left)
  local get_beta __demand(__inline) task get_beta( values : double[6][5], eq : int32 )
    var beta : double[4]
  
    beta[0] = 1.0/3.0*(values[eq][ [v_index(0,is_left)] ]*(4.0*values[eq][ [v_index(0,is_left)] ] - 19.0*values[eq][ [v_index(1,is_left)] ] + 11.0*values[eq][ [v_index(2,is_left)] ])
              + values[eq][ [v_index(1,is_left)] ]*(25.0*values[eq][ [v_index(1,is_left)] ] - 31.0*values[eq][ [v_index(2,is_left)] ]) + 10.0*values[eq][ [v_index(2,is_left)] ]*values[eq][ [v_index(2,is_left)] ])
    
    beta[1] = 1.0/3.0*(values[eq][ [v_index(1,is_left)] ]*(4.0*values[eq][ [v_index(1,is_left)] ] - 13.0*values[eq][ [v_index(2,is_left)] ] + 5.0*values[eq][ [v_index(3,is_left)] ])
              + 13.0*values[eq][ [v_index(2,is_left)] ]*(values[eq][ [v_index(2,is_left)] ] - values[eq][ [v_index(3,is_left)] ]) + 4.0*values[eq][ [v_index(3,is_left)] ]*values[eq][ [v_index(3,is_left)] ])
    
    beta[2] = 1.0/3.0*(values[eq][ [v_index(2,is_left)] ]*(10.0*values[eq][ [v_index(2,is_left)] ] - 31.0*values[eq][ [v_index(3,is_left)] ] + 11.0*values[eq][ [v_index(4,is_left)] ])
              + values[eq][ [v_index(3,is_left)] ]*(25.0*values[eq][ [v_index(3,is_left)] ] - 19.0*values[eq][ [v_index(4,is_left)] ]) + 4.0*values[eq][ [v_index(4,is_left)] ]*values[eq][ [v_index(4,is_left)] ])
    
    beta[3] = (values[eq][ [v_index(0,is_left)] ]*((525910327.0/232243200.0)*values[eq][ [v_index(0,is_left)] ] - (4562164630.0/232243200.0)*values[eq][ [v_index(1,is_left)] ]
              + (7799501420.0/232243200.0)*values[eq][ [v_index(2,is_left)] ] - (6610694540.0/232243200.0)*values[eq][ [v_index(3,is_left)] ] + (2794296070.0/232243200.0)*values[eq][ [v_index(4,is_left)] ]
              - (472758974.0/232243200.0)*values[eq][ [v_index(5,is_left)] ]) + 5.0*values[eq][ [v_index(1,is_left)] ]*((2146987907.0/232243200.0)*values[eq][ [v_index(1,is_left)] ] 
              - (7722406988.0/232243200.0)*values[eq][ [v_index(2,is_left)] ] + (6763559276.0/232243200.0)*values[eq][ [v_index(3,is_left)] ] - (2926461814.0/232243200.0)*values[eq][ [v_index(4,is_left)] ] 
              + (503766638.0/232243200.0)*values[eq][ [v_index(5,is_left)] ]) + 20.0*values[eq][ [v_index(2,is_left)] ]*((1833221603.0/232243200.0)*values[eq][ [v_index(2,is_left)] ] 
              - (3358664662.0/232243200.0)*values[eq][ [v_index(3,is_left)] ] + (1495974539.0/232243200.0)*values[eq][ [v_index(4,is_left)] ] - (263126407.0/232243200.0)*values[eq][ [v_index(5,is_left)] ]) 
              + 20.0*values[eq][ [v_index(3,is_left)] ]*((1607794163.0/232243200.0)*values[eq][ [v_index(3,is_left)] ] - (1486026707.0/232243200.0)*values[eq][ [v_index(4,is_left)] ]
              + (268747951.0/232243200.0)*values[eq][ [v_index(5,is_left)] ]) + 5.0*values[eq][ [v_index(4,is_left)] ]*((1432381427.0/232243200.0)*values[eq][ [v_index(4,is_left)] ] 
              - (536951582.0/232243200.0)*values[eq][ [v_index(5,is_left)] ]) + (263126407.0/232243200.0)*values[eq][ [v_index(5,is_left)] ]*values[eq][ [v_index(5,is_left)] ])
  
    return beta
  end
  return get_beta
end

get_beta_l = make_get_beta(true)
get_beta_r = make_get_beta(false)

local function make_get_nonlinear_weights_LD(get_beta, is_left)
  local get_nonlinear_weights_LD __demand(__inline) task get_nonlinear_weights_LD( values : double[6][5] )
    var nlweights : double[4][5]
  
    var d_central = array((8.0*xi - 5.0)/(16.0*xi + 80.0), 45.0/(16.0*xi + 80.0), 45.0/(16.0*xi + 80.0), (8.0*xi - 5.0)/(16.0*xi + 80.0)) 
    var d_upwind  = array((8.0*xi - 5.0)/(8.0*xi + 40.0), (65.0*xi - 35.0)/(16.0*xi*xi + 72.0*xi - 40.0), (25.0*xi - 10.0)/(16.0*xi*xi + 72.0*xi - 40.0), 0.0)
 
    for eq = 0, 5 do 
      var beta = [get_beta](values, eq)
      
      var C : double = 1.0e10
      -- p = 2
      -- q = 4
      var epsilon   : double = 1.0e-40
      var alpha_tau : double = 55.0
      
      var alpha_2 : double = (values[eq][2] - values[eq][1])
      var alpha_3 : double = (values[eq][3] - values[eq][2])
      var alpha_4 : double = (values[eq][4] - values[eq][3])
  
      var theta_2 : double = cmath.fabs(alpha_2 - alpha_3) / (cmath.fabs(alpha_2) + cmath.fabs(alpha_3) + epsilon)
      var theta_3 : double = cmath.fabs(alpha_3 - alpha_4) / (cmath.fabs(alpha_3) + cmath.fabs(alpha_4) + epsilon)
  
      var sigma : double = cmath.fmax(theta_2, theta_3)
      
      var beta_avg : double = 1.0/8.0*(beta[0] + beta[2] + 6.0*beta[1])
      var tau_6 : double = cmath.fabs(beta[3] - beta_avg)
          
      -- Compute the nonlinear weights
      if tau_6/(beta_avg + epsilon) > alpha_tau then
        
        var omega_central : double[4]
        var sum : double = 0.0
        for i = 0, 4 do
          var dummy : double = (tau_6/(beta[i] + epsilon))
          dummy = dummy*dummy
          dummy = dummy*dummy       -- q = 4
          omega_central[i] = d_central[i]*( C + dummy )
          sum = sum + omega_central[i]
        end
        for i = 0, 4 do
          omega_central[i] = omega_central[i] / sum
        end
  
        var tau_5 : double = cmath.fabs(beta[0] - beta[2])
        var omega_upwind : double[4]
        sum = 0.0
        for i = 0, 4 do
          var dummy : double = (tau_5/(beta[i] + epsilon))
          dummy = dummy*dummy      -- p = 2
          omega_upwind[i] = d_upwind[i]*( 1.0 + dummy )
          sum = sum + omega_upwind[i]
        end
        for i = 0, 4 do
          omega_upwind[i] = omega_upwind[i] / sum
        end
  
        for i = 0, 4 do
          nlweights[eq][ [nl_index(i,is_left)] ] = (sigma)*omega_upwind[i] + (1.0 - sigma)*omega_central[i]
        end
      else
        var omega_central : double[4]
        var sum : double = 0.0
        for i = 0, 4 do
          var dummy : double = (tau_6/(beta[i] + epsilon))
          dummy = dummy*dummy
          dummy = dummy*dummy       -- q = 4
          omega_central[i] = d_central[i]*( C + dummy )
          sum = sum + omega_central[i]
        end
        for i = 0, 4 do
          nlweights[eq][ [nl_index(i,is_left)] ] = omega_central[i] / sum
        end
      end
    end
      
    return nlweights
  end
  return get_nonlinear_weights_LD
end

get_nonlinear_weights_LD_l = make_get_nonlinear_weights_LD(get_beta_l, true )
get_nonlinear_weights_LD_r = make_get_nonlinear_weights_LD(get_beta_r, false)


__demand(__inline)
task get_coefficients_ECI( nlweights : double[4][5] )

  var lcoeff0 = array(0.0,       1.0, 0.0,           3.0/8.0, -5.0/4.0,          15.0/8.0,         0.0,              0.0,               0.0)
  var lcoeff1 = array(-xi + 1.0, xi,  0.0,           0.0,     -xi/2.0 + 3.0/8.0, 3.0/4.0,          xi/2.0 - 1.0/8.0, 0.0,               0.0)
  var lcoeff2 = array(0.0,       xi,  -xi + 1.0,     0.0,     0.0,               xi/2.0 - 1.0/8.0, 3.0/4.0,          -xi/2.0 + 3.0/8.0, 0.0)
  var lcoeff3 = array(0.0,       1.0, 0.0,           0.0,     0.0,               0.0,              15.0/8.0,         -5.0/4.0,          3.0/8.0)

  var coeffs : double[9][5]

  for eq = 0, 5 do
    for i = 0, 9 do
      coeffs[eq][i] = lcoeff0[i]*nlweights[eq][0] + lcoeff1[i]*nlweights[eq][1] + lcoeff2[i]*nlweights[eq][2] + lcoeff3[i]*nlweights[eq][3]
    end
  end

  return coeffs
end

terra wait_for(x : int)
  return x
end

solve_tridiagonal_x_v = make_solve_tridiagonal_x('_2', 'v')
solve_tridiagonal_x_w = make_solve_tridiagonal_x('_3', 'w')

__demand(__inline)
task WCHR_interpolation_x( r_prim_c   : region(ispace(int3d), primitive),
                           r_prim_l   : region(ispace(int3d), primitive),
                           r_prim_r   : region(ispace(int3d), primitive),
                           alpha_l    : region(ispace(int3d), coeffs),
                           beta_l     : region(ispace(int3d), coeffs),
                           gamma_l    : region(ispace(int3d), coeffs),
                           alpha_r    : region(ispace(int3d), coeffs),
                           beta_r     : region(ispace(int3d), coeffs),
                           gamma_r    : region(ispace(int3d), coeffs),
                           rho_avg    : region(ispace(int3d), double),
                           sos_avg    : region(ispace(int3d), double),
                           block_d    : region(ispace(int3d), double[9]),
                           block_Uinv : region(ispace(int3d), double[9]),
                           Nx         : int64,
                           Ny         : int64,
                           Nz         : int64 )
where
  reads(r_prim_c), reads writes(r_prim_l, r_prim_r, alpha_l, beta_l, gamma_l, alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv)
do

  var t_start = c.legion_get_current_time_in_micros()
  var token : int = 0

  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1
  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1
  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1

  var bounds_c = r_prim_c.ispace.bounds
  var bounds_x = r_prim_l.ispace.bounds

  regentlib.assert(bounds_c.lo.x == 0, "Can only perform X interpolation in the X pencil")
  regentlib.assert(bounds_x.lo.x == 0, "Can only perform X interpolation in the X pencil")

  -- var alpha_l = region( ispace(int3d, {nx+1, ny, nz}, bounds_x.lo), coeffs )
  -- var beta_l  = region( ispace(int3d, {nx+1, ny, nz}, bounds_x.lo), coeffs )
  -- var gamma_l = region( ispace(int3d, {nx+1, ny, nz}, bounds_x.lo), coeffs )

  -- var alpha_r = region( ispace(int3d, {nx+1, ny, nz}, bounds_x.lo), coeffs )
  -- var beta_r  = region( ispace(int3d, {nx+1, ny, nz}, bounds_x.lo), coeffs )
  -- var gamma_r = region( ispace(int3d, {nx+1, ny, nz}, bounds_x.lo), coeffs )

  -- var rho_avg = region( ispace(int3d, {nx+1, ny, nz}, bounds_x.lo), double )
  -- var sos_avg = region( ispace(int3d, {nx+1, ny, nz}, bounds_x.lo), double )

  -- var block_d    = region( ispace(int3d, {nx+1, ny, nz}, bounds_x.lo), double[9] )
  -- var block_Uinv = region( ispace(int3d, {nx+1, ny, nz}, bounds_x.lo), double[9] )

  wait_for(token)
  var t_alloc = c.legion_get_current_time_in_micros()

  for i in r_prim_c do
    var rhosos_avg : double[2] = get_rho_sos_avg_x( r_prim_c, i, Nx, Ny, Nz )
    var char_values : double[6][5] = get_char_values_x(r_prim_c, rhosos_avg[0], rhosos_avg[1], i, Nx, Ny, Nz)

    var nlweights_l = get_nonlinear_weights_LD_l(char_values)
    var coeffs_l = get_coefficients_ECI(nlweights_l)
    var nlweights_r = get_nonlinear_weights_LD_r(char_values)
    var coeffs_r = get_coefficients_ECI(nlweights_r)

    alpha_l[i]._0 = coeffs_l[0][0]; beta_l[i]._0 = coeffs_l[0][1]; gamma_l[i]._0 = coeffs_l[0][2];
    alpha_l[i]._1 = coeffs_l[1][0]; beta_l[i]._1 = coeffs_l[1][1]; gamma_l[i]._1 = coeffs_l[1][2];
    alpha_l[i]._2 = coeffs_l[2][0]; beta_l[i]._2 = coeffs_l[2][1]; gamma_l[i]._2 = coeffs_l[2][2];
    alpha_l[i]._3 = coeffs_l[3][0]; beta_l[i]._3 = coeffs_l[3][1]; gamma_l[i]._3 = coeffs_l[3][2];
    alpha_l[i]._4 = coeffs_l[4][0]; beta_l[i]._4 = coeffs_l[4][1]; gamma_l[i]._4 = coeffs_l[4][2];

    alpha_r[i]._0 = coeffs_r[0][0]; beta_r[i]._0 = coeffs_r[0][1]; gamma_r[i]._0 = coeffs_r[0][2];
    alpha_r[i]._1 = coeffs_r[1][0]; beta_r[i]._1 = coeffs_r[1][1]; gamma_r[i]._1 = coeffs_r[1][2];
    alpha_r[i]._2 = coeffs_r[2][0]; beta_r[i]._2 = coeffs_r[2][1]; gamma_r[i]._2 = coeffs_r[2][2];
    alpha_r[i]._3 = coeffs_r[3][0]; beta_r[i]._3 = coeffs_r[3][1]; gamma_r[i]._3 = coeffs_r[3][2];
    alpha_r[i]._4 = coeffs_r[4][0]; beta_r[i]._4 = coeffs_r[4][1]; gamma_r[i]._4 = coeffs_r[4][2];

    rho_avg[i] = rhosos_avg[0]
    sos_avg[i] = rhosos_avg[1]

    -- RHS for left sided interpolation
    r_prim_l[i].rho = coeffs_l[0][3] * char_values[0][0]
                    + coeffs_l[0][4] * char_values[0][1]
                    + coeffs_l[0][5] * char_values[0][2]
                    + coeffs_l[0][6] * char_values[0][3]
                    + coeffs_l[0][7] * char_values[0][4]
                    + coeffs_l[0][8] * char_values[0][5]

    r_prim_l[i].u   = coeffs_l[1][3] * char_values[1][0]
                    + coeffs_l[1][4] * char_values[1][1]
                    + coeffs_l[1][5] * char_values[1][2]
                    + coeffs_l[1][6] * char_values[1][3]
                    + coeffs_l[1][7] * char_values[1][4]
                    + coeffs_l[1][8] * char_values[1][5]

    r_prim_l[i].v   = coeffs_l[2][3] * char_values[2][0]
                    + coeffs_l[2][4] * char_values[2][1]
                    + coeffs_l[2][5] * char_values[2][2]
                    + coeffs_l[2][6] * char_values[2][3]
                    + coeffs_l[2][7] * char_values[2][4]
                    + coeffs_l[2][8] * char_values[2][5]

    r_prim_l[i].w   = coeffs_l[3][3] * char_values[3][0]
                    + coeffs_l[3][4] * char_values[3][1]
                    + coeffs_l[3][5] * char_values[3][2]
                    + coeffs_l[3][6] * char_values[3][3]
                    + coeffs_l[3][7] * char_values[3][4]
                    + coeffs_l[3][8] * char_values[3][5]

    r_prim_l[i].p   = coeffs_l[4][3] * char_values[4][0]
                    + coeffs_l[4][4] * char_values[4][1]
                    + coeffs_l[4][5] * char_values[4][2]
                    + coeffs_l[4][6] * char_values[4][3]
                    + coeffs_l[4][7] * char_values[4][4]
                    + coeffs_l[4][8] * char_values[4][5]

    -- RHS for right sided interpolation
    r_prim_r[i].rho = coeffs_r[0][3] * char_values[0][0]
                    + coeffs_r[0][4] * char_values[0][1]
                    + coeffs_r[0][5] * char_values[0][2]
                    + coeffs_r[0][6] * char_values[0][3]
                    + coeffs_r[0][7] * char_values[0][4]
                    + coeffs_r[0][8] * char_values[0][5]

    r_prim_r[i].u   = coeffs_r[1][3] * char_values[1][0]
                    + coeffs_r[1][4] * char_values[1][1]
                    + coeffs_r[1][5] * char_values[1][2]
                    + coeffs_r[1][6] * char_values[1][3]
                    + coeffs_r[1][7] * char_values[1][4]
                    + coeffs_r[1][8] * char_values[1][5]

    r_prim_r[i].v   = coeffs_r[2][3] * char_values[2][0]
                    + coeffs_r[2][4] * char_values[2][1]
                    + coeffs_r[2][5] * char_values[2][2]
                    + coeffs_r[2][6] * char_values[2][3]
                    + coeffs_r[2][7] * char_values[2][4]
                    + coeffs_r[2][8] * char_values[2][5]

    r_prim_r[i].w   = coeffs_r[3][3] * char_values[3][0]
                    + coeffs_r[3][4] * char_values[3][1]
                    + coeffs_r[3][5] * char_values[3][2]
                    + coeffs_r[3][6] * char_values[3][3]
                    + coeffs_r[3][7] * char_values[3][4]
                    + coeffs_r[3][8] * char_values[3][5]

    r_prim_r[i].p   = coeffs_r[4][3] * char_values[4][0]
                    + coeffs_r[4][4] * char_values[4][1]
                    + coeffs_r[4][5] * char_values[4][2]
                    + coeffs_r[4][6] * char_values[4][3]
                    + coeffs_r[4][7] * char_values[4][4]
                    + coeffs_r[4][8] * char_values[4][5]

  end

  var t_weights = c.legion_get_current_time_in_micros()

  solve_block_tridiagonal_x( alpha_l, beta_l, gamma_l, rho_avg, sos_avg, r_prim_l, block_d, block_Uinv )
  solve_block_tridiagonal_x( alpha_r, beta_r, gamma_r, rho_avg, sos_avg, r_prim_r, block_d, block_Uinv )

  solve_tridiagonal_x_v( alpha_l, beta_l, gamma_l, r_prim_l, rho_avg )
  solve_tridiagonal_x_v( alpha_r, beta_r, gamma_r, r_prim_r, sos_avg )

  solve_tridiagonal_x_w( alpha_l, beta_l, gamma_l, r_prim_l, rho_avg )
  solve_tridiagonal_x_w( alpha_r, beta_r, gamma_r, r_prim_r, sos_avg )

  var t_end = c.legion_get_current_time_in_micros()

  -- c.printf("X: Time to get coefficients and RHS: %12.5e\n", (t_weights-t_alloc)*1e-6)
  -- c.printf("X: Time for block tridiagonal solves: %12.5e\n", (t_end-t_weights)*1e-6)
  -- c.printf("X: Time to get the WCHR interpolation: %12.5e\n", (t_end-t_start)*1e-6)
 return 1
end


solve_tridiagonal_y_u = make_solve_tridiagonal_y('_1', 'u')
solve_tridiagonal_y_w = make_solve_tridiagonal_y('_3', 'w')

__demand(__inline)
task WCHR_interpolation_y( r_prim_c   : region(ispace(int3d), primitive),
                           r_prim_l   : region(ispace(int3d), primitive),
                           r_prim_r   : region(ispace(int3d), primitive),
                           alpha_l    : region(ispace(int3d), coeffs),
                           beta_l     : region(ispace(int3d), coeffs),
                           gamma_l    : region(ispace(int3d), coeffs),
                           alpha_r    : region(ispace(int3d), coeffs),
                           beta_r     : region(ispace(int3d), coeffs),
                           gamma_r    : region(ispace(int3d), coeffs),
                           rho_avg    : region(ispace(int3d), double),
                           sos_avg    : region(ispace(int3d), double),
                           block_d    : region(ispace(int3d), double[9]),
                           block_Uinv : region(ispace(int3d), double[9]),
                           Nx         : int64,
                           Ny         : int64,
                           Nz         : int64 )
where
  reads(r_prim_c), reads writes(r_prim_l, r_prim_r, alpha_l, beta_l, gamma_l, alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv)
do

  var t_start = c.legion_get_current_time_in_micros()

  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1
  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1
  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1

  var bounds_c = r_prim_c.ispace.bounds
  var bounds_y = r_prim_l.ispace.bounds

  regentlib.assert(bounds_c.lo.y == 0, "Can only perform Y interpolation in the Y pencil")
  regentlib.assert(bounds_y.lo.y == 0, "Can only perform Y interpolation in the Y pencil")

  var t_alloc = c.legion_get_current_time_in_micros()

  for i in r_prim_c do
    var rhosos_avg : double[2] = get_rho_sos_avg_y( r_prim_c, i, Nx, Ny, Nz )
    var char_values : double[6][5] = get_char_values_y(r_prim_c, rhosos_avg[0], rhosos_avg[1], i, Nx, Ny, Nz)

    var nlweights_l = get_nonlinear_weights_LD_l(char_values)
    var coeffs_l = get_coefficients_ECI(nlweights_l)
    var nlweights_r = get_nonlinear_weights_LD_r(char_values)
    var coeffs_r = get_coefficients_ECI(nlweights_r)

    alpha_l[i]._0 = coeffs_l[0][0]; beta_l[i]._0 = coeffs_l[0][1]; gamma_l[i]._0 = coeffs_l[0][2];
    alpha_l[i]._1 = coeffs_l[1][0]; beta_l[i]._1 = coeffs_l[1][1]; gamma_l[i]._1 = coeffs_l[1][2];
    alpha_l[i]._2 = coeffs_l[2][0]; beta_l[i]._2 = coeffs_l[2][1]; gamma_l[i]._2 = coeffs_l[2][2];
    alpha_l[i]._3 = coeffs_l[3][0]; beta_l[i]._3 = coeffs_l[3][1]; gamma_l[i]._3 = coeffs_l[3][2];
    alpha_l[i]._4 = coeffs_l[4][0]; beta_l[i]._4 = coeffs_l[4][1]; gamma_l[i]._4 = coeffs_l[4][2];

    alpha_r[i]._0 = coeffs_r[0][0]; beta_r[i]._0 = coeffs_r[0][1]; gamma_r[i]._0 = coeffs_r[0][2];
    alpha_r[i]._1 = coeffs_r[1][0]; beta_r[i]._1 = coeffs_r[1][1]; gamma_r[i]._1 = coeffs_r[1][2];
    alpha_r[i]._2 = coeffs_r[2][0]; beta_r[i]._2 = coeffs_r[2][1]; gamma_r[i]._2 = coeffs_r[2][2];
    alpha_r[i]._3 = coeffs_r[3][0]; beta_r[i]._3 = coeffs_r[3][1]; gamma_r[i]._3 = coeffs_r[3][2];
    alpha_r[i]._4 = coeffs_r[4][0]; beta_r[i]._4 = coeffs_r[4][1]; gamma_r[i]._4 = coeffs_r[4][2];

    rho_avg[i] = rhosos_avg[0]
    sos_avg[i] = rhosos_avg[1]

    -- RHS for left sided interpolation
    r_prim_l[i].rho = coeffs_l[0][3] * char_values[0][0]
                    + coeffs_l[0][4] * char_values[0][1]
                    + coeffs_l[0][5] * char_values[0][2]
                    + coeffs_l[0][6] * char_values[0][3]
                    + coeffs_l[0][7] * char_values[0][4]
                    + coeffs_l[0][8] * char_values[0][5]

    r_prim_l[i].u   = coeffs_l[1][3] * char_values[1][0]
                    + coeffs_l[1][4] * char_values[1][1]
                    + coeffs_l[1][5] * char_values[1][2]
                    + coeffs_l[1][6] * char_values[1][3]
                    + coeffs_l[1][7] * char_values[1][4]
                    + coeffs_l[1][8] * char_values[1][5]

    r_prim_l[i].v   = coeffs_l[2][3] * char_values[2][0]
                    + coeffs_l[2][4] * char_values[2][1]
                    + coeffs_l[2][5] * char_values[2][2]
                    + coeffs_l[2][6] * char_values[2][3]
                    + coeffs_l[2][7] * char_values[2][4]
                    + coeffs_l[2][8] * char_values[2][5]

    r_prim_l[i].w   = coeffs_l[3][3] * char_values[3][0]
                    + coeffs_l[3][4] * char_values[3][1]
                    + coeffs_l[3][5] * char_values[3][2]
                    + coeffs_l[3][6] * char_values[3][3]
                    + coeffs_l[3][7] * char_values[3][4]
                    + coeffs_l[3][8] * char_values[3][5]

    r_prim_l[i].p   = coeffs_l[4][3] * char_values[4][0]
                    + coeffs_l[4][4] * char_values[4][1]
                    + coeffs_l[4][5] * char_values[4][2]
                    + coeffs_l[4][6] * char_values[4][3]
                    + coeffs_l[4][7] * char_values[4][4]
                    + coeffs_l[4][8] * char_values[4][5]

    -- RHS for right sided interpolation
    r_prim_r[i].rho = coeffs_r[0][3] * char_values[0][0]
                    + coeffs_r[0][4] * char_values[0][1]
                    + coeffs_r[0][5] * char_values[0][2]
                    + coeffs_r[0][6] * char_values[0][3]
                    + coeffs_r[0][7] * char_values[0][4]
                    + coeffs_r[0][8] * char_values[0][5]

    r_prim_r[i].u   = coeffs_r[1][3] * char_values[1][0]
                    + coeffs_r[1][4] * char_values[1][1]
                    + coeffs_r[1][5] * char_values[1][2]
                    + coeffs_r[1][6] * char_values[1][3]
                    + coeffs_r[1][7] * char_values[1][4]
                    + coeffs_r[1][8] * char_values[1][5]

    r_prim_r[i].v   = coeffs_r[2][3] * char_values[2][0]
                    + coeffs_r[2][4] * char_values[2][1]
                    + coeffs_r[2][5] * char_values[2][2]
                    + coeffs_r[2][6] * char_values[2][3]
                    + coeffs_r[2][7] * char_values[2][4]
                    + coeffs_r[2][8] * char_values[2][5]

    r_prim_r[i].w   = coeffs_r[3][3] * char_values[3][0]
                    + coeffs_r[3][4] * char_values[3][1]
                    + coeffs_r[3][5] * char_values[3][2]
                    + coeffs_r[3][6] * char_values[3][3]
                    + coeffs_r[3][7] * char_values[3][4]
                    + coeffs_r[3][8] * char_values[3][5]

    r_prim_r[i].p   = coeffs_r[4][3] * char_values[4][0]
                    + coeffs_r[4][4] * char_values[4][1]
                    + coeffs_r[4][5] * char_values[4][2]
                    + coeffs_r[4][6] * char_values[4][3]
                    + coeffs_r[4][7] * char_values[4][4]
                    + coeffs_r[4][8] * char_values[4][5]

  end
      
  var t_weights = c.legion_get_current_time_in_micros()

  solve_block_tridiagonal_y( alpha_l, beta_l, gamma_l, rho_avg, sos_avg, r_prim_l, block_d, block_Uinv )
  solve_block_tridiagonal_y( alpha_r, beta_r, gamma_r, rho_avg, sos_avg, r_prim_r, block_d, block_Uinv )

  solve_tridiagonal_y_u( alpha_l, beta_l, gamma_l, r_prim_l, rho_avg )
  solve_tridiagonal_y_u( alpha_r, beta_r, gamma_r, r_prim_r, sos_avg )

  solve_tridiagonal_y_w( alpha_l, beta_l, gamma_l, r_prim_l, rho_avg )
  solve_tridiagonal_y_w( alpha_r, beta_r, gamma_r, r_prim_r, sos_avg )

  var t_end = c.legion_get_current_time_in_micros()

  -- c.printf("Y: Time to get coefficients and RHS: %12.5e\n", (t_weights-t_alloc)*1e-6)
  -- c.printf("Y: Time for block tridiagonal solves: %12.5e\n", (t_end-t_weights)*1e-6)
  -- c.printf("Y: Time to get the WCHR interpolation: %12.5e\n", (t_end-t_start)*1e-6)
 return 1
end

solve_tridiagonal_z_u = make_solve_tridiagonal_z('_1', 'u')
solve_tridiagonal_z_v = make_solve_tridiagonal_z('_2', 'v')

__demand(__inline)
task WCHR_interpolation_z( r_prim_c   : region(ispace(int3d), primitive),
                           r_prim_l   : region(ispace(int3d), primitive),
                           r_prim_r   : region(ispace(int3d), primitive),
                           alpha_l    : region(ispace(int3d), coeffs),
                           beta_l     : region(ispace(int3d), coeffs),
                           gamma_l    : region(ispace(int3d), coeffs),
                           alpha_r    : region(ispace(int3d), coeffs),
                           beta_r     : region(ispace(int3d), coeffs),
                           gamma_r    : region(ispace(int3d), coeffs),
                           rho_avg    : region(ispace(int3d), double),
                           sos_avg    : region(ispace(int3d), double),
                           block_d    : region(ispace(int3d), double[9]),
                           block_Uinv : region(ispace(int3d), double[9]),
                           Nx         : int64,
                           Ny         : int64,
                           Nz         : int64 )
where
  reads(r_prim_c), reads writes(r_prim_l, r_prim_r, alpha_l, beta_l, gamma_l, alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv)
do

  var t_start = c.legion_get_current_time_in_micros()

  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1
  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1
  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1

  var bounds_c = r_prim_c.ispace.bounds
  var bounds_z = r_prim_l.ispace.bounds

  regentlib.assert(bounds_c.lo.z == 0, "Can only perform Z interpolation in the Z pencil")
  regentlib.assert(bounds_z.lo.z == 0, "Can only perform Z interpolation in the Z pencil")

  var t_alloc = c.legion_get_current_time_in_micros()

  for i in r_prim_c do
    var rhosos_avg : double[2] = get_rho_sos_avg_z( r_prim_c, i, Nx, Ny, Nz )
    var char_values : double[6][5] = get_char_values_z(r_prim_c, rhosos_avg[0], rhosos_avg[1], i, Nx, Ny, Nz)

    var nlweights_l = get_nonlinear_weights_LD_l(char_values)
    var coeffs_l = get_coefficients_ECI(nlweights_l)
    var nlweights_r = get_nonlinear_weights_LD_r(char_values)
    var coeffs_r = get_coefficients_ECI(nlweights_r)

    alpha_l[i]._0 = coeffs_l[0][0]; beta_l[i]._0 = coeffs_l[0][1]; gamma_l[i]._0 = coeffs_l[0][2];
    alpha_l[i]._1 = coeffs_l[1][0]; beta_l[i]._1 = coeffs_l[1][1]; gamma_l[i]._1 = coeffs_l[1][2];
    alpha_l[i]._2 = coeffs_l[2][0]; beta_l[i]._2 = coeffs_l[2][1]; gamma_l[i]._2 = coeffs_l[2][2];
    alpha_l[i]._3 = coeffs_l[3][0]; beta_l[i]._3 = coeffs_l[3][1]; gamma_l[i]._3 = coeffs_l[3][2];
    alpha_l[i]._4 = coeffs_l[4][0]; beta_l[i]._4 = coeffs_l[4][1]; gamma_l[i]._4 = coeffs_l[4][2];

    alpha_r[i]._0 = coeffs_r[0][0]; beta_r[i]._0 = coeffs_r[0][1]; gamma_r[i]._0 = coeffs_r[0][2];
    alpha_r[i]._1 = coeffs_r[1][0]; beta_r[i]._1 = coeffs_r[1][1]; gamma_r[i]._1 = coeffs_r[1][2];
    alpha_r[i]._2 = coeffs_r[2][0]; beta_r[i]._2 = coeffs_r[2][1]; gamma_r[i]._2 = coeffs_r[2][2];
    alpha_r[i]._3 = coeffs_r[3][0]; beta_r[i]._3 = coeffs_r[3][1]; gamma_r[i]._3 = coeffs_r[3][2];
    alpha_r[i]._4 = coeffs_r[4][0]; beta_r[i]._4 = coeffs_r[4][1]; gamma_r[i]._4 = coeffs_r[4][2];

    rho_avg[i] = rhosos_avg[0]
    sos_avg[i] = rhosos_avg[1]

    -- RHS for left sided interpolation
    r_prim_l[i].rho = coeffs_l[0][3] * char_values[0][0]
                    + coeffs_l[0][4] * char_values[0][1]
                    + coeffs_l[0][5] * char_values[0][2]
                    + coeffs_l[0][6] * char_values[0][3]
                    + coeffs_l[0][7] * char_values[0][4]
                    + coeffs_l[0][8] * char_values[0][5]

    r_prim_l[i].u   = coeffs_l[1][3] * char_values[1][0]
                    + coeffs_l[1][4] * char_values[1][1]
                    + coeffs_l[1][5] * char_values[1][2]
                    + coeffs_l[1][6] * char_values[1][3]
                    + coeffs_l[1][7] * char_values[1][4]
                    + coeffs_l[1][8] * char_values[1][5]

    r_prim_l[i].v   = coeffs_l[2][3] * char_values[2][0]
                    + coeffs_l[2][4] * char_values[2][1]
                    + coeffs_l[2][5] * char_values[2][2]
                    + coeffs_l[2][6] * char_values[2][3]
                    + coeffs_l[2][7] * char_values[2][4]
                    + coeffs_l[2][8] * char_values[2][5]

    r_prim_l[i].w   = coeffs_l[3][3] * char_values[3][0]
                    + coeffs_l[3][4] * char_values[3][1]
                    + coeffs_l[3][5] * char_values[3][2]
                    + coeffs_l[3][6] * char_values[3][3]
                    + coeffs_l[3][7] * char_values[3][4]
                    + coeffs_l[3][8] * char_values[3][5]

    r_prim_l[i].p   = coeffs_l[4][3] * char_values[4][0]
                    + coeffs_l[4][4] * char_values[4][1]
                    + coeffs_l[4][5] * char_values[4][2]
                    + coeffs_l[4][6] * char_values[4][3]
                    + coeffs_l[4][7] * char_values[4][4]
                    + coeffs_l[4][8] * char_values[4][5]

    -- RHS for right sided interpolation
    r_prim_r[i].rho = coeffs_r[0][3] * char_values[0][0]
                    + coeffs_r[0][4] * char_values[0][1]
                    + coeffs_r[0][5] * char_values[0][2]
                    + coeffs_r[0][6] * char_values[0][3]
                    + coeffs_r[0][7] * char_values[0][4]
                    + coeffs_r[0][8] * char_values[0][5]

    r_prim_r[i].u   = coeffs_r[1][3] * char_values[1][0]
                    + coeffs_r[1][4] * char_values[1][1]
                    + coeffs_r[1][5] * char_values[1][2]
                    + coeffs_r[1][6] * char_values[1][3]
                    + coeffs_r[1][7] * char_values[1][4]
                    + coeffs_r[1][8] * char_values[1][5]

    r_prim_r[i].v   = coeffs_r[2][3] * char_values[2][0]
                    + coeffs_r[2][4] * char_values[2][1]
                    + coeffs_r[2][5] * char_values[2][2]
                    + coeffs_r[2][6] * char_values[2][3]
                    + coeffs_r[2][7] * char_values[2][4]
                    + coeffs_r[2][8] * char_values[2][5]

    r_prim_r[i].w   = coeffs_r[3][3] * char_values[3][0]
                    + coeffs_r[3][4] * char_values[3][1]
                    + coeffs_r[3][5] * char_values[3][2]
                    + coeffs_r[3][6] * char_values[3][3]
                    + coeffs_r[3][7] * char_values[3][4]
                    + coeffs_r[3][8] * char_values[3][5]

    r_prim_r[i].p   = coeffs_r[4][3] * char_values[4][0]
                    + coeffs_r[4][4] * char_values[4][1]
                    + coeffs_r[4][5] * char_values[4][2]
                    + coeffs_r[4][6] * char_values[4][3]
                    + coeffs_r[4][7] * char_values[4][4]
                    + coeffs_r[4][8] * char_values[4][5]

  end

  var t_weights = c.legion_get_current_time_in_micros()

  solve_block_tridiagonal_z( alpha_l, beta_l, gamma_l, rho_avg, sos_avg, r_prim_l, block_d, block_Uinv )
  solve_block_tridiagonal_z( alpha_r, beta_r, gamma_r, rho_avg, sos_avg, r_prim_r, block_d, block_Uinv )

  solve_tridiagonal_z_u( alpha_l, beta_l, gamma_l, r_prim_l, rho_avg )
  solve_tridiagonal_z_u( alpha_r, beta_r, gamma_r, r_prim_r, sos_avg )

  solve_tridiagonal_z_v( alpha_l, beta_l, gamma_l, r_prim_l, rho_avg )
  solve_tridiagonal_z_v( alpha_r, beta_r, gamma_r, r_prim_r, sos_avg )

  var t_end = c.legion_get_current_time_in_micros()


  -- c.printf("Z: Time to get coefficients and RHS: %12.5e\n", (t_weights-t_alloc)*1e-6)
  -- c.printf("Z: Time for block tridiagonal solves: %12.5e\n", (t_end-t_weights)*1e-6)
  -- c.printf("Z: Time to get the WCHR interpolation: %12.5e\n", (t_end-t_start)*1e-6)
  return 1
end
