import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")

require("fields")
require("IO")
require("EOS")
local superlu = require("superlu_util")

alpha06CI = 3.0/16.0
beta06CI  = 5.0/8.0
gamma06CI = 3.0/16.0

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

local function make_get_nonlinear_weights_JS(get_beta, is_left)
  local get_nonlinear_weights_JS __demand(__inline) task get_nonlinear_weights_JS( values : double[6][5] )
    var nlweights : double[4][5]
  
    for eq = 0, 5 do 
      var beta = [get_beta](values, eq)
  
      var epsilon : double = 1.0e-10
      
      -- Compute the nonlinear weights
      var d = array(1.0/4.0, 1.0/2.0, 1.0/4.0, 0.0)
      -- var d_central = array(1.0/8.0, 3.0/8.0, 3.0/8.0, 1.0/8.0)
      var sum : double = 0.0
      for i = 0, 4 do
        nlweights[eq][ [nl_index(i,is_left)] ] = d[i] / ( (beta[i] + epsilon)*(beta[i] + epsilon) )  -- Hard coded q = 2 here to avoid using power
        sum = sum + nlweights[eq][ [nl_index(i,is_left)] ]
      end
  
      for i = 0, 4 do
        nlweights[eq][ [nl_index(i,is_left)] ] = nlweights[eq][ [nl_index(i,is_left)] ] / sum
      end
    end
  
    return nlweights
  end
  return get_nonlinear_weights_JS
end

local function make_get_nonlinear_weights_Z(get_beta, is_left)
  local get_nonlinear_weights_Z __demand(__inline) task get_nonlinear_weights_Z( values : double[6][5] )
    var nlweights : double[4][5]

    for eq = 0, 5 do 
      var beta = [get_beta](values, eq)

      var epsilon : double = 1.0e-40

      -- Compute the nonlinear weights
      var d = array(1.0/4.0, 1.0/2.0, 1.0/4.0, 0.0)
      var tau_5 = cmath.fabs(beta[0] - beta[2]);

      var sum : double = 0.0
      for i = 0, 4 do
        nlweights[eq][ [nl_index(i,is_left)] ] = d[i]*( 1.0 + (tau_5/(beta[i] + epsilon))*(tau_5/(beta[i] + epsilon)) );
        sum = sum + nlweights[eq][ [nl_index(i,is_left)] ]
      end
      
      for i = 0, 4 do
        nlweights[eq][ [nl_index(i,is_left)] ] = nlweights[eq][ [nl_index(i,is_left)] ] / sum
      end
    end

    return nlweights
  end
  return get_nonlinear_weights_Z
end

local function make_get_nonlinear_weights_LD(get_beta, is_left)
  local get_nonlinear_weights_LD __demand(__inline) task get_nonlinear_weights_LD( values : double[6][5] )
    var nlweights : double[4][5]
  
    var d_central = array(1.0/8.0, 3.0/8.0, 3.0/8.0, 1.0/8.0)
    var d_upwind  = array(1.0/4.0, 1.0/2.0, 1.0/4.0, 0.0)
  
    for eq = 0, 5 do 
      var beta = [get_beta](values, eq)
      
      var C : double = 2.0e3
      -- p = 2
      -- q = 2
      var epsilon : double = 1.0e-6
      var alpha_beta : double = 40.0
      
      var alpha_2 : double = (values[eq][2] - values[eq][1])
      var alpha_3 : double = (values[eq][3] - values[eq][2])
      var alpha_4 : double = (values[eq][4] - values[eq][3])
  
      var theta_2 : double = cmath.fabs(alpha_2 - alpha_3) / (cmath.fabs(alpha_2) + cmath.fabs(alpha_3) + epsilon)
      var theta_3 : double = cmath.fabs(alpha_3 - alpha_4) / (cmath.fabs(alpha_3) + cmath.fabs(alpha_4) + epsilon)
  
      var sigma : double = cmath.fmax(theta_2, theta_3)
      
      var beta_avg : double = 1.0/8.0*(beta[0] + beta[2] + 6.0*beta[1])
      var tau_6 : double = cmath.fabs(beta[3] - beta_avg)
          
      -- Compute the nonlinear weights
      if tau_6/(beta_avg + epsilon) > alpha_beta then
        
        var omega_central : double[4]
        var sum : double = 0.0
        for i = 0, 4 do
          omega_central[i] = d_central[i]*( C + (tau_6/(beta[i] + epsilon))*(tau_6/(beta[i] + epsilon)) )
          sum = sum + omega_central[i]
        end
        for i = 0, 4 do
          omega_central[i] = omega_central[i] / sum
        end
  
        var tau_5 : double = cmath.fabs(beta[0] - beta[2])
        var omega_upwind : double[4]
        sum = 0.0
        for i = 0, 4 do
          omega_upwind[i] = d_upwind[i]*( 1.0 + (tau_5/(beta[i] + epsilon))*(tau_5/(beta[i] + epsilon)) )
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
          omega_central[i] = d_central[i]*( C + (tau_6/(beta[i] + epsilon))*(tau_6/(beta[i] + epsilon)) )
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

get_nonlinear_weights_JS_l = make_get_nonlinear_weights_JS(get_beta_l, true )
get_nonlinear_weights_JS_r = make_get_nonlinear_weights_JS(get_beta_r, false)

get_nonlinear_weights_Z_l  = make_get_nonlinear_weights_Z(get_beta_l, true )
get_nonlinear_weights_Z_r  = make_get_nonlinear_weights_Z(get_beta_r, false)

get_nonlinear_weights_LD_l = make_get_nonlinear_weights_LD(get_beta_l, true )
get_nonlinear_weights_LD_r = make_get_nonlinear_weights_LD(get_beta_r, false)

__demand(__inline)
task get_coefficients( nlweights : double[4][5] )

  var lcoeff0 = array(3.0/4.0, 1.0/4.0, 0.0/4.0, 1.0/4.0, 3.0/4.0, 0.0/4.0, 0.0/4.0)
  var lcoeff1 = array(1.0/4.0, 3.0/4.0, 0.0/4.0, 0.0/4.0, 3.0/4.0, 1.0/4.0, 0.0/4.0)
  var lcoeff2 = array(0.0/4.0, 3.0/4.0, 1.0/4.0, 0.0/4.0, 1.0/4.0, 3.0/4.0, 0.0/4.0)
  var lcoeff3 = array(0.0/4.0, 1.0/4.0, 3.0/4.0, 0.0/4.0, 0.0/4.0, 3.0/4.0, 1.0/4.0)

  var coeffs : double[7][5]

  for eq = 0, 5 do
    for i = 0, 7 do
      coeffs[eq][i] = lcoeff0[i]*nlweights[eq][0] + lcoeff1[i]*nlweights[eq][1] + lcoeff2[i]*nlweights[eq][2] + lcoeff3[i]*nlweights[eq][3]
    end
  end

  return coeffs
end

task WCHR_interpolation_x( r_prim_c   : region(ispace(int3d), primitive),
                           r_prim_l_x : region(ispace(int3d), primitive),
                           r_prim_r_x : region(ispace(int3d), primitive),
                           r_rhs_l    : region(ispace(int3d), primitive),
                           r_rhs_r    : region(ispace(int3d), primitive),
                           matrix_l_x : region(ispace(int2d), superlu.CSR_matrix),
                           matrix_r_x : region(ispace(int2d), superlu.CSR_matrix),
                           slu_x      : region(ispace(int2d), superlu.c.superlu_vars_t) )
where
  reads(r_prim_c), reads writes(r_prim_l_x, r_prim_r_x), reads writes(r_rhs_l, r_rhs_r, matrix_l_x, matrix_r_x, slu_x)
do

  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1
  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1
  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1

  var xdim : int64 = nx+1

  var pr = matrix_l_x.ispace.bounds.hi.x
  var pc = matrix_l_x.ispace.bounds.hi.y
  -- matrix_l_x[{pr,pc}].rowptr[0] = 0

  var bounds_c = r_prim_c.ispace.bounds
  var bounds_x = r_prim_l_x.ispace.bounds

  regentlib.assert(bounds_c.lo.x == 0, "Can only perform X interpolation in the X pencil")
  regentlib.assert(bounds_x.lo.x == 0, "Can only perform X interpolation in the X pencil")

  for i in r_prim_c do
    var rhosos_avg : double[2] = get_rho_sos_avg_x( r_prim_c, i, nx, ny, nz )
    var char_values : double[6][5] = get_char_values_x(r_prim_c, rhosos_avg[0], rhosos_avg[1], i, nx, ny, nz)

    var nlweights_l = get_nonlinear_weights_LD_l(char_values)
    var coeffs_l = get_coefficients(nlweights_l)
    var nlweights_r = get_nonlinear_weights_LD_r(char_values)
    var coeffs_r = get_coefficients(nlweights_r)

    r_rhs_l[i].rho = coeffs_l[0][3] * char_values[0][1]
                   + coeffs_l[0][4] * char_values[0][2]
                   + coeffs_l[0][5] * char_values[0][3]
                   + coeffs_l[0][6] * char_values[0][4]
    r_rhs_l[i].u   = coeffs_l[1][3] * char_values[1][1]
                   + coeffs_l[1][4] * char_values[1][2]
                   + coeffs_l[1][5] * char_values[1][3]
                   + coeffs_l[1][6] * char_values[1][4]
    r_rhs_l[i].v   = coeffs_l[2][3] * char_values[2][1]
                   + coeffs_l[2][4] * char_values[2][2]
                   + coeffs_l[2][5] * char_values[2][3]
                   + coeffs_l[2][6] * char_values[2][4]
    r_rhs_l[i].w   = coeffs_l[3][3] * char_values[3][1]
                   + coeffs_l[3][4] * char_values[3][2]
                   + coeffs_l[3][5] * char_values[3][3]
                   + coeffs_l[3][6] * char_values[3][4]
    r_rhs_l[i].p   = coeffs_l[4][3] * char_values[4][1]
                   + coeffs_l[4][4] * char_values[4][2]
                   + coeffs_l[4][5] * char_values[4][3]
                   + coeffs_l[4][6] * char_values[4][4]

    r_rhs_r[i].rho = coeffs_r[0][3] * char_values[0][1]
                   + coeffs_r[0][4] * char_values[0][2]
                   + coeffs_r[0][5] * char_values[0][3]
                   + coeffs_r[0][6] * char_values[0][4]
    r_rhs_r[i].u   = coeffs_r[1][3] * char_values[1][1]
                   + coeffs_r[1][4] * char_values[1][2]
                   + coeffs_r[1][5] * char_values[1][3]
                   + coeffs_r[1][6] * char_values[1][4]
    r_rhs_r[i].v   = coeffs_r[2][3] * char_values[2][1]
                   + coeffs_r[2][4] * char_values[2][2]
                   + coeffs_r[2][5] * char_values[2][3]
                   + coeffs_r[2][6] * char_values[2][4]
    r_rhs_r[i].w   = coeffs_r[3][3] * char_values[3][1]
                   + coeffs_r[3][4] * char_values[3][2]
                   + coeffs_r[3][5] * char_values[3][3]
                   + coeffs_r[3][6] * char_values[3][4]
    r_rhs_r[i].p   = coeffs_r[4][3] * char_values[4][1]
                   + coeffs_r[4][4] * char_values[4][2]
                   + coeffs_r[4][5] * char_values[4][3]
                   + coeffs_r[4][6] * char_values[4][4]

    var grow : int64 = 5*i.x + i.y*5*xdim + i.z*5*xdim*ny
    var bcounter : int64 = 8*3*i.x + i.y*(8*3*nx+10) + i.z*(8*3*nx+10)*ny

    -- rho
    for j = 0, 3 do
      var bcol : int64 = i.x + j - 1
      var gcol : int64 = ((bcol + nx)%nx)*5 + i.y*5*xdim + i.z*5*xdim*ny -- Top of the block
      var counter : int64 = bcounter + 2*j

      -- matrix_l_x.colind[counter] = gcol+1
      matrix_l_x[{pr,pc}].nzval [counter] = coeffs_l[0][j] * (-0.5*rhosos_avg[0]*rhosos_avg[1])
      matrix_r_x[{pr,pc}].nzval [counter] = coeffs_r[0][j] * (-0.5*rhosos_avg[0]*rhosos_avg[1])

      -- matrix_l_x.colind[counter+1] = gcol+4
      matrix_l_x[{pr,pc}].nzval [counter+1] = coeffs_l[0][j] * (0.5)
      matrix_r_x[{pr,pc}].nzval [counter+1] = coeffs_r[0][j] * (0.5)
    end
    bcounter = bcounter + 3*2
    -- matrix.rowptr[grow+1] = bcounter

    -- u
    for j = 0, 3 do
      var bcol : int64 = i.x + j - 1
      var gcol : int64 = ((bcol + nx)%nx)*5 + i.y*5*xdim + i.z*5*xdim*ny -- Top of the block
      var counter : int64 = bcounter + 2*j

      -- matrix_l_x.colind[counter] = gcol
      matrix_l_x[{pr,pc}].nzval [counter] = coeffs_l[1][j] * (1.0)
      matrix_r_x[{pr,pc}].nzval [counter] = coeffs_r[1][j] * (1.0)

      -- matrix_l_x.colind[counter+1] = gcol+4
      matrix_l_x[{pr,pc}].nzval [counter+1] = coeffs_l[1][j] * (-1.0/(rhosos_avg[1]*rhosos_avg[1]))
      matrix_r_x[{pr,pc}].nzval [counter+1] = coeffs_r[1][j] * (-1.0/(rhosos_avg[1]*rhosos_avg[1]))
    end
    bcounter = bcounter + 3*2
    -- matrix.rowptr[grow+2] = bcounter

    -- v
    for j = 0, 3 do
      var bcol : int64 = i.x + j - 1
      var gcol : int64 = ((bcol + nx)%nx)*5 + i.y*5*xdim + i.z*5*xdim*ny -- Top of the block
      var counter : int64 = bcounter + j

      -- matrix_l_x.colind[counter] = gcol+2
      matrix_l_x[{pr,pc}].nzval [counter] = coeffs_l[2][j] * (1.0)
      matrix_r_x[{pr,pc}].nzval [counter] = coeffs_r[2][j] * (1.0)
    end
    bcounter = bcounter + 3*1
    -- matrix.rowptr[grow+3] = bcounter

    -- w
    for j = 0, 3 do
      var bcol : int64 = i.x + j - 1
      var gcol : int64 = ((bcol + nx)%nx)*5 + i.y*5*xdim + i.z*5*xdim*ny -- Top of the block
      var counter : int64 = bcounter + j

      -- matrix_l_x.colind[counter] = gcol+3
      matrix_l_x[{pr,pc}].nzval [counter] = coeffs_l[3][j] * (1.0)
      matrix_r_x[{pr,pc}].nzval [counter] = coeffs_r[3][j] * (1.0)
    end
    bcounter = bcounter + 3*1
    -- matrix.rowptr[grow+4] = bcounter

    -- p
    for j = 0, 3 do
      var bcol : int64 = i.x + j - 1
      var gcol : int64 = ((bcol + nx)%nx)*5 + i.y*5*xdim + i.z*5*xdim*ny -- Top of the block
      var counter : int64 = bcounter + 2*j

      -- matrix_l_x.colind[counter] = gcol+1
      matrix_l_x[{pr,pc}].nzval [counter] = coeffs_l[4][j] * (0.5*rhosos_avg[0]*rhosos_avg[1])
      matrix_r_x[{pr,pc}].nzval [counter] = coeffs_r[4][j] * (0.5*rhosos_avg[0]*rhosos_avg[1])

      -- matrix_l_x.colind[counter+1] = gcol+4
      matrix_l_x[{pr,pc}].nzval [counter+1] = coeffs_l[4][j] * (0.5)
      matrix_r_x[{pr,pc}].nzval [counter+1] = coeffs_r[4][j] * (0.5)
    end
    bcounter = bcounter + 3*2
    -- matrix.rowptr[grow+5] = bcounter
  end

  -- var perinds : rect3d = { bounds_x.lo + {nx,0,0}, bounds_x.hi }
  -- for i in perinds do
  for iz = bounds_x.lo.z, bounds_x.hi.z+1 do
    for iy = bounds_x.lo.y, bounds_x.hi.y+1 do
      r_rhs_l[{nx,iy,iz}].{rho, u, v, w, p} = 0.0
      r_rhs_r[{nx,iy,iz}].{rho, u, v, w, p} = 0.0
    end
  end

  superlu.MatrixSolve( r_rhs_l, r_prim_l_x, matrix_l_x[{pr,pc}], nx, ny, nz, slu_x )
  superlu.MatrixSolve( r_rhs_r, r_prim_r_x, matrix_r_x[{pr,pc}], nx, ny, nz, slu_x )

  return 1
end
