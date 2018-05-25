import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")

require("fields")
require("IO")
require("SOE")
require("block_tridiagonal")
local problem = require("problem")

local periodic_x = problem.periodic_x
local periodic_y = problem.periodic_y
local periodic_z = problem.periodic_z

local xi  = 2.0/3.0
local xi0 = 9.0/152.0
local xi1 = -14445.0/171608.0
local xi2 = -3182085.0/37433632.0 - (45.0*cmath.sqrt(723535913.0))/37433632.0
local xi3 = -(9.0*cmath.sqrt(723535913.0))/7659176.0 + 96676.0/957397.0

local C         = 1.0e10
local epsilon   = 1.0e-40
local alpha_tau = 55.0

do
  if scheme ~= "WCHR" then
    xi  = 1.0
    xi0 = -1.0/8.0
    xi1 = 0.0
    xi2 = 15.0/8.0
    xi3 = -1.0/8.0
  end
  if scheme == "WCNS-LD" then
    C = 1.0e9
    alpha_tau = 35.0
  end
end
 
alpha06CI = - 45.*( xi - 1. ) / ( 16.*(xi + 5) )
beta06CI  = ( 53.*xi - 5. ) / ( 8.*(xi + 5) )
gamma06CI = - 45.*( xi - 1. ) / ( 16.*(xi + 5) )

local ip = {}
ip.n_ghosts = 3

local function v_index(i,is_left)
  if is_left then
    return i
  else
    return (6 - i - 1)
  end
end

local function v_index_B(i,is_left_boundary)
  if is_left_boundary then
    return i
  else
    return (6 - i)
  end
end

local function v_index_LB(i,is_left_biased)
  if is_left_biased then
    return i
  else
    return (6 - i - 1)
  end
end

local function v_index_RB(i,is_left_biased)
  if is_left_biased then
    return (i + 1)
  else
    return (6 - i)
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

local function make_get_beta_LB(is_left)
  local get_beta_LB __demand(__inline) task get_beta_LB( values : double[7][5], eq : int32 )
    var beta : double[4]
  
    beta[0] = 1.0/3.0*(values[eq][ [v_index_LB(0,is_left)] ]*(4.0*values[eq][ [v_index_LB(0,is_left)] ] - 19.0*values[eq][ [v_index_LB(1,is_left)] ] + 11.0*values[eq][ [v_index_LB(2,is_left)] ])
              + values[eq][ [v_index_LB(1,is_left)] ]*(25.0*values[eq][ [v_index_LB(1,is_left)] ] - 31.0*values[eq][ [v_index_LB(2,is_left)] ]) + 10.0*values[eq][ [v_index_LB(2,is_left)] ]*values[eq][ [v_index_LB(2,is_left)] ])
    
    beta[1] = 1.0/3.0*(values[eq][ [v_index_LB(1,is_left)] ]*(4.0*values[eq][ [v_index_LB(1,is_left)] ] - 13.0*values[eq][ [v_index_LB(2,is_left)] ] + 5.0*values[eq][ [v_index_LB(3,is_left)] ])
              + 13.0*values[eq][ [v_index_LB(2,is_left)] ]*(values[eq][ [v_index_LB(2,is_left)] ] - values[eq][ [v_index_LB(3,is_left)] ]) + 4.0*values[eq][ [v_index_LB(3,is_left)] ]*values[eq][ [v_index_LB(3,is_left)] ])
    
    beta[2] = 1.0/3.0*(values[eq][ [v_index_LB(2,is_left)] ]*(10.0*values[eq][ [v_index_LB(2,is_left)] ] - 31.0*values[eq][ [v_index_LB(3,is_left)] ] + 11.0*values[eq][ [v_index_LB(4,is_left)] ])
              + values[eq][ [v_index_LB(3,is_left)] ]*(25.0*values[eq][ [v_index_LB(3,is_left)] ] - 19.0*values[eq][ [v_index_LB(4,is_left)] ]) + 4.0*values[eq][ [v_index_LB(4,is_left)] ]*values[eq][ [v_index_LB(4,is_left)] ])
    
    beta[3] = (values[eq][ [v_index_LB(0,is_left)] ]*((525910327.0/232243200.0)*values[eq][ [v_index_LB(0,is_left)] ] - (4562164630.0/232243200.0)*values[eq][ [v_index_LB(1,is_left)] ]
              + (7799501420.0/232243200.0)*values[eq][ [v_index_LB(2,is_left)] ] - (6610694540.0/232243200.0)*values[eq][ [v_index_LB(3,is_left)] ] + (2794296070.0/232243200.0)*values[eq][ [v_index_LB(4,is_left)] ]
              - (472758974.0/232243200.0)*values[eq][ [v_index_LB(5,is_left)] ]) + 5.0*values[eq][ [v_index_LB(1,is_left)] ]*((2146987907.0/232243200.0)*values[eq][ [v_index_LB(1,is_left)] ] 
              - (7722406988.0/232243200.0)*values[eq][ [v_index_LB(2,is_left)] ] + (6763559276.0/232243200.0)*values[eq][ [v_index_LB(3,is_left)] ] - (2926461814.0/232243200.0)*values[eq][ [v_index_LB(4,is_left)] ] 
              + (503766638.0/232243200.0)*values[eq][ [v_index_LB(5,is_left)] ]) + 20.0*values[eq][ [v_index_LB(2,is_left)] ]*((1833221603.0/232243200.0)*values[eq][ [v_index_LB(2,is_left)] ] 
              - (3358664662.0/232243200.0)*values[eq][ [v_index_LB(3,is_left)] ] + (1495974539.0/232243200.0)*values[eq][ [v_index_LB(4,is_left)] ] - (263126407.0/232243200.0)*values[eq][ [v_index_LB(5,is_left)] ]) 
              + 20.0*values[eq][ [v_index_LB(3,is_left)] ]*((1607794163.0/232243200.0)*values[eq][ [v_index_LB(3,is_left)] ] - (1486026707.0/232243200.0)*values[eq][ [v_index_LB(4,is_left)] ]
              + (268747951.0/232243200.0)*values[eq][ [v_index_LB(5,is_left)] ]) + 5.0*values[eq][ [v_index_LB(4,is_left)] ]*((1432381427.0/232243200.0)*values[eq][ [v_index_LB(4,is_left)] ] 
              - (536951582.0/232243200.0)*values[eq][ [v_index_LB(5,is_left)] ]) + (263126407.0/232243200.0)*values[eq][ [v_index_LB(5,is_left)] ]*values[eq][ [v_index_LB(5,is_left)] ])
  
    return beta
  end
  return get_beta_LB
end

local function make_get_beta_RB(is_left)
  local get_beta_RB __demand(__inline) task get_beta_RB( values : double[7][5], eq : int32 )
    var beta : double[4]
  
    beta[0] = 1.0/3.0*(values[eq][ [v_index_RB(0,is_left)] ]*(4.0*values[eq][ [v_index_RB(0,is_left)] ] - 19.0*values[eq][ [v_index_RB(1,is_left)] ] + 11.0*values[eq][ [v_index_RB(2,is_left)] ])
              + values[eq][ [v_index_RB(1,is_left)] ]*(25.0*values[eq][ [v_index_RB(1,is_left)] ] - 31.0*values[eq][ [v_index_RB(2,is_left)] ]) + 10.0*values[eq][ [v_index_RB(2,is_left)] ]*values[eq][ [v_index_RB(2,is_left)] ])
    
    beta[1] = 1.0/3.0*(values[eq][ [v_index_RB(1,is_left)] ]*(4.0*values[eq][ [v_index_RB(1,is_left)] ] - 13.0*values[eq][ [v_index_RB(2,is_left)] ] + 5.0*values[eq][ [v_index_RB(3,is_left)] ])
              + 13.0*values[eq][ [v_index_RB(2,is_left)] ]*(values[eq][ [v_index_RB(2,is_left)] ] - values[eq][ [v_index_RB(3,is_left)] ]) + 4.0*values[eq][ [v_index_RB(3,is_left)] ]*values[eq][ [v_index_RB(3,is_left)] ])
    
    beta[2] = 1.0/3.0*(values[eq][ [v_index_RB(2,is_left)] ]*(10.0*values[eq][ [v_index_RB(2,is_left)] ] - 31.0*values[eq][ [v_index_RB(3,is_left)] ] + 11.0*values[eq][ [v_index_RB(4,is_left)] ])
              + values[eq][ [v_index_RB(3,is_left)] ]*(25.0*values[eq][ [v_index_RB(3,is_left)] ] - 19.0*values[eq][ [v_index_RB(4,is_left)] ]) + 4.0*values[eq][ [v_index_RB(4,is_left)] ]*values[eq][ [v_index_RB(4,is_left)] ])
    
    beta[3] = (values[eq][ [v_index_RB(0,is_left)] ]*((525910327.0/232243200.0)*values[eq][ [v_index_RB(0,is_left)] ] - (4562164630.0/232243200.0)*values[eq][ [v_index_RB(1,is_left)] ]
              + (7799501420.0/232243200.0)*values[eq][ [v_index_RB(2,is_left)] ] - (6610694540.0/232243200.0)*values[eq][ [v_index_RB(3,is_left)] ] + (2794296070.0/232243200.0)*values[eq][ [v_index_RB(4,is_left)] ]
              - (472758974.0/232243200.0)*values[eq][ [v_index_RB(5,is_left)] ]) + 5.0*values[eq][ [v_index_RB(1,is_left)] ]*((2146987907.0/232243200.0)*values[eq][ [v_index_RB(1,is_left)] ] 
              - (7722406988.0/232243200.0)*values[eq][ [v_index_RB(2,is_left)] ] + (6763559276.0/232243200.0)*values[eq][ [v_index_RB(3,is_left)] ] - (2926461814.0/232243200.0)*values[eq][ [v_index_RB(4,is_left)] ] 
              + (503766638.0/232243200.0)*values[eq][ [v_index_RB(5,is_left)] ]) + 20.0*values[eq][ [v_index_RB(2,is_left)] ]*((1833221603.0/232243200.0)*values[eq][ [v_index_RB(2,is_left)] ] 
              - (3358664662.0/232243200.0)*values[eq][ [v_index_RB(3,is_left)] ] + (1495974539.0/232243200.0)*values[eq][ [v_index_RB(4,is_left)] ] - (263126407.0/232243200.0)*values[eq][ [v_index_RB(5,is_left)] ]) 
              + 20.0*values[eq][ [v_index_RB(3,is_left)] ]*((1607794163.0/232243200.0)*values[eq][ [v_index_RB(3,is_left)] ] - (1486026707.0/232243200.0)*values[eq][ [v_index_RB(4,is_left)] ]
              + (268747951.0/232243200.0)*values[eq][ [v_index_RB(5,is_left)] ]) + 5.0*values[eq][ [v_index_RB(4,is_left)] ]*((1432381427.0/232243200.0)*values[eq][ [v_index_RB(4,is_left)] ] 
              - (536951582.0/232243200.0)*values[eq][ [v_index_RB(5,is_left)] ]) + (263126407.0/232243200.0)*values[eq][ [v_index_RB(5,is_left)] ]*values[eq][ [v_index_RB(5,is_left)] ])
  
    return beta
  end
  return get_beta_RB
end

get_beta_l = make_get_beta(true)
get_beta_r = make_get_beta(false)

get_beta_LB_l = make_get_beta_LB(true)
get_beta_LB_r = make_get_beta_LB(false)

get_beta_RB_l = make_get_beta_RB(true)
get_beta_RB_r = make_get_beta_RB(false)

local function make_get_nonlinear_weights_LD(get_beta, is_left, scheme)

  local order = 6
  if scheme == "WCNS-JS" or scheme == "WCNS-Z" then
    order = 5
  end
  
  local C_upwind = 1.0
  if scheme == "WCNS-JS" then
    C_upwind = 0.0
  end

  local get_nonlinear_weights_LD __demand(__inline) task get_nonlinear_weights_LD( values : double[6][5] )
    var nlweights : double[4][5]
  
    var d_central = array((8.0*xi - 5.0)/(16.0*xi + 80.0), 45.0/(16.0*xi + 80.0), 45.0/(16.0*xi + 80.0), (8.0*xi - 5.0)/(16.0*xi + 80.0)) 
    var d_upwind  = array((8.0*xi - 5.0)/(8.0*xi + 40.0), (65.0*xi - 35.0)/(16.0*xi*xi + 72.0*xi - 40.0), (25.0*xi - 10.0)/(16.0*xi*xi + 72.0*xi - 40.0), 0.0)
 
    for eq = 0, 5 do 
      var beta = [get_beta](values, eq)
      
      -- p = 2
      -- q = 4
      
      var alpha_2 : double = (values[eq][2] - values[eq][1])
      var alpha_3 : double = (values[eq][3] - values[eq][2])
      var alpha_4 : double = (values[eq][4] - values[eq][3])
  
      var theta_2 : double = cmath.fabs(alpha_2 - alpha_3) / (cmath.fabs(alpha_2) + cmath.fabs(alpha_3) + epsilon)
      var theta_3 : double = cmath.fabs(alpha_3 - alpha_4) / (cmath.fabs(alpha_3) + cmath.fabs(alpha_4) + epsilon)
  
      var sigma : double = cmath.fmax(theta_2, theta_3)
      
      var beta_avg : double = 1.0/8.0*(beta[0] + beta[2] + 6.0*beta[1])
      var tau_6 : double = cmath.fabs(beta[3] - beta_avg)

      if order < 6 then
        sigma = 1.0
      end
          
      -- Compute the nonlinear weights
      if (tau_6/(beta_avg + epsilon) > alpha_tau) or (order < 6) then
        
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
  
        var tau_5 : double = C_upwind*cmath.fabs(beta[0] - beta[2]) + (1.0 - C_upwind)
        var omega_upwind : double[4]
        sum = 0.0
        for i = 0, 4 do
          var dummy : double = (tau_5/(beta[i] + epsilon))
          dummy = dummy*dummy      -- p = 2
          omega_upwind[i] = d_upwind[i]*( C_upwind + dummy )
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

local function make_get_nonlinear_weights_LD_LBLB(get_beta, is_left_biased, is_left_boundary)

  local order = 6
  if scheme == "WCNS-JS" or scheme == "WCNS-Z" then
    order = 5
  end
  
  local C_upwind = 1.0
  if scheme == "WCNS-JS" then
    C_upwind = 0.0
  end

  local get_nonlinear_weights_LD_LBLB __demand(__inline) task get_nonlinear_weights_LD_LBLB( values : double[7][5] )
    var nlweights : double[4][5]

    var d_5 = array( (56.0*xi0 - 5.0)/(24.0*(24.0*xi0 - 5.0)),
                     (5.0*(104.0*xi0 - 11.0))/(24.0*(24.0*xi0 - 5.0)), 
                     -5.0/(48.0*xi0 - 10.0), 
                     0.0 )
    var d_6 = array( (6560.0*xi0*xi1 + 552.0*xi0 - 716.0*xi1 - 75.0)/(24.0*(3648.0*xi0*xi1 + 376.0*xi0 - 1080.0*xi1 - 145.0)),
                     (161984.0*xi0*xi1 + 15480.0*xi0 - 20456.0*xi1 - 2385.0)/(48.0*(3648.0*xi0*xi1 + 376.0*xi0 - 1080.0*xi1 - 145.0)),
                     -(624.0*xi1 + 90.0)/(3648.0*xi0*xi1 + 376.0*xi0 - 1080.0*xi1 - 145.0),
                     (488.0*xi0 - 35.0)/(16.0*(3648.0*xi0*xi1 + 376.0*xi0 - 1080.0*xi1 - 145.0)) )

    for eq = 0, 5 do 
      var beta = [get_beta](values, eq)
      
      -- p = 2
      -- q = 4
  
      var alpha_2 : double = (values[eq][ [v_index_B(2,is_left_boundary)] ] - values[eq][ [v_index_B(1,is_left_boundary)] ])
      var alpha_3 : double = (values[eq][ [v_index_B(3,is_left_boundary)] ] - values[eq][ [v_index_B(2,is_left_boundary)] ])
      var alpha_4 : double = (values[eq][ [v_index_B(4,is_left_boundary)] ] - values[eq][ [v_index_B(3,is_left_boundary)] ])

      var theta_2 : double = cmath.fabs(alpha_2 - alpha_3) / (cmath.fabs(alpha_2) + cmath.fabs(alpha_3) + epsilon)
      var theta_3 : double = cmath.fabs(alpha_3 - alpha_4) / (cmath.fabs(alpha_3) + cmath.fabs(alpha_4) + epsilon)
  
      var sigma : double = cmath.fmax(theta_2, theta_3)
      
      var beta_avg : double = 1.0/8.0*(beta[0] + beta[2] + 6.0*beta[1])
      var tau_6 : double = cmath.fabs(beta[3] - beta_avg)
          
      if order < 6 then
        sigma = 1.0
      end
          
      -- Compute the nonlinear weights
      if (tau_6/(beta_avg + epsilon) > alpha_tau) or (order < 6) then
        
        var omega_6 : double[4]
        var sum : double = 0.0
        for i = 0, 4 do
          var dummy : double = (tau_6/(beta[i] + epsilon))
          dummy = dummy*dummy
          dummy = dummy*dummy       -- q = 4
          omega_6[i] = d_6[i]*( C + dummy )
          sum = sum + omega_6[i]
        end
        for i = 0, 4 do
          omega_6[i] = omega_6[i] / sum
        end
  
        var tau_5 : double = C_upwind*cmath.fabs(beta[0] - beta[2]) + (1.0 - C_upwind)
        var omega_5 : double[4]
        sum = 0.0
        for i = 0, 4 do
          var dummy : double = (tau_5/(beta[i] + epsilon))
          dummy = dummy*dummy      -- p = 2
          omega_5[i] = d_5[i]*( C_upwind + dummy )
          sum = sum + omega_5[i]
        end
        for i = 0, 4 do
          omega_5[i] = omega_5[i] / sum
        end
  
        for i = 0, 4 do
          nlweights[eq][ [nl_index(i,is_left_biased)] ] = (sigma)*omega_5[i] + (1.0 - sigma)*omega_6[i]
        end
      else
        var omega_6 : double[4]
        var sum : double = 0.0
        for i = 0, 4 do
          var dummy : double = (tau_6/(beta[i] + epsilon))
          dummy = dummy*dummy
          dummy = dummy*dummy       -- q = 4
          omega_6[i] = d_6[i]*( C + dummy )
          sum = sum + omega_6[i]
        end
        for i = 0, 4 do
          nlweights[eq][ [nl_index(i,is_left_biased)] ] = omega_6[i] / sum
        end
      end
    end
      
    return nlweights
  end
  return get_nonlinear_weights_LD_LBLB
end

local function make_get_nonlinear_weights_LD_LBRB(get_beta, is_left_biased, is_left_boundary)

  local order = 6
  if scheme == "WCNS-JS" or scheme == "WCNS-Z" then
    order = 5
  end
  
  local C_upwind = 1.0
  if scheme == "WCNS-JS" then
    C_upwind = 0.0
  end

  local get_nonlinear_weights_LD_LBRB __demand(__inline) task get_nonlinear_weights_LD_LBRB( values : double[7][5] )
    var nlweights : double[4][5]

    var d_5 = array( (56.0*xi3 - 5.0)/(32.0*(48.0*xi2*xi3 - 26.0*xi2 + 8.0*xi3 - 5.0)),
                     -(76.0*xi2 + 15.0)/(4.0*(48.0*xi2*xi3 - 26.0*xi2 + 8.0*xi3 - 5.0)),
                     (1536.0*xi2*xi3 - 224.0*xi2 + 200.0*xi3 - 35.0)/(32.0*(48.0*xi2*xi3 - 26.0*xi2 + 8.0*xi3 - 5.0)),
                     0.0 )
    var d_6 = array( (488.0*xi3 - 35.0)/(16.0*(3648.0*xi2*xi3 + 376.0*xi3 - 1080.0*xi2 - 145.0)),
                     -(624.0*xi2 + 90.0)/(3648.0*xi2*xi3 + 376.0*xi3 - 1080.0*xi2 - 145.0),
                     -(161984.0*xi2*xi3 - 20456.0*xi2 + 15480.0*xi3 - 2385.0)/(48.0*(3648.0*xi2*xi3 + 376.0*xi3 - 1080.0*xi2 - 145.0)),
                     (6560.0*xi2*xi3 - 716.0*xi2 + 552.0*xi3 - 75.0)/(24.0*(3648.0*xi2*xi3 + 376.0*xi3 - 1080.0*xi2 - 145.0)) )

    for eq = 0, 5 do 
      var beta = [get_beta](values, eq)
      
      -- p = 2
      -- q = 4
      
      var alpha_2 : double = (values[eq][ [v_index_B(2,is_left_boundary)] ] - values[eq][ [v_index_B(1,is_left_boundary)] ])
      var alpha_3 : double = (values[eq][ [v_index_B(3,is_left_boundary)] ] - values[eq][ [v_index_B(2,is_left_boundary)] ])
      var alpha_4 : double = (values[eq][ [v_index_B(4,is_left_boundary)] ] - values[eq][ [v_index_B(3,is_left_boundary)] ])
  
      var theta_2 : double = cmath.fabs(alpha_2 - alpha_3) / (cmath.fabs(alpha_2) + cmath.fabs(alpha_3) + epsilon)
      var theta_3 : double = cmath.fabs(alpha_3 - alpha_4) / (cmath.fabs(alpha_3) + cmath.fabs(alpha_4) + epsilon)
  
      var sigma : double = cmath.fmax(theta_2, theta_3)
      
      var beta_avg : double = 1.0/8.0*(beta[0] + beta[2] + 6.0*beta[1])
      var tau_6 : double = cmath.fabs(beta[3] - beta_avg)
          
      if order < 6 then
        sigma = 1.0
      end
          
      -- Compute the nonlinear weights
      if (tau_6/(beta_avg + epsilon) > alpha_tau) or (order < 6) then
        
        var omega_6 : double[4]
        var sum : double = 0.0
        for i = 0, 4 do
          var dummy : double = (tau_6/(beta[i] + epsilon))
          dummy = dummy*dummy
          dummy = dummy*dummy       -- q = 4
          omega_6[i] = d_6[i]*( C + dummy )
          sum = sum + omega_6[i]
        end
        for i = 0, 4 do
          omega_6[i] = omega_6[i] / sum
        end
  
        var tau_5 : double = C_upwind*cmath.fabs(beta[0] - beta[2]) + (1.0 - C_upwind)
        var omega_5 : double[4]
        sum = 0.0
        for i = 0, 4 do
          var dummy : double = (tau_5/(beta[i] + epsilon))
          dummy = dummy*dummy      -- p = 2
          omega_5[i] = d_5[i]*( C_upwind + dummy )
          sum = sum + omega_5[i]
        end
        for i = 0, 4 do
          omega_5[i] = omega_5[i] / sum
        end
  
        for i = 0, 4 do
          nlweights[eq][ [nl_index(i,is_left_biased)] ] = (sigma)*omega_5[i] + (1.0 - sigma)*omega_6[i]
        end
      else
        var omega_6 : double[4]
        var sum : double = 0.0
        for i = 0, 4 do
          var dummy : double = (tau_6/(beta[i] + epsilon))
          dummy = dummy*dummy
          dummy = dummy*dummy       -- q = 4
          omega_6[i] = d_6[i]*( C + dummy )
          sum = sum + omega_6[i]
        end
        for i = 0, 4 do
          nlweights[eq][ [nl_index(i,is_left_biased)] ] = omega_6[i] / sum
        end
      end
    end
      
    return nlweights
  end
  return get_nonlinear_weights_LD_LBRB
end

get_nonlinear_weights_LD_l = make_get_nonlinear_weights_LD(get_beta_l, true )
get_nonlinear_weights_LD_r = make_get_nonlinear_weights_LD(get_beta_r, false)

get_nonlinear_weights_LD_LBLB_l = make_get_nonlinear_weights_LD_LBLB(get_beta_LB_l, true,  true)
get_nonlinear_weights_LD_LBLB_r = make_get_nonlinear_weights_LD_LBLB(get_beta_RB_r, false, false)

get_nonlinear_weights_LD_LBRB_l = make_get_nonlinear_weights_LD_LBRB(get_beta_RB_l, true,  false)
get_nonlinear_weights_LD_LBRB_r = make_get_nonlinear_weights_LD_LBRB(get_beta_LB_r, false, true)

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

__demand(__inline)
task get_coefficients_ECI_LBLB( nlweights : double[4][5] )
  var lcoeff0 = array(0.0, 1.0,                0.0,                   3.0/8.0, -5.0/4.0, 15.0/8.0,        0.0,           0.0,               0.0,                0.0)
  var lcoeff1 = array(0.0, 1.0,                0.0,                   0.0,     -1.0/8.0, 3.0/4.0,         3.0/8.0,       0.0,               0.0,                0.0)
  var lcoeff2 = array(0.0, -2.0*xi0 + 3.0/4.0, 2.0*xi0 + 1.0/4.0,     0.0,     0.0,       -xi0 + 1.0/4.0, 3.0/4.0,       xi0,               0.0,                0.0)
  var lcoeff3 = array(0.0, 1.0,                0.0,                   0.0,     0.0,      0.0,             -xi1+15.0/8.0, 3.0*xi1 - 5.0/4.0, -3.0*xi1 + 3.0/8.0, xi1)

  var coeffs : double[10][5]

  for eq = 0, 5 do
    for i = 0, 10 do
      coeffs[eq][i] = lcoeff0[i]*nlweights[eq][0] + lcoeff1[i]*nlweights[eq][1] + lcoeff2[i]*nlweights[eq][2] + lcoeff3[i]*nlweights[eq][3]
    end
  end

  return coeffs
end

__demand(__inline)
task get_coefficients_ECI_LBRB( nlweights : double[4][5] )
  var lcoeff0 = array(0.0,               1.0,                0.0,     xi2, -3.0*xi2 + 3.0/8.0, 3.0*xi2 - 5.0/4.0, -xi2+15.0/8.0,           0.0,            0.0,      0.0)
  var lcoeff1 = array(2.0*xi3 + 1.0/4.0, -2.0*xi3 + 3.0/4.0, 0.0,     0.0, 0.0,                xi3,               3.0/4.0,                 -xi3 + 1.0/4.0, 0.0,      0.0)
  var lcoeff2 = array(0.0,               1.0,                0.0,     0.0, 0.0,                0.0,               3.0/8.0,                 3.0/4.0,        -1.0/8.0, 0.0)
  var lcoeff3 = array(0.0,               1.0,                0.0,     0.0, 0.0,                0.0,               0.0,                     15.0/8.0,       -5.0/4.0, 3.0/8.0)

  var coeffs : double[10][5]

  for eq = 0, 5 do
    for i = 0, 10 do
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
                           block_Uinv : region(ispace(int3d), double[9]) )
where
  reads(r_prim_c), reads writes(r_prim_l, r_prim_r, alpha_l, beta_l, gamma_l, alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv)
do

  -- var t_start = c.legion_get_current_time_in_micros()

  var bounds_c = r_prim_c.ispace.bounds
  var bounds_x = r_prim_l.ispace.bounds

  var Nx = bounds_x.hi.x

  regentlib.assert(bounds_c.lo.x == 0, "Can only perform X interpolation in the X pencil")
  regentlib.assert(bounds_x.lo.x == 0, "Can only perform X interpolation in the X pencil")

  for i in alpha_l do
    var idx_c = int3d { x = i.x + ip.n_ghosts, y = i.y, z = i.z }
    var rhosos_avg : double[2] = get_rho_sos_avg_x( r_prim_c, idx_c )

    if ((not periodic_x) and (i.x == 0)) then
      var char_values : double[7][5] = get_char_values_LB_x( r_prim_c, rhosos_avg[0], rhosos_avg[1], idx_c )

      var nlweights_l = get_nonlinear_weights_LD_LBLB_l(char_values)
      var coeffs_l = get_coefficients_ECI_LBLB(nlweights_l)
      var nlweights_r = get_nonlinear_weights_LD_LBRB_r(char_values)
      var coeffs_r = get_coefficients_ECI_LBLB(nlweights_r)
  
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
                      + coeffs_l[0][9] * char_values[0][6]

      r_prim_l[i].u   = coeffs_l[1][3] * char_values[1][0]
                      + coeffs_l[1][4] * char_values[1][1]
                      + coeffs_l[1][5] * char_values[1][2]
                      + coeffs_l[1][6] * char_values[1][3]
                      + coeffs_l[1][7] * char_values[1][4]
                      + coeffs_l[1][8] * char_values[1][5]
                      + coeffs_l[1][9] * char_values[1][6]

      r_prim_l[i].v   = coeffs_l[2][3] * char_values[2][0]
                      + coeffs_l[2][4] * char_values[2][1]
                      + coeffs_l[2][5] * char_values[2][2]
                      + coeffs_l[2][6] * char_values[2][3]
                      + coeffs_l[2][7] * char_values[2][4]
                      + coeffs_l[2][8] * char_values[2][5]
                      + coeffs_l[2][9] * char_values[2][6]

      r_prim_l[i].w   = coeffs_l[3][3] * char_values[3][0]
                      + coeffs_l[3][4] * char_values[3][1]
                      + coeffs_l[3][5] * char_values[3][2]
                      + coeffs_l[3][6] * char_values[3][3]
                      + coeffs_l[3][7] * char_values[3][4]
                      + coeffs_l[3][8] * char_values[3][5]
                      + coeffs_l[3][9] * char_values[3][6]

      r_prim_l[i].p   = coeffs_l[4][3] * char_values[4][0]
                      + coeffs_l[4][4] * char_values[4][1]
                      + coeffs_l[4][5] * char_values[4][2]
                      + coeffs_l[4][6] * char_values[4][3]
                      + coeffs_l[4][7] * char_values[4][4]
                      + coeffs_l[4][8] * char_values[4][5]
                      + coeffs_l[4][9] * char_values[4][6]

      -- RHS for right sided interpolation
      r_prim_r[i].rho = coeffs_r[0][3] * char_values[0][0]
                      + coeffs_r[0][4] * char_values[0][1]
                      + coeffs_r[0][5] * char_values[0][2]
                      + coeffs_r[0][6] * char_values[0][3]
                      + coeffs_r[0][7] * char_values[0][4]
                      + coeffs_r[0][8] * char_values[0][5]
                      + coeffs_r[0][9] * char_values[0][6]

      r_prim_r[i].u   = coeffs_r[1][3] * char_values[1][0]
                      + coeffs_r[1][4] * char_values[1][1]
                      + coeffs_r[1][5] * char_values[1][2]
                      + coeffs_r[1][6] * char_values[1][3]
                      + coeffs_r[1][7] * char_values[1][4]
                      + coeffs_r[1][8] * char_values[1][5]
                      + coeffs_r[1][9] * char_values[1][6]

      r_prim_r[i].v   = coeffs_r[2][3] * char_values[2][0]
                      + coeffs_r[2][4] * char_values[2][1]
                      + coeffs_r[2][5] * char_values[2][2]
                      + coeffs_r[2][6] * char_values[2][3]
                      + coeffs_r[2][7] * char_values[2][4]
                      + coeffs_r[2][8] * char_values[2][5]
                      + coeffs_r[2][9] * char_values[2][6]

      r_prim_r[i].w   = coeffs_r[3][3] * char_values[3][0]
                      + coeffs_r[3][4] * char_values[3][1]
                      + coeffs_r[3][5] * char_values[3][2]
                      + coeffs_r[3][6] * char_values[3][3]
                      + coeffs_r[3][7] * char_values[3][4]
                      + coeffs_r[3][8] * char_values[3][5]
                      + coeffs_r[3][9] * char_values[3][6]

      r_prim_r[i].p   = coeffs_r[4][3] * char_values[4][0]
                      + coeffs_r[4][4] * char_values[4][1]
                      + coeffs_r[4][5] * char_values[4][2]
                      + coeffs_r[4][6] * char_values[4][3]
                      + coeffs_r[4][7] * char_values[4][4]
                      + coeffs_r[4][8] * char_values[4][5]
                      + coeffs_r[4][9] * char_values[4][6]
      
    elseif ((not periodic_x) and (i.x == Nx)) then
      var char_values : double[7][5] = get_char_values_RB_x( r_prim_c, rhosos_avg[0], rhosos_avg[1], idx_c )

      var nlweights_l = get_nonlinear_weights_LD_LBRB_l(char_values)
      var coeffs_l = get_coefficients_ECI_LBRB(nlweights_l)
      var nlweights_r = get_nonlinear_weights_LD_LBLB_r(char_values)
      var coeffs_r = get_coefficients_ECI_LBRB(nlweights_r)
  
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
                      + coeffs_l[0][9] * char_values[0][6]

      r_prim_l[i].u   = coeffs_l[1][3] * char_values[1][0]
                      + coeffs_l[1][4] * char_values[1][1]
                      + coeffs_l[1][5] * char_values[1][2]
                      + coeffs_l[1][6] * char_values[1][3]
                      + coeffs_l[1][7] * char_values[1][4]
                      + coeffs_l[1][8] * char_values[1][5]
                      + coeffs_l[1][9] * char_values[1][6]

      r_prim_l[i].v   = coeffs_l[2][3] * char_values[2][0]
                      + coeffs_l[2][4] * char_values[2][1]
                      + coeffs_l[2][5] * char_values[2][2]
                      + coeffs_l[2][6] * char_values[2][3]
                      + coeffs_l[2][7] * char_values[2][4]
                      + coeffs_l[2][8] * char_values[2][5]
                      + coeffs_l[2][9] * char_values[2][6]

      r_prim_l[i].w   = coeffs_l[3][3] * char_values[3][0]
                      + coeffs_l[3][4] * char_values[3][1]
                      + coeffs_l[3][5] * char_values[3][2]
                      + coeffs_l[3][6] * char_values[3][3]
                      + coeffs_l[3][7] * char_values[3][4]
                      + coeffs_l[3][8] * char_values[3][5]
                      + coeffs_l[3][9] * char_values[3][6]

      r_prim_l[i].p   = coeffs_l[4][3] * char_values[4][0]
                      + coeffs_l[4][4] * char_values[4][1]
                      + coeffs_l[4][5] * char_values[4][2]
                      + coeffs_l[4][6] * char_values[4][3]
                      + coeffs_l[4][7] * char_values[4][4]
                      + coeffs_l[4][8] * char_values[4][5]
                      + coeffs_l[4][9] * char_values[4][6]

      -- RHS for right sided interpolation
      r_prim_r[i].rho = coeffs_r[0][3] * char_values[0][0]
                      + coeffs_r[0][4] * char_values[0][1]
                      + coeffs_r[0][5] * char_values[0][2]
                      + coeffs_r[0][6] * char_values[0][3]
                      + coeffs_r[0][7] * char_values[0][4]
                      + coeffs_r[0][8] * char_values[0][5]
                      + coeffs_r[0][9] * char_values[0][6]

      r_prim_r[i].u   = coeffs_r[1][3] * char_values[1][0]
                      + coeffs_r[1][4] * char_values[1][1]
                      + coeffs_r[1][5] * char_values[1][2]
                      + coeffs_r[1][6] * char_values[1][3]
                      + coeffs_r[1][7] * char_values[1][4]
                      + coeffs_r[1][8] * char_values[1][5]
                      + coeffs_r[1][9] * char_values[1][6]

      r_prim_r[i].v   = coeffs_r[2][3] * char_values[2][0]
                      + coeffs_r[2][4] * char_values[2][1]
                      + coeffs_r[2][5] * char_values[2][2]
                      + coeffs_r[2][6] * char_values[2][3]
                      + coeffs_r[2][7] * char_values[2][4]
                      + coeffs_r[2][8] * char_values[2][5]
                      + coeffs_r[2][9] * char_values[2][6]

      r_prim_r[i].w   = coeffs_r[3][3] * char_values[3][0]
                      + coeffs_r[3][4] * char_values[3][1]
                      + coeffs_r[3][5] * char_values[3][2]
                      + coeffs_r[3][6] * char_values[3][3]
                      + coeffs_r[3][7] * char_values[3][4]
                      + coeffs_r[3][8] * char_values[3][5]
                      + coeffs_r[3][9] * char_values[3][6]

      r_prim_r[i].p   = coeffs_r[4][3] * char_values[4][0]
                      + coeffs_r[4][4] * char_values[4][1]
                      + coeffs_r[4][5] * char_values[4][2]
                      + coeffs_r[4][6] * char_values[4][3]
                      + coeffs_r[4][7] * char_values[4][4]
                      + coeffs_r[4][8] * char_values[4][5]
                      + coeffs_r[4][9] * char_values[4][6]
  
    else
      var char_values : double[6][5] = get_char_values_x( r_prim_c, rhosos_avg[0], rhosos_avg[1], idx_c )

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
  end

  -- var t_weights = c.legion_get_current_time_in_micros()

  solve_block_tridiagonal_x( alpha_l, beta_l, gamma_l, rho_avg, sos_avg, r_prim_l, block_d, block_Uinv )
  solve_block_tridiagonal_x( alpha_r, beta_r, gamma_r, rho_avg, sos_avg, r_prim_r, block_d, block_Uinv )

  solve_tridiagonal_x_v( alpha_l, beta_l, gamma_l, r_prim_l, rho_avg )
  solve_tridiagonal_x_v( alpha_r, beta_r, gamma_r, r_prim_r, sos_avg )

  solve_tridiagonal_x_w( alpha_l, beta_l, gamma_l, r_prim_l, rho_avg )
  solve_tridiagonal_x_w( alpha_r, beta_r, gamma_r, r_prim_r, sos_avg )

  -- var t_end = c.legion_get_current_time_in_micros()

  -- c.printf("X: Time to get coefficients and RHS: %12.5e\n", (t_weights-t_start)*1e-6)
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
                           block_Uinv : region(ispace(int3d), double[9]) )
where
  reads(r_prim_c), reads writes(r_prim_l, r_prim_r, alpha_l, beta_l, gamma_l, alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv)
do

  -- var t_start = c.legion_get_current_time_in_micros()

  var bounds_c = r_prim_c.ispace.bounds
  var bounds_y = r_prim_l.ispace.bounds

  var Ny = bounds_y.hi.y

  regentlib.assert(bounds_c.lo.y == 0, "Can only perform Y interpolation in the Y pencil")
  regentlib.assert(bounds_y.lo.y == 0, "Can only perform Y interpolation in the Y pencil")

  for i in alpha_l do
    var idx_c = int3d { x = i.x, y = i.y + ip.n_ghosts, z = i.z }
    var rhosos_avg : double[2] = get_rho_sos_avg_y( r_prim_c, idx_c )

    if ((not periodic_y) and (i.y == 0)) then
      var char_values : double[7][5] = get_char_values_LB_y( r_prim_c, rhosos_avg[0], rhosos_avg[1], idx_c )

      var nlweights_l = get_nonlinear_weights_LD_LBLB_l(char_values)
      var coeffs_l = get_coefficients_ECI_LBLB(nlweights_l)
      var nlweights_r = get_nonlinear_weights_LD_LBRB_r(char_values)
      var coeffs_r = get_coefficients_ECI_LBLB(nlweights_r)

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
                      + coeffs_l[0][9] * char_values[0][6]

      r_prim_l[i].u   = coeffs_l[1][3] * char_values[1][0]
                      + coeffs_l[1][4] * char_values[1][1]
                      + coeffs_l[1][5] * char_values[1][2]
                      + coeffs_l[1][6] * char_values[1][3]
                      + coeffs_l[1][7] * char_values[1][4]
                      + coeffs_l[1][8] * char_values[1][5]
                      + coeffs_l[1][9] * char_values[1][6]

      r_prim_l[i].v   = coeffs_l[2][3] * char_values[2][0]
                      + coeffs_l[2][4] * char_values[2][1]
                      + coeffs_l[2][5] * char_values[2][2]
                      + coeffs_l[2][6] * char_values[2][3]
                      + coeffs_l[2][7] * char_values[2][4]
                      + coeffs_l[2][8] * char_values[2][5]
                      + coeffs_l[2][9] * char_values[2][6]

      r_prim_l[i].w   = coeffs_l[3][3] * char_values[3][0]
                      + coeffs_l[3][4] * char_values[3][1]
                      + coeffs_l[3][5] * char_values[3][2]
                      + coeffs_l[3][6] * char_values[3][3]
                      + coeffs_l[3][7] * char_values[3][4]
                      + coeffs_l[3][8] * char_values[3][5]
                      + coeffs_l[3][9] * char_values[3][6]

      r_prim_l[i].p   = coeffs_l[4][3] * char_values[4][0]
                      + coeffs_l[4][4] * char_values[4][1]
                      + coeffs_l[4][5] * char_values[4][2]
                      + coeffs_l[4][6] * char_values[4][3]
                      + coeffs_l[4][7] * char_values[4][4]
                      + coeffs_l[4][8] * char_values[4][5]
                      + coeffs_l[4][9] * char_values[4][6]

      -- RHS for right sided interpolation
      r_prim_r[i].rho = coeffs_r[0][3] * char_values[0][0]
                      + coeffs_r[0][4] * char_values[0][1]
                      + coeffs_r[0][5] * char_values[0][2]
                      + coeffs_r[0][6] * char_values[0][3]
                      + coeffs_r[0][7] * char_values[0][4]
                      + coeffs_r[0][8] * char_values[0][5]
                      + coeffs_r[0][9] * char_values[0][6]

      r_prim_r[i].u   = coeffs_r[1][3] * char_values[1][0]
                      + coeffs_r[1][4] * char_values[1][1]
                      + coeffs_r[1][5] * char_values[1][2]
                      + coeffs_r[1][6] * char_values[1][3]
                      + coeffs_r[1][7] * char_values[1][4]
                      + coeffs_r[1][8] * char_values[1][5]
                      + coeffs_r[1][9] * char_values[1][6]

      r_prim_r[i].v   = coeffs_r[2][3] * char_values[2][0]
                      + coeffs_r[2][4] * char_values[2][1]
                      + coeffs_r[2][5] * char_values[2][2]
                      + coeffs_r[2][6] * char_values[2][3]
                      + coeffs_r[2][7] * char_values[2][4]
                      + coeffs_r[2][8] * char_values[2][5]
                      + coeffs_r[2][9] * char_values[2][6]

      r_prim_r[i].w   = coeffs_r[3][3] * char_values[3][0]
                      + coeffs_r[3][4] * char_values[3][1]
                      + coeffs_r[3][5] * char_values[3][2]
                      + coeffs_r[3][6] * char_values[3][3]
                      + coeffs_r[3][7] * char_values[3][4]
                      + coeffs_r[3][8] * char_values[3][5]
                      + coeffs_r[3][9] * char_values[3][6]

      r_prim_r[i].p   = coeffs_r[4][3] * char_values[4][0]
                      + coeffs_r[4][4] * char_values[4][1]
                      + coeffs_r[4][5] * char_values[4][2]
                      + coeffs_r[4][6] * char_values[4][3]
                      + coeffs_r[4][7] * char_values[4][4]
                      + coeffs_r[4][8] * char_values[4][5]
                      + coeffs_r[4][9] * char_values[4][6]

    elseif ((not periodic_y) and (i.y == Ny)) then
      var char_values : double[7][5] = get_char_values_RB_y( r_prim_c, rhosos_avg[0], rhosos_avg[1], idx_c )

      var nlweights_l = get_nonlinear_weights_LD_LBRB_l(char_values)
      var coeffs_l = get_coefficients_ECI_LBRB(nlweights_l)
      var nlweights_r = get_nonlinear_weights_LD_LBLB_r(char_values)
      var coeffs_r = get_coefficients_ECI_LBRB(nlweights_r)

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
                      + coeffs_l[0][9] * char_values[0][6]

      r_prim_l[i].u   = coeffs_l[1][3] * char_values[1][0]
                      + coeffs_l[1][4] * char_values[1][1]
                      + coeffs_l[1][5] * char_values[1][2]
                      + coeffs_l[1][6] * char_values[1][3]
                      + coeffs_l[1][7] * char_values[1][4]
                      + coeffs_l[1][8] * char_values[1][5]
                      + coeffs_l[1][9] * char_values[1][6]

      r_prim_l[i].v   = coeffs_l[2][3] * char_values[2][0]
                      + coeffs_l[2][4] * char_values[2][1]
                      + coeffs_l[2][5] * char_values[2][2]
                      + coeffs_l[2][6] * char_values[2][3]
                      + coeffs_l[2][7] * char_values[2][4]
                      + coeffs_l[2][8] * char_values[2][5]
                      + coeffs_l[2][9] * char_values[2][6]

      r_prim_l[i].w   = coeffs_l[3][3] * char_values[3][0]
                      + coeffs_l[3][4] * char_values[3][1]
                      + coeffs_l[3][5] * char_values[3][2]
                      + coeffs_l[3][6] * char_values[3][3]
                      + coeffs_l[3][7] * char_values[3][4]
                      + coeffs_l[3][8] * char_values[3][5]
                      + coeffs_l[3][9] * char_values[3][6]

      r_prim_l[i].p   = coeffs_l[4][3] * char_values[4][0]
                      + coeffs_l[4][4] * char_values[4][1]
                      + coeffs_l[4][5] * char_values[4][2]
                      + coeffs_l[4][6] * char_values[4][3]
                      + coeffs_l[4][7] * char_values[4][4]
                      + coeffs_l[4][8] * char_values[4][5]
                      + coeffs_l[4][9] * char_values[4][6]

      -- RHS for right sided interpolation
      r_prim_r[i].rho = coeffs_r[0][3] * char_values[0][0]
                      + coeffs_r[0][4] * char_values[0][1]
                      + coeffs_r[0][5] * char_values[0][2]
                      + coeffs_r[0][6] * char_values[0][3]
                      + coeffs_r[0][7] * char_values[0][4]
                      + coeffs_r[0][8] * char_values[0][5]
                      + coeffs_r[0][9] * char_values[0][6]

      r_prim_r[i].u   = coeffs_r[1][3] * char_values[1][0]
                      + coeffs_r[1][4] * char_values[1][1]
                      + coeffs_r[1][5] * char_values[1][2]
                      + coeffs_r[1][6] * char_values[1][3]
                      + coeffs_r[1][7] * char_values[1][4]
                      + coeffs_r[1][8] * char_values[1][5]
                      + coeffs_r[1][9] * char_values[1][6]

      r_prim_r[i].v   = coeffs_r[2][3] * char_values[2][0]
                      + coeffs_r[2][4] * char_values[2][1]
                      + coeffs_r[2][5] * char_values[2][2]
                      + coeffs_r[2][6] * char_values[2][3]
                      + coeffs_r[2][7] * char_values[2][4]
                      + coeffs_r[2][8] * char_values[2][5]
                      + coeffs_r[2][9] * char_values[2][6]

      r_prim_r[i].w   = coeffs_r[3][3] * char_values[3][0]
                      + coeffs_r[3][4] * char_values[3][1]
                      + coeffs_r[3][5] * char_values[3][2]
                      + coeffs_r[3][6] * char_values[3][3]
                      + coeffs_r[3][7] * char_values[3][4]
                      + coeffs_r[3][8] * char_values[3][5]
                      + coeffs_r[3][9] * char_values[3][6]

      r_prim_r[i].p   = coeffs_r[4][3] * char_values[4][0]
                      + coeffs_r[4][4] * char_values[4][1]
                      + coeffs_r[4][5] * char_values[4][2]
                      + coeffs_r[4][6] * char_values[4][3]
                      + coeffs_r[4][7] * char_values[4][4]
                      + coeffs_r[4][8] * char_values[4][5]
                      + coeffs_r[4][9] * char_values[4][6]

    else
      var char_values : double[6][5] = get_char_values_y( r_prim_c, rhosos_avg[0], rhosos_avg[1], idx_c )
  
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
  end
      
  -- var t_weights = c.legion_get_current_time_in_micros()

  solve_block_tridiagonal_y( alpha_l, beta_l, gamma_l, rho_avg, sos_avg, r_prim_l, block_d, block_Uinv )
  solve_block_tridiagonal_y( alpha_r, beta_r, gamma_r, rho_avg, sos_avg, r_prim_r, block_d, block_Uinv )

  solve_tridiagonal_y_u( alpha_l, beta_l, gamma_l, r_prim_l, rho_avg )
  solve_tridiagonal_y_u( alpha_r, beta_r, gamma_r, r_prim_r, sos_avg )

  solve_tridiagonal_y_w( alpha_l, beta_l, gamma_l, r_prim_l, rho_avg )
  solve_tridiagonal_y_w( alpha_r, beta_r, gamma_r, r_prim_r, sos_avg )

  -- var t_end = c.legion_get_current_time_in_micros()

  -- c.printf("Y: Time to get coefficients and RHS: %12.5e\n", (t_weights-t_start)*1e-6)
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
                           block_Uinv : region(ispace(int3d), double[9]) )
where
  reads(r_prim_c), reads writes(r_prim_l, r_prim_r, alpha_l, beta_l, gamma_l, alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv)
do

  -- var t_start = c.legion_get_current_time_in_micros()

  var bounds_c = r_prim_c.ispace.bounds
  var bounds_z = r_prim_l.ispace.bounds

  var Nz = bounds_z.hi.z

  regentlib.assert(bounds_c.lo.z == 0, "Can only perform Z interpolation in the Z pencil")
  regentlib.assert(bounds_z.lo.z == 0, "Can only perform Z interpolation in the Z pencil")

  for i in alpha_l do
    var idx_c = int3d { x = i.x, y = i.y, z = i.z + ip.n_ghosts }
    var rhosos_avg : double[2] = get_rho_sos_avg_z( r_prim_c, idx_c )

    if ((not periodic_z) and (i.z == 0)) then
      var char_values : double[7][5] = get_char_values_LB_z( r_prim_c, rhosos_avg[0], rhosos_avg[1], idx_c )

      var nlweights_l = get_nonlinear_weights_LD_LBLB_l(char_values)
      var coeffs_l = get_coefficients_ECI_LBLB(nlweights_l)
      var nlweights_r = get_nonlinear_weights_LD_LBRB_r(char_values)
      var coeffs_r = get_coefficients_ECI_LBLB(nlweights_r)

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
                      + coeffs_l[0][9] * char_values[0][6]

      r_prim_l[i].u   = coeffs_l[1][3] * char_values[1][0]
                      + coeffs_l[1][4] * char_values[1][1]
                      + coeffs_l[1][5] * char_values[1][2]
                      + coeffs_l[1][6] * char_values[1][3]
                      + coeffs_l[1][7] * char_values[1][4]
                      + coeffs_l[1][8] * char_values[1][5]
                      + coeffs_l[1][9] * char_values[1][6]

      r_prim_l[i].v   = coeffs_l[2][3] * char_values[2][0]
                      + coeffs_l[2][4] * char_values[2][1]
                      + coeffs_l[2][5] * char_values[2][2]
                      + coeffs_l[2][6] * char_values[2][3]
                      + coeffs_l[2][7] * char_values[2][4]
                      + coeffs_l[2][8] * char_values[2][5]
                      + coeffs_l[2][9] * char_values[2][6]

      r_prim_l[i].w   = coeffs_l[3][3] * char_values[3][0]
                      + coeffs_l[3][4] * char_values[3][1]
                      + coeffs_l[3][5] * char_values[3][2]
                      + coeffs_l[3][6] * char_values[3][3]
                      + coeffs_l[3][7] * char_values[3][4]
                      + coeffs_l[3][8] * char_values[3][5]
                      + coeffs_l[3][9] * char_values[3][6]

      r_prim_l[i].p   = coeffs_l[4][3] * char_values[4][0]
                      + coeffs_l[4][4] * char_values[4][1]
                      + coeffs_l[4][5] * char_values[4][2]
                      + coeffs_l[4][6] * char_values[4][3]
                      + coeffs_l[4][7] * char_values[4][4]
                      + coeffs_l[4][8] * char_values[4][5]
                      + coeffs_l[4][9] * char_values[4][6]

      -- RHS for right sided interpolation
      r_prim_r[i].rho = coeffs_r[0][3] * char_values[0][0]
                      + coeffs_r[0][4] * char_values[0][1]
                      + coeffs_r[0][5] * char_values[0][2]
                      + coeffs_r[0][6] * char_values[0][3]
                      + coeffs_r[0][7] * char_values[0][4]
                      + coeffs_r[0][8] * char_values[0][5]
                      + coeffs_r[0][9] * char_values[0][6]

      r_prim_r[i].u   = coeffs_r[1][3] * char_values[1][0]
                      + coeffs_r[1][4] * char_values[1][1]
                      + coeffs_r[1][5] * char_values[1][2]
                      + coeffs_r[1][6] * char_values[1][3]
                      + coeffs_r[1][7] * char_values[1][4]
                      + coeffs_r[1][8] * char_values[1][5]
                      + coeffs_r[1][9] * char_values[1][6]

      r_prim_r[i].v   = coeffs_r[2][3] * char_values[2][0]
                      + coeffs_r[2][4] * char_values[2][1]
                      + coeffs_r[2][5] * char_values[2][2]
                      + coeffs_r[2][6] * char_values[2][3]
                      + coeffs_r[2][7] * char_values[2][4]
                      + coeffs_r[2][8] * char_values[2][5]
                      + coeffs_r[2][9] * char_values[2][6]

      r_prim_r[i].w   = coeffs_r[3][3] * char_values[3][0]
                      + coeffs_r[3][4] * char_values[3][1]
                      + coeffs_r[3][5] * char_values[3][2]
                      + coeffs_r[3][6] * char_values[3][3]
                      + coeffs_r[3][7] * char_values[3][4]
                      + coeffs_r[3][8] * char_values[3][5]
                      + coeffs_r[3][9] * char_values[3][6]

      r_prim_r[i].p   = coeffs_r[4][3] * char_values[4][0]
                      + coeffs_r[4][4] * char_values[4][1]
                      + coeffs_r[4][5] * char_values[4][2]
                      + coeffs_r[4][6] * char_values[4][3]
                      + coeffs_r[4][7] * char_values[4][4]
                      + coeffs_r[4][8] * char_values[4][5]
                      + coeffs_r[4][9] * char_values[4][6]

    elseif ((not periodic_z) and (i.z == Nz)) then
      var char_values : double[7][5] = get_char_values_RB_z( r_prim_c, rhosos_avg[0], rhosos_avg[1], idx_c )
      
      var nlweights_l = get_nonlinear_weights_LD_LBRB_l(char_values)
      var coeffs_l = get_coefficients_ECI_LBRB(nlweights_l)
      var nlweights_r = get_nonlinear_weights_LD_LBLB_r(char_values)
      var coeffs_r = get_coefficients_ECI_LBRB(nlweights_r)

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
                      + coeffs_l[0][9] * char_values[0][6]

      r_prim_l[i].u   = coeffs_l[1][3] * char_values[1][0]
                      + coeffs_l[1][4] * char_values[1][1]
                      + coeffs_l[1][5] * char_values[1][2]
                      + coeffs_l[1][6] * char_values[1][3]
                      + coeffs_l[1][7] * char_values[1][4]
                      + coeffs_l[1][8] * char_values[1][5]
                      + coeffs_l[1][9] * char_values[1][6]

      r_prim_l[i].v   = coeffs_l[2][3] * char_values[2][0]
                      + coeffs_l[2][4] * char_values[2][1]
                      + coeffs_l[2][5] * char_values[2][2]
                      + coeffs_l[2][6] * char_values[2][3]
                      + coeffs_l[2][7] * char_values[2][4]
                      + coeffs_l[2][8] * char_values[2][5]
                      + coeffs_l[2][9] * char_values[2][6]

      r_prim_l[i].w   = coeffs_l[3][3] * char_values[3][0]
                      + coeffs_l[3][4] * char_values[3][1]
                      + coeffs_l[3][5] * char_values[3][2]
                      + coeffs_l[3][6] * char_values[3][3]
                      + coeffs_l[3][7] * char_values[3][4]
                      + coeffs_l[3][8] * char_values[3][5]
                      + coeffs_l[3][9] * char_values[3][6]

      r_prim_l[i].p   = coeffs_l[4][3] * char_values[4][0]
                      + coeffs_l[4][4] * char_values[4][1]
                      + coeffs_l[4][5] * char_values[4][2]
                      + coeffs_l[4][6] * char_values[4][3]
                      + coeffs_l[4][7] * char_values[4][4]
                      + coeffs_l[4][8] * char_values[4][5]
                      + coeffs_l[4][9] * char_values[4][6]

      -- RHS for right sided interpolation
      r_prim_r[i].rho = coeffs_r[0][3] * char_values[0][0]
                      + coeffs_r[0][4] * char_values[0][1]
                      + coeffs_r[0][5] * char_values[0][2]
                      + coeffs_r[0][6] * char_values[0][3]
                      + coeffs_r[0][7] * char_values[0][4]
                      + coeffs_r[0][8] * char_values[0][5]
                      + coeffs_r[0][9] * char_values[0][6]

      r_prim_r[i].u   = coeffs_r[1][3] * char_values[1][0]
                      + coeffs_r[1][4] * char_values[1][1]
                      + coeffs_r[1][5] * char_values[1][2]
                      + coeffs_r[1][6] * char_values[1][3]
                      + coeffs_r[1][7] * char_values[1][4]
                      + coeffs_r[1][8] * char_values[1][5]
                      + coeffs_r[1][9] * char_values[1][6]

      r_prim_r[i].v   = coeffs_r[2][3] * char_values[2][0]
                      + coeffs_r[2][4] * char_values[2][1]
                      + coeffs_r[2][5] * char_values[2][2]
                      + coeffs_r[2][6] * char_values[2][3]
                      + coeffs_r[2][7] * char_values[2][4]
                      + coeffs_r[2][8] * char_values[2][5]
                      + coeffs_r[2][9] * char_values[2][6]

      r_prim_r[i].w   = coeffs_r[3][3] * char_values[3][0]
                      + coeffs_r[3][4] * char_values[3][1]
                      + coeffs_r[3][5] * char_values[3][2]
                      + coeffs_r[3][6] * char_values[3][3]
                      + coeffs_r[3][7] * char_values[3][4]
                      + coeffs_r[3][8] * char_values[3][5]
                      + coeffs_r[3][9] * char_values[3][6]

      r_prim_r[i].p   = coeffs_r[4][3] * char_values[4][0]
                      + coeffs_r[4][4] * char_values[4][1]
                      + coeffs_r[4][5] * char_values[4][2]
                      + coeffs_r[4][6] * char_values[4][3]
                      + coeffs_r[4][7] * char_values[4][4]
                      + coeffs_r[4][8] * char_values[4][5]
                      + coeffs_r[4][9] * char_values[4][6]

    else
      var char_values : double[6][5] = get_char_values_z( r_prim_c, rhosos_avg[0], rhosos_avg[1], idx_c )

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
  end

  -- var t_weights = c.legion_get_current_time_in_micros()

  solve_block_tridiagonal_z( alpha_l, beta_l, gamma_l, rho_avg, sos_avg, r_prim_l, block_d, block_Uinv )
  solve_block_tridiagonal_z( alpha_r, beta_r, gamma_r, rho_avg, sos_avg, r_prim_r, block_d, block_Uinv )

  solve_tridiagonal_z_u( alpha_l, beta_l, gamma_l, r_prim_l, rho_avg )
  solve_tridiagonal_z_u( alpha_r, beta_r, gamma_r, r_prim_r, sos_avg )

  solve_tridiagonal_z_v( alpha_l, beta_l, gamma_l, r_prim_l, rho_avg )
  solve_tridiagonal_z_v( alpha_r, beta_r, gamma_r, r_prim_r, sos_avg )

  -- var t_end = c.legion_get_current_time_in_micros()


  -- c.printf("Z: Time to get coefficients and RHS: %12.5e\n", (t_weights-t_start)*1e-6)
  -- c.printf("Z: Time for block tridiagonal solves: %12.5e\n", (t_end-t_weights)*1e-6)
  -- c.printf("Z: Time to get the WCHR interpolation: %12.5e\n", (t_end-t_start)*1e-6)
  return 1
end

return ip
