import "regent"

require("fields")
require("derivatives")
require("interpolation")
require("SOE")

local superlu = require("superlu_util")
local problem = require("problem")

local c     = regentlib.c
local cmath = terralib.includec("math.h")


-- Make the node and midpoint-node differencing tasks (Using pentadiagonal solver for this instead of tridiagonal solver)
alpha06d1 = 1.0/3.0
beta06d1  = 0.0
a06d1 = ( 14.0/ 9.0)/2.0
b06d1 = (  1.0/ 9.0)/4.0
c06d1 = (  0.0/100.0)/6.0

-- Compact MND finite difference
-- alpha06MND = -1.0/12.0
-- beta06MND  = 0.0
-- a06MND = 16.0/9.0
-- b06MND = (-17.0/18.0)/2.0
-- c06MND = (0.0)/3.0

-- Compact staggered finite difference
alpha06MND = 9.0/62.0
beta06MND  = 0.0
a06MND = 63.0/62.0
b06MND = (0.0/18.0)/2.0
c06MND = (17.0/62.0)/3.0

-- Explicit MND finite difference
-- alpha06MND = 0.0
-- beta06MND  = 0.0
-- a06MND = 3.0/2.0
-- b06MND = (-3.0/10.0)
-- c06MND = (1.0)/30.0

local r_flux   = regentlib.newsymbol(region(ispace(int3d), conserved), "r_flux")
local r_flux_e = regentlib.newsymbol(region(ispace(int3d), conserved), "r_flux_e")
local r_der    = regentlib.newsymbol(region(ispace(int3d), conserved), "r_der")

local ddx_MND_rho  = make_ddx_MND(r_flux, r_flux_e, "rho",  r_der, "rho",  problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, a06MND, b06MND, c06MND)
local ddx_MND_rhou = make_ddx_MND(r_flux, r_flux_e, "rhou", r_der, "rhou", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, a06MND, b06MND, c06MND)
local ddx_MND_rhov = make_ddx_MND(r_flux, r_flux_e, "rhov", r_der, "rhov", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, a06MND, b06MND, c06MND)
local ddx_MND_rhow = make_ddx_MND(r_flux, r_flux_e, "rhow", r_der, "rhow", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, a06MND, b06MND, c06MND)
local ddx_MND_rhoE = make_ddx_MND(r_flux, r_flux_e, "rhoE", r_der, "rhoE", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, a06MND, b06MND, c06MND)

local ddy_MND_rho  = make_ddy_MND(r_flux, r_flux_e, "rho",  r_der, "rho",  problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, a06MND, b06MND, c06MND)
local ddy_MND_rhou = make_ddy_MND(r_flux, r_flux_e, "rhou", r_der, "rhou", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, a06MND, b06MND, c06MND)
local ddy_MND_rhov = make_ddy_MND(r_flux, r_flux_e, "rhov", r_der, "rhov", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, a06MND, b06MND, c06MND)
local ddy_MND_rhow = make_ddy_MND(r_flux, r_flux_e, "rhow", r_der, "rhow", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, a06MND, b06MND, c06MND)
local ddy_MND_rhoE = make_ddy_MND(r_flux, r_flux_e, "rhoE", r_der, "rhoE", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, a06MND, b06MND, c06MND)

local ddz_MND_rho  = make_ddz_MND(r_flux, r_flux_e, "rho",  r_der, "rho",  problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, a06MND, b06MND, c06MND)
local ddz_MND_rhou = make_ddz_MND(r_flux, r_flux_e, "rhou", r_der, "rhou", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, a06MND, b06MND, c06MND)
local ddz_MND_rhov = make_ddz_MND(r_flux, r_flux_e, "rhov", r_der, "rhov", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, a06MND, b06MND, c06MND)
local ddz_MND_rhow = make_ddz_MND(r_flux, r_flux_e, "rhow", r_der, "rhow", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, a06MND, b06MND, c06MND)
local ddz_MND_rhoE = make_ddz_MND(r_flux, r_flux_e, "rhoE", r_der, "rhoE", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, a06MND, b06MND, c06MND)

local r_prim   = regentlib.newsymbol(region(ispace(int3d), primitive), "r_prim")
local r_prim_e = regentlib.newsymbol(region(ispace(int3d), primitive), "r_prim_e")
local r_der    = regentlib.newsymbol(region(ispace(int3d), tensor2),   "r_der")

local ddx_u   = make_ddx(r_prim, "u", r_der, "_11", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, a06d1, b06d1, c06d1)
local ddy_u   = make_ddy(r_prim, "u", r_der, "_12", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, a06d1, b06d1, c06d1)
local ddz_u   = make_ddz(r_prim, "u", r_der, "_13", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, a06d1, b06d1, c06d1)

local ddx_v   = make_ddx(r_prim, "v", r_der, "_21", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, a06d1, b06d1, c06d1)
local ddy_v   = make_ddy(r_prim, "v", r_der, "_22", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, a06d1, b06d1, c06d1)
local ddz_v   = make_ddz(r_prim, "v", r_der, "_23", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, a06d1, b06d1, c06d1)

local ddx_w   = make_ddx(r_prim, "w", r_der, "_31", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, a06d1, b06d1, c06d1)
local ddy_w   = make_ddy(r_prim, "w", r_der, "_32", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, a06d1, b06d1, c06d1)
local ddz_w   = make_ddz(r_prim, "w", r_der, "_33", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, a06d1, b06d1, c06d1)

-----------------------------------------------------

task add_xflux_der_to_rhs( r_cnsr     : region(ispace(int3d), conserved),
                           r_prim_c   : region(ispace(int3d), primitive),
                           r_prim_l_x : region(ispace(int3d), primitive),
                           r_prim_r_x : region(ispace(int3d), primitive),
                           r_rhs_l_x  : region(ispace(int3d), primitive),
                           r_rhs_r_x  : region(ispace(int3d), primitive),
                           r_flux_c   : region(ispace(int3d), conserved),
                           r_flux_e_x : region(ispace(int3d), conserved),
                           r_fder_c_x : region(ispace(int3d), conserved),
                           r_rhs      : region(ispace(int3d), conserved),
                           LU_x       : region(ispace(int3d), LU_struct),
                           slu_l_x    : region(ispace(int2d), superlu.c.superlu_vars_t),
                           slu_r_x    : region(ispace(int2d), superlu.c.superlu_vars_t),
                           matrix_l_x : region(ispace(int2d), superlu.CSR_matrix),
                           matrix_r_x : region(ispace(int2d), superlu.CSR_matrix),
                           Nx         : int64,
                           Ny         : int64,
                           Nz         : int64 )
where
  reads(r_cnsr, r_prim_c, LU_x),
  reads writes(r_prim_l_x, r_prim_r_x, r_rhs_l_x, r_rhs_r_x, r_flux_c, r_flux_e_x, r_fder_c_x, r_rhs, slu_l_x, slu_r_x, matrix_l_x, matrix_r_x)
do

  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1

  if (nx >= 8) then
    WCHR_interpolation_x( r_prim_c, r_prim_l_x, r_prim_r_x, r_rhs_l_x, r_rhs_r_x, matrix_l_x, matrix_r_x, slu_l_x, slu_r_x, Nx, Ny, Nz )
    positivity_enforcer_x( r_prim_c, r_prim_l_x, r_prim_r_x, Nx, Ny, Nz )
    HLLC_x( r_prim_l_x, r_prim_r_x, r_flux_e_x )
    get_xfluxes_r( r_prim_c, r_cnsr, r_flux_c )
   
    ddx_MND_rho ( r_flux_c, r_flux_e_x, r_fder_c_x, LU_x )
    ddx_MND_rhou( r_flux_c, r_flux_e_x, r_fder_c_x, LU_x )
    ddx_MND_rhov( r_flux_c, r_flux_e_x, r_fder_c_x, LU_x )
    ddx_MND_rhow( r_flux_c, r_flux_e_x, r_fder_c_x, LU_x )
    ddx_MND_rhoE( r_flux_c, r_flux_e_x, r_fder_c_x, LU_x )

    for i in r_rhs do
      r_rhs[i].rho  -= r_fder_c_x[i].rho
      r_rhs[i].rhou -= r_fder_c_x[i].rhou
      r_rhs[i].rhov -= r_fder_c_x[i].rhov
      r_rhs[i].rhow -= r_fder_c_x[i].rhow
      r_rhs[i].rhoE -= r_fder_c_x[i].rhoE
    end
  end
end

task add_yflux_der_to_rhs( r_cnsr     : region(ispace(int3d), conserved),
                           r_prim_c   : region(ispace(int3d), primitive),
                           r_prim_l_y : region(ispace(int3d), primitive),
                           r_prim_r_y : region(ispace(int3d), primitive),
                           r_rhs_l_y  : region(ispace(int3d), primitive),
                           r_rhs_r_y  : region(ispace(int3d), primitive),
                           r_flux_c   : region(ispace(int3d), conserved),
                           r_flux_e_y : region(ispace(int3d), conserved),
                           r_fder_c_y : region(ispace(int3d), conserved),
                           r_rhs      : region(ispace(int3d), conserved),
                           LU_y       : region(ispace(int3d), LU_struct),
                           slu_l_y    : region(ispace(int2d), superlu.c.superlu_vars_t),
                           slu_r_y    : region(ispace(int2d), superlu.c.superlu_vars_t),
                           matrix_l_y : region(ispace(int2d), superlu.CSR_matrix),
                           matrix_r_y : region(ispace(int2d), superlu.CSR_matrix),
                           Nx         : int64,
                           Ny         : int64,
                           Nz         : int64 )
where
  reads(r_cnsr, r_prim_c, LU_y),
  reads writes(r_prim_l_y, r_prim_r_y, r_rhs_l_y, r_rhs_r_y, r_flux_c, r_flux_e_y, r_fder_c_y, r_rhs, slu_l_y, slu_r_y, matrix_l_y, matrix_r_y)
do

  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1

  if (ny >= 8) then
    WCHR_interpolation_y( r_prim_c, r_prim_l_y, r_prim_r_y, r_rhs_l_y, r_rhs_r_y, matrix_l_y, matrix_r_y, slu_l_y, slu_r_y, Nx, Ny, Nz )
    positivity_enforcer_y( r_prim_c, r_prim_l_y, r_prim_r_y, Nx, Ny, Nz )
    HLLC_y( r_prim_l_y, r_prim_r_y, r_flux_e_y )
    get_yfluxes_r( r_prim_c, r_cnsr, r_flux_c )
    
    ddy_MND_rho ( r_flux_c, r_flux_e_y, r_fder_c_y, LU_y )
    ddy_MND_rhou( r_flux_c, r_flux_e_y, r_fder_c_y, LU_y )
    ddy_MND_rhov( r_flux_c, r_flux_e_y, r_fder_c_y, LU_y )
    ddy_MND_rhow( r_flux_c, r_flux_e_y, r_fder_c_y, LU_y )
    ddy_MND_rhoE( r_flux_c, r_flux_e_y, r_fder_c_y, LU_y )
    
    for i in r_rhs do
      r_rhs[i].rho  -= r_fder_c_y[i].rho
      r_rhs[i].rhou -= r_fder_c_y[i].rhou
      r_rhs[i].rhov -= r_fder_c_y[i].rhov
      r_rhs[i].rhow -= r_fder_c_y[i].rhow
      r_rhs[i].rhoE -= r_fder_c_y[i].rhoE
    end
  end
end

task add_zflux_der_to_rhs( r_cnsr     : region(ispace(int3d), conserved),
                           r_prim_c   : region(ispace(int3d), primitive),
                           r_prim_l_z : region(ispace(int3d), primitive),
                           r_prim_r_z : region(ispace(int3d), primitive),
                           r_rhs_l_z  : region(ispace(int3d), primitive),
                           r_rhs_r_z  : region(ispace(int3d), primitive),
                           r_flux_c   : region(ispace(int3d), conserved),
                           r_flux_e_z : region(ispace(int3d), conserved),
                           r_fder_c_z : region(ispace(int3d), conserved),
                           r_rhs      : region(ispace(int3d), conserved),
                           LU_z       : region(ispace(int3d), LU_struct),
                           slu_l_z    : region(ispace(int2d), superlu.c.superlu_vars_t),
                           slu_r_z    : region(ispace(int2d), superlu.c.superlu_vars_t),
                           matrix_l_z : region(ispace(int2d), superlu.CSR_matrix),
                           matrix_r_z : region(ispace(int2d), superlu.CSR_matrix),
                           Nx         : int64,
                           Ny         : int64,
                           Nz         : int64 )
where
  reads(r_cnsr, r_prim_c, LU_z),
  reads writes(r_prim_l_z, r_prim_r_z, r_rhs_l_z, r_rhs_r_z, r_flux_c, r_flux_e_z, r_fder_c_z, r_rhs, slu_l_z, slu_r_z, matrix_l_z, matrix_r_z)
do

  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1

  if (nz >= 8) then
    WCHR_interpolation_z( r_prim_c, r_prim_l_z, r_prim_r_z, r_rhs_l_z, r_rhs_r_z, matrix_l_z, matrix_r_z, slu_l_z, slu_r_z, Nx, Ny, Nz )
    positivity_enforcer_z( r_prim_c, r_prim_l_z, r_prim_r_z, Nx, Ny, Nz )
    HLLC_z( r_prim_l_z, r_prim_r_z, r_flux_e_z )
    get_zfluxes_r( r_prim_c, r_cnsr, r_flux_c )

    ddz_MND_rho ( r_flux_c, r_flux_e_z, r_fder_c_z, LU_z )
    ddz_MND_rhou( r_flux_c, r_flux_e_z, r_fder_c_z, LU_z )
    ddz_MND_rhov( r_flux_c, r_flux_e_z, r_fder_c_z, LU_z )
    ddz_MND_rhow( r_flux_c, r_flux_e_z, r_fder_c_z, LU_z )
    ddz_MND_rhoE( r_flux_c, r_flux_e_z, r_fder_c_z, LU_z )

    for i in r_rhs do
      r_rhs[i].rho  -= r_fder_c_z[i].rho
      r_rhs[i].rhou -= r_fder_c_z[i].rhou
      r_rhs[i].rhov -= r_fder_c_z[i].rhov
      r_rhs[i].rhow -= r_fder_c_z[i].rhow
      r_rhs[i].rhoE -= r_fder_c_z[i].rhoE
    end
  end
end

task update_substep( r_cnsr : region(ispace(int3d), conserved),
                     r_rhs  : region(ispace(int3d), conserved),
                     Q_rhs  : region(ispace(int3d), conserved),
                     dt     : double,
                     A      : double,
                     B      : double )
where
  reads (r_rhs), reads writes(r_cnsr, Q_rhs)
do

  for i in r_rhs do
    Q_rhs[i].rho = dt * r_rhs[i].rho + A*Q_rhs[i].rho
    r_cnsr[i].rho += B*Q_rhs[i].rho

    Q_rhs[i].rhou = dt * r_rhs[i].rhou + A*Q_rhs[i].rhou
    r_cnsr[i].rhou += B*Q_rhs[i].rhou

    Q_rhs[i].rhov = dt * r_rhs[i].rhov + A*Q_rhs[i].rhov
    r_cnsr[i].rhov += B*Q_rhs[i].rhov

    Q_rhs[i].rhow = dt * r_rhs[i].rhow + A*Q_rhs[i].rhow
    r_cnsr[i].rhow += B*Q_rhs[i].rhow

    Q_rhs[i].rhoE = dt * r_rhs[i].rhoE + A*Q_rhs[i].rhoE
    r_cnsr[i].rhoE += B*Q_rhs[i].rhoE
  end
end

task get_velocity_x_derivatives( r_prim_c : region(ispace(int3d), primitive),
                                 r_duidxj : region(ispace(int3d), tensor2),
                                 LU_x     : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c, LU_x), reads writes (r_duidxj.{_11, _21, _31})
do
  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1

  var token = 1
  if (nx < 8) then
    for i in r_duidxj do
      r_duidxj[i]._11 = 0.0
      r_duidxj[i]._21 = 0.0
      r_duidxj[i]._31 = 0.0
    end
    return token
  end

  token += ddx_u(r_prim_c, r_duidxj, LU_x)
  token += ddx_v(r_prim_c, r_duidxj, LU_x)
  token += ddx_w(r_prim_c, r_duidxj, LU_x) 

  return token
end

task get_velocity_y_derivatives( r_prim_c : region(ispace(int3d), primitive),
                                 r_duidxj : region(ispace(int3d), tensor2),
                                 LU_y     : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c, LU_y), reads writes (r_duidxj.{_12, _22, _32})
do
  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1

  var token = 1
  if (ny < 8) then
    for i in r_duidxj do
      r_duidxj[i]._12 = 0.0
      r_duidxj[i]._22 = 0.0
      r_duidxj[i]._32 = 0.0
    end
    return token
  end

  token += ddy_u(r_prim_c, r_duidxj, LU_y)
  token += ddy_v(r_prim_c, r_duidxj, LU_y)
  token += ddy_w(r_prim_c, r_duidxj, LU_y) 

  return token
end

task get_velocity_z_derivatives( r_prim_c : region(ispace(int3d), primitive),
                                 r_duidxj : region(ispace(int3d), tensor2),
                                 LU_z     : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c, LU_z), reads writes (r_duidxj.{_13, _23, _33})
do
  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1

  var token = 1
  if (nz < 8) then
    for i in r_duidxj do
      r_duidxj[i]._13 = 0.0
      r_duidxj[i]._23 = 0.0
      r_duidxj[i]._33 = 0.0
    end
    return token
  end

  token += ddz_u(r_prim_c, r_duidxj, LU_z)
  token += ddz_v(r_prim_c, r_duidxj, LU_z)
  token += ddz_w(r_prim_c, r_duidxj, LU_z) 

  return token
end
