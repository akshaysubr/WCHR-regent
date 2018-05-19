import "regent"

require("fields")
require("derivatives")
interpolation = require("interpolation")
require("SOE")

local problem = require("problem")

local c     = regentlib.c
local cmath = terralib.includec("math.h")

local viscous = problem.viscous

local r_flux   = regentlib.newsymbol(region(ispace(int3d), conserved), "r_flux")
local r_flux_e = regentlib.newsymbol(region(ispace(int3d), conserved), "r_flux_e")
local r_der    = regentlib.newsymbol(region(ispace(int3d), conserved), "r_der")

print("Periodicity in X: ", problem.periodic_x)
print("Periodicity in Y: ", problem.periodic_y)
print("Periodicity in Z: ", problem.periodic_z)

-- Midpoint-and-Node-Differencing tasks to compute flux derivatives
local ddx_MND_rho  = make_ddx_MND(r_flux, r_flux_e, "rho",  r_der, "rho",  problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, problem.periodic_x)
local ddx_MND_rhou = make_ddx_MND(r_flux, r_flux_e, "rhou", r_der, "rhou", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, problem.periodic_x)
local ddx_MND_rhov = make_ddx_MND(r_flux, r_flux_e, "rhov", r_der, "rhov", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, problem.periodic_x)
local ddx_MND_rhow = make_ddx_MND(r_flux, r_flux_e, "rhow", r_der, "rhow", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, problem.periodic_x)
local ddx_MND_rhoE = make_ddx_MND(r_flux, r_flux_e, "rhoE", r_der, "rhoE", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX, problem.periodic_x)

local ddy_MND_rho  = make_ddy_MND(r_flux, r_flux_e, "rho",  r_der, "rho",  problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, problem.periodic_y)
local ddy_MND_rhou = make_ddy_MND(r_flux, r_flux_e, "rhou", r_der, "rhou", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, problem.periodic_y)
local ddy_MND_rhov = make_ddy_MND(r_flux, r_flux_e, "rhov", r_der, "rhov", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, problem.periodic_y)
local ddy_MND_rhow = make_ddy_MND(r_flux, r_flux_e, "rhow", r_der, "rhow", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, problem.periodic_y)
local ddy_MND_rhoE = make_ddy_MND(r_flux, r_flux_e, "rhoE", r_der, "rhoE", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY, problem.periodic_y)

local ddz_MND_rho  = make_ddz_MND(r_flux, r_flux_e, "rho",  r_der, "rho",  problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, problem.periodic_z)
local ddz_MND_rhou = make_ddz_MND(r_flux, r_flux_e, "rhou", r_der, "rhou", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, problem.periodic_z)
local ddz_MND_rhov = make_ddz_MND(r_flux, r_flux_e, "rhov", r_der, "rhov", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, problem.periodic_z)
local ddz_MND_rhow = make_ddz_MND(r_flux, r_flux_e, "rhow", r_der, "rhow", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, problem.periodic_z)
local ddz_MND_rhoE = make_ddz_MND(r_flux, r_flux_e, "rhoE", r_der, "rhoE", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ, problem.periodic_z)

-- Node differencing tasks to compute flux derivatives
local ddx_rho  = make_ddx(r_flux, "rho",  r_der, "rho",  problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddx_rhou = make_ddx(r_flux, "rhou", r_der, "rhou", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddx_rhov = make_ddx(r_flux, "rhov", r_der, "rhov", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddx_rhow = make_ddx(r_flux, "rhow", r_der, "rhow", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddx_rhoE = make_ddx(r_flux, "rhoE", r_der, "rhoE", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)

local ddy_rho  = make_ddy(r_flux, "rho",  r_der, "rho",  problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddy_rhou = make_ddy(r_flux, "rhou", r_der, "rhou", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddy_rhov = make_ddy(r_flux, "rhov", r_der, "rhov", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddy_rhow = make_ddy(r_flux, "rhow", r_der, "rhow", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddy_rhoE = make_ddy(r_flux, "rhoE", r_der, "rhoE", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)

local ddz_rho  = make_ddz(r_flux, "rho",  r_der, "rho",  problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)
local ddz_rhou = make_ddz(r_flux, "rhou", r_der, "rhou", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)
local ddz_rhov = make_ddz(r_flux, "rhov", r_der, "rhov", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)
local ddz_rhow = make_ddz(r_flux, "rhow", r_der, "rhow", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)
local ddz_rhoE = make_ddz(r_flux, "rhoE", r_der, "rhoE", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local r_prim   = regentlib.newsymbol(region(ispace(int3d), primitive), "r_prim")
local r_prim_e = regentlib.newsymbol(region(ispace(int3d), primitive), "r_prim_e")
local r_der    = regentlib.newsymbol(region(ispace(int3d), tensor2),   "r_der")

local ddx_u   = make_ddx(r_prim, "u", r_der, "_11", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_u   = make_ddy(r_prim, "u", r_der, "_12", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddz_u   = make_ddz(r_prim, "u", r_der, "_13", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local ddx_v   = make_ddx(r_prim, "v", r_der, "_21", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_v   = make_ddy(r_prim, "v", r_der, "_22", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddz_v   = make_ddz(r_prim, "v", r_der, "_23", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local ddx_w   = make_ddx(r_prim, "w", r_der, "_31", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_w   = make_ddy(r_prim, "w", r_der, "_32", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddz_w   = make_ddz(r_prim, "w", r_der, "_33", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local r_aux   = regentlib.newsymbol(region(ispace(int3d), auxiliary), "r_aux")
local r_der   = regentlib.newsymbol(region(ispace(int3d), vect),    "r_der")

local ddx_T   = make_ddx(r_aux, "T", r_der, "_1", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_T   = make_ddy(r_aux, "T", r_der, "_2", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddz_T   = make_ddz(r_aux, "T", r_der, "_3", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

-----------------------------------------------------

task get_tauij( r_gradu : region(ispace(int3d), tensor2),
                r_tauij : region(ispace(int3d), tensor2symm),
                r_visc  : region(ispace(int3d), transport_coeffs) )
where
  reads(r_gradu, r_visc.{mu_s, mu_b}), writes(r_tauij)
do

  for i in r_tauij do
    var tau_dil = (r_visc[i].mu_b - (2./3.)*r_visc[i].mu_s) * (r_gradu[i]._11 + r_gradu[i]._22 + r_gradu[i]._33)

    r_tauij[i]._11 = r_visc[i].mu_s * ( r_gradu[i]._11 + r_gradu[i]._11 ) + tau_dil
    r_tauij[i]._12 = r_visc[i].mu_s * ( r_gradu[i]._12 + r_gradu[i]._21 )
    r_tauij[i]._13 = r_visc[i].mu_s * ( r_gradu[i]._13 + r_gradu[i]._31 )
    r_tauij[i]._22 = r_visc[i].mu_s * ( r_gradu[i]._22 + r_gradu[i]._22 ) + tau_dil
    r_tauij[i]._23 = r_visc[i].mu_s * ( r_gradu[i]._23 + r_gradu[i]._32 )
    r_tauij[i]._33 = r_visc[i].mu_s * ( r_gradu[i]._33 + r_gradu[i]._33 ) + tau_dil
  end

  return 1
end




__demand(__inline)
task get_q_x( r_aux_c : region(ispace(int3d), auxiliary),
              r_visc  : region(ispace(int3d), transport_coeffs),
              r_q     : region(ispace(int3d), vect),
              LU_x    : region(ispace(int3d), LU_struct) )
where
  reads(r_aux_c.T, r_visc.kappa, LU_x), reads writes(r_q.{_1})
do

  var token = ddx_T(r_aux_c, r_q, LU_x)

  for i in r_q do
    r_q[i]._1 = - r_visc[i].kappa * r_q[i]._1
  end

  return token
end

__demand(__inline)
task get_q_y( r_aux_c : region(ispace(int3d), auxiliary),
              r_visc  : region(ispace(int3d), transport_coeffs),
              r_q     : region(ispace(int3d), vect),
              LU_y    : region(ispace(int3d), LU_struct) )
where
  reads(r_aux_c.T, r_visc.kappa, LU_y), reads writes(r_q.{_2})
do

  var token = ddy_T(r_aux_c, r_q, LU_y)

  for i in r_q do
    r_q[i]._2 = - r_visc[i].kappa * r_q[i]._2
  end

  return token
end

__demand(__inline)
task get_q_z( r_aux_c : region(ispace(int3d), auxiliary),
              r_visc  : region(ispace(int3d), transport_coeffs),
              r_q     : region(ispace(int3d), vect),
              LU_z    : region(ispace(int3d), LU_struct) )
where
  reads(r_aux_c.T, r_visc.kappa, LU_z), reads writes(r_q.{_3})
do

  var token = ddz_T(r_aux_c, r_q, LU_z)

  for i in r_q do
    r_q[i]._3 = - r_visc[i].kappa * r_q[i]._3
  end

  return token
end




task add_xflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                           r_prim_l_x : region(ispace(int3d), primitive),
                           r_prim_r_x : region(ispace(int3d), primitive),
                           r_flux_c   : region(ispace(int3d), conserved),
                           r_flux_e_x : region(ispace(int3d), conserved),
                           r_fder_c_x : region(ispace(int3d), conserved),
                           r_rhs      : region(ispace(int3d), conserved),
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
                           LU_x       : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c, LU_x ),
  reads writes( r_prim_l_x, r_prim_r_x, r_flux_c, r_flux_e_x, r_fder_c_x, r_rhs),
  reads writes( alpha_l, beta_l, gamma_l, alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv )
do

  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1 - 2*interpolation.n_ghosts

  if (nx >= 8) then
    WCHR_interpolation_x( r_prim_c, r_prim_l_x, r_prim_r_x, alpha_l, beta_l, gamma_l,
                          alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv )
    positivity_enforcer_x( r_prim_c, r_prim_l_x, r_prim_r_x, interpolation.n_ghosts )
    HLLC_x( r_prim_l_x, r_prim_r_x, r_flux_e_x )
    get_xfluxes_r( r_prim_c, r_flux_c )
   
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

task add_viscous_xflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                                   r_aux_c    : region(ispace(int3d), auxiliary),
                                   r_visc     : region(ispace(int3d), transport_coeffs),
                                   r_tauij    : region(ispace(int3d), tensor2symm),
                                   r_q        : region(ispace(int3d), vect),
                                   r_flux_c   : region(ispace(int3d), conserved),
                                   r_fder_c_x : region(ispace(int3d), conserved),
                                   r_rhs      : region(ispace(int3d), conserved),
                                   LU_x       : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c, r_aux_c.T, r_visc.kappa, r_tauij.{_11, _12, _13}, LU_x ),
  reads writes( r_q._1, r_flux_c, r_fder_c_x, r_rhs)
do

  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1 - 2*interpolation.n_ghosts

  if (nx >= 8) then
    if viscous then
      get_q_x(r_aux_c, r_visc, r_q, LU_x)

      for i in r_flux_c do
        r_flux_c[i].rho  = 0.
        r_flux_c[i].rhou = r_tauij[i]._11
        r_flux_c[i].rhov = r_tauij[i]._12
        r_flux_c[i].rhow = r_tauij[i]._13
        r_flux_c[i].rhoE = r_tauij[i]._11 * r_prim_c[i].u + r_tauij[i]._12 * r_prim_c[i].v + r_tauij[i]._13 * r_prim_c[i].w - r_q[i]._1
      end

      ddx_rho ( r_flux_c, r_fder_c_x, LU_x )
      ddx_rhou( r_flux_c, r_fder_c_x, LU_x )
      ddx_rhov( r_flux_c, r_fder_c_x, LU_x )
      ddx_rhow( r_flux_c, r_fder_c_x, LU_x )
      ddx_rhoE( r_flux_c, r_fder_c_x, LU_x )

      for i in r_rhs do
        r_rhs[i].rho  += r_fder_c_x[i].rho
        r_rhs[i].rhou += r_fder_c_x[i].rhou
        r_rhs[i].rhov += r_fder_c_x[i].rhov
        r_rhs[i].rhow += r_fder_c_x[i].rhow
        r_rhs[i].rhoE += r_fder_c_x[i].rhoE
      end
    end

  end
end

task add_yflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                           r_prim_l_y : region(ispace(int3d), primitive),
                           r_prim_r_y : region(ispace(int3d), primitive),
                           r_flux_c   : region(ispace(int3d), conserved),
                           r_flux_e_y : region(ispace(int3d), conserved),
                           r_fder_c_y : region(ispace(int3d), conserved),
                           r_rhs      : region(ispace(int3d), conserved),
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
                           LU_y       : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c, LU_y ),
  reads writes( r_prim_l_y, r_prim_r_y, r_flux_c, r_flux_e_y, r_fder_c_y, r_rhs),
  reads writes( alpha_l, beta_l, gamma_l, alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv )
do

  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1 - 2*interpolation.n_ghosts

  if (ny >= 8) then
    WCHR_interpolation_y( r_prim_c, r_prim_l_y, r_prim_r_y, alpha_l, beta_l, gamma_l,
                          alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv )
    positivity_enforcer_y( r_prim_c, r_prim_l_y, r_prim_r_y, interpolation.n_ghosts )
    HLLC_y( r_prim_l_y, r_prim_r_y, r_flux_e_y )
    get_yfluxes_r( r_prim_c, r_flux_c )
    
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

task add_viscous_yflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                                   r_aux_c    : region(ispace(int3d), auxiliary),
                                   r_visc     : region(ispace(int3d), transport_coeffs),
                                   r_tauij    : region(ispace(int3d), tensor2symm),
                                   r_q        : region(ispace(int3d), vect),
                                   r_flux_c   : region(ispace(int3d), conserved),
                                   r_fder_c_y : region(ispace(int3d), conserved),
                                   r_rhs      : region(ispace(int3d), conserved),
                                   LU_y       : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c, r_aux_c.T, r_visc.kappa, r_tauij.{_12, _22, _23}, LU_y ),
  reads writes( r_q._2, r_flux_c, r_fder_c_y, r_rhs)
do

  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1 - 2*interpolation.n_ghosts

  if (ny >= 8) then
    if viscous then
      get_q_y(r_aux_c, r_visc, r_q, LU_y)

      for i in r_flux_c do
        r_flux_c[i].rho  = 0.
        r_flux_c[i].rhou = r_tauij[i]._12
        r_flux_c[i].rhov = r_tauij[i]._22
        r_flux_c[i].rhow = r_tauij[i]._23
        r_flux_c[i].rhoE = r_tauij[i]._12 * r_prim_c[i].u + r_tauij[i]._22 * r_prim_c[i].v + r_tauij[i]._23 * r_prim_c[i].w - r_q[i]._2
      end

      ddy_rho ( r_flux_c, r_fder_c_y, LU_y )
      ddy_rhou( r_flux_c, r_fder_c_y, LU_y )
      ddy_rhov( r_flux_c, r_fder_c_y, LU_y )
      ddy_rhow( r_flux_c, r_fder_c_y, LU_y )
      ddy_rhoE( r_flux_c, r_fder_c_y, LU_y )

      for i in r_rhs do
        r_rhs[i].rho  += r_fder_c_y[i].rho
        r_rhs[i].rhou += r_fder_c_y[i].rhou
        r_rhs[i].rhov += r_fder_c_y[i].rhov
        r_rhs[i].rhow += r_fder_c_y[i].rhow
        r_rhs[i].rhoE += r_fder_c_y[i].rhoE
      end
    end

  end
end

task add_zflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                           r_prim_l_z : region(ispace(int3d), primitive),
                           r_prim_r_z : region(ispace(int3d), primitive),
                           r_flux_c   : region(ispace(int3d), conserved),
                           r_flux_e_z : region(ispace(int3d), conserved),
                           r_fder_c_z : region(ispace(int3d), conserved),
                           r_rhs      : region(ispace(int3d), conserved),
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
                           LU_z       : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c, LU_z ),
  reads writes( r_prim_l_z, r_prim_r_z, r_flux_c, r_flux_e_z, r_fder_c_z, r_rhs),
  reads writes( alpha_l, beta_l, gamma_l, alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv )
do

  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1 - 2*interpolation.n_ghosts

  if (nz >= 8) then
    WCHR_interpolation_z( r_prim_c, r_prim_l_z, r_prim_r_z, alpha_l, beta_l, gamma_l,
                          alpha_r, beta_r, gamma_r, rho_avg, sos_avg, block_d, block_Uinv )
    positivity_enforcer_z( r_prim_c, r_prim_l_z, r_prim_r_z, interpolation.n_ghosts )
    HLLC_z( r_prim_l_z, r_prim_r_z, r_flux_e_z )
    get_zfluxes_r( r_prim_c, r_flux_c )

    ddz_MND_rho ( r_flux_c, r_flux_e_z, r_fder_c_z, LU_z )
    ddz_MND_rhou( r_flux_c, r_flux_e_z, r_fder_c_z, LU_z )
    ddz_MND_rhov( r_flux_c, r_flux_e_z, r_fder_c_z, LU_z )
    ddz_MND_rhow( r_flux_c, r_flux_e_z, r_fder_c_z, LU_z )
    ddz_MND_rhoE( r_flux_c, r_flux_e_z, r_fder_c_z, LU_z )
    -- var t_derivatives = c.legion_get_current_time_in_micros()

    for i in r_rhs do
      r_rhs[i].rho  -= r_fder_c_z[i].rho
      r_rhs[i].rhou -= r_fder_c_z[i].rhou
      r_rhs[i].rhov -= r_fder_c_z[i].rhov
      r_rhs[i].rhow -= r_fder_c_z[i].rhow
      r_rhs[i].rhoE -= r_fder_c_z[i].rhoE
    end
  end
end

task add_viscous_zflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                                   r_aux_c    : region(ispace(int3d), auxiliary),
                                   r_visc     : region(ispace(int3d), transport_coeffs),
                                   r_tauij    : region(ispace(int3d), tensor2symm),
                                   r_q        : region(ispace(int3d), vect),
                                   r_flux_c   : region(ispace(int3d), conserved),
                                   r_fder_c_z : region(ispace(int3d), conserved),
                                   r_rhs      : region(ispace(int3d), conserved),
                                   LU_z       : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c, r_aux_c.T, r_visc.kappa, r_tauij.{_13, _23, _33}, LU_z ),
  reads writes( r_q._3, r_flux_c, r_fder_c_z, r_rhs)
do

  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1 - 2*interpolation.n_ghosts

  if (nz >= 8) then
    if viscous then
      get_q_z(r_aux_c, r_visc, r_q, LU_z)

      for i in r_flux_c do
        r_flux_c[i].rho  = 0.
        r_flux_c[i].rhou = r_tauij[i]._13
        r_flux_c[i].rhov = r_tauij[i]._23
        r_flux_c[i].rhow = r_tauij[i]._33
        r_flux_c[i].rhoE = r_tauij[i]._13 * r_prim_c[i].u + r_tauij[i]._23 * r_prim_c[i].v + r_tauij[i]._33 * r_prim_c[i].w - r_q[i]._3
      end

      ddz_rho ( r_flux_c, r_fder_c_z, LU_z )
      ddz_rhou( r_flux_c, r_fder_c_z, LU_z )
      ddz_rhov( r_flux_c, r_fder_c_z, LU_z )
      ddz_rhow( r_flux_c, r_fder_c_z, LU_z )
      ddz_rhoE( r_flux_c, r_fder_c_z, LU_z )

      for i in r_rhs do
        r_rhs[i].rho  += r_fder_c_z[i].rho
        r_rhs[i].rhou += r_fder_c_z[i].rhou
        r_rhs[i].rhov += r_fder_c_z[i].rhov
        r_rhs[i].rhow += r_fder_c_z[i].rhow
        r_rhs[i].rhoE += r_fder_c_z[i].rhoE
      end
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
                                 r_gradu  : region(ispace(int3d), tensor2),
                                 LU_x     : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c.{u, v, w}, LU_x), reads writes (r_gradu.{_11, _21, _31})
do
  var nx = r_gradu.ispace.bounds.hi.x - r_gradu.ispace.bounds.lo.x + 1

  var token = 1
  if (nx < 8) then
    for i in r_gradu do
      r_gradu[i]._11 = 0.0
      r_gradu[i]._21 = 0.0
      r_gradu[i]._31 = 0.0
    end
    return token
  end

  token += ddx_u(r_prim_c, r_gradu, LU_x)
  token += ddx_v(r_prim_c, r_gradu, LU_x)
  token += ddx_w(r_prim_c, r_gradu, LU_x) 

  return token
end

task get_velocity_y_derivatives( r_prim_c : region(ispace(int3d), primitive),
                                 r_gradu  : region(ispace(int3d), tensor2),
                                 LU_y     : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c.{u, v, w}, LU_y), reads writes (r_gradu.{_12, _22, _32})
do
  var ny = r_gradu.ispace.bounds.hi.y - r_gradu.ispace.bounds.lo.y + 1

  var token = 1
  if (ny < 8) then
    for i in r_gradu do
      r_gradu[i]._12 = 0.0
      r_gradu[i]._22 = 0.0
      r_gradu[i]._32 = 0.0
    end
    return token
  end

  token += ddy_u(r_prim_c, r_gradu, LU_y)
  token += ddy_v(r_prim_c, r_gradu, LU_y)
  token += ddy_w(r_prim_c, r_gradu, LU_y) 

  return token
end

task get_velocity_z_derivatives( r_prim_c : region(ispace(int3d), primitive),
                                 r_gradu  : region(ispace(int3d), tensor2),
                                 LU_z     : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c.{u, v, w}, LU_z), reads writes (r_gradu.{_13, _23, _33})
do
  var nz = r_gradu.ispace.bounds.hi.z - r_gradu.ispace.bounds.lo.z + 1

  var token = 1
  if (nz < 8) then
    for i in r_gradu do
      r_gradu[i]._13 = 0.0
      r_gradu[i]._23 = 0.0
      r_gradu[i]._33 = 0.0
    end
    return token
  end

  token += ddz_u(r_prim_c, r_gradu, LU_z)
  token += ddz_v(r_prim_c, r_gradu, LU_z)
  token += ddz_w(r_prim_c, r_gradu, LU_z) 

  return token
end




task get_temperature_x_derivatives( r_aux_c : region(ispace(int3d), auxiliary),
                                    r_gradT : region(ispace(int3d), vect),
                                    LU_x    : region(ispace(int3d), LU_struct) )
where
  reads (r_aux_c.{T}, LU_x), reads writes (r_gradT.{_1})
do
  var nx = r_aux_c.ispace.bounds.hi.x - r_aux_c.ispace.bounds.lo.x + 1

  var token = 1
  if (nx < 8) then
    for i in r_gradT do
      r_gradT[i]._1 = 0.0
    end
    return token
  end

  token += ddx_T(r_aux_c, r_gradT, LU_x)

  return token
end

task get_temperature_y_derivatives( r_aux_c : region(ispace(int3d), auxiliary),
                                    r_gradT : region(ispace(int3d), vect),
                                    LU_y    : region(ispace(int3d), LU_struct) )
where
  reads (r_aux_c.{T}, LU_y), reads writes (r_gradT.{_2})
do
  var ny = r_aux_c.ispace.bounds.hi.y - r_aux_c.ispace.bounds.lo.y + 1

  var token = 1
  if (ny < 8) then
    for i in r_gradT do
      r_gradT[i]._2 = 0.0
    end
    return token
  end

  token += ddy_T(r_aux_c, r_gradT, LU_y)

  return token
end

task get_temperature_z_derivatives( r_aux_c : region(ispace(int3d), auxiliary),
                                    r_gradT : region(ispace(int3d), vect),
                                    LU_z    : region(ispace(int3d), LU_struct) )
where
  reads (r_aux_c.{T}, LU_z), reads writes (r_gradT.{_3})
do
  var nz = r_aux_c.ispace.bounds.hi.z - r_aux_c.ispace.bounds.lo.z + 1

  var token = 1
  if (nz < 8) then
    for i in r_gradT do
      r_gradT[i]._3 = 0.0
    end
    return token
  end

  token += ddz_T(r_aux_c, r_gradT, LU_z)

  return token
end

