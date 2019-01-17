import "regent"

require("fields")
require("derivatives")
require("SOE")

local interpolation = require("interpolation")
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

local d2dx2_u = make_d2dx2(r_prim, "u", r_der, "_11", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local d2dy2_u = make_d2dy2(r_prim, "u", r_der, "_12", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local d2dz2_u = make_d2dz2(r_prim, "u", r_der, "_13", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local d2dx2_v = make_d2dx2(r_prim, "v", r_der, "_21", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local d2dy2_v = make_d2dy2(r_prim, "v", r_der, "_22", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local d2dz2_v = make_d2dz2(r_prim, "v", r_der, "_23", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local d2dx2_w = make_d2dx2(r_prim, "w", r_der, "_31", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local d2dy2_w = make_d2dy2(r_prim, "w", r_der, "_32", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local d2dz2_w = make_d2dz2(r_prim, "w", r_der, "_33", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local r_aux   = regentlib.newsymbol(region(ispace(int3d), auxiliary), "r_aux")
local r_der   = regentlib.newsymbol(region(ispace(int3d), scalar),    "r_der")

local ddx_T   = make_ddx(r_aux, "T", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_T   = make_ddy(r_aux, "T", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddz_T   = make_ddz(r_aux, "T", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local d2dx2_T = make_d2dx2(r_aux, "T", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local d2dy2_T = make_d2dy2(r_aux, "T", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local d2dz2_T = make_d2dz2(r_aux, "T", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local r_prim   = regentlib.newsymbol(region(ispace(int3d), primitive), "r_prim")
local r_der   = regentlib.newsymbol(region(ispace(int3d), vect),    "r_der")

local ddx_rho_prim = make_ddx(r_prim, "rho", r_der, "_1", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_rho_prim = make_ddy(r_prim, "rho", r_der, "_2", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddz_rho_prim = make_ddz(r_prim, "rho", r_der, "_3", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local ddx_p   = make_ddx(r_prim, "p", r_der, "_1", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_p   = make_ddy(r_prim, "p", r_der, "_2", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddz_p   = make_ddz(r_prim, "p", r_der, "_3", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local r_visc  = regentlib.newsymbol(region(ispace(int3d), transport_coeffs), "r_visc")
local r_der   = regentlib.newsymbol(region(ispace(int3d), scalar),           "r_der")

local ddx_mu_s   = make_ddx(r_visc, "mu_s", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_mu_s   = make_ddy(r_visc, "mu_s", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddz_mu_s   = make_ddz(r_visc, "mu_s", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local ddx_mu_b   = make_ddx(r_visc, "mu_b", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_mu_b   = make_ddy(r_visc, "mu_b", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddz_mu_b   = make_ddz(r_visc, "mu_b", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local ddx_kappa  = make_ddx(r_visc, "kappa", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_kappa  = make_ddy(r_visc, "kappa", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDY)
local ddz_kappa  = make_ddz(r_visc, "kappa", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDZ)

local r_gradu = regentlib.newsymbol(region(ispace(int3d), tensor2), "r_gradu")
local r_der   = regentlib.newsymbol(region(ispace(int3d), scalar),  "r_der")

local ddx_dudy  = make_ddx(r_gradu, "_12", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddx_dudz  = make_ddx(r_gradu, "_13", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddx_dvdy  = make_ddx(r_gradu, "_22", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddx_dwdz  = make_ddx(r_gradu, "_33", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)

local ddy_dvdx  = make_ddy(r_gradu, "_21", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_dvdz  = make_ddy(r_gradu, "_23", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_dudx  = make_ddy(r_gradu, "_11", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddy_dwdz  = make_ddy(r_gradu, "_33", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)

local ddz_dwdx  = make_ddz(r_gradu, "_31", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddz_dwdy  = make_ddz(r_gradu, "_32", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddz_dudx  = make_ddz(r_gradu, "_11", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)
local ddz_dvdy  = make_ddz(r_gradu, "_22", r_der, "_", problem.NX, problem.NY, problem.NZ, problem.ONEBYDX)

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
task get_tauij_div_x( r_gradu     : region(ispace(int3d), tensor2),
                      r_grad2u    : region(ispace(int3d), tensor2),
                      r_visc      : region(ispace(int3d), transport_coeffs),
                      r_tauij_div : region(ispace(int3d), vect),
                      LU_x        : region(ispace(int3d), LU_struct) )
where
  reads(r_gradu.{_11, _12, _13, _21, _22, _31, _33}, r_grad2u.{_11, _21, _31}, r_visc.{mu_s, mu_b}, LU_x), writes(r_tauij_div)
do

  var nx = r_tauij_div.ispace.bounds.hi.x - r_tauij_div.ispace.bounds.lo.x + 1
  var ny = r_tauij_div.ispace.bounds.hi.y - r_tauij_div.ispace.bounds.lo.y + 1
  var nz = r_tauij_div.ispace.bounds.hi.z - r_tauij_div.ispace.bounds.lo.z + 1

  var bounds_c = r_tauij_div.ispace.bounds

  regentlib.assert(bounds_c.lo.x - interpolation.n_ghosts == 0, "Can only get_tauij_div_x in the X pencil")

  var r_ddx_dudy = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddx_dudz = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddx_dvdy = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddx_dwdz = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddx_mu_s = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddx_mu_b = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )

  var token = 0
  token += ddx_mu_s(r_visc, r_ddx_mu_s, LU_x)
  token += ddx_mu_b(r_visc, r_ddx_mu_b, LU_x)

  token += ddx_dudy(r_gradu, r_ddx_dudy, LU_x)
  token += ddx_dudz(r_gradu, r_ddx_dudz, LU_x)
  token += ddx_dvdy(r_gradu, r_ddx_dvdy, LU_x)
  token += ddx_dwdz(r_gradu, r_ddx_dwdz, LU_x)

  for i in r_tauij_div do
    var dilatation  = r_gradu[i]._11 + r_gradu[i]._22 + r_gradu[i]._33
    var tau_dil_div = (r_visc[i].mu_b - (2./3.)*r_visc[i].mu_s) * (r_grad2u[i]._11 + r_ddx_dvdy[i]._ + r_ddx_dwdz[i]._) 
                    + dilatation * ( r_ddx_mu_b[i]._ - (2./3.)*r_ddx_mu_s[i]._ )

    r_tauij_div[i]._1 = 2.*r_visc[i].mu_s*( r_grad2u[i]._11 ) + 2*r_gradu[i]._11*r_ddx_mu_s[i]._ + tau_dil_div
    r_tauij_div[i]._2 = r_visc[i].mu_s*( r_grad2u[i]._21 + r_ddx_dudy[i]._) + (r_gradu[i]._12 + r_gradu[i]._21)*r_ddx_mu_s[i]._
    r_tauij_div[i]._3 = r_visc[i].mu_s*( r_grad2u[i]._31 + r_ddx_dudz[i]._) + (r_gradu[i]._13 + r_gradu[i]._31)*r_ddx_mu_s[i]._
  end

  regentlib.c.legion_physical_region_destroy(__physical(r_ddx_dudy)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddx_dudz)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddx_dvdy)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddx_dwdz)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddx_mu_s)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddx_mu_b)[0])
  __delete(r_ddx_dudy)
  __delete(r_ddx_dudz)
  __delete(r_ddx_dvdy)
  __delete(r_ddx_dwdz)
  __delete(r_ddx_mu_s)
  __delete(r_ddx_mu_b)

  return 1
end



__demand(__inline)
task get_tauij_div_y( r_gradu     : region(ispace(int3d), tensor2),
                      r_grad2u    : region(ispace(int3d), tensor2),
                      r_visc      : region(ispace(int3d), transport_coeffs),
                      r_tauij_div : region(ispace(int3d), vect),
                      LU_y        : region(ispace(int3d), LU_struct) )
where
  reads(r_gradu.{_11, _12, _21, _22, _23, _32, _33}, r_grad2u.{_12, _22, _32}, r_visc.{mu_s, mu_b}, LU_y), writes(r_tauij_div)
do

  var nx = r_tauij_div.ispace.bounds.hi.x - r_tauij_div.ispace.bounds.lo.x + 1
  var ny = r_tauij_div.ispace.bounds.hi.y - r_tauij_div.ispace.bounds.lo.y + 1
  var nz = r_tauij_div.ispace.bounds.hi.z - r_tauij_div.ispace.bounds.lo.z + 1

  var bounds_c = r_tauij_div.ispace.bounds

  regentlib.assert(bounds_c.lo.y - interpolation.n_ghosts == 0, "Can only get_tauij_div_y in the Y pencil")

  var r_ddy_dvdx = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddy_dvdz = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddy_dudx = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddy_dwdz = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddy_mu_s = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddy_mu_b = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )

  var token = 0
  token += ddy_mu_s(r_visc, r_ddy_mu_s, LU_y)
  token += ddy_mu_b(r_visc, r_ddy_mu_b, LU_y)

  token += ddy_dvdx(r_gradu, r_ddy_dvdx, LU_y)
  token += ddy_dvdz(r_gradu, r_ddy_dvdz, LU_y)
  token += ddy_dudx(r_gradu, r_ddy_dudx, LU_y)
  token += ddy_dwdz(r_gradu, r_ddy_dwdz, LU_y)

  for i in r_tauij_div do
    var dilatation  = r_gradu[i]._11 + r_gradu[i]._22 + r_gradu[i]._33
    var tau_dil_div = (r_visc[i].mu_b - (2./3.)*r_visc[i].mu_s) * (r_ddy_dudx[i]._ + r_grad2u[i]._22 + r_ddy_dwdz[i]._) 
                    + dilatation * ( r_ddy_mu_b[i]._ - (2./3.)*r_ddy_mu_s[i]._ )

    r_tauij_div[i]._1 = r_visc[i].mu_s*( r_grad2u[i]._12 + r_ddy_dvdx[i]._) + (r_gradu[i]._12 + r_gradu[i]._21)*r_ddy_mu_s[i]._
    r_tauij_div[i]._2 = 2.*r_visc[i].mu_s*( r_grad2u[i]._22 ) + 2*r_gradu[i]._22*r_ddy_mu_s[i]._ + tau_dil_div
    r_tauij_div[i]._3 = r_visc[i].mu_s*( r_grad2u[i]._32 + r_ddy_dvdz[i]._) + (r_gradu[i]._23 + r_gradu[i]._32)*r_ddy_mu_s[i]._
  end

  regentlib.c.legion_physical_region_destroy(__physical(r_ddy_dvdx)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddy_dvdz)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddy_dudx)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddy_dwdz)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddy_mu_s)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddy_mu_b)[0])
  __delete(r_ddy_dvdx)
  __delete(r_ddy_dvdz)
  __delete(r_ddy_dudx)
  __delete(r_ddy_dwdz)
  __delete(r_ddy_mu_s)
  __delete(r_ddy_mu_b)

  return 1
end



__demand(__inline)
task get_tauij_div_z( r_gradu     : region(ispace(int3d), tensor2),
                      r_grad2u    : region(ispace(int3d), tensor2),
                      r_visc      : region(ispace(int3d), transport_coeffs),
                      r_tauij_div : region(ispace(int3d), vect),
                      LU_z        : region(ispace(int3d), LU_struct) )
where
  reads(r_gradu.{_11, _13, _22, _23, _31, _32, _33}, r_grad2u.{_13, _23, _33}, r_visc.{mu_s, mu_b}, LU_z), writes(r_tauij_div)
do

  var nx = r_tauij_div.ispace.bounds.hi.x - r_tauij_div.ispace.bounds.lo.x + 1
  var ny = r_tauij_div.ispace.bounds.hi.y - r_tauij_div.ispace.bounds.lo.y + 1
  var nz = r_tauij_div.ispace.bounds.hi.z - r_tauij_div.ispace.bounds.lo.z + 1

  var bounds_c = r_tauij_div.ispace.bounds

  regentlib.assert(bounds_c.lo.z - interpolation.n_ghosts == 0, "Can only get_tauij_div_z in the Z pencil")

  var r_ddz_dwdx = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddz_dwdy = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddz_dudx = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddz_dvdy = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddz_mu_s = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddz_mu_b = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )

  var token = 0
  token += ddz_mu_s(r_visc, r_ddz_mu_s, LU_z)
  token += ddz_mu_b(r_visc, r_ddz_mu_b, LU_z)

  token += ddz_dwdx(r_gradu, r_ddz_dwdx, LU_z)
  token += ddz_dwdy(r_gradu, r_ddz_dwdy, LU_z)
  token += ddz_dudx(r_gradu, r_ddz_dudx, LU_z)
  token += ddz_dvdy(r_gradu, r_ddz_dvdy, LU_z)

  for i in r_tauij_div do
    var dilatation  = r_gradu[i]._11 + r_gradu[i]._22 + r_gradu[i]._33
    var tau_dil_div = (r_visc[i].mu_b - (2./3.)*r_visc[i].mu_s) * (r_ddz_dudx[i]._ + r_ddz_dvdy[i]._ + r_grad2u[i]._33) 
                    + dilatation * ( r_ddz_mu_b[i]._ - (2./3.)*r_ddz_mu_s[i]._ )

    r_tauij_div[i]._1 = r_visc[i].mu_s*( r_grad2u[i]._13 + r_ddz_dwdx[i]._) + (r_gradu[i]._13 + r_gradu[i]._31)*r_ddz_mu_s[i]._
    r_tauij_div[i]._2 = r_visc[i].mu_s*( r_grad2u[i]._23 + r_ddz_dwdy[i]._) + (r_gradu[i]._23 + r_gradu[i]._32)*r_ddz_mu_s[i]._
    r_tauij_div[i]._3 = 2.*r_visc[i].mu_s*( r_grad2u[i]._33 ) + 2*r_gradu[i]._33*r_ddz_mu_s[i]._ + tau_dil_div
  end

  regentlib.c.legion_physical_region_destroy(__physical(r_ddz_dwdx)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddz_dwdy)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddz_dudx)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddz_dvdy)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddz_mu_s)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddz_mu_b)[0])
  __delete(r_ddz_dwdx)
  __delete(r_ddz_dwdy)
  __delete(r_ddz_dudx)
  __delete(r_ddz_dvdy)
  __delete(r_ddz_mu_s)
  __delete(r_ddz_mu_b)

  return 1
end



__demand(__inline)
task get_q_x( r_aux_c : region(ispace(int3d), auxiliary),
              r_visc  : region(ispace(int3d), transport_coeffs),
              r_q     : region(ispace(int3d), vect),
              LU_x    : region(ispace(int3d), LU_struct) )
where
  reads(r_aux_c.T, r_visc.kappa, LU_x), writes(r_q.{_1})
do

  var nx = r_q.ispace.bounds.hi.x - r_q.ispace.bounds.lo.x + 1
  var ny = r_q.ispace.bounds.hi.y - r_q.ispace.bounds.lo.y + 1
  var nz = r_q.ispace.bounds.hi.z - r_q.ispace.bounds.lo.z + 1

  var bounds_c = r_q.ispace.bounds

  regentlib.assert(bounds_c.lo.x - interpolation.n_ghosts == 0, "Can only get_q_x in the X pencil")

  var r_ddx_T = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )

  var token = ddx_T(r_aux_c, r_ddx_T, LU_x)

  for i in r_q do
    r_q[i]._1 = - r_visc[i].kappa * r_ddx_T[i]._
  end

  regentlib.c.legion_physical_region_destroy(__physical(r_ddx_T)[0])
  __delete(r_ddx_T)

  return token
end



__demand(__inline)
task get_q_y( r_aux_c : region(ispace(int3d), auxiliary),
              r_visc  : region(ispace(int3d), transport_coeffs),
              r_q     : region(ispace(int3d), vect),
              LU_y    : region(ispace(int3d), LU_struct) )
where
  reads(r_aux_c.T, r_visc.kappa, LU_y), writes(r_q.{_2})
do

  var nx = r_q.ispace.bounds.hi.x - r_q.ispace.bounds.lo.x + 1
  var ny = r_q.ispace.bounds.hi.y - r_q.ispace.bounds.lo.y + 1
  var nz = r_q.ispace.bounds.hi.z - r_q.ispace.bounds.lo.z + 1

  var bounds_c = r_q.ispace.bounds

  regentlib.assert(bounds_c.lo.y - interpolation.n_ghosts == 0, "Can only get_q_y in the Y pencil")

  var r_ddy_T = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )

  var token = ddy_T(r_aux_c, r_ddy_T, LU_y)

  for i in r_q do
    r_q[i]._2 = - r_visc[i].kappa * r_ddy_T[i]._
  end

  regentlib.c.legion_physical_region_destroy(__physical(r_ddy_T)[0])
  __delete(r_ddy_T)

  return token
end



__demand(__inline)
task get_q_z( r_aux_c : region(ispace(int3d), auxiliary),
              r_visc  : region(ispace(int3d), transport_coeffs),
              r_q     : region(ispace(int3d), vect),
              LU_z    : region(ispace(int3d), LU_struct) )
where
  reads(r_aux_c.T, r_visc.kappa, LU_z), writes(r_q.{_3})
do

  var nx = r_q.ispace.bounds.hi.x - r_q.ispace.bounds.lo.x + 1
  var ny = r_q.ispace.bounds.hi.y - r_q.ispace.bounds.lo.y + 1
  var nz = r_q.ispace.bounds.hi.z - r_q.ispace.bounds.lo.z + 1

  var bounds_c = r_q.ispace.bounds

  regentlib.assert(bounds_c.lo.z - interpolation.n_ghosts == 0, "Can only get_q_z in the Z pencil")

  var r_ddz_T = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )

  var token = ddz_T(r_aux_c, r_ddz_T, LU_z)

  for i in r_q do
    r_q[i]._3 = - r_visc[i].kappa * r_ddz_T[i]._
  end

  regentlib.c.legion_physical_region_destroy(__physical(r_ddz_T)[0])
  __delete(r_ddz_T)

  return token
end



__demand(__inline)
task get_q_div_x( r_aux_c : region(ispace(int3d), auxiliary),
                  r_visc  : region(ispace(int3d), transport_coeffs),
                  r_q_div : region(ispace(int3d), scalar),
                  LU_x    : region(ispace(int3d), LU_struct),
                  LU2_x   : region(ispace(int3d), LU_struct) )
where
  reads(r_aux_c.T, r_visc.kappa, LU_x, LU2_x), writes(r_q_div)
do

  var nx = r_q_div.ispace.bounds.hi.x - r_q_div.ispace.bounds.lo.x + 1
  var ny = r_q_div.ispace.bounds.hi.y - r_q_div.ispace.bounds.lo.y + 1
  var nz = r_q_div.ispace.bounds.hi.z - r_q_div.ispace.bounds.lo.z + 1

  var bounds_c = r_q_div.ispace.bounds

  regentlib.assert(bounds_c.lo.x - interpolation.n_ghosts == 0, "Can only get_q_div_x in the X pencil")

  var r_ddx_T     = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_d2dx2_T   = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddx_kappa = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )

  var token = 0
  token += ddx_T(r_aux_c, r_ddx_T, LU_x)
  token += d2dx2_T(r_aux_c, r_d2dx2_T, LU2_x)
  token += ddx_kappa(r_visc, r_ddx_kappa, LU_x)

  for i in r_q_div do
    r_q_div[i]._ = - (r_ddx_kappa[i]._ * r_ddx_T[i]._ + r_visc[i].kappa * r_d2dx2_T[i]._)
  end

  regentlib.c.legion_physical_region_destroy(__physical(r_ddx_T)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_d2dx2_T)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddx_kappa)[0])
  __delete(r_ddx_T)
  __delete(r_d2dx2_T)
  __delete(r_ddx_kappa)

  return token
end



__demand(__inline)
task get_q_div_y( r_aux_c : region(ispace(int3d), auxiliary),
                  r_visc  : region(ispace(int3d), transport_coeffs),
                  r_q_div : region(ispace(int3d), scalar),
                  LU_y    : region(ispace(int3d), LU_struct),
                  LU2_y   : region(ispace(int3d), LU_struct) )
where
  reads(r_aux_c.T, r_visc.kappa, LU_y, LU2_y), writes(r_q_div)
do

  var nx = r_q_div.ispace.bounds.hi.x - r_q_div.ispace.bounds.lo.x + 1
  var ny = r_q_div.ispace.bounds.hi.y - r_q_div.ispace.bounds.lo.y + 1
  var nz = r_q_div.ispace.bounds.hi.z - r_q_div.ispace.bounds.lo.z + 1

  var bounds_c = r_q_div.ispace.bounds

  regentlib.assert(bounds_c.lo.y - interpolation.n_ghosts == 0, "Can only get_q_div_y in the Y pencil")

  var r_ddy_T     = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_d2dy2_T   = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddy_kappa = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )

  var token = 0
  token += ddy_T(r_aux_c, r_ddy_T, LU_y)
  token += d2dy2_T(r_aux_c, r_d2dy2_T, LU2_y)
  token += ddy_kappa(r_visc, r_ddy_kappa, LU_y)

  for i in r_q_div do
    r_q_div[i]._ = - (r_ddy_kappa[i]._ * r_ddy_T[i]._ + r_visc[i].kappa * r_d2dy2_T[i]._)
  end

  regentlib.c.legion_physical_region_destroy(__physical(r_ddy_T)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_d2dy2_T)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddy_kappa)[0])
  __delete(r_ddy_T)
  __delete(r_d2dy2_T)
  __delete(r_ddy_kappa)

  return token
end



__demand(__inline)
task get_q_div_z( r_aux_c : region(ispace(int3d), auxiliary),
                  r_visc  : region(ispace(int3d), transport_coeffs),
                  r_q_div : region(ispace(int3d), scalar),
                  LU_z    : region(ispace(int3d), LU_struct),
                  LU2_z   : region(ispace(int3d), LU_struct) )
where
  reads(r_aux_c.T, r_visc.kappa, LU_z, LU2_z), writes(r_q_div)
do

  var nx = r_q_div.ispace.bounds.hi.x - r_q_div.ispace.bounds.lo.x + 1
  var ny = r_q_div.ispace.bounds.hi.y - r_q_div.ispace.bounds.lo.y + 1
  var nz = r_q_div.ispace.bounds.hi.z - r_q_div.ispace.bounds.lo.z + 1

  var bounds_c = r_q_div.ispace.bounds

  regentlib.assert(bounds_c.lo.z - interpolation.n_ghosts == 0, "Can only get_q_div_z in the Z pencil")

  var r_ddz_T     = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_d2dz2_T   = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
  var r_ddz_kappa = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )

  var token = 0
  token += ddz_T(r_aux_c, r_ddz_T, LU_z)
  token += d2dz2_T(r_aux_c, r_d2dz2_T, LU2_z)
  token += ddz_kappa(r_visc, r_ddz_kappa, LU_z)

  for i in r_q_div do
    r_q_div[i]._ = - (r_ddz_kappa[i]._ * r_ddz_T[i]._ + r_visc[i].kappa * r_d2dz2_T[i]._)
  end

  regentlib.c.legion_physical_region_destroy(__physical(r_ddz_T)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_d2dz2_T)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_ddz_kappa)[0])
  __delete(r_ddz_T)
  __delete(r_d2dz2_T)
  __delete(r_ddz_kappa)

  return token
end



task add_xflux_der_to_rhs( r_prim_c    : region(ispace(int3d), primitive),
                           r_prim_c_wo : region(ispace(int3d), primitive),
                           r_rhs       : region(ispace(int3d), conserved),
                           LU_x        : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c, r_prim_c_wo, LU_x ),
  reads writes( r_rhs )
do

  var bounds_c = r_prim_c.ispace.bounds

  var Nx   = problem.NX
  var Nx_g = Nx + 2*interpolation.n_ghosts

  regentlib.assert(bounds_c.lo.x == 0,      "Can only add X flux derivative in the X pencil")
  regentlib.assert(bounds_c.hi.x == Nx_g-1, "Can only add X flux derivative in the X pencil")

  var bounds_der_lo = {bounds_c.lo.x + interpolation.n_ghosts, bounds_c.lo.y, bounds_c.lo.z}

  var Nx_e = Nx + 1

  var ny = bounds_c.hi.y - bounds_c.lo.y + 1
  var nz = bounds_c.hi.z - bounds_c.lo.z + 1

  var r_prim_l_x = region( ispace(int3d, {Nx_e, ny, nz}, bounds_c.lo), primitive )
  var r_prim_r_x = region( ispace(int3d, {Nx_e, ny, nz}, bounds_c.lo), primitive )

  var r_flux_c   = region( ispace(int3d, {Nx_g, ny, nz}, bounds_c.lo),   conserved )
  var r_flux_e_x = region( ispace(int3d, {Nx_e, ny, nz}, bounds_c.lo),   conserved )
  var r_fder_c_x = region( ispace(int3d, {Nx,   ny, nz}, bounds_der_lo), conserved )

  if (Nx >= 8) then
    WCHR_interpolation_x( r_prim_c, r_prim_l_x, r_prim_r_x )
    positivity_enforcer_x( r_prim_c, r_prim_l_x, r_prim_r_x, interpolation.n_ghosts )

    if (problem.Riemann_solver == "HLL") then
      HLL_x( r_prim_l_x, r_prim_r_x, r_flux_e_x )
    elseif (problem.Riemann_solver == "HLLC") then
      HLLC_x( r_prim_l_x, r_prim_r_x, r_flux_e_x )
    else
      var r_gradu_l_x       = region( ispace(int3d, {Nx_e, ny, nz}, bounds_c.lo), tensor2 )
      var r_gradu_r_x       = region( ispace(int3d, {Nx_e, ny, nz}, bounds_c.lo), tensor2 )
      var r_theta_avg_x     = region( ispace(int3d, {Nx_e, ny, nz}, bounds_c.lo), double )
      var r_omega_mag_avg_x = region( ispace(int3d, {Nx_e, ny, nz}, bounds_c.lo), double )

      get_velocity_derivatives_x(r_prim_c_wo, r_gradu_l_x, r_gradu_r_x, interpolation.n_ghosts)
      compute_theta_avg(r_gradu_l_x, r_gradu_r_x, r_theta_avg_x)
      compute_omega_mag_avg(r_gradu_l_x, r_gradu_r_x, r_omega_mag_avg_x)

      HLLC_HLL_x( r_prim_l_x, r_prim_r_x, r_theta_avg_x, r_omega_mag_avg_x, r_flux_e_x )

      regentlib.c.legion_physical_region_destroy(__physical(r_gradu_l_x)[0])
      regentlib.c.legion_physical_region_destroy(__physical(r_gradu_r_x)[0])
      regentlib.c.legion_physical_region_destroy(__physical(r_theta_avg_x)[0])
      regentlib.c.legion_physical_region_destroy(__physical(r_omega_mag_avg_x)[0])

      __delete(r_gradu_l_x)
      __delete(r_gradu_r_x)
      __delete(r_theta_avg_x)
      __delete(r_omega_mag_avg_x)
    end

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

  regentlib.c.legion_physical_region_destroy(__physical(r_prim_l_x)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_prim_r_x)[0])

  regentlib.c.legion_physical_region_destroy(__physical(r_flux_c)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_flux_e_x)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_fder_c_x)[0])

  __delete(r_prim_l_x)
  __delete(r_prim_r_x)

  __delete(r_flux_c)
  __delete(r_flux_e_x)
  __delete(r_fder_c_x)
end



task add_viscous_xflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                                   r_aux_c    : region(ispace(int3d), auxiliary),
                                   r_visc     : region(ispace(int3d), transport_coeffs),
                                   r_tauij    : region(ispace(int3d), tensor2symm),
                                   r_rhs      : region(ispace(int3d), conserved),
                                   LU_x       : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c, r_aux_c.T, r_visc.kappa, r_tauij.{_11, _12, _13}, LU_x ),
  reads writes( r_rhs )
do

  var bounds_c = r_prim_c.ispace.bounds

  var Nx   = problem.NX

  regentlib.assert(bounds_c.lo.x == interpolation.n_ghosts,          "Can only add X flux derivative in the X pencil")
  regentlib.assert(bounds_c.hi.x == Nx + interpolation.n_ghosts - 1, "Can only add X flux derivative in the X pencil")

  var ny = bounds_c.hi.y - bounds_c.lo.y + 1
  var nz = bounds_c.hi.z - bounds_c.lo.z + 1

  var r_q = region( ispace(int3d, {Nx, ny, nz}, bounds_c.lo), vect )

  var r_flux_c   = region( ispace(int3d, {Nx, ny, nz}, bounds_c.lo), conserved )
  var r_fder_c_x = region( ispace(int3d, {Nx, ny, nz}, bounds_c.lo), conserved )

  if (Nx >= 8) then
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

  regentlib.c.legion_physical_region_destroy(__physical(r_q)[0])

  regentlib.c.legion_physical_region_destroy(__physical(r_flux_c)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_fder_c_x)[0])

  __delete(r_q)

  __delete(r_flux_c)
  __delete(r_fder_c_x)
end



task add_nonconservative_viscous_xflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                                                   r_aux_c    : region(ispace(int3d), auxiliary),
                                                   r_visc     : region(ispace(int3d), transport_coeffs),
                                                   r_gradu    : region(ispace(int3d), tensor2),
                                                   r_grad2u   : region(ispace(int3d), tensor2),
                                                   r_tauij    : region(ispace(int3d), tensor2symm),
                                                   r_rhs      : region(ispace(int3d), conserved),
                                                   LU_x       : region(ispace(int3d), LU_struct),
                                                   LU2_x      : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c.{u, v, w}, r_aux_c.T, r_visc, r_tauij.{_11, _12, _13}, LU_x, LU2_x ),
  reads( r_gradu.{_11, _12, _13, _21, _22, _31, _33}, r_grad2u.{_11, _21, _31} ),
  reads writes( r_rhs.{rhou, rhov, rhow, rhoE} )
do

  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1
  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1
  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1

  var bounds_c = r_prim_c.ispace.bounds

  regentlib.assert(bounds_c.lo.x - interpolation.n_ghosts == 0, "Can only add_nonconservative_viscous_xflux_der_to_rhs in the X pencil")

  if (nx >= 8) then
    if viscous then
      var r_q_div     = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
      var r_tauij_div = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), vect )

      get_tauij_div_x(r_gradu, r_grad2u, r_visc, r_tauij_div, LU_x)
      get_q_div_x(r_aux_c, r_visc, r_q_div, LU_x, LU2_x)

      for i in r_rhs do
        var viscpower = r_tauij[i]._11 * r_gradu[i]._11 + r_tauij[i]._12 * r_gradu[i]._21 + r_tauij[i]._13 * r_gradu[i]._31
                      + r_prim_c[i].u * r_tauij_div[i]._1 + r_prim_c[i].v * r_tauij_div[i]._2+ r_prim_c[i].w * r_tauij_div[i]._3

        r_rhs[i].rhou += r_tauij_div[i]._1
        r_rhs[i].rhov += r_tauij_div[i]._2
        r_rhs[i].rhow += r_tauij_div[i]._3
        r_rhs[i].rhoE += - r_q_div[i]._ + viscpower
      end

      regentlib.c.legion_physical_region_destroy(__physical(r_q_div)[0])
      regentlib.c.legion_physical_region_destroy(__physical(r_tauij_div)[0])
      __delete(r_q_div)
      __delete(r_tauij_div)
    end

  end
end



task add_yflux_der_to_rhs( r_prim_c    : region(ispace(int3d), primitive),
                           r_prim_c_wo : region(ispace(int3d), primitive),
                           r_rhs       : region(ispace(int3d), conserved),
                           LU_y        : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c, r_prim_c_wo, LU_y ),
  reads writes( r_rhs )
do

  var bounds_c = r_prim_c.ispace.bounds

  var Ny   = problem.NY
  var Ny_g = Ny + 2*interpolation.n_ghosts

  regentlib.assert(bounds_c.lo.y == 0,      "Can only add Y flux derivative in the Y pencil")
  regentlib.assert(bounds_c.hi.y == Ny_g-1, "Can only add Y flux derivative in the Y pencil")

  var bounds_der_lo = {bounds_c.lo.x, bounds_c.lo.y + interpolation.n_ghosts, bounds_c.lo.z}

  var Ny_e = Ny + 1

  var nx = bounds_c.hi.x - bounds_c.lo.x + 1
  var nz = bounds_c.hi.z - bounds_c.lo.z + 1

  var r_prim_l_y = region( ispace(int3d, {nx, Ny_e, nz}, bounds_c.lo), primitive )
  var r_prim_r_y = region( ispace(int3d, {nx, Ny_e, nz}, bounds_c.lo), primitive )

  var r_flux_c   = region( ispace(int3d, {nx, Ny_g, nz}, bounds_c.lo),   conserved )
  var r_flux_e_y = region( ispace(int3d, {nx, Ny_e, nz}, bounds_c.lo),   conserved )
  var r_fder_c_y = region( ispace(int3d, {nx, Ny,   nz}, bounds_der_lo), conserved )

  if (Ny >= 8) then
    WCHR_interpolation_y( r_prim_c, r_prim_l_y, r_prim_r_y )
    positivity_enforcer_y( r_prim_c, r_prim_l_y, r_prim_r_y, interpolation.n_ghosts )

    if (problem.Riemann_solver == "HLL") then
      HLL_y( r_prim_l_y, r_prim_r_y, r_flux_e_y )
    elseif (problem.Riemann_solver == "HLLC") then
      HLLC_y( r_prim_l_y, r_prim_r_y, r_flux_e_y )
    else
      var r_gradu_l_y       = region( ispace(int3d, {nx, Ny_e, nz}, bounds_c.lo), tensor2 )
      var r_gradu_r_y       = region( ispace(int3d, {nx, Ny_e, nz}, bounds_c.lo), tensor2 )
      var r_theta_avg_y     = region( ispace(int3d, {nx, Ny_e, nz}, bounds_c.lo), double )
      var r_omega_mag_avg_y = region( ispace(int3d, {nx, Ny_e, nz}, bounds_c.lo), double )

      get_velocity_derivatives_y(r_prim_c_wo, r_gradu_l_y, r_gradu_r_y, interpolation.n_ghosts)
      compute_theta_avg(r_gradu_l_y, r_gradu_r_y, r_theta_avg_y)
      compute_omega_mag_avg(r_gradu_l_y, r_gradu_r_y, r_omega_mag_avg_y)

      HLLC_HLL_y( r_prim_l_y, r_prim_r_y, r_theta_avg_y, r_omega_mag_avg_y, r_flux_e_y )

      regentlib.c.legion_physical_region_destroy(__physical(r_gradu_l_y)[0])
      regentlib.c.legion_physical_region_destroy(__physical(r_gradu_r_y)[0])
      regentlib.c.legion_physical_region_destroy(__physical(r_theta_avg_y)[0])
      regentlib.c.legion_physical_region_destroy(__physical(r_omega_mag_avg_y)[0])

      __delete(r_gradu_l_y)
      __delete(r_gradu_r_y)
      __delete(r_theta_avg_y)
      __delete(r_omega_mag_avg_y)
    end

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

  regentlib.c.legion_physical_region_destroy(__physical(r_prim_l_y)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_prim_r_y)[0])

  regentlib.c.legion_physical_region_destroy(__physical(r_flux_c)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_flux_e_y)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_fder_c_y)[0])

  __delete(r_prim_l_y)
  __delete(r_prim_r_y)

  __delete(r_flux_c)
  __delete(r_flux_e_y)
  __delete(r_fder_c_y)
end



task add_viscous_yflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                                   r_aux_c    : region(ispace(int3d), auxiliary),
                                   r_visc     : region(ispace(int3d), transport_coeffs),
                                   r_tauij    : region(ispace(int3d), tensor2symm),
                                   r_rhs      : region(ispace(int3d), conserved),
                                   LU_y       : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c, r_aux_c.T, r_visc.kappa, r_tauij.{_12, _22, _23}, LU_y ),
  reads writes( r_rhs )
do

  var bounds_c = r_prim_c.ispace.bounds

  var Ny   = problem.NY

  regentlib.assert(bounds_c.lo.y == interpolation.n_ghosts,          "Can only add Y flux derivative in the Y pencil")
  regentlib.assert(bounds_c.hi.y == Ny + interpolation.n_ghosts - 1, "Can only add Y flux derivative in the Y pencil")

  var nx = bounds_c.hi.x - bounds_c.lo.x + 1
  var nz = bounds_c.hi.z - bounds_c.lo.z + 1

  var r_q = region( ispace(int3d, {nx, Ny, nz}, bounds_c.lo), vect )

  var r_flux_c   = region( ispace(int3d, {nx, Ny, nz}, bounds_c.lo), conserved )
  var r_fder_c_y = region( ispace(int3d, {nx, Ny, nz}, bounds_c.lo), conserved )

  if (Ny >= 8) then
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

  regentlib.c.legion_physical_region_destroy(__physical(r_q)[0])

  regentlib.c.legion_physical_region_destroy(__physical(r_flux_c)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_fder_c_y)[0])

  __delete(r_q)

  __delete(r_flux_c)
  __delete(r_fder_c_y)
end



task add_nonconservative_viscous_yflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                                                   r_aux_c    : region(ispace(int3d), auxiliary),
                                                   r_visc     : region(ispace(int3d), transport_coeffs),
                                                   r_gradu    : region(ispace(int3d), tensor2),
                                                   r_grad2u   : region(ispace(int3d), tensor2),
                                                   r_tauij    : region(ispace(int3d), tensor2symm),
                                                   r_rhs      : region(ispace(int3d), conserved),
                                                   LU_y       : region(ispace(int3d), LU_struct),
                                                   LU2_y      : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c.{u, v, w}, r_aux_c.T, r_visc, r_tauij.{_12, _22, _23}, LU_y, LU2_y ),
  reads( r_gradu.{_11, _12, _21, _22, _23, _32, _33}, r_grad2u.{_12, _22, _32} ),
  reads writes( r_rhs.{rhou, rhov, rhow, rhoE} )
do

  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1
  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1
  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1

  var bounds_c = r_prim_c.ispace.bounds

  regentlib.assert(bounds_c.lo.y - interpolation.n_ghosts == 0, "Can only add_nonconservative_viscous_yflux_der_to_rhs in the Y pencil")

  if (ny >= 8) then
    if viscous then
      var r_q_div     = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
      var r_tauij_div = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), vect )

      get_tauij_div_y(r_gradu, r_grad2u, r_visc, r_tauij_div, LU_y)
      get_q_div_y(r_aux_c, r_visc, r_q_div, LU_y, LU2_y)

      for i in r_rhs do
        var viscpower = r_tauij[i]._12 * r_gradu[i]._12 + r_tauij[i]._22 * r_gradu[i]._22 + r_tauij[i]._23 * r_gradu[i]._32
                      + r_prim_c[i].u * r_tauij_div[i]._1 + r_prim_c[i].v * r_tauij_div[i]._2+ r_prim_c[i].w * r_tauij_div[i]._3

        r_rhs[i].rhou += r_tauij_div[i]._1
        r_rhs[i].rhov += r_tauij_div[i]._2
        r_rhs[i].rhow += r_tauij_div[i]._3
        r_rhs[i].rhoE += - r_q_div[i]._ + viscpower
      end

      regentlib.c.legion_physical_region_destroy(__physical(r_q_div)[0])
      regentlib.c.legion_physical_region_destroy(__physical(r_tauij_div)[0])
      __delete(r_q_div)
      __delete(r_tauij_div)
    end

  end
end



task add_zflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                           r_prim_c_wo : region(ispace(int3d), primitive),
                           r_rhs      : region(ispace(int3d), conserved),
                           LU_z       : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c, r_prim_c_wo, LU_z ),
  reads writes( r_rhs )
do

  var bounds_c = r_prim_c.ispace.bounds

  var Nz   = problem.NZ
  var Nz_g = Nz + 2*interpolation.n_ghosts

  regentlib.assert(bounds_c.lo.z == 0,      "Can only add Z flux derivative in the Z pencil")
  regentlib.assert(bounds_c.hi.z == Nz_g-1, "Can only add Z flux derivative in the Z pencil")

  var bounds_der_lo = {bounds_c.lo.x, bounds_c.lo.y, bounds_c.lo.z + interpolation.n_ghosts}

  var Nz_e = Nz + 1

  var nx = bounds_c.hi.x - bounds_c.lo.x + 1
  var ny = bounds_c.hi.y - bounds_c.lo.y + 1

  var r_prim_l_z = region( ispace(int3d, {nx, ny, Nz_e}, bounds_c.lo), primitive )
  var r_prim_r_z = region( ispace(int3d, {nx, ny, Nz_e}, bounds_c.lo), primitive )

  var r_flux_c   = region( ispace(int3d, {nx, ny, Nz_g}, bounds_c.lo),   conserved )
  var r_flux_e_z = region( ispace(int3d, {nx, ny, Nz_e}, bounds_c.lo),   conserved )
  var r_fder_c_z = region( ispace(int3d, {nx, ny, Nz},   bounds_der_lo), conserved )

  if (Nz >= 8) then
    WCHR_interpolation_z( r_prim_c, r_prim_l_z, r_prim_r_z )
    positivity_enforcer_z( r_prim_c, r_prim_l_z, r_prim_r_z, interpolation.n_ghosts )

    if (problem.Riemann_solver == "HLL") then
      HLL_z( r_prim_l_z, r_prim_r_z, r_flux_e_z )
    elseif (problem.Riemann_solver == "HLLC") then
      HLLC_z( r_prim_l_z, r_prim_r_z, r_flux_e_z )
    else
      var r_gradu_l_z       = region( ispace(int3d, {nx, ny, Nz_e}, bounds_c.lo), tensor2 )
      var r_gradu_r_z       = region( ispace(int3d, {nx, ny, Nz_e}, bounds_c.lo), tensor2 )
      var r_theta_avg_z     = region( ispace(int3d, {nx, ny, Nz_e}, bounds_c.lo), double )
      var r_omega_mag_avg_z = region( ispace(int3d, {nx, ny, Nz_e}, bounds_c.lo), double )

      get_velocity_derivatives_z(r_prim_c_wo, r_gradu_l_z, r_gradu_r_z, interpolation.n_ghosts)
      compute_theta_avg(r_gradu_l_z, r_gradu_r_z, r_theta_avg_z)
      compute_omega_mag_avg(r_gradu_l_z, r_gradu_r_z, r_omega_mag_avg_z)

      HLLC_HLL_z( r_prim_l_z, r_prim_r_z, r_theta_avg_z, r_omega_mag_avg_z, r_flux_e_z )

      regentlib.c.legion_physical_region_destroy(__physical(r_gradu_l_z)[0])
      regentlib.c.legion_physical_region_destroy(__physical(r_gradu_r_z)[0])
      regentlib.c.legion_physical_region_destroy(__physical(r_theta_avg_z)[0])
      regentlib.c.legion_physical_region_destroy(__physical(r_omega_mag_avg_z)[0])

      __delete(r_gradu_l_z)
      __delete(r_gradu_r_z)
      __delete(r_theta_avg_z)
      __delete(r_omega_mag_avg_z)
    end

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

  regentlib.c.legion_physical_region_destroy(__physical(r_prim_l_z)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_prim_r_z)[0])

  regentlib.c.legion_physical_region_destroy(__physical(r_flux_c)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_flux_e_z)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_fder_c_z)[0])

  __delete(r_prim_l_z)
  __delete(r_prim_r_z)

  __delete(r_flux_c)
  __delete(r_flux_e_z)
  __delete(r_fder_c_z)
end



task add_viscous_zflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                                   r_aux_c    : region(ispace(int3d), auxiliary),
                                   r_visc     : region(ispace(int3d), transport_coeffs),
                                   r_tauij    : region(ispace(int3d), tensor2symm),
                                   r_rhs      : region(ispace(int3d), conserved),
                                   LU_z       : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c, r_aux_c.T, r_visc.kappa, r_tauij.{_13, _23, _33}, LU_z ),
  reads writes( r_rhs )
do

  var bounds_c = r_prim_c.ispace.bounds

  var Nz   = problem.NZ

  regentlib.assert(bounds_c.lo.z == interpolation.n_ghosts,          "Can only add Z flux derivative in the Z pencil")
  regentlib.assert(bounds_c.hi.z == Nz + interpolation.n_ghosts - 1, "Can only add Z flux derivative in the Z pencil")

  var nx = bounds_c.hi.x - bounds_c.lo.x + 1
  var ny = bounds_c.hi.y - bounds_c.lo.y + 1

  var r_q = region( ispace(int3d, {nx, ny, Nz}, bounds_c.lo), vect )

  var r_flux_c   = region( ispace(int3d, {nx, ny, Nz}, bounds_c.lo), conserved )
  var r_fder_c_z = region( ispace(int3d, {nx, ny, Nz}, bounds_c.lo), conserved )

  if (Nz >= 8) then
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

  regentlib.c.legion_physical_region_destroy(__physical(r_q)[0])

  regentlib.c.legion_physical_region_destroy(__physical(r_flux_c)[0])
  regentlib.c.legion_physical_region_destroy(__physical(r_fder_c_z)[0])

  __delete(r_q)

  __delete(r_flux_c)
  __delete(r_fder_c_z)
end



task add_nonconservative_viscous_zflux_der_to_rhs( r_prim_c   : region(ispace(int3d), primitive),
                                                   r_aux_c    : region(ispace(int3d), auxiliary),
                                                   r_visc     : region(ispace(int3d), transport_coeffs),
                                                   r_gradu    : region(ispace(int3d), tensor2),
                                                   r_grad2u   : region(ispace(int3d), tensor2),
                                                   r_tauij    : region(ispace(int3d), tensor2symm),
                                                   r_rhs      : region(ispace(int3d), conserved),
                                                   LU_z       : region(ispace(int3d), LU_struct),
                                                   LU2_z      : region(ispace(int3d), LU_struct) )
where
  reads( r_prim_c.{u, v, w}, r_aux_c.T, r_visc, r_tauij.{_13, _23, _33}, LU_z, LU2_z ),
  reads( r_gradu.{_11, _13, _22, _23, _31, _32, _33}, r_grad2u.{_13, _23, _33} ),
  reads writes( r_rhs.{rhou, rhov, rhow, rhoE} )
do

  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1
  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1
  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1

  var bounds_c = r_prim_c.ispace.bounds

  regentlib.assert(bounds_c.lo.z - interpolation.n_ghosts == 0, "Can only add_nonconservative_viscous_zflux_der_to_rhs in the Z pencil")

  if (nz >= 8) then
    if viscous then
      var r_q_div     = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), scalar )
      var r_tauij_div = region( ispace(int3d, {nx, ny, nz}, bounds_c.lo), vect )

      get_tauij_div_z(r_gradu, r_grad2u, r_visc, r_tauij_div, LU_z)
      get_q_div_z(r_aux_c, r_visc, r_q_div, LU_z, LU2_z)

      for i in r_rhs do
        var viscpower = r_tauij[i]._13 * r_gradu[i]._13 + r_tauij[i]._23 * r_gradu[i]._23 + r_tauij[i]._33 * r_gradu[i]._33
                      + r_prim_c[i].u * r_tauij_div[i]._1 + r_prim_c[i].v * r_tauij_div[i]._2+ r_prim_c[i].w * r_tauij_div[i]._3

        r_rhs[i].rhou += r_tauij_div[i]._1
        r_rhs[i].rhov += r_tauij_div[i]._2
        r_rhs[i].rhow += r_tauij_div[i]._3
        r_rhs[i].rhoE += - r_q_div[i]._ + viscpower
      end

      regentlib.c.legion_physical_region_destroy(__physical(r_q_div)[0])
      regentlib.c.legion_physical_region_destroy(__physical(r_tauij_div)[0])
      __delete(r_q_div)
      __delete(r_tauij_div)
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
                                 r_grad2u : region(ispace(int3d), tensor2),
                                 LU_x     : region(ispace(int3d), LU_struct),
                                 LU2_x    : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c.{u, v, w}, LU_x, LU2_x), reads writes (r_gradu.{_11, _21, _31}, r_grad2u.{_11, _21, _31})
do
  var nx = r_gradu.ispace.bounds.hi.x - r_gradu.ispace.bounds.lo.x + 1

  var token = 1
  if (nx < 8) then
    for i in r_gradu do
      r_gradu[i]._11 = 0.0
      r_gradu[i]._21 = 0.0
      r_gradu[i]._31 = 0.0
      r_grad2u[i]._11 = 0.0
      r_grad2u[i]._21 = 0.0
      r_grad2u[i]._31 = 0.0
    end
    return token
  end

  token += ddx_u(r_prim_c, r_gradu, LU_x)
  token += ddx_v(r_prim_c, r_gradu, LU_x)
  token += ddx_w(r_prim_c, r_gradu, LU_x) 

  token += d2dx2_u(r_prim_c, r_grad2u, LU2_x)
  token += d2dx2_v(r_prim_c, r_grad2u, LU2_x)
  token += d2dx2_w(r_prim_c, r_grad2u, LU2_x) 

  return token
end



task get_velocity_y_derivatives( r_prim_c : region(ispace(int3d), primitive),
                                 r_gradu  : region(ispace(int3d), tensor2),
                                 r_grad2u : region(ispace(int3d), tensor2),
                                 LU_y     : region(ispace(int3d), LU_struct),
                                 LU2_y    : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c.{u, v, w}, LU_y, LU2_y), reads writes (r_gradu.{_12, _22, _32}, r_grad2u.{_12, _22, _32})
do
  var ny = r_gradu.ispace.bounds.hi.y - r_gradu.ispace.bounds.lo.y + 1

  var token = 1
  if (ny < 8) then
    for i in r_gradu do
      r_gradu[i]._12 = 0.0
      r_gradu[i]._22 = 0.0
      r_gradu[i]._32 = 0.0
      r_grad2u[i]._12 = 0.0
      r_grad2u[i]._22 = 0.0
      r_grad2u[i]._32 = 0.0
    end
    return token
  end

  token += ddy_u(r_prim_c, r_gradu, LU_y)
  token += ddy_v(r_prim_c, r_gradu, LU_y)
  token += ddy_w(r_prim_c, r_gradu, LU_y) 

  token += d2dy2_u(r_prim_c, r_grad2u, LU2_y)
  token += d2dy2_v(r_prim_c, r_grad2u, LU2_y)
  token += d2dy2_w(r_prim_c, r_grad2u, LU2_y) 

  return token
end



task get_velocity_z_derivatives( r_prim_c : region(ispace(int3d), primitive),
                                 r_gradu  : region(ispace(int3d), tensor2),
                                 r_grad2u : region(ispace(int3d), tensor2),
                                 LU_z     : region(ispace(int3d), LU_struct),
                                 LU2_z    : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c.{u, v, w}, LU_z, LU2_z), reads writes (r_gradu.{_13, _23, _33}, r_grad2u.{_13, _23, _33})
do
  var nz = r_gradu.ispace.bounds.hi.z - r_gradu.ispace.bounds.lo.z + 1

  var token = 1
  if (nz < 8) then
    for i in r_gradu do
      r_gradu[i]._13 = 0.0
      r_gradu[i]._23 = 0.0
      r_gradu[i]._33 = 0.0
      r_grad2u[i]._13 = 0.0
      r_grad2u[i]._23 = 0.0
      r_grad2u[i]._33 = 0.0
    end
    return token
  end

  token += ddz_u(r_prim_c, r_gradu, LU_z)
  token += ddz_v(r_prim_c, r_gradu, LU_z)
  token += ddz_w(r_prim_c, r_gradu, LU_z) 

  token += d2dz2_u(r_prim_c, r_grad2u, LU2_z)
  token += d2dz2_v(r_prim_c, r_grad2u, LU2_z)
  token += d2dz2_w(r_prim_c, r_grad2u, LU2_z) 

  return token
end



-- task get_temperature_x_derivatives( r_aux_c : region(ispace(int3d), auxiliary),
--                                     r_gradT : region(ispace(int3d), vect),
--                                     LU_x    : region(ispace(int3d), LU_struct) )
-- where
--   reads (r_aux_c.{T}, LU_x), reads writes (r_gradT.{_1})
-- do
--   var nx = r_aux_c.ispace.bounds.hi.x - r_aux_c.ispace.bounds.lo.x + 1
-- 
--   var token = 1
--   if (nx < 8) then
--     for i in r_gradT do
--       r_gradT[i]._1 = 0.0
--     end
--     return token
--   end
-- 
--   token += ddx_T(r_aux_c, r_gradT, LU_x)
-- 
--   return token
-- end
-- 
-- task get_temperature_y_derivatives( r_aux_c : region(ispace(int3d), auxiliary),
--                                     r_gradT : region(ispace(int3d), vect),
--                                     LU_y    : region(ispace(int3d), LU_struct) )
-- where
--   reads (r_aux_c.{T}, LU_y), reads writes (r_gradT.{_2})
-- do
--   var ny = r_aux_c.ispace.bounds.hi.y - r_aux_c.ispace.bounds.lo.y + 1
-- 
--   var token = 1
--   if (ny < 8) then
--     for i in r_gradT do
--       r_gradT[i]._2 = 0.0
--     end
--     return token
--   end
-- 
--   token += ddy_T(r_aux_c, r_gradT, LU_y)
-- 
--   return token
-- end
-- 
-- task get_temperature_z_derivatives( r_aux_c : region(ispace(int3d), auxiliary),
--                                     r_gradT : region(ispace(int3d), vect),
--                                     LU_z    : region(ispace(int3d), LU_struct) )
-- where
--   reads (r_aux_c.{T}, LU_z), reads writes (r_gradT.{_3})
-- do
--   var nz = r_aux_c.ispace.bounds.hi.z - r_aux_c.ispace.bounds.lo.z + 1
-- 
--   var token = 1
--   if (nz < 8) then
--     for i in r_gradT do
--       r_gradT[i]._3 = 0.0
--     end
--     return token
--   end
-- 
--   token += ddz_T(r_aux_c, r_gradT, LU_z)
-- 
--   return token
-- end



task get_density_x_derivatives( r_prim_c  : region(ispace(int3d), primitive),
                                r_gradrho : region(ispace(int3d), vect),
                                LU_x      : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c.{rho}, LU_x), reads writes (r_gradrho.{_1})
do
  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1

  var token = 1
  if (nx < 8) then
    for i in r_gradrho do
      r_gradrho[i]._1 = 0.0
    end
    return token
  end

  token += ddx_rho_prim(r_prim_c, r_gradrho, LU_x)

  return token
end



task get_density_y_derivatives( r_prim_c  : region(ispace(int3d), primitive),
                                r_gradrho : region(ispace(int3d), vect),
                                LU_y      : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c.{rho}, LU_y), reads writes (r_gradrho.{_2})
do
  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1

  var token = 1
  if (ny < 8) then
    for i in r_gradrho do
      r_gradrho[i]._2 = 0.0
    end
    return token
  end

  token += ddy_rho_prim(r_prim_c, r_gradrho, LU_y)

  return token
end



task get_density_z_derivatives( r_prim_c  : region(ispace(int3d), primitive),
                                r_gradrho : region(ispace(int3d), vect),
                                LU_z      : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c.{rho}, LU_z), reads writes (r_gradrho.{_3})
do
  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1

  var token = 1
  if (nz < 8) then
    for i in r_gradrho do
      r_gradrho[i]._3 = 0.0
    end
    return token
  end

  token += ddz_rho_prim(r_prim_c, r_gradrho, LU_z)

  return token
end



task get_pressure_x_derivatives( r_prim_c : region(ispace(int3d), primitive),
                                 r_gradp  : region(ispace(int3d), vect),
                                 LU_x     : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c.{p}, LU_x), reads writes (r_gradp.{_1})
do
  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1

  var token = 1
  if (nx < 8) then
    for i in r_gradp do
      r_gradp[i]._1 = 0.0
    end
    return token
  end

  token += ddx_p(r_prim_c, r_gradp, LU_x)

  return token
end



task get_pressure_y_derivatives( r_prim_c : region(ispace(int3d), primitive),
                                 r_gradp  : region(ispace(int3d), vect),
                                 LU_y     : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c.{p}, LU_y), reads writes (r_gradp.{_2})
do
  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1

  var token = 1
  if (ny < 8) then
    for i in r_gradp do
      r_gradp[i]._2 = 0.0
    end
    return token
  end

  token += ddy_p(r_prim_c, r_gradp, LU_y)

  return token
end



task get_pressure_z_derivatives( r_prim_c : region(ispace(int3d), primitive),
                                 r_gradp  : region(ispace(int3d), vect),
                                 LU_z     : region(ispace(int3d), LU_struct) )
where
  reads (r_prim_c.{p}, LU_z), reads writes (r_gradp.{_3})
do
  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1

  var token = 1
  if (nz < 8) then
    for i in r_gradp do
      r_gradp[i]._3 = 0.0
    end
    return token
  end

  token += ddz_p(r_prim_c, r_gradp, LU_z)

  return token
end


-- __demand(__inline)
-- task pre_substep()
-- where
-- do
--       -- if problem.viscous then
--       --   -- Get the transport coefficients.
--       --   __demand(__parallel)
--       --   for i in pencil_interior do
--       --     problem.get_transport_coeffs( p_prim_c_y[i], p_aux_c_y[i], p_visc_y[i] )
--       --   end
-- 
--       --   -- Get the viscous stress tensor.
--       --   __demand(__parallel)
--       --   for i in pencil_interior do
--       --     get_tauij( p_gradu_y[i], p_tauij_y[i], p_visc_y[i] )
--       --   end
--       -- end
-- 
-- end


-- __demand(__inline)
-- task get_RHS( r_rhs : region(ispace(int3d), conserved),
--               p_rhs_x : partition(disjoint, r_rhs, ispace(int2d)),
--               p_rhs_y : partition(disjoint, r_rhs, ispace(int2d)),
--               p_rhs_z : partition(disjoint, r_rhs, ispace(int2d)),
-- 
--               r_prim_c : region(ispace(int3d), primitive),
--               p_prim_c_x : partition(disjoint, r_prim_c, ispace(int2d)),
--               p_prim_c_y : partition(disjoint, r_prim_c, ispace(int2d)),
--               p_prim_c_z : partition(disjoint, r_prim_c, ispace(int2d)),
--               p_prim_c_x_wg : partition(disjoint, r_prim_c, ispace(int2d)),
--               p_prim_c_y_wg : partition(disjoint, r_prim_c, ispace(int2d)),
--               p_prim_c_z_wg : partition(disjoint, r_prim_c, ispace(int2d)),
-- 
--               r_aux_c : region(ispace(int3d), auxiliary),
--               p_aux_c_x : partition(disjoint, r_aux_c, ispace(int2d)),
--               p_aux_c_y : partition(disjoint, r_aux_c, ispace(int2d)),
--               p_aux_c_z : partition(disjoint, r_aux_c, ispace(int2d)),
-- 
--               r_visc : region(ispace(int3d), transport_coeffs),
--               p_visc_x : partition(disjoint, r_visc, ispace(int2d)),
--               p_visc_y : partition(disjoint, r_visc, ispace(int2d)),
--               p_visc_z : partition(disjoint, r_visc, ispace(int2d)),
-- 
--               r_gradu : region(ispace(int3d), tensor2),
--               p_gradu_x : partition(disjoint, r_gradu, ispace(int2d)),
--               p_gradu_y : partition(disjoint, r_gradu, ispace(int2d)),
--               p_gradu_z : partition(disjoint, r_gradu, ispace(int2d)),
-- 
--               r_grad2u : region(ispace(int3d), tensor2),
--               p_grad2u_x : partition(disjoint, r_grad2u, ispace(int2d)),
--               p_grad2u_y : partition(disjoint, r_grad2u, ispace(int2d)),
--               p_grad2u_z : partition(disjoint, r_grad2u, ispace(int2d)),
-- 
--               r_tauij : region(ispace(int3d), tensor2symm),
--               p_tauij_x : partition(disjoint, r_tauij, ispace(int2d)),
--               p_tauij_y : partition(disjoint, r_tauij, ispace(int2d)),
--               p_tauij_z : partition(disjoint, r_tauij, ispace(int2d)),
-- 
--               LU_x : region(ispace(int3d), LU_struct),
--               p_LU_x : partition(disjoint, LU_x, ispace(int2d)),
--               LU_y : region(ispace(int3d), LU_struct),
--               p_LU_y : partition(disjoint, LU_y, ispace(int2d)),
--               LU_z : region(ispace(int3d), LU_struct),
--               p_LU_z : partition(disjoint, LU_z, ispace(int2d)),
-- 
--               LU_N_x : region(ispace(int3d), LU_struct),
--               p_LU_N_x : partition(disjoint, LU_N_x, ispace(int2d)),
--               LU_N_y : region(ispace(int3d), LU_struct),
--               p_LU_N_y : partition(disjoint, LU_N_y, ispace(int2d)),
--               LU_N_z : region(ispace(int3d), LU_struct),
--               p_LU_N_z : partition(disjoint, LU_N_z, ispace(int2d)),
-- 
--               LU2_N_x : region(ispace(int3d), LU_struct),
--               p_LU2_N_x : partition(disjoint, LU2_N_x, ispace(int2d)),
--               LU2_N_y : region(ispace(int3d), LU_struct),
--               p_LU2_N_y : partition(disjoint, LU2_N_y, ispace(int2d)),
--               LU2_N_z : region(ispace(int3d), LU_struct),
--               p_LU2_N_z : partition(disjoint, LU2_N_z, ispace(int2d)),
-- 
--               pencil_interior : ispace(int2d) )
-- where
--   reads writes(r_rhs),
--   reads(r_prim_c, r_aux_c.T, r_visc, r_tauij, r_gradu, r_grad2u),
--   reads(LU_x, LU_y, LU_z, LU_N_x, LU_N_y, LU_N_z, LU2_N_x, LU2_N_y, LU2_N_z)
-- do
-- 
--       -- Set RHS to zero.
--       __demand(__parallel)
--       for i in pencil_interior do
--         set_zero_cnsr( p_rhs_y[i] )
--       end
--       
--       -- Add x-direction convective flux derivative to RHS.
--       __demand(__parallel)
--       for i in pencil_interior do
--         add_xflux_der_to_rhs( p_prim_c_x_wg[i], p_rhs_x[i], p_LU_x[i] )
--       end
-- 
--       -- Add x-direction viscous flux derivative to RHS.
--       if problem.conservative_viscous_terms then
--         __demand(__parallel)
--         for i in pencil_interior do
--           add_viscous_xflux_der_to_rhs( p_prim_c_x[i], p_aux_c_x[i], p_visc_x[i], p_tauij_x[i], p_rhs_x[i], p_LU_N_x[i] )
--         end
--       else
--         __demand(__parallel)
--         for i in pencil_interior do
--           add_nonconservative_viscous_xflux_der_to_rhs( p_prim_c_x[i], p_aux_c_x[i], p_visc_x[i], p_gradu_x[i], p_grad2u_x[i], 
--                                                         p_tauij_x[i], p_rhs_x[i], p_LU_N_x[i], p_LU2_N_x[i] )
--         end
--       end
-- 
--       -- Add y-direction convective flux derivative to RHS.
--       __demand(__parallel)
--       for i in pencil_interior do
--         add_yflux_der_to_rhs( p_prim_c_y_wg[i], p_rhs_y[i], p_LU_y[i] )
--       end
-- 
--       -- Add y-direction viscous flux derivative to RHS.
--       if problem.conservative_viscous_terms then
--         __demand(__parallel)
--         for i in pencil_interior do
--           add_viscous_yflux_der_to_rhs( p_prim_c_y[i], p_aux_c_y[i], p_visc_y[i], p_tauij_y[i], p_rhs_y[i], p_LU_N_y[i] )
--         end
--       else
--         __demand(__parallel)
--         for i in pencil_interior do
--           add_nonconservative_viscous_yflux_der_to_rhs( p_prim_c_y[i], p_aux_c_y[i], p_visc_y[i], p_gradu_y[i], p_grad2u_y[i], 
--                                                         p_tauij_y[i], p_rhs_y[i], p_LU_N_y[i], p_LU2_N_y[i] )
--         end
--       end
-- 
--       -- Add z-direction convective flux derivative to RHS.
--       __demand(__parallel)
--       for i in pencil_interior do
--         add_zflux_der_to_rhs( p_prim_c_z_wg[i], p_rhs_z[i], p_LU_z[i] )
--       end
-- 
--       -- Add z-direction viscous flux derivative to RHS.
--       if problem.conservative_viscous_terms then
--         __demand(__parallel)
--         for i in pencil_interior do
--           add_viscous_zflux_der_to_rhs( p_prim_c_z[i], p_aux_c_z[i], p_visc_z[i], p_tauij_z[i], p_rhs_z[i], p_LU_N_z[i] )
--         end
--       else
--         __demand(__parallel)
--         for i in pencil_interior do
--           add_nonconservative_viscous_zflux_der_to_rhs( p_prim_c_z[i], p_aux_c_z[i], p_visc_z[i], p_gradu_z[i], p_grad2u_z[i], 
--                                                         p_tauij_z[i], p_rhs_z[i], p_LU_N_z[i], p_LU2_N_z[i] )
--         end
--       end
-- 
-- end

-- __demand(__inline)
-- task post_substep()
-- where
-- do
--       -- Update the primitive variables.
--       __demand(__parallel)
--       for i in pencil_interior do
--         token += get_primitive_r(p_cnsr_y[i], p_prim_c_y[i])
--       end
-- 
--       -- Update temperature.
--       __demand(__parallel)
--       for i in pencil_interior do
--         token += get_temperature_r( p_prim_c_y[i], p_aux_c_y[i] )
--       end
-- 
--       -- Update velocity gradient tensor.
--       __demand(__parallel)
--       for i in pencil_interior do
--         token += get_velocity_x_derivatives( p_prim_c_x[i], p_gradu_x[i], p_grad2u_x[i], p_LU_N_x[i], p_LU2_N_x[i] )
--       end
--       __demand(__parallel)
--       for i in pencil_interior do
--         token += get_velocity_y_derivatives( p_prim_c_y[i], p_gradu_y[i], p_grad2u_y[i], p_LU_N_y[i], p_LU2_N_y[i] )
--       end
--       __demand(__parallel)
--       for i in pencil_interior do
--         token += get_velocity_z_derivatives( p_prim_c_z[i], p_gradu_z[i], p_grad2u_z[i], p_LU_N_z[i], p_LU2_N_z[i] )
--       end
-- 
--       -- Get the density derivatives
--       __demand(__parallel)
--       for i in pencil_interior do
--         token += get_density_x_derivatives( p_prim_c_x[i], p_gradrho_x[i], p_LU_N_x[i] )
--       end
--       __demand(__parallel)
--       for i in pencil_interior do
--         token += get_density_y_derivatives( p_prim_c_y[i], p_gradrho_y[i], p_LU_N_y[i] )
--       end
--       __demand(__parallel)
--       for i in pencil_interior do
--         token += get_density_z_derivatives( p_prim_c_z[i], p_gradrho_z[i], p_LU_N_z[i] )
--       end
--  
--       -- Get the pressure derivatives
--       __demand(__parallel)
--       for i in pencil_interior do
--         token += get_pressure_x_derivatives( p_prim_c_x[i], p_gradp_x[i], p_LU_N_x[i] )
--       end
--       __demand(__parallel)
--       for i in pencil_interior do
--         token += get_pressure_y_derivatives( p_prim_c_y[i], p_gradp_y[i], p_LU_N_y[i] )
--       end
--       __demand(__parallel)
--       for i in pencil_interior do
--         token += get_pressure_z_derivatives( p_prim_c_z[i], p_gradp_z[i], p_LU_N_z[i] )
--       end
--  
--       -- Fill ghost cells in non-periodic directions first
--       if not problem.periodic_x then
--         __demand(__parallel)
--         for i in pencil_interior do
--           -- Fill in ghost cells
--           nonperiodic_ghost_cells_x(p_coords_x[i], p_prim_c_x_wg[i], p_gradrho_x[i], p_gradu_x[i], p_gradp_x[i], tsim, n_ghosts)
--         end
--       end
--       if not problem.periodic_y then
--         __demand(__parallel)
--         for i in pencil_interior do
--           -- Fill in ghost cells
--           nonperiodic_ghost_cells_y(p_coords_y[i], p_prim_c_y_wg[i], tsim, n_ghosts)
--         end
--       end
--       if not problem.periodic_z then
--         __demand(__parallel)
--         for i in pencil_interior do
--           -- Fill in ghost cells
--           nonperiodic_ghost_cells_z(p_coords_z[i], p_prim_c_z_wg[i], tsim, n_ghosts)
--         end
--       end
-- 
--       -- Fill ghost cells in periodic directions next
--       if problem.periodic_x then
--         __demand(__parallel)
--         for i in pencil_interior do
--           -- Fill in ghost cells
--           periodic_ghost_cells_x(p_prim_c_x_wg[i], n_ghosts)
--         end
--       end
--       if problem.periodic_y then
--         __demand(__parallel)
--         for i in pencil_interior do
--           -- Fill in ghost cells
--           periodic_ghost_cells_y(p_prim_c_y_wg[i], n_ghosts)
--         end
--       end
--       if problem.periodic_z then
--         __demand(__parallel)
--         for i in pencil_interior do
--           -- Fill in ghost cells
--           periodic_ghost_cells_z(p_prim_c_z_wg[i], n_ghosts)
--         end
--       end
-- end
