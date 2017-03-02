import "regent"

require("fields")

local c     = regentlib.c
local cmath = terralib.includec("math.h")

local gamma = 1.4
local gamma_m1 = gamma - 1.0
local onebygm1 = 1.0 / gamma_m1

__demand(__inline)
task get_primitive( rho  : double,
                    rhou : double,
                    rhov : double,
                    rhow : double,
                    rhoE : double )

  var prim : double[5]

  var onebyrho : double = 1.0 / rho
  prim[0] = rho
  prim[1] = rhou * onebyrho
  prim[2] = rhov * onebyrho
  prim[3] = rhow * onebyrho
  prim[4] = gamma_m1 * ( rhoE - 0.5 * rho * (prim[1]*prim[1] + prim[2]*prim[2] + prim[3]*prim[3]) )

  return prim
end

__demand(__inline)
task get_conserved( rho : double,
                    u   : double,
                    v   : double,
                    w   : double,
                    p   : double )

  var cnsr : double[5]

  cnsr[0] = rho
  cnsr[1] = rho * u
  cnsr[2] = rho * v
  cnsr[3] = rho * w
  cnsr[4] = onebygm1 * p + 0.5 * rho * (u*u + v*v + w*w)

  return cnsr
end

__demand(__inline)
task get_xfluxes( rho  : double,
                  u    : double,
                  v    : double,
                  w    : double,
                  p    : double, 
                  rhou : double,
                  rhov : double,
                  rhow : double,
                  rhoE : double )

  var flux : double[5]

  flux[0] =  rhou
  flux[1] =  rhou * u + p
  flux[2] =  rhou * v
  flux[3] =  rhou * w
  flux[4] = (rhoE + p) * u

  return flux
end

__demand(__inline)
task get_sos( rho : double, p : double )
  return cmath.sqrt(gamma*p/rho)
end

__demand(__inline)
task get_rho_sos_avg_x( r_prim_c : region(ispace(int3d), primitive),
                        idx      : int3d,
                        Nx       : int64,
                        Ny       : int64,
                        Nz       : int64)
where
  reads(r_prim_c)
do
  var rhosos : double[2]

  var idxm1 = [poff(idx, -1, 0, 0, Nx, Ny, Nz)]
  rhosos[0] = 0.5*( r_prim_c[ idxm1 ].rho + r_prim_c[ idx ].rho )
  rhosos[1] = 0.5*( get_sos(r_prim_c[idxm1].rho, r_prim_c[idxm1].p) + get_sos(r_prim_c[idx].rho, r_prim_c[idx].p))

  return rhosos
end

__demand(__inline)
task get_rho_sos_avg_y( r_prim_c : region(ispace(int3d), primitive),
                        idx      : int3d,
                        Nx       : int64,
                        Ny       : int64,
                        Nz       : int64)
where
  reads(r_prim_c)
do
  var rhosos : double[2]

  var idxm1 = [poff(idx, 0, -1, 0, Nx, Ny, Nz)]
  rhosos[0] = 0.5*( r_prim_c[ idxm1 ].rho + r_prim_c[ idx ].rho )
  rhosos[1] = 0.5*( get_sos(r_prim_c[idxm1].rho, r_prim_c[idxm1].p) + get_sos(r_prim_c[idx].rho, r_prim_c[idx].p))

  return rhosos
end

__demand(__inline)
task get_rho_sos_avg_z( r_prim_c : region(ispace(int3d), primitive),
                        idx      : int3d,
                        Nx       : int64,
                        Ny       : int64,
                        Nz       : int64)
where
  reads(r_prim_c)
do
  var rhosos : double[2]

  var idxm1 = [poff(idx, 0, 0, -1, Nx, Ny, Nz)]
  rhosos[0] = 0.5*( r_prim_c[ idxm1 ].rho + r_prim_c[ idx ].rho )
  rhosos[1] = 0.5*( get_sos(r_prim_c[idxm1].rho, r_prim_c[idxm1].p) + get_sos(r_prim_c[idx].rho, r_prim_c[idx].p))

  return rhosos
end

task get_primitive_r( r_cnsr : region(ispace(int3d), conserved),
                      r_prim : region(ispace(int3d), primitive) )
where
  reads(r_cnsr), writes(r_prim)
do

  for i in r_prim do
    var prim : double[5] = get_primitive( r_cnsr[i].rho, r_cnsr[i].rhou, r_cnsr[i].rhov, r_cnsr[i].rhow, r_cnsr[i].rhoE )

    r_prim[i].rho = prim[0]
    r_prim[i].u   = prim[1]
    r_prim[i].v   = prim[2]
    r_prim[i].w   = prim[3]
    r_prim[i].p   = prim[4]
  end

  return 1
end

task get_conserved_r( r_prim : region(ispace(int3d), primitive),
                      r_cnsr : region(ispace(int3d), conserved) )
where
  reads(r_prim), writes(r_cnsr)
do

  for i in r_cnsr do
    var cnsr : double[5] = get_conserved( r_prim[i].rho, r_prim[i].u, r_prim[i].v, r_prim[i].w, r_prim[i].p )
    r_cnsr[i].rho  = cnsr[0]
    r_cnsr[i].rhou = cnsr[1]
    r_cnsr[i].rhov = cnsr[2]
    r_cnsr[i].rhow = cnsr[3]
    r_cnsr[i].rhoE = cnsr[4]
  end

  return 1
end

__demand(__inline)
task get_char_values_x( r_prim_c : region(ispace(int3d), primitive),
                        rho_avg  : double,
                        sos_avg  : double,
                        idx      : int3d,
                        Nx       : int64,
                        Ny       : int64,
                        Nz       : int64)
where
  reads(r_prim_c)
do
  var char_values : double[6][5]

  for i = -3, 3 do
    var p = [poff(idx, i, 0, 0, Nx, Ny, Nz)]
    char_values[0][i+3] = -0.5*rho_avg*sos_avg * r_prim_c[p].u + 0.5*r_prim_c[p].p
    char_values[1][i+3] = r_prim_c[p].rho - r_prim_c[p].p/(sos_avg*sos_avg)
    char_values[2][i+3] = r_prim_c[p].v
    char_values[3][i+3] = r_prim_c[p].w
    char_values[4][i+3] = 0.5*rho_avg*sos_avg * r_prim_c[p].u + 0.5*r_prim_c[p].p
  end

  return char_values
end

__demand(__inline)
task get_char_values_y( r_prim_c : region(ispace(int3d), primitive),
                        rho_avg  : double,
                        sos_avg  : double,
                        idx      : int3d,
                        Nx       : int64,
                        Ny       : int64,
                        Nz       : int64)
where
  reads(r_prim_c)
do
  var char_values : double[6][5]

  for i = -3, 3 do
    var p = [poff(idx, 0, i, 0, Nx, Ny, Nz)]
    char_values[0][i+3] = -0.5*rho_avg*sos_avg * r_prim_c[p].v + 0.5*r_prim_c[p].p
    char_values[1][i+3] = r_prim_c[p].u
    char_values[2][i+3] = r_prim_c[p].rho - r_prim_c[p].p/(sos_avg*sos_avg)
    char_values[3][i+3] = r_prim_c[p].w
    char_values[4][i+3] = 0.5*rho_avg*sos_avg * r_prim_c[p].v + 0.5*r_prim_c[p].p
  end

  return char_values
end

__demand(__inline)
task get_char_values_z( r_prim_c : region(ispace(int3d), primitive),
                        rho_avg  : double,
                        sos_avg  : double,
                        idx      : int3d,
                        Nx       : int64,
                        Ny       : int64,
                        Nz       : int64)
where
  reads(r_prim_c)
do
  var char_values : double[6][5]

  for i = -3, 3 do
    var p = [poff(idx, 0, 0, i, Nx, Ny, Nz)]
    char_values[0][i+3] = -0.5*rho_avg*sos_avg * r_prim_c[p].w + 0.5*r_prim_c[p].p
    char_values[1][i+3] = r_prim_c[p].u
    char_values[2][i+3] = r_prim_c[p].v
    char_values[3][i+3] = r_prim_c[p].rho - r_prim_c[p].p/(sos_avg*sos_avg)
    char_values[4][i+3] = 0.5*rho_avg*sos_avg * r_prim_c[p].w + 0.5*r_prim_c[p].p
  end

  return char_values
end

task get_xfluxes_r( r_prim : region(ispace(int3d), primitive),
                    r_cnsr : region(ispace(int3d), conserved),
                    r_flux : region(ispace(int3d), conserved) )
where
  reads (r_prim, r_cnsr), writes (r_flux)
do
  for i in r_prim do
    var flux : double[5] =  get_xfluxes( r_prim[i].rho ,
                                         r_prim[i].u   ,
                                         r_prim[i].v   ,
                                         r_prim[i].w   ,
                                         r_prim[i].p   , 
                                         r_cnsr[i].rhou,
                                         r_cnsr[i].rhov,
                                         r_cnsr[i].rhow,
                                         r_cnsr[i].rhoE )
    r_flux[i].rho  = flux[0]
    r_flux[i].rhou = flux[1]
    r_flux[i].rhov = flux[2]
    r_flux[i].rhow = flux[3]
    r_flux[i].rhoE = flux[4]
  end
end

task get_yfluxes_r( r_prim : region(ispace(int3d), primitive),
                    r_cnsr : region(ispace(int3d), conserved),
                    r_flux : region(ispace(int3d), conserved) )
where
  reads (r_prim, r_cnsr), writes (r_flux)
do
  for i in r_prim do
    r_flux[i].rho  = r_cnsr[i].rhov
    r_flux[i].rhou = r_cnsr[i].rhov * r_prim[i].u
    r_flux[i].rhov = r_cnsr[i].rhov * r_prim[i].v + r_prim[i].p
    r_flux[i].rhow = r_cnsr[i].rhov * r_prim[i].w
    r_flux[i].rhoE =(r_cnsr[i].rhoE + r_prim[i].p) * r_prim[i].v
  end
end

task get_zfluxes_r( r_prim : region(ispace(int3d), primitive),
                    r_cnsr : region(ispace(int3d), conserved),
                    r_flux : region(ispace(int3d), conserved) )
where
  reads (r_prim, r_cnsr), writes (r_flux)
do
  for i in r_prim do
    r_flux[i].rho  = r_cnsr[i].rhow
    r_flux[i].rhou = r_cnsr[i].rhow * r_prim[i].u
    r_flux[i].rhov = r_cnsr[i].rhow * r_prim[i].v
    r_flux[i].rhow = r_cnsr[i].rhow * r_prim[i].w + r_prim[i].p
    r_flux[i].rhoE =(r_cnsr[i].rhoE + r_prim[i].p) * r_prim[i].w
  end
end

terra sign(x : double)
  return 2*([int](x >= 0)) - 1
end

task HLLC_x( r_prim_l_x : region(ispace(int3d), primitive),
             r_prim_r_x : region(ispace(int3d), primitive),
             r_flux_e_x : region(ispace(int3d), conserved) )
where
  reads(r_prim_l_x, r_prim_r_x), writes(r_flux_e_x)
do
  for i in r_flux_e_x do
    var u_avg   : double = 0.5*( r_prim_l_x[i].u + r_prim_r_x[i].u )
    var sos_l   : double = get_sos(r_prim_l_x[i].rho, r_prim_l_x[i].p)
    var sos_r   : double = get_sos(r_prim_r_x[i].rho, r_prim_r_x[i].p)
    var sos_avg : double = 0.5*( sos_l + sos_r )

    var s_L : double = cmath.fmin(u_avg - sos_avg, r_prim_l_x[i].u - sos_l)
    var s_R : double = cmath.fmax(u_avg + sos_avg, r_prim_r_x[i].u + sos_r)

    var s_m : double = cmath.fmin(0.0, s_L)
    var s_p : double = cmath.fmax(0.0, s_R)

    var s_star : double = ( r_prim_r_x[i].p - r_prim_l_x[i].p 
                            + r_prim_l_x[i].rho * r_prim_l_x[i].u * (s_L - r_prim_l_x[i].u) 
                            - r_prim_r_x[i].rho * r_prim_r_x[i].u * (s_R - r_prim_r_x[i].u) ) 
                          / ( r_prim_l_x[i].rho * (s_L - r_prim_l_x[i].u) 
                            - r_prim_r_x[i].rho * (s_R - r_prim_r_x[i].u))

    var Q_L : double[5] = get_conserved( r_prim_l_x[i].rho, r_prim_l_x[i].u, r_prim_l_x[i].v, r_prim_l_x[i].w, r_prim_l_x[i].p )
    var F_L : double[5] = get_xfluxes( r_prim_l_x[i].rho, r_prim_l_x[i].u, r_prim_l_x[i].v, r_prim_l_x[i].w, r_prim_l_x[i].p,
                                       Q_L[1], Q_L[2], Q_L[3], Q_L[4] )

    var Q_R : double[5] = get_conserved( r_prim_r_x[i].rho, r_prim_r_x[i].u, r_prim_r_x[i].v, r_prim_r_x[i].w, r_prim_r_x[i].p )
    var F_R : double[5] = get_xfluxes( r_prim_r_x[i].rho, r_prim_r_x[i].u, r_prim_r_x[i].v, r_prim_r_x[i].w, r_prim_r_x[i].p,
                                       Q_R[1], Q_R[2], Q_R[3], Q_R[4] )

    var chi_star_L : double = ( s_L - r_prim_l_x[i].u ) / ( s_L - s_star )
    var chi_star_R : double = ( s_R - r_prim_r_x[i].u ) / ( s_R - s_star )

    var Q_star_L : double[5]
    Q_star_L[0] = chi_star_L * ( r_prim_l_x[i].rho )
    Q_star_L[1] = chi_star_L * ( r_prim_l_x[i].rho * s_star )
    Q_star_L[2] = chi_star_L * ( Q_L[2] )
    Q_star_L[3] = chi_star_L * ( Q_L[3] )
    Q_star_L[4] = chi_star_L * ( Q_L[4] + (s_star - r_prim_l_x[i].u)*( r_prim_l_x[i].rho * s_star + r_prim_l_x[i].p/(s_L - r_prim_l_x[i].u) ) )

    var Q_star_R : double[5]
    Q_star_R[0] = chi_star_R * ( r_prim_r_x[i].rho )
    Q_star_R[1] = chi_star_R * ( r_prim_r_x[i].rho * s_star )
    Q_star_R[2] = chi_star_R * ( Q_R[2] )
    Q_star_R[3] = chi_star_R * ( Q_R[3] )
    Q_star_R[4] = chi_star_R * ( Q_R[4] + (s_star - r_prim_r_x[i].u)*( r_prim_r_x[i].rho * s_star + r_prim_r_x[i].p/(s_R - r_prim_r_x[i].u) ) )

    var switch = 0.5*(1 + sign(s_star))
    r_flux_e_x[i].rho  = switch * (F_L[0] + s_m*(Q_star_L[0] - Q_L[0])) + (1-switch) * (F_R[0] + s_p*(Q_star_R[0] - Q_R[0]))
    r_flux_e_x[i].rhou = switch * (F_L[1] + s_m*(Q_star_L[1] - Q_L[1])) + (1-switch) * (F_R[1] + s_p*(Q_star_R[1] - Q_R[1]))
    r_flux_e_x[i].rhov = switch * (F_L[2] + s_m*(Q_star_L[2] - Q_L[2])) + (1-switch) * (F_R[2] + s_p*(Q_star_R[2] - Q_R[2]))
    r_flux_e_x[i].rhow = switch * (F_L[3] + s_m*(Q_star_L[3] - Q_L[3])) + (1-switch) * (F_R[3] + s_p*(Q_star_R[3] - Q_R[3]))
    r_flux_e_x[i].rhoE = switch * (F_L[4] + s_m*(Q_star_L[4] - Q_L[4])) + (1-switch) * (F_R[4] + s_p*(Q_star_R[4] - Q_R[4]))
  end
end
