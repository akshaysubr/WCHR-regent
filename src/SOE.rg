import "regent"

require("fields")
require("EOS")

local c     = regentlib.c
local cmath = terralib.includec("math.h")

__demand(__inline)
task get_primitive( rho  : double,
                    rhou : double,
                    rhov : double,
                    rhow : double,
                    rhoE : double )

  var prim : double[5]
  
  var onebyrho = 1.0 / rho

  prim[0] = rho
  prim[1] = rhou * onebyrho
  prim[2] = rhov * onebyrho
  prim[3] = rhow * onebyrho

  var rhoe : double = rhoE - 0.5 * rho * (prim[1]*prim[1] + prim[2]*prim[2] + prim[3]*prim[3])
  prim[4] = get_pressure(rho, rhoe)

  return prim
end

__demand(__inline)
task get_conserved( rho : double,
                    u   : double,
                    v   : double,
                    w   : double,
                    p   : double )

  var cnsr : double[5]
  
  var rhoe : double = get_internal_energy(rho, p)
  
  cnsr[0] = rho
  cnsr[1] = rho * u
  cnsr[2] = rho * v
  cnsr[3] = rho * w
  cnsr[4] = rhoe + 0.5 * rho * (u*u + v*v + w*w)

  return cnsr
end

terra get_Rinv( rho : double, sos : double, Rinv : &double )
  Rinv[0 + 3*0] = 0.; Rinv[0 + 3*1] = -0.5*rho*sos; Rinv[0 + 3*2] = 0.5;
  Rinv[1 + 3*0] = 1.; Rinv[1 + 3*1] = 0.;           Rinv[1 + 3*2] = -1./(sos*sos);
  Rinv[2 + 3*0] = 0.; Rinv[2 + 3*1] =  0.5*rho*sos; Rinv[2 + 3*2] = 0.5;
end
get_Rinv:setinlined(true)

__demand(__inline)
task get_Rinv_r( rho : double, sos : double, r : region(ispace(int3d), double[9]), i : int3d(double[9], r) )
where
  reads writes (r)
do
  (@i)[0 + 3*0] = 0.; (@i)[0 + 3*1] = -0.5*rho*sos; (@i)[0 + 3*2] = 0.5;
  (@i)[1 + 3*0] = 1.; (@i)[1 + 3*1] = 0.;           (@i)[1 + 3*2] = -1./(sos*sos);
  (@i)[2 + 3*0] = 0.; (@i)[2 + 3*1] =  0.5*rho*sos; (@i)[2 + 3*2] = 0.5;
end

__demand(__inline)
task get_xfluxes( rho  : double,
                  u    : double,
                  v    : double,
                  w    : double,
                  p    : double )

  var flux : double[5]

  var rhoe : double = get_internal_energy(rho, p)
  var rhoE : double = rhoe + 0.5 * rho * (u*u + v*v + w*w)

  flux[0] =  rho * u
  flux[1] =  rho * u * u + p
  flux[2] =  rho * u * v
  flux[3] =  rho * u * w
  flux[4] = (rhoE + p) * u

  return flux
end

__demand(__inline)
task get_yfluxes( rho  : double,
                  u    : double,
                  v    : double,
                  w    : double,
                  p    : double )

  var flux : double[5]

  var rhoe : double = get_internal_energy(rho, p)
  var rhoE : double = rhoe + 0.5 * rho * (u*u + v*v + w*w)

  flux[0] =  rho * v
  flux[1] =  rho * v * u
  flux[2] =  rho * v * v + p
  flux[3] =  rho * v * w
  flux[4] = (rhoE + p) * v

  return flux
end

__demand(__inline)
task get_zfluxes( rho  : double,
                  u    : double,
                  v    : double,
                  w    : double,
                  p    : double )

  var flux : double[5]

  var rhoe : double = get_internal_energy(rho, p)
  var rhoE : double = rhoe + 0.5 * rho * (u*u + v*v + w*w)

  flux[0] =  rho * w
  flux[1] =  rho * w * u
  flux[2] =  rho * w * v
  flux[3] =  rho * w * w + p
  flux[4] = (rhoE + p) * w

  return flux
end

__demand(__inline)
task get_rho_sos_avg_x( r_prim_c : region(ispace(int3d), primitive),
                        idx      : int3d)
where
  reads(r_prim_c)
do
  var rhosos : double[2]

  var idxm1 = int3d {x = idx.x-1, y = idx.y, z = idx.z}
  rhosos[0] = 0.5*( r_prim_c[ idxm1 ].rho + r_prim_c[ idx ].rho )
  rhosos[1] = 0.5*( get_sos(r_prim_c[idxm1].rho, r_prim_c[idxm1].p) + get_sos(r_prim_c[idx].rho, r_prim_c[idx].p))

  return rhosos
end

__demand(__inline)
task get_rho_sos_avg_y( r_prim_c : region(ispace(int3d), primitive),
                        idx      : int3d)
where
  reads(r_prim_c)
do
  var rhosos : double[2]

  var idxm1 = int3d {x = idx.x, y = idx.y-1, z = idx.z}
  rhosos[0] = 0.5*( r_prim_c[ idxm1 ].rho + r_prim_c[ idx ].rho )
  rhosos[1] = 0.5*( get_sos(r_prim_c[idxm1].rho, r_prim_c[idxm1].p) + get_sos(r_prim_c[idx].rho, r_prim_c[idx].p))

  return rhosos
end

__demand(__inline)
task get_rho_sos_avg_z( r_prim_c : region(ispace(int3d), primitive),
                        idx      : int3d)
where
  reads(r_prim_c)
do
  var rhosos : double[2]

  var idxm1 = int3d {x = idx.x, y = idx.y, z = idx.z-1}
  rhosos[0] = 0.5*( r_prim_c[ idxm1 ].rho + r_prim_c[ idx ].rho )
  rhosos[1] = 0.5*( get_sos(r_prim_c[idxm1].rho, r_prim_c[idxm1].p) + get_sos(r_prim_c[idx].rho, r_prim_c[idx].p))

  return rhosos
end

task get_max_stable_dt_1d( r_prim_c : region(ispace(int3d), primitive),
                           dx       : double)
where
  reads(r_prim_c)
do
  var max_spectral_radius : double = 0.0

  for i in r_prim_c do
    var sos = get_sos(r_prim_c[i].rho, r_prim_c[i].p)
    var local_spectral_radius_x = (sos + cmath.fabs(r_prim_c[i].u))/dx

    max_spectral_radius = cmath.fmax(max_spectral_radius, local_spectral_radius_x)
  end

  return 1.0/max_spectral_radius
end

task get_max_stable_dt_2d( r_prim_c : region(ispace(int3d), primitive),
                           dx       : double,
                           dy       : double)
where
  reads(r_prim_c)
do
  var max_spectral_radius : double = 0.0

  for i in r_prim_c do
    var sos = get_sos(r_prim_c[i].rho, r_prim_c[i].p)
    var local_spectral_radius_x = (sos + cmath.fabs(r_prim_c[i].u))/dx
    var local_spectral_radius_y = (sos + cmath.fabs(r_prim_c[i].v))/dy
    
    var local_spectral_radius = local_spectral_radius_x + local_spectral_radius_y

    max_spectral_radius = cmath.fmax(max_spectral_radius, local_spectral_radius)
  end

  return 1.0/max_spectral_radius
end

task get_max_stable_dt_3d( r_prim_c : region(ispace(int3d), primitive),
                           dx       : double,
                           dy       : double,
                           dz       : double)
where
  reads(r_prim_c)
do
  var max_spectral_radius : double = 0.0

  for i in r_prim_c do
    var sos = get_sos(r_prim_c[i].rho, r_prim_c[i].p)
    var local_spectral_radius_x = (sos + cmath.fabs(r_prim_c[i].u))/dx
    var local_spectral_radius_y = (sos + cmath.fabs(r_prim_c[i].v))/dy
    var local_spectral_radius_z = (sos + cmath.fabs(r_prim_c[i].w))/dz
    
    var local_spectral_radius = local_spectral_radius_x + local_spectral_radius_y + local_spectral_radius_z

    max_spectral_radius = cmath.fmax(max_spectral_radius, local_spectral_radius)
  end
 
  return 1.0/max_spectral_radius
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

task get_temperature_r( r_prim : region(ispace(int3d), primitive),
                       r_aux  : region(ispace(int3d), auxiliary) )
where
  reads(r_prim.{rho, p}), writes(r_aux.T)
do
  for i in r_aux do
    r_aux[i].T = get_temperature( r_prim[i].rho, r_prim[i].p )
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
                        idx      : int3d )
where
  reads(r_prim_c)
do
  var char_values : double[6][5]

  for i = -3, 3 do
    var p = int3d {x = idx.x+i, y = idx.y, z = idx.z}
    char_values[0][i+3] = -0.5*rho_avg*sos_avg * r_prim_c[p].u + 0.5*r_prim_c[p].p
    char_values[1][i+3] = r_prim_c[p].rho - r_prim_c[p].p/(sos_avg*sos_avg)
    char_values[2][i+3] = r_prim_c[p].v
    char_values[3][i+3] = r_prim_c[p].w
    char_values[4][i+3] = 0.5*rho_avg*sos_avg * r_prim_c[p].u + 0.5*r_prim_c[p].p
  end

  return char_values
end

__demand(__inline)
task get_char_values_LB_x( r_prim_c : region(ispace(int3d), primitive),
                           rho_avg  : double,
                           sos_avg  : double,
                           idx      : int3d)
where
  reads(r_prim_c)
do
  var char_values : double[7][5]

  for i = -3, 4 do
    var p = int3d {x = idx.x+i, y = idx.y, z = idx.z}
    char_values[0][i+3] = -0.5*rho_avg*sos_avg * r_prim_c[p].u + 0.5*r_prim_c[p].p
    char_values[1][i+3] = r_prim_c[p].rho - r_prim_c[p].p/(sos_avg*sos_avg)
    char_values[2][i+3] = r_prim_c[p].v
    char_values[3][i+3] = r_prim_c[p].w
    char_values[4][i+3] = 0.5*rho_avg*sos_avg * r_prim_c[p].u + 0.5*r_prim_c[p].p
  end

  return char_values
end

__demand(__inline)
task get_char_values_RB_x( r_prim_c : region(ispace(int3d), primitive),
                           rho_avg  : double,
                           sos_avg  : double,
                           idx      : int3d)
where
  reads(r_prim_c)
do
  var char_values : double[7][5]

  for i = -4, 3 do
    var p = int3d {x = idx.x+i, y = idx.y, z = idx.z}
    char_values[0][i+4] = -0.5*rho_avg*sos_avg * r_prim_c[p].u + 0.5*r_prim_c[p].p
    char_values[1][i+4] = r_prim_c[p].rho - r_prim_c[p].p/(sos_avg*sos_avg)
    char_values[2][i+4] = r_prim_c[p].v
    char_values[3][i+4] = r_prim_c[p].w
    char_values[4][i+4] = 0.5*rho_avg*sos_avg * r_prim_c[p].u + 0.5*r_prim_c[p].p
  end

  return char_values
end

__demand(__inline)
task get_char_values_y( r_prim_c : region(ispace(int3d), primitive),
                        rho_avg  : double,
                        sos_avg  : double,
                        idx      : int3d )
where
  reads(r_prim_c)
do
  var char_values : double[6][5]

  for i = -3, 3 do
    var p = int3d {x = idx.x, y = idx.y+i, z = idx.z}
    char_values[0][i+3] = -0.5*rho_avg*sos_avg * r_prim_c[p].v + 0.5*r_prim_c[p].p
    char_values[1][i+3] = r_prim_c[p].u
    char_values[2][i+3] = r_prim_c[p].rho - r_prim_c[p].p/(sos_avg*sos_avg)
    char_values[3][i+3] = r_prim_c[p].w
    char_values[4][i+3] = 0.5*rho_avg*sos_avg * r_prim_c[p].v + 0.5*r_prim_c[p].p
  end

  return char_values
end

__demand(__inline)
task get_char_values_LB_y( r_prim_c : region(ispace(int3d), primitive),
                           rho_avg  : double,
                           sos_avg  : double,
                           idx      : int3d)
where
  reads(r_prim_c)
do
  var char_values : double[7][5]

  for i = -3, 4 do
    var p = int3d {x = idx.x, y = idx.y+i, z = idx.z}
    char_values[0][i+3] = -0.5*rho_avg*sos_avg * r_prim_c[p].v + 0.5*r_prim_c[p].p
    char_values[1][i+3] = r_prim_c[p].u
    char_values[2][i+3] = r_prim_c[p].rho - r_prim_c[p].p/(sos_avg*sos_avg)
    char_values[3][i+3] = r_prim_c[p].w
    char_values[4][i+3] = 0.5*rho_avg*sos_avg * r_prim_c[p].v + 0.5*r_prim_c[p].p
  end

  return char_values
end

__demand(__inline)
task get_char_values_RB_y( r_prim_c : region(ispace(int3d), primitive),
                           rho_avg  : double,
                           sos_avg  : double,
                           idx      : int3d)
where
  reads(r_prim_c)
do
  var char_values : double[7][5]

  for i = -4, 3 do
    var p = int3d {x = idx.x, y = idx.y+i, z = idx.z}
    char_values[0][i+4] = -0.5*rho_avg*sos_avg * r_prim_c[p].v + 0.5*r_prim_c[p].p
    char_values[1][i+4] = r_prim_c[p].u
    char_values[2][i+4] = r_prim_c[p].rho - r_prim_c[p].p/(sos_avg*sos_avg)
    char_values[3][i+4] = r_prim_c[p].w
    char_values[4][i+4] = 0.5*rho_avg*sos_avg * r_prim_c[p].v + 0.5*r_prim_c[p].p
  end

  return char_values
end

__demand(__inline)
task get_char_values_z( r_prim_c : region(ispace(int3d), primitive),
                        rho_avg  : double,
                        sos_avg  : double,
                        idx      : int3d)
where
  reads(r_prim_c)
do
  var char_values : double[6][5]

  for i = -3, 3 do
    var p = int3d {x = idx.x, y = idx.y, z = idx.z+i}
    char_values[0][i+3] = -0.5*rho_avg*sos_avg * r_prim_c[p].w + 0.5*r_prim_c[p].p
    char_values[1][i+3] = r_prim_c[p].u
    char_values[2][i+3] = r_prim_c[p].v
    char_values[3][i+3] = r_prim_c[p].rho - r_prim_c[p].p/(sos_avg*sos_avg)
    char_values[4][i+3] = 0.5*rho_avg*sos_avg * r_prim_c[p].w + 0.5*r_prim_c[p].p
  end

  return char_values
end

__demand(__inline)
task get_char_values_LB_z( r_prim_c : region(ispace(int3d), primitive),
                           rho_avg  : double,
                           sos_avg  : double,
                           idx      : int3d)
where
  reads(r_prim_c)
do
  var char_values : double[7][5]

  for i = -3, 4 do
    var p = int3d {x = idx.x, y = idx.y, z = idx.z+i}
    char_values[0][i+3] = -0.5*rho_avg*sos_avg * r_prim_c[p].w + 0.5*r_prim_c[p].p
    char_values[1][i+3] = r_prim_c[p].u
    char_values[2][i+3] = r_prim_c[p].v
    char_values[3][i+3] = r_prim_c[p].rho - r_prim_c[p].p/(sos_avg*sos_avg)
    char_values[4][i+3] = 0.5*rho_avg*sos_avg * r_prim_c[p].w + 0.5*r_prim_c[p].p
  end

  return char_values
end

__demand(__inline)
task get_char_values_RB_z( r_prim_c : region(ispace(int3d), primitive),
                           rho_avg  : double,
                           sos_avg  : double,
                           idx      : int3d)
where
  reads(r_prim_c)
do
  var char_values : double[7][5]

  for i = -4, 3 do
    var p = int3d {x = idx.x, y = idx.y, z = idx.z+i}
    char_values[0][i+4] = -0.5*rho_avg*sos_avg * r_prim_c[p].w + 0.5*r_prim_c[p].p
    char_values[1][i+4] = r_prim_c[p].u
    char_values[2][i+4] = r_prim_c[p].v
    char_values[3][i+4] = r_prim_c[p].rho - r_prim_c[p].p/(sos_avg*sos_avg)
    char_values[4][i+4] = 0.5*rho_avg*sos_avg * r_prim_c[p].w + 0.5*r_prim_c[p].p
  end

  return char_values
end

__demand(__inline)
task get_xfluxes_r( r_prim : region(ispace(int3d), primitive),
                    r_flux : region(ispace(int3d), conserved) )
where
  reads (r_prim), writes (r_flux)
do
  for i in r_prim do
    var flux : double[5] =  get_xfluxes( r_prim[i].rho ,
                                         r_prim[i].u   ,
                                         r_prim[i].v   ,
                                         r_prim[i].w   ,
                                         r_prim[i].p   )
    r_flux[i].rho  = flux[0]
    r_flux[i].rhou = flux[1]
    r_flux[i].rhov = flux[2]
    r_flux[i].rhow = flux[3]
    r_flux[i].rhoE = flux[4]
  end
end

__demand(__inline)
task get_yfluxes_r( r_prim : region(ispace(int3d), primitive),
                    r_flux : region(ispace(int3d), conserved) )
where
  reads (r_prim), writes (r_flux)
do
  for i in r_prim do
    var flux : double[5] = get_yfluxes( r_prim[i].rho ,
                                        r_prim[i].u   ,
                                        r_prim[i].v   ,
                                        r_prim[i].w   ,
                                        r_prim[i].p   )
    
    r_flux[i].rho  = flux[0]
    r_flux[i].rhou = flux[1]
    r_flux[i].rhov = flux[2]
    r_flux[i].rhow = flux[3]
    r_flux[i].rhoE = flux[4]
  end
end

__demand(__inline)
task get_zfluxes_r( r_prim : region(ispace(int3d), primitive),
                    r_flux : region(ispace(int3d), conserved) )
where
  reads (r_prim), writes (r_flux)
do
  for i in r_prim do
    var flux : double[5] = get_zfluxes( r_prim[i].rho ,
                                        r_prim[i].u   ,
                                        r_prim[i].v   ,
                                        r_prim[i].w   ,
                                        r_prim[i].p   )
    r_flux[i].rho  = flux[0]
    r_flux[i].rhou = flux[1]
    r_flux[i].rhov = flux[2]
    r_flux[i].rhow = flux[3]
    r_flux[i].rhoE = flux[4]
  end
end

terra sign(x : double)
  return 2*([int](x >= 0)) - 1
end

__demand(__inline)
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

    -- var s_star : double = ( r_prim_r_x[i].p - r_prim_l_x[i].p 
    --                         + r_prim_l_x[i].rho * r_prim_l_x[i].u * (s_L - r_prim_l_x[i].u) 
    --                         - r_prim_r_x[i].rho * r_prim_r_x[i].u * (s_R - r_prim_r_x[i].u) ) 
    --                       / ( r_prim_l_x[i].rho * (s_L - r_prim_l_x[i].u) 
    --                         - r_prim_r_x[i].rho * (s_R - r_prim_r_x[i].u))

    var Q_L : double[5] = get_conserved( r_prim_l_x[i].rho, r_prim_l_x[i].u, r_prim_l_x[i].v, r_prim_l_x[i].w, r_prim_l_x[i].p )
    var F_L : double[5] = get_xfluxes( r_prim_l_x[i].rho, r_prim_l_x[i].u, r_prim_l_x[i].v, r_prim_l_x[i].w, r_prim_l_x[i].p )

    var Q_R : double[5] = get_conserved( r_prim_r_x[i].rho, r_prim_r_x[i].u, r_prim_r_x[i].v, r_prim_r_x[i].w, r_prim_r_x[i].p )
    var F_R : double[5] = get_xfluxes( r_prim_r_x[i].rho, r_prim_r_x[i].u, r_prim_r_x[i].v, r_prim_r_x[i].w, r_prim_r_x[i].p )

    -- var chi_star_L : double = ( s_L - r_prim_l_x[i].u ) / ( s_L - s_star )
    -- var chi_star_R : double = ( s_R - r_prim_r_x[i].u ) / ( s_R - s_star )

    -- var Q_star_L : double[5]
    -- Q_star_L[0] = chi_star_L * ( r_prim_l_x[i].rho )
    -- Q_star_L[1] = chi_star_L * ( r_prim_l_x[i].rho * s_star )
    -- Q_star_L[2] = chi_star_L * ( Q_L[2] )
    -- Q_star_L[3] = chi_star_L * ( Q_L[3] )
    -- Q_star_L[4] = chi_star_L * ( Q_L[4] + (s_star - r_prim_l_x[i].u)*( r_prim_l_x[i].rho * s_star + r_prim_l_x[i].p/(s_L - r_prim_l_x[i].u) ) )

    -- var Q_star_R : double[5]
    -- Q_star_R[0] = chi_star_R * ( r_prim_r_x[i].rho )
    -- Q_star_R[1] = chi_star_R * ( r_prim_r_x[i].rho * s_star )
    -- Q_star_R[2] = chi_star_R * ( Q_R[2] )
    -- Q_star_R[3] = chi_star_R * ( Q_R[3] )
    -- Q_star_R[4] = chi_star_R * ( Q_R[4] + (s_star - r_prim_r_x[i].u)*( r_prim_r_x[i].rho * s_star + r_prim_r_x[i].p/(s_R - r_prim_r_x[i].u) ) )

    -- HLLC Riemann Fluxes
    -- var switch = 0.5*(1 + sign(s_star))
    -- r_flux_e_x[i].rho  = switch * (F_L[0] + s_m*(Q_star_L[0] - Q_L[0])) + (1-switch) * (F_R[0] + s_p*(Q_star_R[0] - Q_R[0]))
    -- r_flux_e_x[i].rhou = switch * (F_L[1] + s_m*(Q_star_L[1] - Q_L[1])) + (1-switch) * (F_R[1] + s_p*(Q_star_R[1] - Q_R[1]))
    -- r_flux_e_x[i].rhov = switch * (F_L[2] + s_m*(Q_star_L[2] - Q_L[2])) + (1-switch) * (F_R[2] + s_p*(Q_star_R[2] - Q_R[2]))
    -- r_flux_e_x[i].rhow = switch * (F_L[3] + s_m*(Q_star_L[3] - Q_L[3])) + (1-switch) * (F_R[3] + s_p*(Q_star_R[3] - Q_R[3]))
    -- r_flux_e_x[i].rhoE = switch * (F_L[4] + s_m*(Q_star_L[4] - Q_L[4])) + (1-switch) * (F_R[4] + s_p*(Q_star_R[4] - Q_R[4]))

    -- HLL Riemann Fluxes
    var switch_L = 0.5*(1 + sign(s_L))
    var switch_R = 0.5*(1 - sign(s_R))
    r_flux_e_x[i].rho  = switch_L * F_L[0] + switch_R * F_R[0] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[0] - s_L*F_R[0] + s_R*s_L * ( Q_R[0] - Q_L[0] ) ) / (s_R - s_L)
    r_flux_e_x[i].rhou = switch_L * F_L[1] + switch_R * F_R[1] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[1] - s_L*F_R[1] + s_R*s_L * ( Q_R[1] - Q_L[1] ) ) / (s_R - s_L)
    r_flux_e_x[i].rhov = switch_L * F_L[2] + switch_R * F_R[2] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[2] - s_L*F_R[2] + s_R*s_L * ( Q_R[2] - Q_L[2] ) ) / (s_R - s_L)
    r_flux_e_x[i].rhow = switch_L * F_L[3] + switch_R * F_R[3] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[3] - s_L*F_R[3] + s_R*s_L * ( Q_R[3] - Q_L[3] ) ) / (s_R - s_L)
    r_flux_e_x[i].rhoE = switch_L * F_L[4] + switch_R * F_R[4] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[4] - s_L*F_R[4] + s_R*s_L * ( Q_R[4] - Q_L[4] ) ) / (s_R - s_L)
  end
end

__demand(__inline)
task HLLC_y( r_prim_l_y : region(ispace(int3d), primitive),
             r_prim_r_y : region(ispace(int3d), primitive),
             r_flux_e_y : region(ispace(int3d), conserved) )
where
  reads(r_prim_l_y, r_prim_r_y), writes(r_flux_e_y)
do
  for i in r_flux_e_y do
    var v_avg   : double = 0.5*( r_prim_l_y[i].v + r_prim_r_y[i].v )
    var sos_l   : double = get_sos(r_prim_l_y[i].rho, r_prim_l_y[i].p)
    var sos_r   : double = get_sos(r_prim_r_y[i].rho, r_prim_r_y[i].p)
    var sos_avg : double = 0.5*( sos_l + sos_r )

    var s_L : double = cmath.fmin(v_avg - sos_avg, r_prim_l_y[i].v - sos_l)
    var s_R : double = cmath.fmax(v_avg + sos_avg, r_prim_r_y[i].v + sos_r)

    var s_m : double = cmath.fmin(0.0, s_L)
    var s_p : double = cmath.fmax(0.0, s_R)

    -- var s_star : double = ( r_prim_r_y[i].p - r_prim_l_y[i].p 
    --                         + r_prim_l_y[i].rho * r_prim_l_y[i].v * (s_L - r_prim_l_y[i].v) 
    --                         - r_prim_r_y[i].rho * r_prim_r_y[i].v * (s_R - r_prim_r_y[i].v) ) 
    --                       / ( r_prim_l_y[i].rho * (s_L - r_prim_l_y[i].v) 
    --                         - r_prim_r_y[i].rho * (s_R - r_prim_r_y[i].v))

    var Q_L : double[5] = get_conserved( r_prim_l_y[i].rho, r_prim_l_y[i].u, r_prim_l_y[i].v, r_prim_l_y[i].w, r_prim_l_y[i].p )
    var F_L : double[5] = get_yfluxes( r_prim_l_y[i].rho, r_prim_l_y[i].u, r_prim_l_y[i].v, r_prim_l_y[i].w, r_prim_l_y[i].p )

    var Q_R : double[5] = get_conserved( r_prim_r_y[i].rho, r_prim_r_y[i].u, r_prim_r_y[i].v, r_prim_r_y[i].w, r_prim_r_y[i].p )
    var F_R : double[5] = get_yfluxes( r_prim_r_y[i].rho, r_prim_r_y[i].u, r_prim_r_y[i].v, r_prim_r_y[i].w, r_prim_r_y[i].p )

    -- var chi_star_L : double = ( s_L - r_prim_l_y[i].v ) / ( s_L - s_star )
    -- var chi_star_R : double = ( s_R - r_prim_r_y[i].v ) / ( s_R - s_star )

    -- var Q_star_L : double[5]
    -- Q_star_L[0] = chi_star_L * ( r_prim_l_y[i].rho )
    -- Q_star_L[1] = chi_star_L * ( Q_L[1] )
    -- Q_star_L[2] = chi_star_L * ( r_prim_l_y[i].rho * s_star ) 
    -- Q_star_L[3] = chi_star_L * ( Q_L[3] )
    -- Q_star_L[4] = chi_star_L * ( Q_L[4] + (s_star - r_prim_l_y[i].v)*( r_prim_l_y[i].rho * s_star + r_prim_l_y[i].p/(s_L - r_prim_l_y[i].v) ) )

    -- var Q_star_R : double[5]
    -- Q_star_R[0] = chi_star_R * ( r_prim_r_y[i].rho )
    -- Q_star_R[1] = chi_star_R * ( Q_R[1] )
    -- Q_star_R[2] = chi_star_R * ( r_prim_r_y[i].rho * s_star )
    -- Q_star_R[3] = chi_star_R * ( Q_R[3] )
    -- Q_star_R[4] = chi_star_R * ( Q_R[4] + (s_star - r_prim_r_y[i].v)*( r_prim_r_y[i].rho * s_star + r_prim_r_y[i].p/(s_R - r_prim_r_y[i].v) ) )

    -- HLLC Riemann Fluxes
    -- var switch = 0.5*(1 + sign(s_star))
    -- r_flux_e_y[i].rho  = switch * (F_L[0] + s_m*(Q_star_L[0] - Q_L[0])) + (1-switch) * (F_R[0] + s_p*(Q_star_R[0] - Q_R[0]))
    -- r_flux_e_y[i].rhou = switch * (F_L[1] + s_m*(Q_star_L[1] - Q_L[1])) + (1-switch) * (F_R[1] + s_p*(Q_star_R[1] - Q_R[1]))
    -- r_flux_e_y[i].rhov = switch * (F_L[2] + s_m*(Q_star_L[2] - Q_L[2])) + (1-switch) * (F_R[2] + s_p*(Q_star_R[2] - Q_R[2]))
    -- r_flux_e_y[i].rhow = switch * (F_L[3] + s_m*(Q_star_L[3] - Q_L[3])) + (1-switch) * (F_R[3] + s_p*(Q_star_R[3] - Q_R[3]))
    -- r_flux_e_y[i].rhoE = switch * (F_L[4] + s_m*(Q_star_L[4] - Q_L[4])) + (1-switch) * (F_R[4] + s_p*(Q_star_R[4] - Q_R[4]))

    -- HLL Riemann Fluxes
    var switch_L = 0.5*(1 + sign(s_L))
    var switch_R = 0.5*(1 - sign(s_R))
    r_flux_e_y[i].rho  = switch_L * F_L[0] + switch_R * F_R[0] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[0] - s_L*F_R[0] + s_R*s_L * ( Q_R[0] - Q_L[0] ) ) / (s_R - s_L)
    r_flux_e_y[i].rhou = switch_L * F_L[1] + switch_R * F_R[1] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[1] - s_L*F_R[1] + s_R*s_L * ( Q_R[1] - Q_L[1] ) ) / (s_R - s_L)
    r_flux_e_y[i].rhov = switch_L * F_L[2] + switch_R * F_R[2] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[2] - s_L*F_R[2] + s_R*s_L * ( Q_R[2] - Q_L[2] ) ) / (s_R - s_L)
    r_flux_e_y[i].rhow = switch_L * F_L[3] + switch_R * F_R[3] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[3] - s_L*F_R[3] + s_R*s_L * ( Q_R[3] - Q_L[3] ) ) / (s_R - s_L)
    r_flux_e_y[i].rhoE = switch_L * F_L[4] + switch_R * F_R[4] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[4] - s_L*F_R[4] + s_R*s_L * ( Q_R[4] - Q_L[4] ) ) / (s_R - s_L)
  end
end

__demand(__inline)
task HLLC_z( r_prim_l_z : region(ispace(int3d), primitive),
             r_prim_r_z : region(ispace(int3d), primitive),
             r_flux_e_z : region(ispace(int3d), conserved) )
where
  reads(r_prim_l_z, r_prim_r_z), writes(r_flux_e_z)
do
  for i in r_flux_e_z do
    var w_avg   : double = 0.5*( r_prim_l_z[i].w + r_prim_r_z[i].w )
    var sos_l   : double = get_sos(r_prim_l_z[i].rho, r_prim_l_z[i].p)
    var sos_r   : double = get_sos(r_prim_r_z[i].rho, r_prim_r_z[i].p)
    var sos_avg : double = 0.5*( sos_l + sos_r )

    var s_L : double = cmath.fmin(w_avg - sos_avg, r_prim_l_z[i].w - sos_l)
    var s_R : double = cmath.fmax(w_avg + sos_avg, r_prim_r_z[i].w + sos_r)

    var s_m : double = cmath.fmin(0.0, s_L)
    var s_p : double = cmath.fmax(0.0, s_R)

    -- var s_star : double = ( r_prim_r_z[i].p - r_prim_l_z[i].p 
    --                         + r_prim_l_z[i].rho * r_prim_l_z[i].w * (s_L - r_prim_l_z[i].w) 
    --                         - r_prim_r_z[i].rho * r_prim_r_z[i].w * (s_R - r_prim_r_z[i].w) ) 
    --                       / ( r_prim_l_z[i].rho * (s_L - r_prim_l_z[i].w) 
    --                         - r_prim_r_z[i].rho * (s_R - r_prim_r_z[i].w))

    var Q_L : double[5] = get_conserved( r_prim_l_z[i].rho, r_prim_l_z[i].u, r_prim_l_z[i].v, r_prim_l_z[i].w, r_prim_l_z[i].p )
    var F_L : double[5] = get_zfluxes( r_prim_l_z[i].rho, r_prim_l_z[i].u, r_prim_l_z[i].v, r_prim_l_z[i].w, r_prim_l_z[i].p )

    var Q_R : double[5] = get_conserved( r_prim_r_z[i].rho, r_prim_r_z[i].u, r_prim_r_z[i].v, r_prim_r_z[i].w, r_prim_r_z[i].p )
    var F_R : double[5] = get_zfluxes( r_prim_r_z[i].rho, r_prim_r_z[i].u, r_prim_r_z[i].v, r_prim_r_z[i].w, r_prim_r_z[i].p )

    -- var chi_star_L : double = ( s_L - r_prim_l_z[i].w ) / ( s_L - s_star )
    -- var chi_star_R : double = ( s_R - r_prim_r_z[i].w ) / ( s_R - s_star )

    -- var Q_star_L : double[5]
    -- Q_star_L[0] = chi_star_L * ( r_prim_l_z[i].rho )
    -- Q_star_L[1] = chi_star_L * ( Q_L[1] )
    -- Q_star_L[2] = chi_star_L * ( Q_L[2] ) 
    -- Q_star_L[3] = chi_star_L * ( r_prim_l_z[i].rho * s_star )
    -- Q_star_L[4] = chi_star_L * ( Q_L[4] + (s_star - r_prim_l_z[i].w)*( r_prim_l_z[i].rho * s_star + r_prim_l_z[i].p/(s_L - r_prim_l_z[i].w) ) )

    -- var Q_star_R : double[5]
    -- Q_star_R[0] = chi_star_R * ( r_prim_r_z[i].rho )
    -- Q_star_R[1] = chi_star_R * ( Q_R[1] )
    -- Q_star_R[2] = chi_star_R * ( Q_R[2] )
    -- Q_star_R[3] = chi_star_R * ( r_prim_r_z[i].rho * s_star ) 
    -- Q_star_R[4] = chi_star_R * ( Q_R[4] + (s_star - r_prim_r_z[i].w)*( r_prim_r_z[i].rho * s_star + r_prim_r_z[i].p/(s_R - r_prim_r_z[i].w) ) )

    -- HLLC Riemann Fluxes
    -- var switch = 0.5*(1 + sign(s_star))
    -- r_flux_e_z[i].rho  = switch * (F_L[0] + s_m*(Q_star_L[0] - Q_L[0])) + (1-switch) * (F_R[0] + s_p*(Q_star_R[0] - Q_R[0]))
    -- r_flux_e_z[i].rhou = switch * (F_L[1] + s_m*(Q_star_L[1] - Q_L[1])) + (1-switch) * (F_R[1] + s_p*(Q_star_R[1] - Q_R[1]))
    -- r_flux_e_z[i].rhov = switch * (F_L[2] + s_m*(Q_star_L[2] - Q_L[2])) + (1-switch) * (F_R[2] + s_p*(Q_star_R[2] - Q_R[2]))
    -- r_flux_e_z[i].rhow = switch * (F_L[3] + s_m*(Q_star_L[3] - Q_L[3])) + (1-switch) * (F_R[3] + s_p*(Q_star_R[3] - Q_R[3]))
    -- r_flux_e_z[i].rhoE = switch * (F_L[4] + s_m*(Q_star_L[4] - Q_L[4])) + (1-switch) * (F_R[4] + s_p*(Q_star_R[4] - Q_R[4]))

    -- HLL Riemann Fluxes
    var switch_L = 0.5*(1 + sign(s_L))
    var switch_R = 0.5*(1 - sign(s_R))
    r_flux_e_z[i].rho  = switch_L * F_L[0] + switch_R * F_R[0] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[0] - s_L*F_R[0] + s_R*s_L * ( Q_R[0] - Q_L[0] ) ) / (s_R - s_L)
    r_flux_e_z[i].rhou = switch_L * F_L[1] + switch_R * F_R[1] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[1] - s_L*F_R[1] + s_R*s_L * ( Q_R[1] - Q_L[1] ) ) / (s_R - s_L)
    r_flux_e_z[i].rhov = switch_L * F_L[2] + switch_R * F_R[2] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[2] - s_L*F_R[2] + s_R*s_L * ( Q_R[2] - Q_L[2] ) ) / (s_R - s_L)
    r_flux_e_z[i].rhow = switch_L * F_L[3] + switch_R * F_R[3] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[3] - s_L*F_R[3] + s_R*s_L * ( Q_R[3] - Q_L[3] ) ) / (s_R - s_L)
    r_flux_e_z[i].rhoE = switch_L * F_L[4] + switch_R * F_R[4] + (1 - switch_L)*(1 - switch_R) * ( s_R*F_L[4] - s_L*F_R[4] + s_R*s_L * ( Q_R[4] - Q_L[4] ) ) / (s_R - s_L)
  end
end

__demand(__inline)
task positivity_enforcer_x( r_prim_c   : region(ispace(int3d), primitive),
                            r_prim_l_x : region(ispace(int3d), primitive),
                            r_prim_r_x : region(ispace(int3d), primitive),
                            n_ghosts   : int64 )
where
  reads(r_prim_c), reads writes(r_prim_l_x, r_prim_r_x)
do
  var counter : int64 = 0
  for i in r_prim_l_x do
    
    if (r_prim_l_x[i].rho <= 0) or (r_prim_r_x[i].rho <= 0) or (r_prim_l_x[i].p <= 0) or (r_prim_r_x[i].p <= 0) then
      var idxm1 = int3d { x = i.x + n_ghosts - 1, y = i.y, z = i.z}
      var idxp1 = int3d { x = i.x + n_ghosts + 0, y = i.y, z = i.z}

      r_prim_l_x[i].rho = r_prim_c[idxm1].rho
      r_prim_l_x[i].u   = r_prim_c[idxm1].u
      r_prim_l_x[i].v   = r_prim_c[idxm1].v
      r_prim_l_x[i].w   = r_prim_c[idxm1].w
      r_prim_l_x[i].p   = r_prim_c[idxm1].p

      r_prim_r_x[i].rho = r_prim_c[idxp1].rho
      r_prim_r_x[i].u   = r_prim_c[idxp1].u
      r_prim_r_x[i].v   = r_prim_c[idxp1].v
      r_prim_r_x[i].w   = r_prim_c[idxp1].w
      r_prim_r_x[i].p   = r_prim_c[idxp1].p

      counter += 1
    end

  end

  if counter > 0 then
    c.printf("WARNING: Positivity enforcer was used in X %d times!\n", counter)
  end
end

__demand(__inline)
task positivity_enforcer_y( r_prim_c   : region(ispace(int3d), primitive),
                            r_prim_l_y : region(ispace(int3d), primitive),
                            r_prim_r_y : region(ispace(int3d), primitive),
                            n_ghosts   : int64 )
where
  reads(r_prim_c), reads writes(r_prim_l_y, r_prim_r_y)
do
  var counter : int64 = 0
  for i in r_prim_l_y do
    
    if (r_prim_l_y[i].rho <= 0) or (r_prim_r_y[i].rho <= 0) or (r_prim_l_y[i].p <= 0) or (r_prim_r_y[i].p <= 0) then
      var idxm1 = int3d { x = i.x, y = i.y + n_ghosts - 1, z = i.z}
      var idxp1 = int3d { x = i.x, y = i.y + n_ghosts + 0, z = i.z}

      r_prim_l_y[i].rho = r_prim_c[idxm1].rho
      r_prim_l_y[i].u   = r_prim_c[idxm1].u
      r_prim_l_y[i].v   = r_prim_c[idxm1].v
      r_prim_l_y[i].w   = r_prim_c[idxm1].w
      r_prim_l_y[i].p   = r_prim_c[idxm1].p

      r_prim_r_y[i].rho = r_prim_c[idxp1].rho
      r_prim_r_y[i].u   = r_prim_c[idxp1].u
      r_prim_r_y[i].v   = r_prim_c[idxp1].v
      r_prim_r_y[i].w   = r_prim_c[idxp1].w
      r_prim_r_y[i].p   = r_prim_c[idxp1].p

      counter += 1
    end

  end

  if counter > 0 then
    c.printf("WARNING: Positivity enforcer was used in Y %d times!\n", counter)
  end
end

__demand(__inline)
task positivity_enforcer_z( r_prim_c   : region(ispace(int3d), primitive),
                            r_prim_l_z : region(ispace(int3d), primitive),
                            r_prim_r_z : region(ispace(int3d), primitive),
                            n_ghosts   : int64 )
where
  reads(r_prim_c), reads writes(r_prim_l_z, r_prim_r_z)
do
  var counter : int64 = 0
  for i in r_prim_l_z do
    
    if (r_prim_l_z[i].rho <= 0) or (r_prim_r_z[i].rho <= 0) or (r_prim_l_z[i].p <= 0) or (r_prim_r_z[i].p <= 0) then
      var idxm1 = int3d { x = i.x, y = i.y, z = i.z + n_ghosts - 1}
      var idxp1 = int3d { x = i.x, y = i.y, z = i.z + n_ghosts + 0}

      r_prim_l_z[i].rho = r_prim_c[idxm1].rho
      r_prim_l_z[i].u   = r_prim_c[idxm1].u
      r_prim_l_z[i].v   = r_prim_c[idxm1].v
      r_prim_l_z[i].w   = r_prim_c[idxm1].w
      r_prim_l_z[i].p   = r_prim_c[idxm1].p

      r_prim_r_z[i].rho = r_prim_c[idxp1].rho
      r_prim_r_z[i].u   = r_prim_c[idxp1].u
      r_prim_r_z[i].v   = r_prim_c[idxp1].v
      r_prim_r_z[i].w   = r_prim_c[idxp1].w
      r_prim_r_z[i].p   = r_prim_c[idxp1].p

      counter += 1
    end

  end

  if counter > 0 then
    c.printf("WARNING: Positivity enforcer was used in Z %d times!\n", counter)
  end
end
