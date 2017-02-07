import "regent"

require("fields")

local gamma = 1.4
local gamma_m1 = gamma - 1.0
local onebygm1 = 1.0 / gamma_m1

task get_primitive( r_cnsr : region(ispace(int3d), conserved),
                    r_prim : region(ispace(int3d), primitive) )
where
  reads(r_cnsr), reads writes(r_prim)
do

  for i in r_prim do
    r_prim[i].rho = r_cnsr[i].rho
    var onebyrho : double = 1.0 / r_cnsr[i].rho
    r_prim[i].u   = r_cnsr[i].rhou * onebyrho
    r_prim[i].v   = r_cnsr[i].rhov * onebyrho
    r_prim[i].w   = r_cnsr[i].rhow * onebyrho
    r_prim[i].p   = gamma_m1 * ( r_cnsr[i].rhoE - 0.5 * r_cnsr[i].rho * (r_prim[i].u*r_prim[i].u + r_prim[i].v*r_prim[i].v + r_prim[i].w*r_prim[i].w) )
  end
end

task get_conserved( r_prim : region(ispace(int3d), primitive),
                    r_cnsr : region(ispace(int3d), conserved) )
where
  reads(r_prim), reads writes(r_cnsr)
do

  get_e_from_p(r_prim)

  for i in r_cnsr do
    r_cnsr[i].rho  = r_prim[i].rho
    r_cnsr[i].rhou = r_prim[i].rho * r_prim[i].u
    r_cnsr[i].rhov = r_prim[i].rho * r_prim[i].v
    r_cnsr[i].rhow = r_prim[i].rho * r_prim[i].w
    r_cnsr[i].rhoE = onebygm1 * r_prim[i].p + 0.5 * r_prim[i].rho * (r_prim[i].u*r_prim[i].u + r_prim[i].v*r_prim[i].v + r_prim[i].w*r_prim[i].w)
  end
end

task get_xfluxes( r_prim : region(ispace(int3d), primitive),
                  r_cnsr : region(ispace(int3d), conserved),
                  r_flux : region(ispace(int3d), conserved) )
where
  reads (r_prim, r_cnsr), reads writes (r_flux)
do
  for i in r_prim do
    r_flux[i].rho  = r_cnsr[i].rhou
    r_flux[i].rhou = r_cnsr[i].rhou * r_prim[i].u + r_prim[i].p
    r_flux[i].rhov = r_cnsr[i].rhou * r_prim[i].v
    r_flux[i].rhow = r_cnsr[i].rhou * r_prim[i].w
    r_flux[i].rhoE =(r_cnsr[i].rhoE + r_prim[i].p) * r_prim[i].u
  end
end

task get_yfluxes( r_prim : region(ispace(int3d), primitive),
                  r_cnsr : region(ispace(int3d), conserved),
                  r_flux : region(ispace(int3d), conserved) )
where
  reads (r_prim, r_cnsr), reads writes (r_flux)
do
  for i in r_prim do
    r_flux[i].rho  = r_cnsr[i].rhov
    r_flux[i].rhou = r_cnsr[i].rhov * r_prim[i].u
    r_flux[i].rhov = r_cnsr[i].rhov * r_prim[i].v + r_prim[i].p
    r_flux[i].rhow = r_cnsr[i].rhov * r_prim[i].w
    r_flux[i].rhoE =(r_cnsr[i].rhoE + r_prim[i].p) * r_prim[i].v
  end
end

task get_zfluxes( r_prim : region(ispace(int3d), primitive),
                  r_cnsr : region(ispace(int3d), conserved),
                  r_flux : region(ispace(int3d), conserved) )
where
  reads (r_prim, r_cnsr), reads writes (r_flux)
do
  for i in r_prim do
    r_flux[i].rho  = r_cnsr[i].rhow
    r_flux[i].rhou = r_cnsr[i].rhow * r_prim[i].u
    r_flux[i].rhov = r_cnsr[i].rhow * r_prim[i].v
    r_flux[i].rhow = r_cnsr[i].rhow * r_prim[i].w + r_prim[i].p
    r_flux[i].rhoE =(r_cnsr[i].rhoE + r_prim[i].p) * r_prim[i].w
  end
end

