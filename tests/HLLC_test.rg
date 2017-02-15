import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")
require("EOS")

task initialize( r_prim_l_x : region(ispace(int3d), primitive),
                 r_prim_r_x : region(ispace(int3d), primitive) )
where
  reads writes(r_prim_l_x, r_prim_r_x)
do
  for i in r_prim_l_x.ispace do
    r_prim_l_x[i].rho = 1.0
    r_prim_l_x[i].u   = 0.0
    r_prim_l_x[i].v   = 0.0
    r_prim_l_x[i].w   = 0.0
    r_prim_l_x[i].p   = 1.0

    r_prim_r_x[i].rho = 0.125
    r_prim_r_x[i].u   = 0.0
    r_prim_r_x[i].v   = 0.0
    r_prim_r_x[i].w   = 0.0
    r_prim_r_x[i].p   = 0.1
  end

  return 1
end

terra wait_for(x : int)
  return x
end

task main()

  --------------------------------------------------------------------------------------------
  --                       DATA STUCTURES
  --------------------------------------------------------------------------------------------
  var grid_e_x   = ispace(int3d, {x = 1, y = 1, z = 1})  -- x cell edge index space

  var r_prim_l_x = region(grid_e_x, primitive)  -- Primitive variables at left x cell edge
  var r_prim_r_x = region(grid_e_x, primitive)  -- Primitive variables at right x cell edge
  
  var r_flux_e_x   = region(grid_e_x, conserved)  -- Flux at the x edges
  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  var token = initialize(r_prim_l_x, r_prim_r_x)
  wait_for(token)

  c.printf("Left  state: rho = %g, u = %g, v = %g, w = %g, p = %g\n", r_prim_l_x[{0,0,0}].rho, r_prim_l_x[{0,0,0}].u, r_prim_l_x[{0,0,0}].v, r_prim_l_x[{0,0,0}].w, r_prim_l_x[{0,0,0}].p )
  c.printf("Right state: rho = %g, u = %g, v = %g, w = %g, p = %g\n", r_prim_r_x[{0,0,0}].rho, r_prim_r_x[{0,0,0}].u, r_prim_r_x[{0,0,0}].v, r_prim_r_x[{0,0,0}].w, r_prim_r_x[{0,0,0}].p )
  
  HLLC_x(r_prim_l_x, r_prim_r_x, r_flux_e_x)

  c.printf("HLLC  Flux: rho = %g, rhou = %g, rhov = %g, rhow = %g, rhoE = %g\n", r_flux_e_x[{0,0,0}].rho, r_flux_e_x[{0,0,0}].rhou, r_flux_e_x[{0,0,0}].rhov, r_flux_e_x[{0,0,0}].rhow, r_flux_e_x[{0,0,0}].rhoE )
  c.printf("Exact Flux: rho = %g, rhou = %g, rhov = %g, rhow = %g, rhoE = %g\n", 4.3187225579215521e-01, 4.8900185572528088e-01, 0.0, 0.0, 1.1640167368346308e+00 )

  regentlib.assert( cmath.fabs(r_flux_e_x[{0,0,0}].rho  - 4.3187225579215521e-01) < 1.0e-14, "rho  flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_x[{0,0,0}].rhou - 4.8900185572528088e-01) < 1.0e-14, "rhou flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_x[{0,0,0}].rhov - 0.0)                    < 1.0e-14, "rhov flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_x[{0,0,0}].rhow - 0.0)                    < 1.0e-14, "rhow flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_x[{0,0,0}].rhoE - 1.1640167368346308e+00) < 1.0e-14, "rhoE flux does not match" )

end

regentlib.start(main)
