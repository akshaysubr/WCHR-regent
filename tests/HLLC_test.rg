import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

require("fields")
require("EOS")

task initialize_x( r_prim_l_x : region(ispace(int3d), primitive),
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

task initialize_y( r_prim_l_y : region(ispace(int3d), primitive),
                 r_prim_r_y : region(ispace(int3d), primitive) )
where
  reads writes(r_prim_l_y, r_prim_r_y)
do
  for i in r_prim_l_y.ispace do
    r_prim_l_y[i].rho = 1.0
    r_prim_l_y[i].u   = 0.0
    r_prim_l_y[i].v   = 0.0
    r_prim_l_y[i].w   = 0.0
    r_prim_l_y[i].p   = 1.0

    r_prim_r_y[i].rho = 0.125
    r_prim_r_y[i].u   = 0.0
    r_prim_r_y[i].v   = 0.0
    r_prim_r_y[i].w   = 0.0
    r_prim_r_y[i].p   = 0.1
  end

  return 1
end

task initialize_z( r_prim_l_z : region(ispace(int3d), primitive),
                 r_prim_r_z : region(ispace(int3d), primitive) )
where
  reads writes(r_prim_l_z, r_prim_r_z)
do
  for i in r_prim_l_z.ispace do
    r_prim_l_z[i].rho = 1.0
    r_prim_l_z[i].u   = 0.0
    r_prim_l_z[i].v   = 0.0
    r_prim_l_z[i].w   = 0.0
    r_prim_l_z[i].p   = 1.0

    r_prim_r_z[i].rho = 0.125
    r_prim_r_z[i].u   = 0.0
    r_prim_r_z[i].v   = 0.0
    r_prim_r_z[i].w   = 0.0
    r_prim_r_z[i].p   = 0.1
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
  var grid_e_y   = ispace(int3d, {x = 1, y = 1, z = 1})  -- y cell edge index space
  var grid_e_z   = ispace(int3d, {x = 1, y = 1, z = 1})  -- z cell edge index space

  var r_prim_l_x = region(grid_e_x, primitive)  -- Primitive variables at left x cell edge
  var r_prim_r_x = region(grid_e_x, primitive)  -- Primitive variables at right x cell edge
  
  var r_prim_l_y = region(grid_e_y, primitive)  -- Primitive variables at left y cell edge
  var r_prim_r_y = region(grid_e_y, primitive)  -- Primitive variables at right y cell edge
  
  var r_prim_l_z = region(grid_e_z, primitive)  -- Primitive variables at left z cell edge
  var r_prim_r_z = region(grid_e_z, primitive)  -- Primitive variables at right z cell edge
  
  var r_flux_e_x   = region(grid_e_x, conserved)  -- Flux at the x edges
  var r_flux_e_y   = region(grid_e_y, conserved)  -- Flux at the y edges
  var r_flux_e_z   = region(grid_e_z, conserved)  -- Flux at the z edges
  --------------------------------------------------------------------------------------------
  --------------------------------------------------------------------------------------------

  var token = initialize_x(r_prim_l_x, r_prim_r_x)
  wait_for(token)

  token = initialize_y(r_prim_l_y, r_prim_r_y)
  wait_for(token)
  
  token = initialize_z(r_prim_l_z, r_prim_r_z)
  wait_for(token)

  c.printf("Left  state: rho = %g, u = %g, v = %g, w = %g, p = %g\n", r_prim_l_x[{0,0,0}].rho, r_prim_l_x[{0,0,0}].u, r_prim_l_x[{0,0,0}].v, r_prim_l_x[{0,0,0}].w, r_prim_l_x[{0,0,0}].p )
  c.printf("Right state: rho = %g, u = %g, v = %g, w = %g, p = %g\n", r_prim_r_x[{0,0,0}].rho, r_prim_r_x[{0,0,0}].u, r_prim_r_x[{0,0,0}].v, r_prim_r_x[{0,0,0}].w, r_prim_r_x[{0,0,0}].p )
  
  HLLC_x(r_prim_l_x, r_prim_r_x, r_flux_e_x)

  c.printf("HLLC  Flux: rho = %g, rhou = %g, rhov = %g, rhow = %g, rhoE = %g\n", r_flux_e_x[{0,0,0}].rho, r_flux_e_x[{0,0,0}].rhou, r_flux_e_x[{0,0,0}].rhov, r_flux_e_x[{0,0,0}].rhow, r_flux_e_x[{0,0,0}].rhoE )
  c.printf("Exact Flux: rho = %g, rhou = %g, rhov = %g, rhow = %g, rhoE = %g\n\n", 4.3187225579215521e-01, 4.8900185572528088e-01, 0.0, 0.0, 1.1640167368346308e+00 )

  regentlib.assert( cmath.fabs(r_flux_e_x[{0,0,0}].rho  - 4.3187225579215521e-01) < 1.0e-14, "rho  flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_x[{0,0,0}].rhou - 4.8900185572528088e-01) < 1.0e-14, "rhou flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_x[{0,0,0}].rhov - 0.0)                    < 1.0e-14, "rhov flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_x[{0,0,0}].rhow - 0.0)                    < 1.0e-14, "rhow flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_x[{0,0,0}].rhoE - 1.1640167368346308e+00) < 1.0e-14, "rhoE flux does not match" )

  c.printf("Bottom  state: rho = %g, u = %g, v = %g, w = %g, p = %g\n", r_prim_l_y[{0,0,0}].rho, r_prim_l_y[{0,0,0}].u, r_prim_l_y[{0,0,0}].v, r_prim_l_y[{0,0,0}].w, r_prim_l_y[{0,0,0}].p )
  c.printf("Top state: rho = %g, u = %g, v = %g, w = %g, p = %g\n", r_prim_r_y[{0,0,0}].rho, r_prim_r_y[{0,0,0}].u, r_prim_r_y[{0,0,0}].v, r_prim_r_y[{0,0,0}].w, r_prim_r_y[{0,0,0}].p )
  
  HLLC_y(r_prim_l_y, r_prim_r_y, r_flux_e_y)

  c.printf("HLLC  Flux: rho = %g, rhou = %g, rhov = %g, rhow = %g, rhoE = %g\n", r_flux_e_y[{0,0,0}].rho, r_flux_e_y[{0,0,0}].rhou, r_flux_e_y[{0,0,0}].rhov, r_flux_e_y[{0,0,0}].rhow, r_flux_e_y[{0,0,0}].rhoE )
  c.printf("Exact Flux: rho = %g, rhou = %g, rhov = %g, rhow = %g, rhoE = %g\n\n", 4.3187225579215521e-01, 0.0, 4.8900185572528088e-01, 0.0, 1.1640167368346308e+00 )

  regentlib.assert( cmath.fabs(r_flux_e_y[{0,0,0}].rho  - 4.3187225579215521e-01) < 1.0e-14, "rho  flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_y[{0,0,0}].rhou - 0.0)                    < 1.0e-14, "rhou flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_y[{0,0,0}].rhov - 4.8900185572528088e-01) < 1.0e-14, "rhov flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_y[{0,0,0}].rhow - 0.0)                    < 1.0e-14, "rhow flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_y[{0,0,0}].rhoE - 1.1640167368346308e+00) < 1.0e-14, "rhoE flux does not match" )

  c.printf("Back  state: rho = %g, u = %g, v = %g, w = %g, p = %g\n", r_prim_l_z[{0,0,0}].rho, r_prim_l_z[{0,0,0}].u, r_prim_l_z[{0,0,0}].v, r_prim_l_z[{0,0,0}].w, r_prim_l_z[{0,0,0}].p )
  c.printf("Front state: rho = %g, u = %g, v = %g, w = %g, p = %g\n", r_prim_r_z[{0,0,0}].rho, r_prim_r_z[{0,0,0}].u, r_prim_r_z[{0,0,0}].v, r_prim_r_z[{0,0,0}].w, r_prim_r_z[{0,0,0}].p )
  
  HLLC_z(r_prim_l_z, r_prim_r_z, r_flux_e_z)

  c.printf("HLLC  Flux: rho = %g, rhou = %g, rhov = %g, rhow = %g, rhoE = %g\n", r_flux_e_z[{0,0,0}].rho, r_flux_e_z[{0,0,0}].rhou, r_flux_e_z[{0,0,0}].rhov, r_flux_e_z[{0,0,0}].rhow, r_flux_e_z[{0,0,0}].rhoE )
  c.printf("Exact Flux: rho = %g, rhou = %g, rhov = %g, rhow = %g, rhoE = %g\n\n", 4.3187225579215521e-01, 0.0, 0.0, 4.8900185572528088e-01, 1.1640167368346308e+00 )

  regentlib.assert( cmath.fabs(r_flux_e_z[{0,0,0}].rho  - 4.3187225579215521e-01) < 1.0e-14, "rho  flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_z[{0,0,0}].rhou - 0.0)                    < 1.0e-14, "rhou flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_z[{0,0,0}].rhov - 0.0)                    < 1.0e-14, "rhov flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_z[{0,0,0}].rhow - 4.8900185572528088e-01) < 1.0e-14, "rhow flux does not match" )
  regentlib.assert( cmath.fabs(r_flux_e_z[{0,0,0}].rhoE - 1.1640167368346308e+00) < 1.0e-14, "rhoE flux does not match" )


end

regentlib.start(main)
