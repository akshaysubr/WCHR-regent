import "regent"

local c = regentlib.c

-- Field spaces
fspace coordinates {
  x_c : double,
  y_c : double,
  z_c : double,
}

fspace primitive {
  rho : double,
  u   : double,
  v   : double,
  w   : double,
  p   : double,
}

fspace coeffs {
  _0 : double,
  _1 : double,
  _2 : double,
  _3 : double,
  _4 : double,
}

fspace auxiliary {
  T   : double,
}

fspace conserved {
  rho  : double,
  rhou : double,
  rhov : double,
  rhow : double,
  rhoE : double,
}

fspace transport_coeffs {
  mu_s  : double,
  mu_b  : double,
  kappa : double,
}

fspace vect {
  _1  : double,
  _2  : double,
  _3  : double,
}

fspace tensor2 {
  _11 : double,
  _12 : double,
  _13 : double,
  _21 : double,
  _22 : double,
  _23 : double,
  _31 : double,
  _32 : double,
  _33 : double,
}

fspace tensor2symm {
  _11 : double,
  _12 : double,
  _13 : double,
  _22 : double,
  _23 : double,
  _33 : double,
}

fspace fluxvar {
  f       : double,
  dfdx    : double,
  dfdy    : double,
  dfdz    : double,
  f_rhs_x : double,
  f_rhs_y : double,
  f_rhs_z : double,
}

fspace soe_vector {
  _0 : double,
  _1 : double,
  _2 : double,
  _3 : double,
  _4 : double,
}

fspace LU_struct {
  b  : double,
  eg : double,
  k  : double,
  l  : double,
  g  : double,
  h  : double,
  ff : double,
  v  : double,
  w  : double,
}

fspace LU_coeffs {
  e  : double,
  a  : double,
  d  : double,
  c  : double,
  f  : double,
}

function poff(i, x, y, z, Nx, Ny, Nz)
  return rexpr int3d { x = (i.x + x + Nx)%Nx, y = (i.y + y + Ny)%Ny, z = (i.z + z + Nz)%Nz } end
end

task set_zero_cnsr( r_cnsr : region(ispace(int3d), conserved) )
where
  writes(r_cnsr)
do
  for i in r_cnsr do
    r_cnsr[i].{rho, rhou, rhov, rhow, rhoE} = 0.0
  end

  return 1
end

task set_zero_prim( r_prim : region(ispace(int3d), primitive) )
where
  writes(r_prim)
do
  for i in r_prim do
    r_prim[i].{rho, u, v, w, p} = 0.0
  end

  return 1
end

-- task add_value_cnsr( r_cnsr : region(ispace(int3d), conserved),
--                      r_rhs  : region(ispace(int3d), conserved),
--                      coeff  : double )
-- where
--   reads (r_rhs), reads writes(r_cnsr)
-- do
-- 
--   for i in r_rhs do
--     r_cnsr[i].rho  += coeff*r_rhs[i].rho
--     r_cnsr[i].rhou += coeff*r_rhs[i].rhou
--     r_cnsr[i].rhov += coeff*r_rhs[i].rhov
--     r_cnsr[i].rhow += coeff*r_rhs[i].rhow
--     r_cnsr[i].rhoE += coeff*r_rhs[i].rhoE
--   end
-- end
-- 
-- task self_multiply_cnsr( r_cnsr : region(ispace(int3d), conserved),
--                          coeff  : double )
-- where
--   reads writes(r_cnsr)
-- do
-- 
--   for i in r_cnsr do
--     r_cnsr[i].rho  = coeff*r_cnsr[i].rho
--     r_cnsr[i].rhou = coeff*r_cnsr[i].rhou
--     r_cnsr[i].rhov = coeff*r_cnsr[i].rhov
--     r_cnsr[i].rhow = coeff*r_cnsr[i].rhow
--     r_cnsr[i].rhoE = coeff*r_cnsr[i].rhoE
--   end
-- end

-- task min_rho_p( r_prim : region(ispace(int3d), primitive) )
-- where
--   reads (r_prim)
-- do
--   var minrho : double = r_prim[ r_prim.ispace.bounds.lo ].rho
--   for i in r_prim do
--     var rho : double = r_prim[i].rho
--     if ( rho < minrho ) then
--       minrho = rho
--     end
--   end
--   return minrho
-- end
-- 
-- task max_rho_p( r_prim : region(ispace(int3d), primitive) )
-- where
--   reads (r_prim)
-- do
--   var maxrho : double = r_prim[ r_prim.ispace.bounds.lo ].rho
--   for i in r_prim do
--     var rho : double = r_prim[i].rho
--     if ( rho > maxrho ) then
--       maxrho = rho
--     end
--   end
--   return maxrho
-- end
-- 
-- task min_p_p( r_prim : region(ispace(int3d), primitive) )
-- where
--   reads (r_prim)
-- do
--   var minp : double = r_prim[ r_prim.ispace.bounds.lo ].p
--   for i in r_prim do
--     var p : double = r_prim[i].p
--     if ( p < minp ) then
--       minp = p
--     end
--   end
--   return minp
-- end
-- 
-- task max_p_p( r_prim : region(ispace(int3d), primitive) )
-- where
--   reads (r_prim)
-- do
--   var maxp : double = r_prim[ r_prim.ispace.bounds.lo ].p
--   for i in r_prim do
--     var p : double = r_prim[i].p
--     if ( p > maxp ) then
--       maxp = p
--     end
--   end
--   return maxp
-- end
-- 
-- task min_rho_c( r_prim : region(ispace(int3d), conserved) )
-- where
--   reads (r_prim)
-- do
--   var minrho : double = r_prim[ r_prim.ispace.bounds.lo ].rho
--   for i in r_prim do
--     var rho : double = r_prim[i].rho
--     if ( rho < minrho ) then
--       minrho = rho
--     end
--   end
--   return minrho
-- end
