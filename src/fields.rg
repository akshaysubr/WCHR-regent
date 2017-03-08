import "regent"

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

fspace conserved {
  rho  : double,
  rhou : double,
  rhov : double,
  rhow : double,
  rhoE : double,
}

fspace auxiliary {
  e    : double,
  sos  : double,
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

function poff(i, x, y, z, Nx, Ny, Nz)
  return rexpr int3d { x = (i.x + x + Nx)%Nx, y = (i.y + y + Ny)%Ny, z = (i.z + z + Nz)%Nz } end
end

task min_rho_p( r_prim : region(ispace(int3d), primitive) )
where
  reads (r_prim)
do
  var minrho : double = r_prim[ r_prim.ispace.bounds.lo ].rho
  for i in r_prim do
    var rho : double = r_prim[i].rho
    if ( rho < minrho ) then
      minrho = rho
    end
  end
  return minrho
end

task min_rho_c( r_prim : region(ispace(int3d), conserved) )
where
  reads (r_prim)
do
  var minrho : double = r_prim[ r_prim.ispace.bounds.lo ].rho
  for i in r_prim do
    var rho : double = r_prim[i].rho
    if ( rho < minrho ) then
      minrho = rho
    end
  end
  return minrho
end
