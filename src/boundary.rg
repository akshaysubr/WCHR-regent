import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")

require("fields")

task periodic_ghost_cells_x( r_prim_c   : region(ispace(int3d), primitive),
                             n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Nx_g = bounds_c.hi.x + 1
  var Nx   = Nx_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z+1 do
    for j = bounds_c.lo.y, bounds_c.hi.y+1 do
      for i = 0,n_ghosts do
        var ghost_l = int3d {i,j,k}
        var int_r   = int3d {Nx+i,j,k}
        r_prim_c[ghost_l].rho = r_prim_c[int_r].rho
        r_prim_c[ghost_l].u   = r_prim_c[int_r].u
        r_prim_c[ghost_l].v   = r_prim_c[int_r].v
        r_prim_c[ghost_l].w   = r_prim_c[int_r].w
        r_prim_c[ghost_l].p   = r_prim_c[int_r].p

        var ghost_r = int3d {Nx+n_ghosts+i,j,k}
        var int_l   = int3d {n_ghosts+i,j,k}
        r_prim_c[ghost_r].rho = r_prim_c[int_l].rho
        r_prim_c[ghost_r].u   = r_prim_c[int_l].u
        r_prim_c[ghost_r].v   = r_prim_c[int_l].v
        r_prim_c[ghost_r].w   = r_prim_c[int_l].w
        r_prim_c[ghost_r].p   = r_prim_c[int_l].p
      end
    end
  end

  return 1
end

task periodic_ghost_cells_y( r_prim_c   : region(ispace(int3d), primitive),
                             n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Ny_g = bounds_c.hi.y + 1
  var Ny   = Ny_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z+1 do
    for i = bounds_c.lo.x, bounds_c.hi.x+1 do
      for j = 0,n_ghosts do
        var ghost_l = int3d {i,j,k}
        var int_r   = int3d {i,Ny+j,k}
        r_prim_c[ghost_l].rho = r_prim_c[int_r].rho
        r_prim_c[ghost_l].u   = r_prim_c[int_r].u
        r_prim_c[ghost_l].v   = r_prim_c[int_r].v
        r_prim_c[ghost_l].w   = r_prim_c[int_r].w
        r_prim_c[ghost_l].p   = r_prim_c[int_r].p

        var ghost_r = int3d {i,Ny+n_ghosts+j,k}
        var int_l   = int3d {i,n_ghosts+j,k}
        r_prim_c[ghost_r].rho = r_prim_c[int_l].rho
        r_prim_c[ghost_r].u   = r_prim_c[int_l].u
        r_prim_c[ghost_r].v   = r_prim_c[int_l].v
        r_prim_c[ghost_r].w   = r_prim_c[int_l].w
        r_prim_c[ghost_r].p   = r_prim_c[int_l].p
      end
    end
  end

  return 1
end

task periodic_ghost_cells_z( r_prim_c   : region(ispace(int3d), primitive),
                             n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Nz_g = bounds_c.hi.z + 1
  var Nz   = Nz_g - 2*n_ghosts

  for j = bounds_c.lo.y, bounds_c.hi.y+1 do
    for i = bounds_c.lo.x, bounds_c.hi.x+1 do
      for k = 0,n_ghosts do
        var ghost_l = int3d {i,j,k}
        var int_r   = int3d {i,j,Nz+k}
        r_prim_c[ghost_l].rho = r_prim_c[int_r].rho
        r_prim_c[ghost_l].u   = r_prim_c[int_r].u
        r_prim_c[ghost_l].v   = r_prim_c[int_r].v
        r_prim_c[ghost_l].w   = r_prim_c[int_r].w
        r_prim_c[ghost_l].p   = r_prim_c[int_r].p

        var ghost_r = int3d {i,j,Nz+n_ghosts+k}
        var int_l   = int3d {i,j,n_ghosts+k}
        r_prim_c[ghost_r].rho = r_prim_c[int_l].rho
        r_prim_c[ghost_r].u   = r_prim_c[int_l].u
        r_prim_c[ghost_r].v   = r_prim_c[int_l].v
        r_prim_c[ghost_r].w   = r_prim_c[int_l].w
        r_prim_c[ghost_r].p   = r_prim_c[int_l].p
      end
    end
  end

  return 1
end

