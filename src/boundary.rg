import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")

require("fields")
require("EOS")

local problem = require("problem")

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

task dirchlet_ghost_cells_l_x( r_prim_c   : region(ispace(int3d), primitive),
                               n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Nx_g = bounds_c.hi.x + 1
  var Nx   = Nx_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z + 1 do
    for j = bounds_c.lo.y, bounds_c.hi.y + 1 do
      for i = 0, n_ghosts do
        var ghost_l = int3d {i, j, k}
        r_prim_c[ghost_l].rho = problem.boundary_l_x.rho
        r_prim_c[ghost_l].u   = problem.boundary_l_x.u
        r_prim_c[ghost_l].v   = problem.boundary_l_x.v
        r_prim_c[ghost_l].w   = problem.boundary_l_x.w
        r_prim_c[ghost_l].p   = problem.boundary_l_x.p
      end
    end
  end

  return 1
end

task dirchlet_ghost_cells_r_x( r_prim_c   : region(ispace(int3d), primitive),
                               n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Nx_g = bounds_c.hi.x + 1
  var Nx   = Nx_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z + 1 do
    for j = bounds_c.lo.y, bounds_c.hi.y + 1 do
      for i = 0, n_ghosts do
        var ghost_r = int3d {Nx + n_ghosts + i, j, k}
        r_prim_c[ghost_r].rho = problem.boundary_r_x.rho
        r_prim_c[ghost_r].u   = problem.boundary_r_x.u
        r_prim_c[ghost_r].v   = problem.boundary_r_x.v
        r_prim_c[ghost_r].w   = problem.boundary_r_x.w
        r_prim_c[ghost_r].p   = problem.boundary_r_x.p
      end
    end
  end

  return 1
end

task dirchlet_ghost_cells_l_y( r_prim_c   : region(ispace(int3d), primitive),
                               n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Ny_g = bounds_c.hi.y + 1
  var Ny   = Ny_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z + 1 do
    for j = 0, n_ghosts do
      for i = bounds_c.lo.x, bounds_c.hi.x + 1 do
        var ghost_l = int3d {i, j, k}
        r_prim_c[ghost_l].rho = problem.boundary_l_y.rho
        r_prim_c[ghost_l].u   = problem.boundary_l_y.u
        r_prim_c[ghost_l].v   = problem.boundary_l_y.v  
        r_prim_c[ghost_l].w   = problem.boundary_l_y.w  
        r_prim_c[ghost_l].p   = problem.boundary_l_y.p  
      end
    end
  end

  return 1
end

task dirchlet_ghost_cells_r_y( r_prim_c   : region(ispace(int3d), primitive),
                               n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Ny_g = bounds_c.hi.y + 1
  var Ny   = Ny_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z + 1 do
    for j = 0, n_ghosts do
      for i = bounds_c.lo.x, bounds_c.hi.x + 1 do
        var ghost_r = int3d {i, Ny + n_ghosts + j, k}
        r_prim_c[ghost_r].rho = problem.boundary_r_y.rho
        r_prim_c[ghost_r].u   = problem.boundary_r_y.u
        r_prim_c[ghost_r].v   = problem.boundary_r_y.v  
        r_prim_c[ghost_r].w   = problem.boundary_r_y.w  
        r_prim_c[ghost_r].p   = problem.boundary_r_y.p  
      end
    end
  end

  return 1
end

task dirchlet_ghost_cells_l_z( r_prim_c   : region(ispace(int3d), primitive),
                               n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Nz_g = bounds_c.hi.z + 1
  var Nz   = Nz_g - 2*n_ghosts

  for k = 0, n_ghosts do
    for j = bounds_c.lo.y, bounds_c.hi.y + 1 do
      for i = bounds_c.lo.x, bounds_c.hi.x + 1 do
        var ghost_l = int3d {i, j, k}
        r_prim_c[ghost_l].rho = problem.boundary_l_z.rho
        r_prim_c[ghost_l].u   = problem.boundary_l_z.u
        r_prim_c[ghost_l].v   = problem.boundary_l_z.v   
        r_prim_c[ghost_l].w   = problem.boundary_l_z.w   
        r_prim_c[ghost_l].p   = problem.boundary_l_z.p   
      end
    end
  end

  return 1
end

task dirchlet_ghost_cells_r_z( r_prim_c   : region(ispace(int3d), primitive),
                               n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Nz_g = bounds_c.hi.z + 1
  var Nz   = Nz_g - 2*n_ghosts

  for k = 0, n_ghosts do
    for j = bounds_c.lo.y, bounds_c.hi.y + 1 do
      for i = bounds_c.lo.x, bounds_c.hi.x + 1 do
        var ghost_r = int3d {i, j, Nz + n_ghosts + k}
        r_prim_c[ghost_r].rho = problem.boundary_r_z.rho
        r_prim_c[ghost_r].u   = problem.boundary_r_z.u
        r_prim_c[ghost_r].v   = problem.boundary_r_z.v  
        r_prim_c[ghost_r].w   = problem.boundary_r_z.w  
        r_prim_c[ghost_r].p   = problem.boundary_r_z.p  
      end
    end
  end

  return 1
end

task extrapolation_ghost_cells_l_x( r_prim_c   : region(ispace(int3d), primitive),
                                    n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Nx_g = bounds_c.hi.x + 1
  var Nx   = Nx_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z + 1 do
    for j = bounds_c.lo.y, bounds_c.hi.y + 1 do
      for i = 0, n_ghosts do
        var ghost_l = int3d {i, j, k}
        var int_idx = int3d {n_ghosts, j, k}
        r_prim_c[ghost_l].rho = r_prim_c[int_idx].rho
        r_prim_c[ghost_l].u   = r_prim_c[int_idx].u  
        r_prim_c[ghost_l].v   = r_prim_c[int_idx].v   
        r_prim_c[ghost_l].w   = r_prim_c[int_idx].w   
        r_prim_c[ghost_l].p   = r_prim_c[int_idx].p   
      end
    end
  end

  return 1
end

task extrapolation_ghost_cells_r_x( r_prim_c   : region(ispace(int3d), primitive),
                                    n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Nx_g = bounds_c.hi.x + 1
  var Nx   = Nx_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z + 1 do
    for j = bounds_c.lo.y, bounds_c.hi.y + 1 do
      for i = 0, n_ghosts do
        var ghost_r = int3d {Nx + n_ghosts + i, j, k}
        var int_idx = int3d {Nx + n_ghosts - 1, j, k}
        r_prim_c[ghost_r].rho = r_prim_c[int_idx].rho
        r_prim_c[ghost_r].u   = r_prim_c[int_idx].u  
        r_prim_c[ghost_r].v   = r_prim_c[int_idx].v  
        r_prim_c[ghost_r].w   = r_prim_c[int_idx].w  
        r_prim_c[ghost_r].p   = r_prim_c[int_idx].p  
      end
    end
  end

  return 1
end

task extrapolation_ghost_cells_l_y( r_prim_c   : region(ispace(int3d), primitive),
                                    n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Ny_g = bounds_c.hi.y + 1
  var Ny   = Ny_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z + 1 do
    for j = 0, n_ghosts do
      for i = bounds_c.lo.x, bounds_c.hi.x + 1 do
        var ghost_l = int3d {i, j, k}
        var int_idx = int3d {i, n_ghosts, k}
        r_prim_c[ghost_l].rho = r_prim_c[int_idx].rho
        r_prim_c[ghost_l].u   = r_prim_c[int_idx].u  
        r_prim_c[ghost_l].v   = r_prim_c[int_idx].v   
        r_prim_c[ghost_l].w   = r_prim_c[int_idx].w   
        r_prim_c[ghost_l].p   = r_prim_c[int_idx].p   
      end
    end
  end

  return 1
end

task extrapolation_ghost_cells_r_y( r_prim_c   : region(ispace(int3d), primitive),
                                    n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Ny_g = bounds_c.hi.y + 1
  var Ny   = Ny_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z + 1 do
    for j = 0, n_ghosts do
      for i = bounds_c.lo.x, bounds_c.hi.x + 1 do
        var ghost_r = int3d {i, Ny + n_ghosts + j, k}
        var int_idx = int3d {i, Ny + n_ghosts - 1, k}
        r_prim_c[ghost_r].rho = r_prim_c[int_idx].rho
        r_prim_c[ghost_r].u   = r_prim_c[int_idx].u  
        r_prim_c[ghost_r].v   = r_prim_c[int_idx].v  
        r_prim_c[ghost_r].w   = r_prim_c[int_idx].w  
        r_prim_c[ghost_r].p   = r_prim_c[int_idx].p  
      end
    end
  end

  return 1
end

task extrapolation_ghost_cells_l_z( r_prim_c   : region(ispace(int3d), primitive),
                                    n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Nz_g = bounds_c.hi.z + 1
  var Nz   = Nz_g - 2*n_ghosts

  for k = 0, n_ghosts do
    for j = bounds_c.lo.y, bounds_c.hi.y + 1 do
      for i = bounds_c.lo.x, bounds_c.hi.x + 1 do
        var ghost_l = int3d {i, j, k}
        var int_idx = int3d {i, j, n_ghosts}
        r_prim_c[ghost_l].rho = r_prim_c[int_idx].rho
        r_prim_c[ghost_l].u   = r_prim_c[int_idx].u  
        r_prim_c[ghost_l].v   = r_prim_c[int_idx].v   
        r_prim_c[ghost_l].w   = r_prim_c[int_idx].w   
        r_prim_c[ghost_l].p   = r_prim_c[int_idx].p   
      end
    end
  end

  return 1
end

task extrapolation_ghost_cells_r_z( r_prim_c   : region(ispace(int3d), primitive),
                                    n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var bounds_c = r_prim_c.ispace.bounds
  var Nz_g = bounds_c.hi.z + 1
  var Nz   = Nz_g - 2*n_ghosts

  for k = 0, n_ghosts do
    for j = bounds_c.lo.y, bounds_c.hi.y + 1 do
      for i = bounds_c.lo.x, bounds_c.hi.x + 1 do
        var ghost_r = int3d {i, j, Nz + n_ghosts + k}
        var int_idx = int3d {i, j, Nz + n_ghosts - 1}
        r_prim_c[ghost_r].rho = r_prim_c[int_idx].rho
        r_prim_c[ghost_r].u   = r_prim_c[int_idx].u  
        r_prim_c[ghost_r].v   = r_prim_c[int_idx].v  
        r_prim_c[ghost_r].w   = r_prim_c[int_idx].w  
        r_prim_c[ghost_r].p   = r_prim_c[int_idx].p  
      end
    end
  end

  return 1
end

function make_ddx_left_sided(r_func, f_func, ONEBYDX)
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func, f_func)
  
  local ddx_left_sided __demand(__inline) task ddx_left_sided( [r_func],
                                                               idx : int3d )
  where
    [privileges_r_func]
  do
    var dQdx = 0.5*[ONEBYDX]* (   3.0*[r_func][ int3d { idx.x,     idx.y, idx.z } ].[f_func]
                                - 4.0*[r_func][ int3d { idx.x - 1, idx.y, idx.z } ].[f_func]
                                + 1.0*[r_func][ int3d { idx.x - 2, idx.y, idx.z } ].[f_func] )
    return dQdx
  end
  return ddx_left_sided
end

function make_ddx_right_sided(r_func, f_func, ONEBYDX)
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func, f_func)
  
  local ddx_right_sided __demand(__inline) task ddx_right_sided( [r_func],
                                                                 idx : int3d )
  where
    [privileges_r_func]
  do
    var dQdx = 0.5*[ONEBYDX]* ( - 3.0*[r_func][ int3d { idx.x,     idx.y, idx.z } ].[f_func]
                                + 4.0*[r_func][ int3d { idx.x + 1, idx.y, idx.z } ].[f_func]
                                - 1.0*[r_func][ int3d { idx.x + 2, idx.y, idx.z } ].[f_func] )
    return dQdx
  end
  return ddx_right_sided
end

function make_ddy_left_sided(r_func, f_func, ONEBYDY)
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func, f_func)
  
  local ddy_left_sided __demand(__inline) task ddy_left_sided( [r_func],
                                                               idx : int3d )
  where
    [privileges_r_func]
  do
    var dQdy = 0.5*[ONEBYDY]* (   3.0*[r_func][ int3d { idx.x, idx.y,     idx.z } ].[f_func]
                                - 4.0*[r_func][ int3d { idx.x, idx.y - 1, idx.z } ].[f_func]
                                + 1.0*[r_func][ int3d { idx.x, idx.y - 2, idx.z } ].[f_func] )
    return dQdy
  end
  return ddy_left_sided
end

function make_ddy_right_sided(r_func, f_func, ONEBYDY)
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func, f_func)
  
  local ddy_right_sided __demand(__inline) task ddy_right_sided( [r_func],
                                                                 idx : int3d )
  where
    [privileges_r_func]
  do
    var dQdy = 0.5*[ONEBYDY]* ( - 3.0*[r_func][ int3d { idx.x, idx.y,     idx.z } ].[f_func]
                                + 4.0*[r_func][ int3d { idx.x, idx.y + 1, idx.z } ].[f_func]
                                - 1.0*[r_func][ int3d { idx.x, idx.y + 2, idx.z } ].[f_func] )
    return dQdy
  end
  return ddy_right_sided
end

function make_ddz_left_sided(r_func, f_func, ONEBYDZ)
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func, f_func)
  
  local ddz_left_sided __demand(__inline) task ddz_left_sided( [r_func],
                                                               idx : int3d )
  where
    [privileges_r_func]
  do
    var dQdz = 0.5*[ONEBYDZ]* (   3.0*[r_func][ int3d { idx.x, idx.y, idx.z     } ].[f_func]
                                - 4.0*[r_func][ int3d { idx.x, idx.y, idx.z - 1 } ].[f_func]
                                + 1.0*[r_func][ int3d { idx.x, idx.y, idx.z - 2 } ].[f_func] )
    return dQdz
  end
  return ddz_left_sided
end

function make_ddz_right_sided(r_func, f_func, ONEBYDZ)
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func, f_func)
  
  local ddz_right_sided __demand(__inline) task ddz_right_sided( [r_func],
                                                                 idx : int3d )
  where
    [privileges_r_func]
  do
    var dQdz = 0.5*[ONEBYDZ]* ( - 3.0*[r_func][ int3d { idx.x, idx.y, idx.z     } ].[f_func]
                                + 4.0*[r_func][ int3d { idx.x, idx.y, idx.z + 1 } ].[f_func]
                                - 1.0*[r_func][ int3d { idx.x, idx.y, idx.z + 2 } ].[f_func] )
    return dQdz
  end
  return ddz_right_sided
end

function make_extrapolate_l_x(r_func, f_func, DX)
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func, f_func)
  
  local extrapolate_l_x __demand(__inline) task extrapolate_l_x( [r_func],
                                                                 idx      : int3d,
                                                                 dQdx     : double )
  where
    [privileges_r_func]
  do
    var Q_g : double[4]

    Q_g[3] = [r_func][ int3d { idx.x + 1, idx.y, idx.z } ].[f_func]
             - 2.0*DX*dQdx

    Q_g[2] = - 2.0*[r_func][ int3d { idx.x + 1, idx.y, idx.z } ].[f_func]
             - 3.0*[r_func][ int3d { idx.x,     idx.y, idx.z } ].[f_func]
             + 6.0*Q_g[3]
             + 6.0*DX*dQdx

    Q_g[1] =    3.0*[r_func][ int3d { idx.x + 1, idx.y, idx.z } ].[f_func]
             + 10.0*[r_func][ int3d { idx.x,     idx.y, idx.z } ].[f_func]
             - 18.0*Q_g[3]
             +  6.0*Q_g[2]
             - 12.0*DX*dQdx

    Q_g[0] = -      4.0*[r_func][ int3d { idx.x + 1, idx.y, idx.z } ].[f_func]
             - 65.0/3.0*[r_func][ int3d { idx.x,     idx.y, idx.z } ].[f_func]
             +     40.0*Q_g[3]
             -     20.0*Q_g[2]
             + 20.0/3.0*Q_g[1]
             + 20.0*DX*dQdx

    return Q_g
  end
  return extrapolate_l_x
end

function make_extrapolate_r_x(r_func, f_func, DX)
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func, f_func)
  
  local extrapolate_r_x __demand(__inline) task extrapolate_r_x( [r_func],
                                                                 idx      : int3d,
                                                                 dQdx     : double )
  where
    [privileges_r_func]
  do
    var Q_g : double[4]

    Q_g[0] = [r_func][ int3d { idx.x - 1, idx.y, idx.z } ].[f_func]
             + 2.0*DX*dQdx

    Q_g[1] = - 2.0*[r_func][ int3d { idx.x - 1, idx.y, idx.z } ].[f_func]
             - 3.0*[r_func][ int3d { idx.x,     idx.y, idx.z } ].[f_func]
             + 6.0*Q_g[0]
             - 6.0*DX*dQdx

    Q_g[2] =    3.0*[r_func][ int3d { idx.x - 1, idx.y, idx.z } ].[f_func]
             + 10.0*[r_func][ int3d { idx.x,     idx.y, idx.z } ].[f_func]
             - 18.0*Q_g[0]
             +  6.0*Q_g[1]
             + 12.0*DX*dQdx

    Q_g[3] = -      4.0*[r_func][ int3d { idx.x - 1, idx.y, idx.z } ].[f_func]
             - 65.0/3.0*[r_func][ int3d { idx.x,     idx.y, idx.z } ].[f_func]
             +     40.0*Q_g[0]
             -     20.0*Q_g[1]
             + 20.0/3.0*Q_g[2]
             - 20.0*DX*dQdx

    return Q_g
  end
  return extrapolate_r_x
end

function make_extrapolate_l_y(r_func, f_func, DY)
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func, f_func)
  
  local extrapolate_l_y __demand(__inline) task extrapolate_l_y( [r_func],
                                                                 idx      : int3d,
                                                                 dQdy     : double )
  where
    [privileges_r_func]
  do
    var Q_g : double[4]

    Q_g[3] = [r_func][ int3d { idx.x, idx.y + 1, idx.z } ].[f_func]
             - 2.0*DY*dQdy

    Q_g[2] = - 2.0*[r_func][ int3d { idx.x, idx.y + 1, idx.z } ].[f_func]
             - 3.0*[r_func][ int3d { idx.x, idx.y    , idx.z } ].[f_func]
             + 6.0*Q_g[3]
             + 6.0*DY*dQdy

    Q_g[1] =    3.0*[r_func][ int3d { idx.x, idx.y + 1, idx.z } ].[f_func]
             + 10.0*[r_func][ int3d { idx.x, idx.y    , idx.z } ].[f_func]
             - 18.0*Q_g[3]
             +  6.0*Q_g[2]
             - 12.0*DY*dQdy

    Q_g[0] = -      4.0*[r_func][ int3d { idx.x, idx.y + 1, idx.z } ].[f_func]
             - 65.0/3.0*[r_func][ int3d { idx.x, idx.y    , idx.z } ].[f_func]
             +     40.0*Q_g[3]
             -     20.0*Q_g[2]
             + 20.0/3.0*Q_g[1]
             + 20.0*DY*dQdy

    return Q_g
  end
  return extrapolate_l_y
end

function make_extrapolate_r_y(r_func, f_func, DY)
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func, f_func)
  
  local extrapolate_r_y __demand(__inline) task extrapolate_r_y( [r_func],
                                                                 idx      : int3d,
                                                                 dQdy     : double )
  where
    [privileges_r_func]
  do
    var Q_g : double[4]

    Q_g[0] = [r_func][ int3d { idx.x, idx.y - 1, idx.z } ].[f_func]
             + 2.0*DY*dQdy

    Q_g[1] = - 2.0*[r_func][ int3d { idx.x, idx.y - 1, idx.z } ].[f_func]
             - 3.0*[r_func][ int3d { idx.x, idx.y    , idx.z } ].[f_func]
             + 6.0*Q_g[0]
             - 6.0*DY*dQdy

    Q_g[2] =    3.0*[r_func][ int3d { idx.x, idx.y - 1, idx.z } ].[f_func]
             + 10.0*[r_func][ int3d { idx.x, idx.y    , idx.z } ].[f_func]
             - 18.0*Q_g[0]
             +  6.0*Q_g[1]
             + 12.0*DY*dQdy

    Q_g[3] = -      4.0*[r_func][ int3d { idx.x, idx.y - 1, idx.z } ].[f_func]
             - 65.0/3.0*[r_func][ int3d { idx.x, idx.y    , idx.z } ].[f_func]
             +     40.0*Q_g[0]
             -     20.0*Q_g[1]
             + 20.0/3.0*Q_g[2]
             - 20.0*DY*dQdy

    return Q_g
  end
  return extrapolate_r_y
end

function make_extrapolate_l_z(r_func, f_func, DZ)
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func, f_func)
  
  local extrapolate_l_z __demand(__inline) task extrapolate_l_z( [r_func],
                                                                 idx      : int3d,
                                                                 dQdz     : double )
  where
    [privileges_r_func]
  do
    var Q_g : double[4]

    Q_g[3] = [r_func][ int3d { idx.x, idx.y, idx.z + 1 } ].[f_func]
             - 2.0*DZ*dQdz

    Q_g[2] = - 2.0*[r_func][ int3d { idx.x, idx.y, idx.z + 1 } ].[f_func]
             - 3.0*[r_func][ int3d { idx.x, idx.y, idx.z     } ].[f_func]
             + 6.0*Q_g[3]
             + 6.0*DZ*dQdz

    Q_g[1] =    3.0*[r_func][ int3d { idx.x, idx.y, idx.z + 1 } ].[f_func]
             + 10.0*[r_func][ int3d { idx.x, idx.y, idx.z     } ].[f_func]
             - 18.0*Q_g[3]
             +  6.0*Q_g[2]
             - 12.0*DZ*dQdz

    Q_g[0] = -      4.0*[r_func][ int3d { idx.x, idx.y, idx.z + 1 } ].[f_func]
             - 65.0/3.0*[r_func][ int3d { idx.x, idx.y, idx.z     } ].[f_func]
             +     40.0*Q_g[3]
             -     20.0*Q_g[2]
             + 20.0/3.0*Q_g[1]
             + 20.0*DZ*dQdz

    return Q_g
  end
  return extrapolate_l_z
end

function make_extrapolate_r_z(r_func, f_func, DZ)
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func, f_func)
  
  local extrapolate_r_z __demand(__inline) task extrapolate_r_z( [r_func],
                                                                 idx      : int3d,
                                                                 dQdz     : double )
  where
    [privileges_r_func]
  do
    var Q_g : double[4]

    Q_g[0] = [r_func][ int3d { idx.x, idx.y, idx.z - 1 } ].[f_func]
             + 2.0*DZ*dQdz

    Q_g[1] = - 2.0*[r_func][ int3d { idx.x, idx.y, idx.z - 1 } ].[f_func]
             - 3.0*[r_func][ int3d { idx.x, idx.y, idx.z     } ].[f_func]
             + 6.0*Q_g[0]
             - 6.0*DZ*dQdz

    Q_g[2] =    3.0*[r_func][ int3d { idx.x, idx.y, idx.z - 1 } ].[f_func]
             + 10.0*[r_func][ int3d { idx.x, idx.y, idx.z     } ].[f_func]
             - 18.0*Q_g[0]
             +  6.0*Q_g[1]
             + 12.0*DZ*dQdz

    Q_g[3] = -      4.0*[r_func][ int3d { idx.x, idx.y, idx.z - 1 } ].[f_func]
             - 65.0/3.0*[r_func][ int3d { idx.x, idx.y, idx.z     } ].[f_func]
             +     40.0*Q_g[0]
             -     20.0*Q_g[1]
             + 20.0/3.0*Q_g[2]
             - 20.0*DZ*dQdz

    return Q_g
  end
  return extrapolate_r_z
end



local r_prim = regentlib.newsymbol(region(ispace(int3d), primitive), "r_primitive")

local ddx_left_sided_rho  = make_ddx_left_sided(r_prim,  "rho", problem.ONEBYDX)
local ddx_right_sided_rho = make_ddx_right_sided(r_prim, "rho", problem.ONEBYDX)
local ddx_left_sided_u    = make_ddx_left_sided(r_prim,  "u",   problem.ONEBYDX)
local ddx_right_sided_u   = make_ddx_right_sided(r_prim, "u",   problem.ONEBYDX)
local ddx_left_sided_v    = make_ddx_left_sided(r_prim,  "v",   problem.ONEBYDX)
local ddx_right_sided_v   = make_ddx_right_sided(r_prim, "v",   problem.ONEBYDX)
local ddx_left_sided_w    = make_ddx_left_sided(r_prim,  "w",   problem.ONEBYDX)
local ddx_right_sided_w   = make_ddx_right_sided(r_prim, "w",   problem.ONEBYDX)
local ddx_left_sided_p    = make_ddx_left_sided(r_prim,  "p",   problem.ONEBYDX)
local ddx_right_sided_p   = make_ddx_right_sided(r_prim, "p",   problem.ONEBYDX)

local ddy_left_sided_rho  = make_ddy_left_sided(r_prim,  "rho", problem.ONEBYDY)
local ddy_right_sided_rho = make_ddy_right_sided(r_prim, "rho", problem.ONEBYDY)
local ddy_left_sided_u    = make_ddy_left_sided(r_prim,  "u",   problem.ONEBYDY)
local ddy_right_sided_u   = make_ddy_right_sided(r_prim, "u",   problem.ONEBYDY)
local ddy_left_sided_v    = make_ddy_left_sided(r_prim,  "v",   problem.ONEBYDY)
local ddy_right_sided_v   = make_ddy_right_sided(r_prim, "v",   problem.ONEBYDY)
local ddy_left_sided_w    = make_ddy_left_sided(r_prim,  "w",   problem.ONEBYDY)
local ddy_right_sided_w   = make_ddy_right_sided(r_prim, "w",   problem.ONEBYDY)
local ddy_left_sided_p    = make_ddy_left_sided(r_prim,  "p",   problem.ONEBYDY)
local ddy_right_sided_p   = make_ddy_right_sided(r_prim, "p",   problem.ONEBYDY)

local ddz_left_sided_rho  = make_ddz_left_sided(r_prim,  "rho", problem.ONEBYDZ)
local ddz_right_sided_rho = make_ddz_right_sided(r_prim, "rho", problem.ONEBYDZ)
local ddz_left_sided_u    = make_ddz_left_sided(r_prim,  "u",   problem.ONEBYDZ)
local ddz_right_sided_u   = make_ddz_right_sided(r_prim, "u",   problem.ONEBYDZ)
local ddz_left_sided_v    = make_ddz_left_sided(r_prim,  "v",   problem.ONEBYDZ)
local ddz_right_sided_v   = make_ddz_right_sided(r_prim, "v",   problem.ONEBYDZ)
local ddz_left_sided_w    = make_ddz_left_sided(r_prim,  "w",   problem.ONEBYDZ)
local ddz_right_sided_w   = make_ddz_right_sided(r_prim, "w",   problem.ONEBYDZ)
local ddz_left_sided_p    = make_ddz_left_sided(r_prim,  "p",   problem.ONEBYDZ)
local ddz_right_sided_p   = make_ddz_right_sided(r_prim, "p",   problem.ONEBYDZ)

local extrapolate_l_x_rho = make_extrapolate_l_x(r_prim, "rho", problem.DX)
local extrapolate_r_x_rho = make_extrapolate_r_x(r_prim, "rho", problem.DX)
local extrapolate_l_x_u   = make_extrapolate_l_x(r_prim, "u",   problem.DX)
local extrapolate_r_x_u   = make_extrapolate_r_x(r_prim, "u",   problem.DX)
local extrapolate_l_x_v   = make_extrapolate_l_x(r_prim, "v",   problem.DX)
local extrapolate_r_x_v   = make_extrapolate_r_x(r_prim, "v",   problem.DX)
local extrapolate_l_x_w   = make_extrapolate_l_x(r_prim, "w",   problem.DX)
local extrapolate_r_x_w   = make_extrapolate_r_x(r_prim, "w",   problem.DX)
local extrapolate_l_x_p   = make_extrapolate_l_x(r_prim, "p",   problem.DX)
local extrapolate_r_x_p   = make_extrapolate_r_x(r_prim, "p",   problem.DX)

local extrapolate_l_y_rho = make_extrapolate_l_y(r_prim, "rho", problem.DY)
local extrapolate_r_y_rho = make_extrapolate_r_y(r_prim, "rho", problem.DY)
local extrapolate_l_y_u   = make_extrapolate_l_y(r_prim, "u",   problem.DY)
local extrapolate_r_y_u   = make_extrapolate_r_y(r_prim, "u",   problem.DY)
local extrapolate_l_y_v   = make_extrapolate_l_y(r_prim, "v",   problem.DY)
local extrapolate_r_y_v   = make_extrapolate_r_y(r_prim, "v",   problem.DY)
local extrapolate_l_y_w   = make_extrapolate_l_y(r_prim, "w",   problem.DY)
local extrapolate_r_y_w   = make_extrapolate_r_y(r_prim, "w",   problem.DY)
local extrapolate_l_y_p   = make_extrapolate_l_y(r_prim, "p",   problem.DY)
local extrapolate_r_y_p   = make_extrapolate_r_y(r_prim, "p",   problem.DY)

local extrapolate_l_z_rho = make_extrapolate_l_z(r_prim, "rho", problem.DZ)
local extrapolate_r_z_rho = make_extrapolate_r_z(r_prim, "rho", problem.DZ)
local extrapolate_l_z_u   = make_extrapolate_l_z(r_prim, "u",   problem.DZ)
local extrapolate_r_z_u   = make_extrapolate_r_z(r_prim, "u",   problem.DZ)
local extrapolate_l_z_v   = make_extrapolate_l_z(r_prim, "v",   problem.DZ)
local extrapolate_r_z_v   = make_extrapolate_r_z(r_prim, "v",   problem.DZ)
local extrapolate_l_z_w   = make_extrapolate_l_z(r_prim, "w",   problem.DZ)
local extrapolate_r_z_w   = make_extrapolate_r_z(r_prim, "w",   problem.DZ)
local extrapolate_l_z_p   = make_extrapolate_l_z(r_prim, "p",   problem.DZ)
local extrapolate_r_z_p   = make_extrapolate_r_z(r_prim, "p",   problem.DZ)

task subsonic_inflow_ghost_cells_l_x( r_prim_c   : region(ispace(int3d), primitive),
                                      n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do

  var rho_inflow = problem.boundary_l_x.rho
  var u_inflow = problem.boundary_l_x.u
  var v_inflow = problem.boundary_l_x.v
  var w_inflow = problem.boundary_l_x.w
  var p_inflow = problem.boundary_l_x.p
  
  var T_inflow = get_temperature( rho_inflow, p_inflow )
  var RT_inflow = problem.Rgas*T_inflow

  var l_x = problem.boundary_l_x.L_x
  var eta_2 = problem.boundary_l_x.eta_2
  var eta_3 = problem.boundary_l_x.eta_3
  var eta_4 = problem.boundary_l_x.eta_4
  var eta_5 = problem.boundary_l_x.eta_5

  var bounds_c = r_prim_c.ispace.bounds
  var Nx_g = bounds_c.hi.x + 1
  var Nx   = Nx_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z + 1 do
    for j = bounds_c.lo.y, bounds_c.hi.y + 1 do
      var int_idx = int3d {n_ghosts, j, k}

      var rho_x = 0.0
      var u_x   = 0.0
      var v_x   = 0.0
      var w_x   = 0.0
      var p_x   = 0.0

      u_x = ddx_right_sided_u( r_prim_c, int_idx )
      p_x = ddx_right_sided_p( r_prim_c, int_idx )
      var sos = get_sos( r_prim_c[int_idx].rho, r_prim_c[int_idx].p )
      var Mach = r_prim_c[int_idx].u / sos

      var T_2 = 0.0
      var T_3 = 0.0
      var T_4 = 0.0
      var T_5 = 0.0

      var L_1 = (p_x - r_prim_c[int_idx].rho*sos*u_x)
      var L_2 = eta_2 * (r_prim_c[int_idx].rho*sos/l_x) * (r_prim_c[int_idx].p/r_prim_c[int_idx].rho - RT_inflow) - T_2
      var L_3 = eta_3 * sos/l_x*( r_prim_c[int_idx].v - v_inflow ) - T_3
      var L_4 = eta_4 * sos/l_x*( r_prim_c[int_idx].w - w_inflow ) - T_4
      var L_5 = eta_5 * (r_prim_c[int_idx].rho*sos*sos*(1-Mach*Mach)/l_x) * (r_prim_c[int_idx].u - u_inflow) - T_5

      L_2 = L_2 / (r_prim_c[int_idx].u)
      L_3 = L_3 / (r_prim_c[int_idx].u)
      L_4 = L_4 / (r_prim_c[int_idx].u)
      L_5 = L_5 / (r_prim_c[int_idx].u + sos)

      rho_x = 1./(sos*sos) * (0.5*L_1 + L_2 + 0.5*L_5)
      u_x = (1./(2.*r_prim_c[int_idx].rho*sos)) * (-L_1 + L_5)
      v_x = L_3
      w_x = L_4
      p_x = 0.5*(L_1 + L_5)

      var rho_g_L = extrapolate_l_x_rho( r_prim_c, int_idx, rho_x )
      var u_g_L   = extrapolate_l_x_u  ( r_prim_c, int_idx, u_x   )
      var v_g_L   = extrapolate_l_x_v  ( r_prim_c, int_idx, v_x   )
      var w_g_L   = extrapolate_l_x_w  ( r_prim_c, int_idx, w_x   )
      var p_g_L   = extrapolate_l_x_p  ( r_prim_c, int_idx, p_x   )

      for i = 0, n_ghosts do
        var ghost_l = int3d {i, j, k}
        r_prim_c[ghost_l].rho = rho_g_L[i]
        r_prim_c[ghost_l].u   = u_g_L[i]  
        r_prim_c[ghost_l].v   = v_g_L[i]   
        r_prim_c[ghost_l].w   = w_g_L[i]   
        r_prim_c[ghost_l].p   = p_g_L[i]   
      end
    end
  end

  return 1
end

task subsonic_inflow_ghost_cells_r_x( r_prim_c   : region(ispace(int3d), primitive),
                                      n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var rho_inflow = problem.boundary_r_x.rho
  var u_inflow   = problem.boundary_r_x.u
  var v_inflow   = problem.boundary_r_x.v
  var w_inflow   = problem.boundary_r_x.w
  var p_inflow   = problem.boundary_r_x.p
  
  var T_inflow = get_temperature( rho_inflow, p_inflow )
  var RT_inflow = problem.Rgas*T_inflow

  var l_x = problem.boundary_r_x.L_x
  var eta_1 = problem.boundary_r_x.eta_1
  var eta_2 = problem.boundary_r_x.eta_2
  var eta_3 = problem.boundary_r_x.eta_3
  var eta_4 = problem.boundary_r_x.eta_4

  var bounds_c = r_prim_c.ispace.bounds
  var Nx_g = bounds_c.hi.x + 1
  var Nx   = Nx_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z + 1 do
    for j = bounds_c.lo.y, bounds_c.hi.y + 1 do
      var int_idx = int3d {Nx + n_ghosts - 1, j, k}

      var rho_x = 0.0
      var u_x   = 0.0
      var v_x   = 0.0
      var w_x   = 0.0
      var p_x   = 0.0

      u_x = ddx_left_sided_u( r_prim_c, int_idx )
      p_x = ddx_left_sided_p( r_prim_c, int_idx )
      var sos = get_sos( r_prim_c[int_idx].rho, r_prim_c[int_idx].p )
      var Mach = r_prim_c[int_idx].u / sos

      var T_1 = 0.0
      var T_2 = 0.0
      var T_3 = 0.0
      var T_4 = 0.0

      var L_1 = eta_1 * (r_prim_c[int_idx].rho*sos*sos*(1-Mach*Mach)/l_x) * (r_prim_c[int_idx].u - u_inflow) - T_1
      var L_2 = eta_2 * (r_prim_c[int_idx].rho*sos/l_x) * (r_prim_c[int_idx].p/r_prim_c[int_idx].rho - RT_inflow) - T_2
      var L_3 = eta_3 * sos/l_x*( r_prim_c[int_idx].v - v_inflow ) - T_3
      var L_4 = eta_4 * sos/l_x*( r_prim_c[int_idx].w - w_inflow ) - T_4
      var L_5 = (p_x + r_prim_c[int_idx].rho*sos*u_x)

      L_1 = L_1 / (r_prim_c[int_idx].u - sos)
      L_2 = L_2 / (r_prim_c[int_idx].u)
      L_3 = L_3 / (r_prim_c[int_idx].u)
      L_4 = L_4 / (r_prim_c[int_idx].u)

      rho_x = 1./(sos*sos) * (0.5*L_1 + L_2 + 0.5*L_5)
      u_x = (1./(2.*r_prim_c[int_idx].rho*sos)) * (-L_1 + L_5)
      v_x = L_3
      w_x = L_4
      p_x = 0.5*(L_1 + L_5)

      var rho_g_R = extrapolate_r_x_rho( r_prim_c, int_idx, rho_x )
      var u_g_R   = extrapolate_r_x_u  ( r_prim_c, int_idx, u_x   )
      var v_g_R   = extrapolate_r_x_v  ( r_prim_c, int_idx, v_x   )
      var w_g_R   = extrapolate_r_x_w  ( r_prim_c, int_idx, w_x   )
      var p_g_R   = extrapolate_r_x_p  ( r_prim_c, int_idx, p_x   )

      for i = 0, n_ghosts do
        var ghost_r = int3d {Nx + n_ghosts + i, j, k}
        r_prim_c[ghost_r].rho = rho_g_R[i]
        r_prim_c[ghost_r].u   = u_g_R[i]
        r_prim_c[ghost_r].v   = v_g_R[i]
        r_prim_c[ghost_r].w   = w_g_R[i]
        r_prim_c[ghost_r].p   = p_g_R[i]
      end
    end
  end

  return 1
end

task subsonic_outflow_ghost_cells_l_x( r_prim_c   : region(ispace(int3d), primitive),
                                       n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var rho_outflow = problem.boundary_l_x.rho
  var u_outflow   = problem.boundary_l_x.u
  var v_outflow   = problem.boundary_l_x.v
  var w_outflow   = problem.boundary_l_x.w
  var p_outflow   = problem.boundary_l_x.p
  
  var l_x = problem.boundary_l_x.L_x
  var sigma = problem.boundary_l_x.sigma

  var bounds_c = r_prim_c.ispace.bounds
  var Nx_g = bounds_c.hi.x + 1
  var Nx   = Nx_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z + 1 do
    for j = bounds_c.lo.y, bounds_c.hi.y + 1 do
      var int_idx = int3d {n_ghosts, j, k}

      var rho_x = 0.0
      var u_x   = 0.0
      var v_x   = 0.0
      var w_x   = 0.0
      var p_x   = 0.0

      rho_x = ddx_right_sided_rho( r_prim_c, int_idx )
      u_x   = ddx_right_sided_u( r_prim_c, int_idx )
      v_x   = ddx_right_sided_v( r_prim_c, int_idx )
      w_x   = ddx_right_sided_w( r_prim_c, int_idx )
      p_x   = ddx_right_sided_p( r_prim_c, int_idx )

      var sos = get_sos( r_prim_c[int_idx].rho, r_prim_c[int_idx].p )
      var Mach = r_prim_c[int_idx].u / sos
      var K = sigma * sos * (1.0 - Mach*Mach) / l_x

      var T_5 = 0.0

      var L_1 = (p_x - r_prim_c[int_idx].rho*sos*u_x)
      var L_2 = (sos*sos * rho_x - p_x)
      var L_3 = v_x
      var L_4 = w_x
      var L_5 = K * (r_prim_c[int_idx].p - p_outflow) - T_5

      L_5 = L_5 / (r_prim_c[int_idx].u + sos)

      rho_x = 1./(sos*sos) * (0.5*L_1 + L_2 + 0.5*L_5)
      u_x = (1./(2.*r_prim_c[int_idx].rho*sos)) * (-L_1 + L_5)
      v_x = L_3
      w_x = L_4
      p_x = 0.5*(L_1 + L_5)

      var rho_g_L = extrapolate_l_x_rho( r_prim_c, int_idx, rho_x )
      var u_g_L   = extrapolate_l_x_u  ( r_prim_c, int_idx, u_x   )
      var v_g_L   = extrapolate_l_x_v  ( r_prim_c, int_idx, v_x   )
      var w_g_L   = extrapolate_l_x_w  ( r_prim_c, int_idx, w_x   )
      var p_g_L   = extrapolate_l_x_p  ( r_prim_c, int_idx, p_x   )

      for i = 0, n_ghosts do
        var ghost_l = int3d {i, j, k}
        r_prim_c[ghost_l].rho = rho_g_L[i]
        r_prim_c[ghost_l].u   = u_g_L[i]
        r_prim_c[ghost_l].v   = v_g_L[i]
        r_prim_c[ghost_l].w   = w_g_L[i]
        r_prim_c[ghost_l].p   = p_g_L[i]
      end
    end
  end

  return 1
end

task subsonic_outflow_ghost_cells_r_x( r_prim_c   : region(ispace(int3d), primitive),
                                       n_ghosts   : int64 )
where
  reads writes(r_prim_c)
do
  var rho_outflow = problem.boundary_r_x.rho
  var u_outflow   = problem.boundary_r_x.u
  var v_outflow   = problem.boundary_r_x.v
  var w_outflow   = problem.boundary_r_x.w
  var p_outflow   = problem.boundary_r_x.p
  
  var l_x = problem.boundary_r_x.L_x
  var sigma = problem.boundary_r_x.sigma

  var bounds_c = r_prim_c.ispace.bounds
  var Nx_g = bounds_c.hi.x + 1
  var Nx   = Nx_g - 2*n_ghosts

  for k = bounds_c.lo.z, bounds_c.hi.z + 1 do
    for j = bounds_c.lo.y, bounds_c.hi.y + 1 do
      var int_idx = int3d {Nx + n_ghosts - 1, j, k}

      var rho_x = 0.0
      var u_x   = 0.0
      var v_x   = 0.0
      var w_x   = 0.0
      var p_x   = 0.0

      rho_x = ddx_left_sided_rho( r_prim_c, int_idx )
      u_x   = ddx_left_sided_u( r_prim_c, int_idx )
      v_x   = ddx_left_sided_v( r_prim_c, int_idx )
      w_x   = ddx_left_sided_w( r_prim_c, int_idx )
      p_x   = ddx_left_sided_p( r_prim_c, int_idx )

      var sos = get_sos( r_prim_c[int_idx].rho, r_prim_c[int_idx].p )
      var Mach = r_prim_c[int_idx].u / sos
      var K = sigma * sos * (1.0 - Mach*Mach) / l_x

      var T_1 = 0.0

      var L_1 = K * (r_prim_c[int_idx].p - p_outflow) - T_1
      var L_2 = (sos*sos * rho_x - p_x)
      var L_3 = v_x
      var L_4 = w_x
      var L_5 = (p_x + r_prim_c[int_idx].rho*sos*u_x)

      L_1 = L_1 / (r_prim_c[int_idx].u - sos)

      rho_x = 1./(sos*sos) * (0.5*L_1 + L_2 + 0.5*L_5)
      u_x = (1./(2.*r_prim_c[int_idx].rho*sos)) * (-L_1 + L_5)
      v_x = L_3
      w_x = L_4
      p_x = 0.5*(L_1 + L_5)

      var rho_g_R = extrapolate_r_x_rho( r_prim_c, int_idx, rho_x )
      var u_g_R   = extrapolate_r_x_u  ( r_prim_c, int_idx, u_x   )
      var v_g_R   = extrapolate_r_x_v  ( r_prim_c, int_idx, v_x   )
      var w_g_R   = extrapolate_r_x_w  ( r_prim_c, int_idx, w_x   )
      var p_g_R   = extrapolate_r_x_p  ( r_prim_c, int_idx, p_x   )

      for i = 0, n_ghosts do
        var ghost_r = int3d {Nx + n_ghosts + i, j, k}
        r_prim_c[ghost_r].rho = rho_g_R[i]
        r_prim_c[ghost_r].u   = u_g_R[i]
        r_prim_c[ghost_r].v   = v_g_R[i]
        r_prim_c[ghost_r].w   = w_g_R[i]
        r_prim_c[ghost_r].p   = p_g_R[i]
      end
    end
  end

  return 1
end



task nonperiodic_ghost_cells_x( r_prim_c  : region(ispace(int3d), primitive),
                                n_ghosts  : int64 )
where
  reads writes(r_prim_c)
do
  -- Left ghost cells
  if (problem.boundary_l_x.condition == "DIRICHLET") then
    dirchlet_ghost_cells_l_x( r_prim_c, n_ghosts )
  elseif (problem.boundary_l_x.condition == "EXTRAPOLATION") then
    extrapolation_ghost_cells_l_x( r_prim_c, n_ghosts )
  elseif (problem.boundary_l_x.condition == "SUBSONIC_INFLOW") then
    subsonic_inflow_ghost_cells_l_x( r_prim_c, n_ghosts )
  elseif (problem.boundary_l_x.condition == "SUBSONIC_OUTFLOW") then
    subsonic_outflow_ghost_cells_l_x( r_prim_c, n_ghosts )
  else
    regentlib.assert(false, "Unknown left boundary condition in the x-direction")
  end

  -- Right ghost cells
  if (problem.boundary_r_x.condition == "DIRICHLET") then
    dirchlet_ghost_cells_r_x( r_prim_c, n_ghosts )
  elseif (problem.boundary_r_x.condition == "EXTRAPOLATION") then
    extrapolation_ghost_cells_r_x( r_prim_c, n_ghosts )
  elseif (problem.boundary_r_x.condition == "SUBSONIC_INFLOW") then
    subsonic_inflow_ghost_cells_r_x( r_prim_c, n_ghosts )
  elseif (problem.boundary_r_x.condition == "SUBSONIC_OUTFLOW") then
    subsonic_outflow_ghost_cells_r_x( r_prim_c, n_ghosts )
  else
    regentlib.assert(false, "Unknown right boundary condition in the x-direction")
  end
end

task nonperiodic_ghost_cells_y( r_prim_c  : region(ispace(int3d), primitive),
                                n_ghosts  : int64 )
where
  reads writes(r_prim_c)
do
  -- Left ghost cells
  if (problem.boundary_l_y.condition == "DIRICHLET") then
    dirchlet_ghost_cells_l_y( r_prim_c, n_ghosts )
  elseif (problem.boundary_l_y.condition == "EXTRAPOLATION") then
    extrapolation_ghost_cells_l_y( r_prim_c, n_ghosts )
  else
    regentlib.assert(false, "Unknown left boundary condition in the y-direction")
  end

  -- Right ghost cells
  if (problem.boundary_r_y.condition == "DIRICHLET") then
    dirchlet_ghost_cells_r_y( r_prim_c, n_ghosts )
  elseif (problem.boundary_r_y.condition == "EXTRAPOLATION") then
    extrapolation_ghost_cells_r_y( r_prim_c, n_ghosts )
  else
    regentlib.assert(false, "Unknown right boundary condition in the y-direction")
  end
end

task nonperiodic_ghost_cells_z( r_prim_c  : region(ispace(int3d), primitive),
                                n_ghosts  : int64 )
where
  reads writes(r_prim_c)
do
  -- Left ghost cells
  if (problem.boundary_l_z.condition == "DIRICHLET") then
    dirchlet_ghost_cells_l_z( r_prim_c, n_ghosts )
  elseif (problem.boundary_l_z.condition == "EXTRAPOLATION") then
    extrapolation_ghost_cells_l_z( r_prim_c, n_ghosts )
  else
    regentlib.assert(false, "Unknown left boundary condition in the z-direction")
  end

  -- Right ghost cells
  if (problem.boundary_r_z.condition == "DIRICHLET") then
    dirchlet_ghost_cells_r_z( r_prim_c, n_ghosts )
  elseif (problem.boundary_r_z.condition == "EXTRAPOLATION") then
    extrapolation_ghost_cells_r_z( r_prim_c, n_ghosts )
  else
    regentlib.assert(false, "Unknown right boundary condition in the z-direction")
  end
end
