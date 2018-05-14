import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")

require("fields")

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
  else
    regentlib.assert(false, "Unknown left boundary condition in the x-direction")
  end

  -- Right ghost cells
  if (problem.boundary_r_x.condition == "DIRICHLET") then
    dirchlet_ghost_cells_r_x( r_prim_c, n_ghosts )
  elseif (problem.boundary_r_x.condition == "EXTRAPOLATION") then
    extrapolation_ghost_cells_r_x( r_prim_c, n_ghosts )
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
