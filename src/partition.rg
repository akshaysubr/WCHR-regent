import "regent"

local c     = regentlib.c

task partitionLU( LU     : region(ispace(int3d), LU_struct),
                  pencil : ispace(int2d) )
  var coloring = c.legion_domain_point_coloring_create()

  var prow = pencil.bounds.hi.x + 1
  var pcol = pencil.bounds.hi.y + 1

  var bounds = LU.ispace.bounds
  var N = bounds.hi.x + 1

  for i in pencil do
    var lo = int3d { x = 0,   y = i.x, z = i.y }
    var hi = int3d { x = N-1, y = i.x, z = i.y }
    var rect = rect3d { lo = lo, hi = hi }
    c.legion_domain_point_coloring_color_domain(coloring, i, rect)
  end
  var p = partition(disjoint, LU, coloring, pencil)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

local function make_partition_xpencil(r)
  local task partition_xpencil( [r],
                                xpencil : ispace(int2d) )
    var coloring = c.legion_domain_point_coloring_create()
  
    var prow = xpencil.bounds.hi.x + 1
    var pcol = xpencil.bounds.hi.y + 1
  
    var bounds = [r].ispace.bounds
    var Nx = bounds.hi.x + 1
    var Ny = bounds.hi.y + 1
    var Nz = bounds.hi.z + 1
  
    for i in xpencil do
      var lo = int3d { x = 0, y = i.x*(Ny/prow), z = i.y*(Nz/pcol) }
      var hi = int3d { x = Nx-1, y = (i.x+1)*(Ny/prow)-1, z = (i.y+1)*(Nz/pcol)-1 }
      var rect = rect3d { lo = lo, hi = hi }
      c.legion_domain_point_coloring_color_domain(coloring, i, rect)
    end
    var p = partition(disjoint, [r], coloring, xpencil)
    c.legion_domain_point_coloring_destroy(coloring)
    return p
  end
  return partition_xpencil
end

local function make_partition_ypencil(r)
  local task partition_ypencil( [r],
                                ypencil : ispace(int2d) )
    var coloring = c.legion_domain_point_coloring_create()
  
    var prow = ypencil.bounds.hi.x + 1
    var pcol = ypencil.bounds.hi.y + 1
  
    var bounds = [r].ispace.bounds
    var Nx = bounds.hi.x + 1
    var Ny = bounds.hi.y + 1
    var Nz = bounds.hi.z + 1
  
    for i in ypencil do
      var lo = int3d { x = i.x*(Nx/prow), y = 0, z = i.y*(Nz/pcol) }
      var hi = int3d { x = (i.x+1)*(Nx/prow)-1, y = Ny-1, z = (i.y+1)*(Nz/pcol)-1 }
      var rect = rect3d { lo = lo, hi = hi }
      c.legion_domain_point_coloring_color_domain(coloring, i, rect)
    end
    var p = partition(disjoint, [r], coloring, ypencil)
    c.legion_domain_point_coloring_destroy(coloring)
    return p
  end
  return partition_ypencil
end

local function make_partition_zpencil(r)
  local task partition_zpencil( [r],
                                zpencil : ispace(int2d) )
    var coloring = c.legion_domain_point_coloring_create()
  
    var prow = zpencil.bounds.hi.x + 1
    var pcol = zpencil.bounds.hi.y + 1
  
    var bounds = [r].ispace.bounds
    var Nx = bounds.hi.x + 1
    var Ny = bounds.hi.y + 1
    var Nz = bounds.hi.z + 1
  
    for i in zpencil do
      var lo = int3d { x = i.x*(Nx/prow), y = i.y*(Ny/pcol), z = 0 }
      var hi = int3d { x = (i.x+1)*(Nx/prow)-1, y = (i.y+1)*(Ny/pcol)-1, z = Nz-1 }
      var rect = rect3d { lo = lo, hi = hi }
      c.legion_domain_point_coloring_color_domain(coloring, i, rect)
    end
    var p = partition(disjoint, [r], coloring, zpencil)
    c.legion_domain_point_coloring_destroy(coloring)
    return p
  end
  return partition_zpencil
end
