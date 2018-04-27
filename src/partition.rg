import "regent"

local c     = regentlib.c

require("fields")

task partition_LU( LU     : region(ispace(int3d), LU_struct),
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

local function make_partition2D(r)
  local task partition2D( [r],
                          pencil : ispace(int2d) )
    var coloring = c.legion_domain_point_coloring_create()
  
    var prow = pencil.bounds.hi.x + 1
    var pcol = pencil.bounds.hi.y + 1
  
    for i in pencil do
      var lo = int2d { x = i.x, y = i.y }
      var hi = int2d { x = i.x, y = i.y }
      var rect = rect2d { lo = lo, hi = hi }
      c.legion_domain_point_coloring_color_domain(coloring, i, rect)
    end
    var p = partition(disjoint, [r], coloring, pencil)
    c.legion_domain_point_coloring_destroy(coloring)
    return p
  end
  return partition2D
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

local coords   = regentlib.newsymbol(region(ispace(int3d), coordinates), "coords")
partition_xpencil_coords = make_partition_xpencil(coords)
partition_ypencil_coords = make_partition_ypencil(coords)
partition_zpencil_coords = make_partition_zpencil(coords)

local r_prim   = regentlib.newsymbol(region(ispace(int3d), primitive), "r_prim")
partition_xpencil_prim = make_partition_xpencil(r_prim)
partition_ypencil_prim = make_partition_ypencil(r_prim)
partition_zpencil_prim = make_partition_zpencil(r_prim)

local r_cnsr   = regentlib.newsymbol(region(ispace(int3d), conserved), "r_cnsr")
partition_xpencil_cnsr = make_partition_xpencil(r_cnsr)
partition_ypencil_cnsr = make_partition_ypencil(r_cnsr)
partition_zpencil_cnsr = make_partition_zpencil(r_cnsr)

local r_aux   = regentlib.newsymbol(region(ispace(int3d), auxiliary), "r_aux")
partition_xpencil_aux = make_partition_xpencil(r_aux)
partition_ypencil_aux = make_partition_ypencil(r_aux)
partition_zpencil_aux = make_partition_zpencil(r_aux)

local r_visc   = regentlib.newsymbol(region(ispace(int3d), transport_coeffs), "r_visc")
partition_xpencil_visc = make_partition_xpencil(r_visc)
partition_ypencil_visc = make_partition_ypencil(r_visc)
partition_zpencil_visc = make_partition_zpencil(r_visc)

local r_vect   = regentlib.newsymbol(region(ispace(int3d), vect), "r_vect")
partition_xpencil_vect = make_partition_xpencil(r_vect)
partition_ypencil_vect = make_partition_ypencil(r_vect)
partition_zpencil_vect = make_partition_zpencil(r_vect)

local r_tnsr2  = regentlib.newsymbol(region(ispace(int3d), tensor2), "r_tnsr2")
partition_xpencil_tnsr2 = make_partition_xpencil(r_tnsr2)
partition_ypencil_tnsr2 = make_partition_ypencil(r_tnsr2)
partition_zpencil_tnsr2 = make_partition_zpencil(r_tnsr2)

local r_tnsr2symm  = regentlib.newsymbol(region(ispace(int3d), tensor2symm), "r_tnsr2symm")
partition_xpencil_tnsr2symm = make_partition_xpencil(r_tnsr2symm)
partition_ypencil_tnsr2symm = make_partition_ypencil(r_tnsr2symm)
partition_zpencil_tnsr2symm = make_partition_zpencil(r_tnsr2symm)

local r_coeffs  = regentlib.newsymbol(region(ispace(int3d), coeffs), "r_coeffs")
partition_xpencil_coeffs = make_partition_xpencil(r_coeffs)
partition_ypencil_coeffs = make_partition_ypencil(r_coeffs)
partition_zpencil_coeffs = make_partition_zpencil(r_coeffs)

local r_double  = regentlib.newsymbol(region(ispace(int3d), double), "r_double")
partition_xpencil_double = make_partition_xpencil(r_double)
partition_ypencil_double = make_partition_ypencil(r_double)
partition_zpencil_double = make_partition_zpencil(r_double)

local r_double9  = regentlib.newsymbol(region(ispace(int3d), double[9]), "r_double9")
partition_xpencil_double9 = make_partition_xpencil(r_double9)
partition_ypencil_double9 = make_partition_ypencil(r_double9)
partition_zpencil_double9 = make_partition_zpencil(r_double9)

