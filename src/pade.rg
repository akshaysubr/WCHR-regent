import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

local max = regentlib.fmax

-- Some problem parameters
local NX = 32
local NY = 32
local NZ = 32
local LL = 2.0*math.pi
local DX = LL / NX
local DY = LL / NY
local DZ = LL / NZ
local ONEBYDX = 1.0 / (DX)
local ONEBYDY = 1.0 / (DY)
local ONEBYDZ = 1.0 / (DZ)

local parallelism = 1

local a10d1 = ( 17.0/ 12.0)/2.0
local b10d1 = (101.0/150.0)/4.0
local c10d1 = (  1.0/100.0)/6.0

local a10d2 = (1065.0/1798.0)/1.0
local b10d2 = (1038.0/ 899.0)/4.0
local c10d2 = (  79.0/1798.0)/9.0

fspace coordinates {
  x   : double,
  y   : double,
  z   : double,
}

fspace point {
  f   : double,
  dfx : double,
  dfy : double,
  dfz : double,
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

task factorize(parallelism : int) : int2d
  var limit = [int](cmath.sqrt([double](parallelism)))
  var size_x = 1
  var size_y = parallelism
  for i = 1, limit + 1 do
    if parallelism % i == 0 then
      size_x, size_y = i, parallelism / i
      if size_x > size_y then
        size_x, size_y = size_y, size_x
      end
    end
  end
  return int2d { size_x, size_y }
end

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

task make_xpencil( points  : region(ispace(int3d), point),
                   xpencil : ispace(int2d) )
  var coloring = c.legion_domain_point_coloring_create()

  var prow = xpencil.bounds.hi.x + 1
  var pcol = xpencil.bounds.hi.y + 1

  var bounds = points.ispace.bounds
  var Nx = bounds.hi.x + 1
  var Ny = bounds.hi.y + 1
  var Nz = bounds.hi.z + 1

  --c.printf("make_xpencil:\n")
  for i in xpencil do
    var lo = int3d { x = 0, y = i.x*(Ny/prow), z = i.y*(Nz/pcol) }
    var hi = int3d { x = Nx-1, y = (i.x+1)*(Ny/prow)-1, z = (i.y+1)*(Nz/pcol)-1 }
    var rect = rect3d { lo = lo, hi = hi }
    c.legion_domain_point_coloring_color_domain(coloring, i, rect)
    --c.printf("    (%d,%d), lo = %d, %d, %d, hi = %d, %d, %d\n",i.x,i.y,lo.x,lo.y,lo.z,hi.x,hi.y,hi.z)
  end
  var p = partition(disjoint, points, coloring, xpencil)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

task make_xpencil_c( points  : region(ispace(int3d), coordinates),
                     xpencil : ispace(int2d) )
  var coloring = c.legion_domain_point_coloring_create()

  var prow = xpencil.bounds.hi.x + 1
  var pcol = xpencil.bounds.hi.y + 1

  var bounds = points.ispace.bounds
  var Nx = bounds.hi.x + 1
  var Ny = bounds.hi.y + 1
  var Nz = bounds.hi.z + 1

  for i in xpencil do
    var lo = int3d { x = 0, y = i.x*(Ny/prow), z = i.y*(Nz/pcol) }
    var hi = int3d { x = Nx-1, y = (i.x+1)*(Ny/prow)-1, z = (i.y+1)*(Nz/pcol)-1 }
    var rect = rect3d { lo = lo, hi = hi }
    c.legion_domain_point_coloring_color_domain(coloring, i, rect)
  end
  var p = partition(disjoint, points, coloring, xpencil)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

task make_ypencil( points  : region(ispace(int3d), point),
                   ypencil : ispace(int2d) )
  var coloring = c.legion_domain_point_coloring_create()

  var prow = ypencil.bounds.hi.x + 1
  var pcol = ypencil.bounds.hi.y + 1

  var bounds = points.ispace.bounds
  var Nx = bounds.hi.x + 1
  var Ny = bounds.hi.y + 1
  var Nz = bounds.hi.z + 1

  --c.printf("make_ypencil:\n")
  for i in ypencil do
    var lo = int3d { x = i.x*(Nx/prow), y = 0, z = i.y*(Nz/pcol) }
    var hi = int3d { x = (i.x+1)*(Nx/prow)-1, y = Ny-1, z = (i.y+1)*(Nz/pcol)-1 }
    var rect = rect3d { lo = lo, hi = hi }
    c.legion_domain_point_coloring_color_domain(coloring, i, rect)
    --c.printf("    (%d,%d), lo = %d, %d, %d, hi = %d, %d, %d\n",i.x,i.y,lo.x,lo.y,lo.z,hi.x,hi.y,hi.z)
  end
  var p = partition(disjoint, points, coloring, ypencil)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

task make_ypencil_c( points  : region(ispace(int3d), coordinates),
                     ypencil : ispace(int2d) )
  var coloring = c.legion_domain_point_coloring_create()

  var prow = ypencil.bounds.hi.x + 1
  var pcol = ypencil.bounds.hi.y + 1

  var bounds = points.ispace.bounds
  var Nx = bounds.hi.x + 1
  var Ny = bounds.hi.y + 1
  var Nz = bounds.hi.z + 1

  for i in ypencil do
    var lo = int3d { x = i.x*(Nx/prow), y = 0, z = i.y*(Nz/pcol) }
    var hi = int3d { x = (i.x+1)*(Nx/prow)-1, y = Ny-1, z = (i.y+1)*(Nz/pcol)-1 }
    var rect = rect3d { lo = lo, hi = hi }
    c.legion_domain_point_coloring_color_domain(coloring, i, rect)
  end
  var p = partition(disjoint, points, coloring, ypencil)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

task make_zpencil( points  : region(ispace(int3d), point),
                   zpencil : ispace(int2d) )
  var coloring = c.legion_domain_point_coloring_create()

  var prow = zpencil.bounds.hi.x + 1
  var pcol = zpencil.bounds.hi.y + 1

  var bounds = points.ispace.bounds
  var Nx = bounds.hi.x + 1
  var Ny = bounds.hi.y + 1
  var Nz = bounds.hi.z + 1

  --c.printf("make_zpencil:\n")
  for i in zpencil do
    var lo = int3d { x = i.x*(Nx/prow), y = i.y*(Ny/pcol), z = 0 }
    var hi = int3d { x = (i.x+1)*(Nx/prow)-1, y = (i.y+1)*(Ny/pcol)-1, z = Nz-1 }
    var rect = rect3d { lo = lo, hi = hi }
    c.legion_domain_point_coloring_color_domain(coloring, i, rect)
    --c.printf("    (%d,%d), lo = %d, %d, %d, hi = %d, %d, %d\n",i.x,i.y,lo.x,lo.y,lo.z,hi.x,hi.y,hi.z)
  end
  var p = partition(disjoint, points, coloring, zpencil)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

task make_zpencil_c( points  : region(ispace(int3d), coordinates),
                     zpencil : ispace(int2d) )
  var coloring = c.legion_domain_point_coloring_create()

  var prow = zpencil.bounds.hi.x + 1
  var pcol = zpencil.bounds.hi.y + 1

  var bounds = points.ispace.bounds
  var Nx = bounds.hi.x + 1
  var Ny = bounds.hi.y + 1
  var Nz = bounds.hi.z + 1

  for i in zpencil do
    var lo = int3d { x = i.x*(Nx/prow), y = i.y*(Ny/pcol), z = 0 }
    var hi = int3d { x = (i.x+1)*(Nx/prow)-1, y = (i.y+1)*(Ny/pcol)-1, z = Nz-1 }
    var rect = rect3d { lo = lo, hi = hi }
    c.legion_domain_point_coloring_color_domain(coloring, i, rect)
  end
  var p = partition(disjoint, points, coloring, zpencil)
  c.legion_domain_point_coloring_destroy(coloring)
  return p
end

task get_LU_decomposition(LU : region(ispace(int3d), LU_struct),
                          e  : double,
                          a  : double,
                          d  : double,
                          cc : double,
                          f  : double)
where
  reads writes( LU )
do

  var N : int64 = LU.ispace.bounds.hi.x + 1
  var pr = LU.ispace.bounds.hi.y
  var pc = LU.ispace.bounds.hi.z

  -- Step 1
  LU[{0,pr,pc}].g = d
  LU[{1,pr,pc}].b = a/LU[{0,pr,pc}].g
  LU[{0,pr,pc}].h = cc
  LU[{0,pr,pc}].k = f/LU[{0,pr,pc}].g
  LU[{0,pr,pc}].w = a
  LU[{0,pr,pc}].v = e
  LU[{0,pr,pc}].l = cc/LU[{0,pr,pc}].g
  
  LU[{1,pr,pc}].g = d - LU[{1,pr,pc}].b*LU[{0,pr,pc}].h
  LU[{1,pr,pc}].k = -LU[{0,pr,pc}].k*LU[{0,pr,pc}].h/LU[{1,pr,pc}].g
  LU[{1,pr,pc}].w = e - LU[{1,pr,pc}].b*LU[{0,pr,pc}].w
  LU[{1,pr,pc}].v = -LU[{1,pr,pc}].b*LU[{0,pr,pc}].v
  LU[{1,pr,pc}].l = (f - LU[{0,pr,pc}].l*LU[{0,pr,pc}].h) / LU[{1,pr,pc}].g
  LU[{1,pr,pc}].h = cc - LU[{1,pr,pc}].b*f

  -- Step 2
  for i = 2,N-3 do
    LU[{i,pr,pc}].b = ( a - ( e/LU[{i-2,pr,pc}].g )*LU[{i-2,pr,pc}].h ) / LU[{i-1,pr,pc}].g
    LU[{i,pr,pc}].h = cc - LU[{i,pr,pc}].b*f
    LU[{i,pr,pc}].g = d - ( e/LU[{i-2,pr,pc}].g )*f - LU[{i,pr,pc}].b*LU[{i-1,pr,pc}].h
  end

  -- Step 3
  LU[{N-3,pr,pc}].b = ( a - ( e/LU[{N-5,pr,pc}].g )*LU[{N-5,pr,pc}].h ) / LU[{N-4,pr,pc}].g
  LU[{N-3,pr,pc}].g = d - ( e/LU[{N-5,pr,pc}].g )*f - LU[{N-3,pr,pc}].b*LU[{N-4,pr,pc}].h

  -- Step 4
  for i = 2,N-4 do
    LU[{i,pr,pc}].k = -( LU[{i-2,pr,pc}].k*f + LU[{i-1,pr,pc}].k*LU[{i-1,pr,pc}].h ) / LU[{i,pr,pc}].g
    LU[{i,pr,pc}].v = -( e/LU[{i-2,pr,pc}].g )*LU[{i-2,pr,pc}].v - LU[{i,pr,pc}].b*LU[{i-1,pr,pc}].v
  end

  -- Step 5
  LU[{N-4,pr,pc}].k = ( e - LU[{N-6,pr,pc}].k*f - LU[{N-5,pr,pc}].k*LU[{N-5,pr,pc}].h ) / LU[{N-4,pr,pc}].g
  LU[{N-3,pr,pc}].k = ( a - LU[{N-5,pr,pc}].k*f - LU[{N-4,pr,pc}].k*LU[{N-4,pr,pc}].h ) / LU[{N-3,pr,pc}].g
  LU[{N-4,pr,pc}].v = f  - ( e/LU[{N-6,pr,pc}].g )*LU[{N-6,pr,pc}].v - LU[{N-4,pr,pc}].b*LU[{N-5,pr,pc}].v
  LU[{N-3,pr,pc}].v = cc - ( e/LU[{N-5,pr,pc}].g )*LU[{N-5,pr,pc}].v - LU[{N-3,pr,pc}].b*LU[{N-4,pr,pc}].v
  LU[{N-2,pr,pc}].g = d
  for i = 0,N-2 do
    LU[{N-2,pr,pc}].g -= LU[{i,pr,pc}].k*LU[{i,pr,pc}].v
  end

  -- Step 6
  for i = 2,N-3 do
    LU[{i,pr,pc}].w = -( e/LU[{i-2,pr,pc}].g )*LU[{i-2,pr,pc}].w - LU[{i,pr,pc}].b*LU[{i-1,pr,pc}].w
    LU[{i,pr,pc}].l = -( LU[{i-2,pr,pc}].l*f + LU[{i-1,pr,pc}].l*LU[{i-1,pr,pc}].h ) / LU[{i,pr,pc}].g
  end

  -- Step 7
  LU[{N-3,pr,pc}].w = f - ( e/LU[{N-5,pr,pc}].g )*LU[{N-5,pr,pc}].w - LU[{N-3,pr,pc}].b*LU[{N-4,pr,pc}].w
  LU[{N-2,pr,pc}].w = cc
  for i = 0,N-2 do
    LU[{N-2,pr,pc}].w -= LU[{i,pr,pc}].k*LU[{i,pr,pc}].w
  end
  LU[{N-3,pr,pc}].l = ( e - LU[{N-5,pr,pc}].l*f - LU[{N-4,pr,pc}].l*LU[{N-4,pr,pc}].h ) / LU[{N-3,pr,pc}].g
  LU[{N-2,pr,pc}].l = a
  for i = 0,N-2 do
    LU[{N-2,pr,pc}].l -= LU[{i,pr,pc}].l*LU[{i,pr,pc}].v
  end
  LU[{N-2,pr,pc}].l = LU[{N-2,pr,pc}].l / LU[{N-2,pr,pc}].g
  LU[{N-1,pr,pc}].g = d
  for i = 0,N-1 do
    LU[{N-1,pr,pc}].g -= LU[{i,pr,pc}].l*LU[{i,pr,pc}].w
  end

  -- Set eg = e/g
  for i = 2,N-2 do
    LU[{i,pr,pc}].eg = e/LU[{i-2,pr,pc}].g
  end

  -- Set ff = f
  for i = 0,N-4 do
    LU[{i,pr,pc}].ff = f
  end

  -- Set g = 1/g
  for i = 0,N do
    LU[{i,pr,pc}].g = 1.0/LU[{i,pr,pc}].g
  end

  -- c.printf("LU decomposition:\n")
  -- for i = 0,N do
  --   c.printf("%8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f\n",LU[i].b,LU[i].eg,LU[i].k,LU[i].l,LU[i].g,LU[i].h,LU[i].ff,LU[i].v,LU[i].w)
  -- end

end

local function make_SolveXLU(field)
  local points = regentlib.newsymbol(region(ispace(int3d), point), "points")
  local reads_points_dfx = regentlib.privilege(regentlib.reads, points, field)
  local writes_points_dfx = regentlib.privilege(regentlib.writes, points, field)
  local reads_writes_points_dfx = terralib.newlist({reads_points_dfx, writes_points_dfx})

  local task SolveXLU_meta( [points],
                            LU : region(ispace(int3d), LU_struct) )
  where
    [reads_writes_points_dfx], reads(LU)
  do
    var bounds = points.ispace.bounds
    var N = bounds.hi.x + 1
    var pr = LU.ispace.bounds.hi.y
    var pc = LU.ispace.bounds.hi.z
  
    for j = bounds.lo.y, bounds.hi.y+1 do
      for k = bounds.lo.z, bounds.hi.z+1 do
  
        -- Step 8
        points[{1,j,k}].dfx = points[{1,j,k}].dfx - LU[{1,pr,pc}].b*points[{0,j,k}].dfx
        var sum1 : double = LU[{0,pr,pc}].k*points[{0,j,k}].dfx + LU[{1,pr,pc}].k*points[{1,j,k}].dfx
        var sum2 : double = LU[{0,pr,pc}].l*points[{0,j,k}].dfx + LU[{1,pr,pc}].l*points[{1,j,k}].dfx
  
        -- Step 9
        for i = 2,N-2 do
          points[{i,j,k}].dfx = points[{i,j,k}].dfx - LU[{i,pr,pc}].b*points[{i-1,j,k}].dfx - LU[{i,pr,pc}].eg*points[{i-2,j,k}].dfx
          sum1 += LU[{i,pr,pc}].k*points[{i,j,k}].dfx
          sum2 += LU[{i,pr,pc}].l*points[{i,j,k}].dfx
        end
  
        -- Step 10
        points[{N-2,j,k}].dfx = points[{N-2,j,k}].dfx - sum1
        points[{N-1,j,k}].dfx = ( points[{N-1,j,k}].dfx - sum2 - LU[{N-2,pr,pc}].l*points[{N-2,j,k}].dfx )*LU[{N-1,pr,pc}].g
  
        -- Step 11
        points[{N-2,j,k}].dfx = ( points[{N-2,j,k}].dfx - LU[{N-2,pr,pc}].w*points[{N-1,j,k}].dfx )*LU[{N-2,pr,pc}].g
        points[{N-3,j,k}].dfx = ( points[{N-3,j,k}].dfx - LU[{N-3,pr,pc}].v*points[{N-2,j,k}].dfx - LU[{N-3,pr,pc}].w*points[{N-1,j,k}].dfx )*LU[{N-3,pr,pc}].g
        points[{N-4,j,k}].dfx = ( points[{N-4,j,k}].dfx - LU[{N-4,pr,pc}].h*points[{N-3,j,k}].dfx - LU[{N-4,pr,pc}].v*points[{N-2,j,k}].dfx - LU[{N-4,pr,pc}].w*points[{N-1,j,k}].dfx )*LU[{N-4,pr,pc}].g
        for i = N-5,-1,-1 do
          points[{i,j,k}].dfx = ( points[{i,j,k}].dfx - LU[{i,pr,pc}].h*points[{i+1,j,k}].dfx - LU[{i,pr,pc}].ff*points[{i+2,j,k}].dfx - LU[{i,pr,pc}].v*points[{N-2,j,k}].dfx - LU[{i,pr,pc}].w*points[{N-1,j,k}].dfx )*LU[{i,pr,pc}].g
        end
  
      end
    end
    return 1
  end
  return SolveXLU_meta
end

SolveXLU = make_SolveXLU("dfx")

task SolveYLU( points : region(ispace(int3d), point),
               LU     : region(ispace(int3d), LU_struct) )
where
  reads writes(points.dfy), reads(LU)
do
  var bounds = points.ispace.bounds
  var N = bounds.hi.y + 1
  var pr = LU.ispace.bounds.hi.y
  var pc = LU.ispace.bounds.hi.z

  for i = bounds.lo.x, bounds.hi.x+1 do
    for k = bounds.lo.z, bounds.hi.z+1 do

      -- Step 8
      points[{i,1,k}].dfy = points[{i,1,k}].dfy - LU[{1,pr,pc}].b*points[{i,0,k}].dfy
      var sum1 : double = LU[{0,pr,pc}].k*points[{i,0,k}].dfy + LU[{1,pr,pc}].k*points[{i,1,k}].dfy
      var sum2 : double = LU[{0,pr,pc}].l*points[{i,0,k}].dfy + LU[{1,pr,pc}].l*points[{i,1,k}].dfy

      -- Step 9
      for j = 2,N-2 do
        points[{i,j,k}].dfy = points[{i,j,k}].dfy - LU[{j,pr,pc}].b*points[{i,j-1,k}].dfy - LU[{j,pr,pc}].eg*points[{i,j-2,k}].dfy
        sum1 += LU[{j,pr,pc}].k*points[{i,j,k}].dfy
        sum2 += LU[{j,pr,pc}].l*points[{i,j,k}].dfy
      end

      -- Step 10
      points[{i,N-2,k}].dfy = points[{i,N-2,k}].dfy - sum1
      points[{i,N-1,k}].dfy = ( points[{i,N-1,k}].dfy - sum2 - LU[{N-2,pr,pc}].l*points[{i,N-2,k}].dfy )*LU[{N-1,pr,pc}].g

      -- Step 11
      points[{i,N-2,k}].dfy = ( points[{i,N-2,k}].dfy - LU[{N-2,pr,pc}].w*points[{i,N-1,k}].dfy )*LU[{N-2,pr,pc}].g
      points[{i,N-3,k}].dfy = ( points[{i,N-3,k}].dfy - LU[{N-3,pr,pc}].v*points[{i,N-2,k}].dfy - LU[{N-3,pr,pc}].w*points[{i,N-1,k}].dfy )*LU[{N-3,pr,pc}].g
      points[{i,N-4,k}].dfy = ( points[{i,N-4,k}].dfy - LU[{N-4,pr,pc}].h*points[{i,N-3,k}].dfy - LU[{N-4,pr,pc}].v*points[{i,N-2,k}].dfy - LU[{N-4,pr,pc}].w*points[{i,N-1,k}].dfy )*LU[{N-4,pr,pc}].g
      for j = N-5,-1,-1 do
        points[{i,j,k}].dfy = ( points[{i,j,k}].dfy - LU[{j,pr,pc}].h*points[{i,j+1,k}].dfy - LU[{j,pr,pc}].ff*points[{i,j+2,k}].dfy - LU[{j,pr,pc}].v*points[{i,N-2,k}].dfy - LU[{j,pr,pc}].w*points[{i,N-1,k}].dfy )*LU[{j,pr,pc}].g
      end

    end
  end

  return 1
end

task SolveZLU( points : region(ispace(int3d), point),
               LU     : region(ispace(int3d), LU_struct) )
where
  reads writes(points.dfz), reads(LU)
do
  var bounds = points.ispace.bounds
  var N = bounds.hi.z + 1
  var pr = LU.ispace.bounds.hi.y
  var pc = LU.ispace.bounds.hi.z

  for i = bounds.lo.x, bounds.hi.x+1 do
    for j = bounds.lo.y, bounds.hi.y+1 do

      -- Step 8
      points[{i,j,1}].dfz = points[{i,j,1}].dfz - LU[{1,pr,pc}].b*points[{i,j,0}].dfz
      var sum1 : double = LU[{0,pr,pc}].k*points[{i,j,0}].dfz + LU[{1,pr,pc}].k*points[{i,j,1}].dfz
      var sum2 : double = LU[{0,pr,pc}].l*points[{i,j,0}].dfz + LU[{1,pr,pc}].l*points[{i,j,1}].dfz

      -- Step 9
      for k = 2,N-2 do
        points[{i,j,k}].dfz = points[{i,j,k}].dfz - LU[{k,pr,pc}].b*points[{i,j,k-1}].dfz - LU[{k,pr,pc}].eg*points[{i,j,k-2}].dfz
        sum1 += LU[{k,pr,pc}].k*points[{i,j,k}].dfz
        sum2 += LU[{k,pr,pc}].l*points[{i,j,k}].dfz
      end

      -- Step 10
      points[{i,j,N-2}].dfz = points[{i,j,N-2}].dfz - sum1
      points[{i,j,N-1}].dfz = ( points[{i,j,N-1}].dfz - sum2 - LU[{N-2,pr,pc}].l*points[{i,j,N-2}].dfz )*LU[{N-1,pr,pc}].g

      -- Step 11
      points[{i,j,N-2}].dfz = ( points[{i,j,N-2}].dfz - LU[{N-2,pr,pc}].w*points[{i,j,N-1}].dfz )*LU[{N-2,pr,pc}].g
      points[{i,j,N-3}].dfz = ( points[{i,j,N-3}].dfz - LU[{N-3,pr,pc}].v*points[{i,j,N-2}].dfz - LU[{N-3,pr,pc}].w*points[{i,j,N-1}].dfz )*LU[{N-3,pr,pc}].g
      points[{i,j,N-4}].dfz = ( points[{i,j,N-4}].dfz - LU[{N-4,pr,pc}].h*points[{i,j,N-3}].dfz - LU[{N-4,pr,pc}].v*points[{i,j,N-2}].dfz - LU[{N-4,pr,pc}].w*points[{i,j,N-1}].dfz )*LU[{N-4,pr,pc}].g
      for k = N-5,-1,-1 do
        points[{i,j,k}].dfz = ( points[{i,j,k}].dfz - LU[{k,pr,pc}].h*points[{i,j,k+1}].dfz - LU[{k,pr,pc}].ff*points[{i,j,k+2}].dfz - LU[{k,pr,pc}].v*points[{i,j,N-2}].dfz - LU[{k,pr,pc}].w*points[{i,j,N-1}].dfz )*LU[{k,pr,pc}].g
      end

    end
  end
  return 1
end

local function poff(i, x, y, z, Nx, Ny, Nz)
  return rexpr int3d { x = (i.x + x + Nx)%Nx, y = (i.y + y + Ny)%Ny, z = (i.z + z + Nz)%Nz } end
end

local function make_stencil_pattern(points, index, a10, b10, c10, Nx, Ny, Nz, onebydx, dir, der)
  local value

  if dir == 0 then      -- x direction stencil
    if der == 1 then
      value = rexpr       - c10*points[ [poff(index, -3, 0, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value - b10*points[ [poff(index, -2, 0, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value - a10*points[ [poff(index, -1, 0, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value + a10*points[ [poff(index,  1, 0, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value + b10*points[ [poff(index,  2, 0, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value + c10*points[ [poff(index,  3, 0, 0, Nx, Ny, Nz)] ].f end
      value = rexpr onebydx * ( value ) end
    elseif der == 2 then
      value = rexpr         a10*( points[ [poff(index, -1, 0, 0, Nx, Ny, Nz)] ].f - 2.0*points[ index ].f + points[ [poff(index, 1, 0, 0, Nx, Ny, Nz)] ].f ) end
      value = rexpr value + b10*( points[ [poff(index, -2, 0, 0, Nx, Ny, Nz)] ].f - 2.0*points[ index ].f + points[ [poff(index, 2, 0, 0, Nx, Ny, Nz)] ].f ) end
      value = rexpr value + c10*( points[ [poff(index, -3, 0, 0, Nx, Ny, Nz)] ].f - 2.0*points[ index ].f + points[ [poff(index, 3, 0, 0, Nx, Ny, Nz)] ].f ) end
      value = rexpr onebydx*onebydx * (value) end
    end
  elseif dir == 1 then  -- y direction stencil
    if der == 1 then
      value = rexpr       - c10*points[ [poff(index, 0, -3, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value - b10*points[ [poff(index, 0, -2, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value - a10*points[ [poff(index, 0, -1, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value + a10*points[ [poff(index, 0,  1, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value + b10*points[ [poff(index, 0,  2, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value + c10*points[ [poff(index, 0,  3, 0, Nx, Ny, Nz)] ].f end
      value = rexpr onebydx * ( value ) end
    elseif der == 2 then
      value = rexpr         a10*( points[ [poff(index, 0, -1, 0, Nx, Ny, Nz)] ].f - 2.0*points[ index ].f + points[ [poff(index, 0, 1, 0, Nx, Ny, Nz)] ].f ) end
      value = rexpr value + b10*( points[ [poff(index, 0, -2, 0, Nx, Ny, Nz)] ].f - 2.0*points[ index ].f + points[ [poff(index, 0, 2, 0, Nx, Ny, Nz)] ].f ) end
      value = rexpr value + c10*( points[ [poff(index, 0, -3, 0, Nx, Ny, Nz)] ].f - 2.0*points[ index ].f + points[ [poff(index, 0, 3, 0, Nx, Ny, Nz)] ].f ) end
      value = rexpr onebydx*onebydx * (value) end
    end
  elseif dir == 2 then  -- z direction stencil
    if der == 1 then
      value = rexpr       - c10*points[ [poff(index, 0, 0, -3, Nx, Ny, Nz)] ].f end
      value = rexpr value - b10*points[ [poff(index, 0, 0, -2, Nx, Ny, Nz)] ].f end
      value = rexpr value - a10*points[ [poff(index, 0, 0, -1, Nx, Ny, Nz)] ].f end
      value = rexpr value + a10*points[ [poff(index, 0, 0,  1, Nx, Ny, Nz)] ].f end
      value = rexpr value + b10*points[ [poff(index, 0, 0,  2, Nx, Ny, Nz)] ].f end
      value = rexpr value + c10*points[ [poff(index, 0, 0,  3, Nx, Ny, Nz)] ].f end
      value = rexpr onebydx * ( value ) end
    elseif der == 2 then
      value = rexpr         a10*( points[ [poff(index, 0, 0, -1, Nx, Ny, Nz)] ].f - 2.0*points[ index ].f + points[ [poff(index, 0, 0, 1, Nx, Ny, Nz)] ].f ) end
      value = rexpr value + b10*( points[ [poff(index, 0, 0, -2, Nx, Ny, Nz)] ].f - 2.0*points[ index ].f + points[ [poff(index, 0, 0, 2, Nx, Ny, Nz)] ].f ) end
      value = rexpr value + c10*( points[ [poff(index, 0, 0, -3, Nx, Ny, Nz)] ].f - 2.0*points[ index ].f + points[ [poff(index, 0, 0, 3, Nx, Ny, Nz)] ].f ) end
      value = rexpr onebydx*onebydx * (value) end
    end
  end
  return value
end

local function make_stencil_x(Nx, Ny, Nz, onebydx, a10, b10, c10, der)
  local task rhs_x( points : region(ispace(int3d), point) )
  where
    reads(points.f), writes(points.dfx)
  do
    for i in points do
      points[i].dfx = [make_stencil_pattern(points, i, a10, b10, c10, Nx, Ny, Nz, onebydx, 0, der)]
    end
  end
  return rhs_x
end

local function make_stencil_y(Nx, Ny, Nz, onebydy, a10, b10, c10, der)
  local task rhs_y( points : region(ispace(int3d), point) )
  where
    reads(points.f), writes(points.dfy)
  do
    for i in points do
      points[i].dfy = [make_stencil_pattern(points, i, a10, b10, c10, Nx, Ny, Nz, onebydy, 1, der)]
    end
  end
  return rhs_y
end

local function make_stencil_z(Nx, Ny, Nz, onebydz, a10, b10, c10, der)
  local task rhs_z( points : region(ispace(int3d), point) )
  where
    reads(points.f), writes(points.dfz)
  do
    for i in points do
      points[i].dfz = [make_stencil_pattern(points, i, a10, b10, c10, Nx, Ny, Nz, onebydz, 2, der)]
    end
  end
  return rhs_z
end

local ComputeXRHS  = make_stencil_x(NX, NY, NZ, ONEBYDX, a10d1, b10d1, c10d1, 1)
local ComputeYRHS  = make_stencil_y(NX, NY, NZ, ONEBYDY, a10d1, b10d1, c10d1, 1)
local ComputeZRHS  = make_stencil_z(NX, NY, NZ, ONEBYDZ, a10d1, b10d1, c10d1, 1)

local ComputeX2RHS = make_stencil_x(NX, NY, NZ, ONEBYDX, a10d2, b10d2, c10d2, 2)
local ComputeY2RHS = make_stencil_y(NX, NY, NZ, ONEBYDY, a10d2, b10d2, c10d2, 2)
local ComputeZ2RHS = make_stencil_z(NX, NY, NZ, ONEBYDZ, a10d2, b10d2, c10d2, 2)

task ddx( points : region(ispace(int3d), point),
          LU     : region(ispace(int3d), LU_struct) )
where
  reads(LU, points.f), reads writes(points.dfx)
do
  ComputeXRHS(points)
  var token = SolveXLU(points,LU)
  --c.printf("In ddx\n")
  --for p in points do
  --  if p.x == 0 and p.y == 0 and p.z == 0 then
  --    c.printf("{%d,%d,%d}: %8.5f, %8.5f\n",p.x,p.y,p.z,points[p].f,points[p].dfx)
  --  end
  --end
  token = points[points.ispace.bounds.lo].dfx
  return token
end

task d2dx2( points : region(ispace(int3d), point),
            LU     : region(ispace(int3d), LU_struct) )
where
  reads(LU, points.f), reads writes(points.dfx)
do
  ComputeX2RHS(points)
  var token = SolveXLU(points,LU)
  --c.printf("In d2dx2\n")
  --for p in points do
  --  if p.x == 0 and p.y == 0 and p.z == 0 then
  --    c.printf("{%d,%d,%d}: %8.5f, %8.5f\n",p.x,p.y,p.z,points[p].f,points[p].dfx)
  --  end
  --end
  token = points[points.ispace.bounds.lo].dfx
  return token
end

task ddy( points : region(ispace(int3d), point),
          LU     : region(ispace(int3d), LU_struct) )
where
  reads(LU, points.f), reads writes(points.dfy)
do
  ComputeYRHS(points)
  var token = SolveYLU(points,LU)
  --c.printf("In ddy\n")
  --for p in points do
  --  if p.x == 0 and p.y == 0 and p.z == 0 then
  --    c.printf("{%d,%d,%d}: %8.5f, %8.5f\n",p.x,p.y,p.z,points[p].f,points[p].dfy)
  --  end
  --end
  token = points[points.ispace.bounds.lo].dfy
  return token
end

task d2dy2( points : region(ispace(int3d), point),
            LU     : region(ispace(int3d), LU_struct) )
where
  reads(LU, points.f), reads writes(points.dfy)
do
  ComputeY2RHS(points)
  var token = SolveYLU(points,LU)
  --c.printf("In d2dy2\n")
  --for p in points do
  --  if p.x == 0 and p.y == 0 and p.z == 0 then
  --    c.printf("{%d,%d,%d}: %8.5f, %8.5f\n",p.x,p.y,p.z,points[p].f,points[p].dfy)
  --  end
  --end
  token = points[points.ispace.bounds.lo].dfy
  return token
end

task ddz( points : region(ispace(int3d), point),
          LU     : region(ispace(int3d), LU_struct) )
where
  reads(LU, points.f), reads writes(points.dfz)
do
  ComputeZRHS(points)
  var token = SolveZLU(points,LU)
  --c.printf("In ddz\n")
  --for p in points do
  --  if p.x == 0 and p.y == 0 and p.z == 0 then
  --    c.printf("{%d,%d,%d}: %8.5f, %8.5f\n",p.x,p.y,p.z,points[p].f,points[p].dfz)
  --  end
  --end
  token = points[points.ispace.bounds.lo].dfz
  return token
end

task d2dz2( points : region(ispace(int3d), point),
            LU     : region(ispace(int3d), LU_struct) )
where
  reads(LU, points.f), reads writes(points.dfz)
do
  ComputeZ2RHS(points)
  var token = SolveZLU(points,LU)
  --c.printf("In d2dz2\n")
  --for p in points do
  --  if p.x == 0 and p.y == 0 and p.z == 0 then
  --    c.printf("{%d,%d,%d}: %8.5f, %8.5f\n",p.x,p.y,p.z,points[p].f,points[p].dfz)
  --  end
  --end
  token = points[points.ispace.bounds.lo].dfz
  return token
end

task gradient( points_x : region(ispace(int3d), point),
               points_y : region(ispace(int3d), point),
               points_z : region(ispace(int3d), point),
               LU_x     : region(ispace(int3d), LU_struct),
               LU_y     : region(ispace(int3d), LU_struct),
               LU_z     : region(ispace(int3d), LU_struct) )
where
  reads(LU_x, LU_y, LU_z), reads(points_x.f, points_y.f, points_z.f),
  reads writes(points_x.dfx, points_y.dfy, points_z.dfz)
do 
  var token = 0
  token += ddx(points_x,LU_x)
  token += ddy(points_y,LU_y)
  token += ddz(points_z,LU_z)
  return token
end

task initialize( points : region(ispace(int3d), point),
                 exact  : region(ispace(int3d), point),
                 coords : region(ispace(int3d), coordinates),
                 dx     : double,
                 dy     : double,
                 dz     : double )
where
  reads writes(coords.x, coords.y, coords.z, points.f, exact.f, exact.dfx, exact.dfy, exact.dfz)
do
  var bounds = points.ispace.bounds

  for i = bounds.lo.x, bounds.hi.x+1 do
    for j = bounds.lo.y, bounds.hi.y+1 do
      for k = bounds.lo.z, bounds.hi.z+1 do
        var e : int3d = { x = i, y = j, z = k }
        coords[e].x   = i*dx
        coords[e].y   = j*dy
        coords[e].z   = k*dz
        points[e].f   = cmath.sin(coords[e].x) + cmath.sin(coords[e].y) + cmath.sin(coords[e].z)
        exact [e].f   = cmath.sin(coords[e].x) + cmath.sin(coords[e].y) + cmath.sin(coords[e].z)
        exact [e].dfx = cmath.cos(coords[e].x)
        exact [e].dfy = cmath.cos(coords[e].y)
        exact [e].dfz = cmath.cos(coords[e].z)
      end
    end
  end
 
  --c.printf("In initialize\n")
  --for p in points do
  --  if p.x == 0 and p.z == 0 then
  --    c.printf("{%d,%d,%d}: %8.5f\n",p.x,p.y,p.z,points[p].f)
  --  end
  --end
  return 0
end

task get_error_x( points : region(ispace(int3d), point),
                  exact  : region(ispace(int3d), point) )
where
  reads(points.dfx, exact.dfx)
do
  var err : double = 0.0
  for i in points do
    err = max(err, cmath.fabs(points[i].dfx - exact[i].dfx))
  end
  return err
end

task get_error_y( points : region(ispace(int3d), point),
                  exact  : region(ispace(int3d), point) )
where
  reads(points.dfy, exact.dfy)
do
  var err : double = 0.0
  for i in points do
    err = max(err, cmath.fabs(points[i].dfy - exact[i].dfy))
    --if i.x == 0 and i.z == 0 then
    --  c.printf("{%d,%d,%d}: %8.5f, %8.5f\n",i.x,i.y,i.z,points[i].dfy,exact[i].dfy)
    --end
  end
  return err
end

task get_error_z( points : region(ispace(int3d), point),
                  exact  : region(ispace(int3d), point) )
where
  reads(points.dfz, exact.dfz)
do
  var err : double = 0.0
  for i in points do
    err = max(err, cmath.fabs(points[i].dfz - exact[i].dfz))
  end
  return err
end

task get_error_d2( points : region(ispace(int3d), point) )
where
  reads(points.f, points.dfx, points.dfy, points.dfz)
do
  var err : double = 0.0
  for i in points do
    err = max(err, cmath.fabs(points[i].dfx + points[i].dfy + points[i].dfz + points[i].f))
  end
  return err
end

terra wait_for(x : int)
  return x
end

task run_main( points   : region(ispace(int3d), point),
               exact    : region(ispace(int3d), point),
               coords   : region(ispace(int3d), coordinates),
               LU_x     : region(ispace(int3d), LU_struct),
               LU_x2    : region(ispace(int3d), LU_struct),
               LU_y     : region(ispace(int3d), LU_struct),
               LU_y2    : region(ispace(int3d), LU_struct),
               LU_z     : region(ispace(int3d), LU_struct),
               LU_z2    : region(ispace(int3d), LU_struct),
               pencil   : ispace(int2d),
               points_x : partition(disjoint, points, ispace(int2d)),
               points_y : partition(disjoint, points, ispace(int2d)),
               points_z : partition(disjoint, points, ispace(int2d)),
               exact_x  : partition(disjoint, exact,  ispace(int2d)),
               exact_y  : partition(disjoint, exact,  ispace(int2d)),
               exact_z  : partition(disjoint, exact,  ispace(int2d)),
               coords_x : partition(disjoint, coords, ispace(int2d)),
               coords_y : partition(disjoint, coords, ispace(int2d)),
               coords_z : partition(disjoint, coords, ispace(int2d)),
               pLU_x    : partition(disjoint, LU_x,   ispace(int2d)),
               pLU_x2   : partition(disjoint, LU_x2,  ispace(int2d)),
               pLU_y    : partition(disjoint, LU_y,   ispace(int2d)),
               pLU_y2   : partition(disjoint, LU_y2,  ispace(int2d)),
               pLU_z    : partition(disjoint, LU_z,   ispace(int2d)),
               pLU_z2   : partition(disjoint, LU_z2,  ispace(int2d)),
               dx       : double,
               dy       : double,
               dz       : double )
where
  reads(LU_x, LU_x2, LU_y, LU_y2, LU_z, LU_z2), reads(coords, exact, points.f), reads writes (points.dfx, points.dfy, points.dfz)
do

  var token = 0
  wait_for(token)
  var ts_start = c.legion_get_current_time_in_micros()
  
  -- Get df/dx, df/dy, df/dz
  --must_epoch
    __demand(__parallel)
    for i in pencil do
      token += ddx(points_x[i],pLU_x[i])
    end
  --end

  --must_epoch
    __demand(__parallel)
    for i in pencil do
      token += ddy(points_y[i],pLU_y[i])
    end
  --end

  --must_epoch
    __demand(__parallel)
    for i in pencil do
      token += ddz(points_z[i],pLU_z[i])
    end
  --end
  
  wait_for(token)
  var ts_d1 = c.legion_get_current_time_in_micros() - ts_start
  
  var err_x = 0.0
  __demand(__parallel)
  for i in pencil do
    err_x += get_error_x(points_x[i],exact_x[i])
  end
  err_x = err_x / parallelism

  var err_y = 0.0
  __demand(__parallel)
  for i in pencil do
    err_y += get_error_y(points_x[i],exact_x[i]) 
  end
  err_y = err_y / parallelism

  var err_z = 0.0
  __demand(__parallel)
  for i in pencil do
    err_z += get_error_z(points_x[i],exact_x[i]) 
  end
  err_z = err_z / parallelism
  
  wait_for(err_x)
  wait_for(err_y)
  wait_for(err_z)
  ts_start = c.legion_get_current_time_in_micros()
  
  -- Get d2f/dx2, d2f/dy2, d2f/dz2
  --must_epoch
    __demand(__parallel)
    for i in pencil do
      token += d2dx2(points_x[i],pLU_x2[i])
    end
  --end

  --must_epoch
    __demand(__parallel)
    for i in pencil do
      token += d2dy2(points_y[i],pLU_y2[i])
    end
  --end

  --must_epoch
    __demand(__parallel)
    for i in pencil do
      token += d2dz2(points_z[i],pLU_z2[i])
    end
  --end
  
  wait_for(token)
  var ts_d2 = c.legion_get_current_time_in_micros() - ts_start
  
  var err_d2 = 0.0
  __demand(__parallel)
  for i in pencil do
    err_d2 += get_error_d2(points_x[i])
  end
  err_d2 = err_d2 / parallelism

  c.printf("Time to get the 1st derivatives: %12.5e\n", (ts_d1)*1e-6)
  c.printf("  Maximum error in x = %12.5e\n", err_x)
  c.printf("  Maximum error in y = %12.5e\n", err_y)
  c.printf("  Maximum error in z = %12.5e\n", err_z)
  c.printf("Time to get the 2nd derivatives: %12.5e\n", (ts_d2)*1e-6)
  c.printf("  Maximum error in laplacian = %12.5e\n", err_d2)
  
end

task main()
  var L  : double = LL      -- Domain length
  var Nx : int64  = NX      -- Grid size
  var Ny : int64  = NY      -- Grid size
  var Nz : int64  = NZ      -- Grid size
  var dx : double = DX      -- Grid spacing
  var dy : double = DY      -- Grid spacing
  var dz : double = DZ      -- Grid spacing

  c.printf("================ Problem parameters ================\n")
  c.printf("           grid size  = %d x %d x %d\n", Nx, Ny, Nz )
  c.printf("                   L  = %f\n", L )
  c.printf("           dx, dy, dz = %f, %f, %f\n", dx, dy, dz)
  c.printf("          parallelism = %d\n", parallelism)
  c.printf("====================================================\n")

  -- Coefficients for the 10th order 1st derivative
  var alpha10d1 : double = 1.0/2.0
  var beta10d1  : double = 1.0/20.0
  
  -- Coefficients for the 10th order 2nd derivative
  var alpha10d2 : double = 334.0/899.0
  var beta10d2  : double = 43.0/1798.0

  var prowcol = factorize(parallelism)
  var pencil = ispace(int2d, prowcol)
  
  var grid_x = ispace(int3d, { x = Nx, y = prowcol.x, z = prowcol.y } )
  var LU_x   = region(grid_x, LU_struct)
  var LU_x2  = region(grid_x, LU_struct)
  
  var pLU_x  = partitionLU(LU_x,  pencil)
  var pLU_x2 = partitionLU(LU_x2, pencil)
  must_epoch
    __demand(__parallel)
    for i in pencil do
      get_LU_decomposition(pLU_x [i], beta10d1, alpha10d1, 1.0, alpha10d1, beta10d1)
    end
  end
  must_epoch
    __demand(__parallel)
    for i in pencil do
      get_LU_decomposition(pLU_x2[i], beta10d2, alpha10d2, 1.0, alpha10d2, beta10d2)
    end
  end

  var grid_y = ispace(int3d, { x = Ny, y = prowcol.x, z = prowcol.y } )
  var LU_y   = region(grid_y, LU_struct)
  var LU_y2  = region(grid_y, LU_struct)

  var pLU_y  = partitionLU(LU_y,  pencil)
  var pLU_y2 = partitionLU(LU_y2, pencil)
  must_epoch
    __demand(__parallel)
    for i in pencil do
      get_LU_decomposition(pLU_y [i], beta10d1, alpha10d1, 1.0, alpha10d1, beta10d1)
    end
  end
  must_epoch
    __demand(__parallel)
    for i in pencil do
      get_LU_decomposition(pLU_y2[i], beta10d2, alpha10d2, 1.0, alpha10d2, beta10d2)
    end
  end

  var grid_z = ispace(int3d, { x = Nz, y = prowcol.x, z = prowcol.y } )
  var LU_z   = region(grid_z, LU_struct)
  var LU_z2  = region(grid_z, LU_struct)

  var pLU_z  = partitionLU(LU_z,  pencil)
  var pLU_z2 = partitionLU(LU_z2, pencil)  
  must_epoch
    __demand(__parallel)
    for i in pencil do
      get_LU_decomposition(pLU_z [i], beta10d1, alpha10d1, 1.0, alpha10d1, beta10d1)
    end
  end
  must_epoch
    __demand(__parallel)
    for i in pencil do
      get_LU_decomposition(pLU_z2[i], beta10d2, alpha10d2, 1.0, alpha10d2, beta10d2)
    end
  end

  var grid   = ispace(int3d, { x = Nx, y = Ny, z = Nz })
  var coords = region(grid, coordinates)
  var points = region(grid, point)
  var exact  = region(grid, point)

  var points_x = make_xpencil(points, pencil) -- Partition of x-pencils
  var points_y = make_ypencil(points, pencil) -- Partition of y-pencils
  var points_z = make_zpencil(points, pencil) -- Partition of z-pencils

  var exact_x  = make_xpencil(exact,  pencil) -- Partition of x-pencils
  var exact_y  = make_ypencil(exact,  pencil) -- Partition of y-pencils
  var exact_z  = make_zpencil(exact,  pencil) -- Partition of z-pencils

  var coords_x = make_xpencil_c(coords, pencil) -- Partition of x-pencils
  var coords_y = make_ypencil_c(coords, pencil) -- Partition of y-pencils
  var coords_z = make_zpencil_c(coords, pencil) -- Partition of z-pencils
 
  var token = 0 
  -- Initialize function f
  --must_epoch
    __demand(__parallel)
    for i in pencil do
      token += initialize(points_x[i], exact_x[i], coords_x[i], dx, dy, dz)
    end
  --end

  -- run_main( points, exact, coords, LU_x, LU_x2, LU_y, LU_y2, LU_z, LU_z2,
  --           pencil, points_x, points_y, points_z, exact_x, exact_y, exact_z,
  --           coords_x, coords_y, coords_z, pLU_x, pLU_x2, pLU_y, pLU_y2, pLU_z, pLU_z2,
  --           dx, dy, dz )

  for iter = 0,10 do  
    wait_for(token)
    var ts_start = c.legion_get_current_time_in_micros()
    
    -- Get df/dx, df/dy, df/dz
    --must_epoch
      __demand(__parallel)
      for i in pencil do
        token += ddx(points_x[i],pLU_x[i])
      end
    --end

    --must_epoch
      __demand(__parallel)
      for i in pencil do
        token += ddy(points_y[i],pLU_y[i])
      end
    --end

    --must_epoch
      __demand(__parallel)
      for i in pencil do
        token += ddz(points_z[i],pLU_z[i])
      end
    --end
    
    wait_for(token)
    var ts_d1 = c.legion_get_current_time_in_micros() - ts_start
    
    var err_x = 0.0
    __demand(__parallel)
    for i in pencil do
      err_x += get_error_x(points_x[i],exact_x[i])
    end
    err_x = err_x / parallelism

    var err_y = 0.0
    __demand(__parallel)
    for i in pencil do
      err_y += get_error_y(points_x[i],exact_x[i]) 
    end
    err_y = err_y / parallelism

    var err_z = 0.0
    __demand(__parallel)
    for i in pencil do
      err_z += get_error_z(points_x[i],exact_x[i]) 
    end
    err_z = err_z / parallelism
    
    wait_for(err_x)
    wait_for(err_y)
    wait_for(err_z)
    ts_start = c.legion_get_current_time_in_micros()
    
    -- Get d2f/dx2, d2f/dy2, d2f/dz2
    --must_epoch
      __demand(__parallel)
      for i in pencil do
        token += d2dx2(points_x[i],pLU_x2[i])
      end
    --end

    --must_epoch
      __demand(__parallel)
      for i in pencil do
        token += d2dy2(points_y[i],pLU_y2[i])
      end
    --end

    --must_epoch
      __demand(__parallel)
      for i in pencil do
        token += d2dz2(points_z[i],pLU_z2[i])
      end
    --end
    
    wait_for(token)
    var ts_d2 = c.legion_get_current_time_in_micros() - ts_start
    
    var err_d2 = 0.0
    __demand(__parallel)
    for i in pencil do
      err_d2 += get_error_d2(points_x[i])
    end
    err_d2 = err_d2 / parallelism

    c.printf("Time to get the 1st derivatives: %12.5e\n", (ts_d1)*1e-6)
    c.printf("  Maximum error in x = %12.5e\n", err_x)
    c.printf("  Maximum error in y = %12.5e\n", err_y)
    c.printf("  Maximum error in z = %12.5e\n", err_z)
    c.printf("Time to get the 2nd derivatives: %12.5e\n", (ts_d2)*1e-6)
    c.printf("  Maximum error in laplacian = %12.5e\n", err_d2)
  end
  
end

regentlib.start(main)
