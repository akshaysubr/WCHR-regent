import "regent"

local c     = regentlib.c
local cmath = terralib.includec("math.h")
local PI    = cmath.M_PI

local superlu = require("superlu_util")

-- Grid dimensions
local NX = 16
local NY = 16
local NZ = 1

-- Domain size
local LX = 2.0*math.pi
local LY = 2.0*math.pi
local LZ = 2.0*math.pi

-- Grid spacing
local DX = LX / NX
local DY = LY / NY
local DZ = LZ / NZ

local ONEBYDX = 1.0 / DX
local ONEBYDY = 1.0 / DY
local ONEBYDZ = 1.0 / DZ

-- 10th order derivative parameters
local a10d1 = (17./12.)/2.
local b10d1 = (101./150.)/4.
local c10d1 = (1./100.)/6.
local alpha10d1 = 1./2.
local beta10d1  = 1./20.

fspace point {
  x_c : double,
  y_c : double,
  z_c : double,
  f   : double,
  df  : double,
  dX  : double
}

local function poff(i, x, y, z, Nx, Ny, Nz)
  return rexpr int3d { x = (i.x + x + Nx)%Nx, y = (i.y + y + Ny)%Ny, z = (i.z + z + Nz)%Nz } end
end

local function make_stencil_pattern(points, index, a10, b10, c10, Nx, Ny, Nz, onebydx, dir)
  local value

  if dir == 0 then      -- x direction stencil
      value = rexpr       - c10*points[ [poff(index, -3, 0, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value - b10*points[ [poff(index, -2, 0, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value - a10*points[ [poff(index, -1, 0, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value + a10*points[ [poff(index,  1, 0, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value + b10*points[ [poff(index,  2, 0, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value + c10*points[ [poff(index,  3, 0, 0, Nx, Ny, Nz)] ].f end
      value = rexpr onebydx * ( value ) end
  elseif dir == 1 then  -- y direction stencil
      value = rexpr       - c10*points[ [poff(index, 0, -3, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value - b10*points[ [poff(index, 0, -2, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value - a10*points[ [poff(index, 0, -1, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value + a10*points[ [poff(index, 0,  1, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value + b10*points[ [poff(index, 0,  2, 0, Nx, Ny, Nz)] ].f end
      value = rexpr value + c10*points[ [poff(index, 0,  3, 0, Nx, Ny, Nz)] ].f end
      value = rexpr onebydx * ( value ) end
  elseif dir == 2 then  -- z direction stencil
      value = rexpr       - c10*points[ [poff(index, 0, 0, -3, Nx, Ny, Nz)] ].f end
      value = rexpr value - b10*points[ [poff(index, 0, 0, -2, Nx, Ny, Nz)] ].f end
      value = rexpr value - a10*points[ [poff(index, 0, 0, -1, Nx, Ny, Nz)] ].f end
      value = rexpr value + a10*points[ [poff(index, 0, 0,  1, Nx, Ny, Nz)] ].f end
      value = rexpr value + b10*points[ [poff(index, 0, 0,  2, Nx, Ny, Nz)] ].f end
      value = rexpr value + c10*points[ [poff(index, 0, 0,  3, Nx, Ny, Nz)] ].f end
      value = rexpr onebydx * ( value ) end
  end
  return value
end

local function make_stencil_x(Nx, Ny, Nz, onebydx, a10, b10, c10)
  local rhs_x __demand(__inline) task rhs_x( points : region(ispace(int3d), point) )
  where
    reads(points.f), writes(points.dX)
  do
    for i in points do
      points[i].dX = [make_stencil_pattern(points, i, a10, b10, c10, Nx, Ny, Nz, onebydx, 0)]
    end
  end
  return rhs_x
end

local ComputeXRHS  = make_stencil_x(NX, NY, NZ, ONEBYDX, a10d1, b10d1, c10d1)

task ddx( points : region(ispace(int3d), point),
          matrix : region(ispace(int1d), superlu.CSR_matrix),
          slu    : region(ispace(int1d), superlu.c.superlu_vars_t) )
where
  reads(points.f), reads writes(points.dX, points.df, matrix, slu)
do
  var bounds = points.ispace.bounds
  var nx = bounds.hi.x + 1 - bounds.lo.x
  var ny = bounds.hi.y + 1 - bounds.lo.y
  var nz = bounds.hi.z + 1 - bounds.lo.z
  ComputeXRHS(points)
  superlu.MatrixSolve(__physical(points.dX), __fields(points.dX), __physical(points.df), __fields(points.df), points.bounds, matrix[0], nx, ny, nz, __physical(slu)[0], __fields(slu)[0], slu.bounds ) 
  var token = points[points.ispace.bounds.lo].df
  return token
end

task initialize_fields( points : region(ispace(int3d), point),
                        dx     : double,
                        dy     : double,
                        dz     : double )
where
  reads writes(points)
do
  for i in points do
    points[i].x_c = i.x * dx
    points[i].y_c = i.y * dy
    points[i].z_c = i.z * dz
    points[i].f = cmath.sin(points[i].x_c) * cmath.cos(points[i].y_c) * cmath.cos(points[i].z_c)
  end
  return 1
end

task check_error( points : region(ispace(int3d), point) )
where
  reads(points)
do
  var err : double = 0.0
  for i in points do
    var this_err = cmath.fabs( points[i].df - (cmath.cos(points[i].x_c) * cmath.cos(points[i].y_c) * cmath.cos(points[i].z_c)) )
    if this_err >= err then
      err = this_err
    end
  end
  return err
end

terra wait_for(x : int)
  return x
end

task main()

  var Nx : int64 = NX
  var Ny : int64 = NY
  var Nz : int64 = NZ

  var Lx : double = LX
  var Ly : double = LY
  var Lz : double = LZ

  var dx : double = DX
  var dy : double = DY
  var dz : double = DZ

  c.printf("================ Problem parameters ================\n")
  c.printf("           grid size  = %d x %d x %d\n", Nx, Ny, Nz )
  c.printf("         domain size  = %f x %f x %f\n", Lx, Ly, Lz )
  c.printf("           dx, dy, dz = %f, %f, %f\n", dx, dy, dz)
  c.printf("====================================================\n")

  var grid   = ispace(int3d, { x = Nx, y = Ny, z = Nz })
  var points = region(grid, point)

  var slu    = region(ispace(int1d, 1), superlu.c.superlu_vars_t)
  var matrix = region(ispace(int1d, 1), superlu.CSR_matrix)

  var token = 0 
  token += initialize_fields(points, dx, dy, dz)

  -- Initialize SuperLU stuff
  matrix[0] = superlu.initialize_matrix(alpha10d1, beta10d1, Nx, Ny, Nz)
  slu[0]    = superlu.initialize_superlu_vars( matrix[0], Nx*Ny*Nz, __physical(points.dX), __fields(points.dX),
                                               __physical(points.df), __fields(points.df), points.bounds)

  wait_for(token)
  var ts_start = c.legion_get_current_time_in_micros()

  -- Compute x derivative
  token += ddx(points, matrix, slu)

  wait_for(token)
  var ts_d1 = c.legion_get_current_time_in_micros() - ts_start
  c.printf("Time to get the 1st derivatives: %12.5e\n", (ts_d1)*1e-6)
  c.printf("Maximum error in 1st derivative: %12.5e\n", check_error(points))
    
end

regentlib.start(main)
