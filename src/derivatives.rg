import "regent"

require("fields")
local superlu = require("superlu_util")

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

end

local function poff(i, x, y, z, Nx, Ny, Nz)
  return rexpr int3d { x = (i.x + x + Nx)%Nx, y = (i.y + y + Ny)%Ny, z = (i.z + z + Nz)%Nz } end
end

local function make_stencil_pattern(points, f, index, a, b, c, Nx, Ny, Nz, onebydx, dir, der)
  local value

  if dir == 0 then      -- x direction stencil
    if der == 1 then
      value = rexpr       - c*points[ [poff(index, -3, 0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value - b*points[ [poff(index, -2, 0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value - a*points[ [poff(index, -1, 0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + a*points[ [poff(index,  1, 0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + b*points[ [poff(index,  2, 0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + c*points[ [poff(index,  3, 0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr onebydx * ( value ) end
    elseif der == 2 then
      value = rexpr         a*( points[ [poff(index, -1, 0, 0, Nx, Ny, Nz)] ].[f] - 2.0*points[ index ].[f] + points[ [poff(index, 1, 0, 0, Nx, Ny, Nz)] ].[f] ) end
      value = rexpr value + b*( points[ [poff(index, -2, 0, 0, Nx, Ny, Nz)] ].[f] - 2.0*points[ index ].[f] + points[ [poff(index, 2, 0, 0, Nx, Ny, Nz)] ].[f] ) end
      value = rexpr value + c*( points[ [poff(index, -3, 0, 0, Nx, Ny, Nz)] ].[f] - 2.0*points[ index ].[f] + points[ [poff(index, 3, 0, 0, Nx, Ny, Nz)] ].[f] ) end
      value = rexpr onebydx*onebydx * (value) end
    end
  elseif dir == 1 then  -- y direction stencil
    if der == 1 then
      value = rexpr       - c*points[ [poff(index, 0, -3, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value - b*points[ [poff(index, 0, -2, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value - a*points[ [poff(index, 0, -1, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + a*points[ [poff(index, 0,  1, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + b*points[ [poff(index, 0,  2, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + c*points[ [poff(index, 0,  3, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr onebydx * ( value ) end
    elseif der == 2 then
      value = rexpr         a*( points[ [poff(index, 0, -1, 0, Nx, Ny, Nz)] ].[f] - 2.0*points[ index ].[f] + points[ [poff(index, 0, 1, 0, Nx, Ny, Nz)] ].[f] ) end
      value = rexpr value + b*( points[ [poff(index, 0, -2, 0, Nx, Ny, Nz)] ].[f] - 2.0*points[ index ].[f] + points[ [poff(index, 0, 2, 0, Nx, Ny, Nz)] ].[f] ) end
      value = rexpr value + c*( points[ [poff(index, 0, -3, 0, Nx, Ny, Nz)] ].[f] - 2.0*points[ index ].[f] + points[ [poff(index, 0, 3, 0, Nx, Ny, Nz)] ].[f] ) end
      value = rexpr onebydx*onebydx * (value) end
    end
  elseif dir == 2 then  -- z direction stencil
    if der == 1 then
      value = rexpr       - c*points[ [poff(index, 0, 0, -3, Nx, Ny, Nz)] ].[f] end
      value = rexpr value - b*points[ [poff(index, 0, 0, -2, Nx, Ny, Nz)] ].[f] end
      value = rexpr value - a*points[ [poff(index, 0, 0, -1, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + a*points[ [poff(index, 0, 0,  1, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + b*points[ [poff(index, 0, 0,  2, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + c*points[ [poff(index, 0, 0,  3, Nx, Ny, Nz)] ].[f] end
      value = rexpr onebydx * ( value ) end
    elseif der == 2 then
      value = rexpr         a*( points[ [poff(index, 0, 0, -1, Nx, Ny, Nz)] ].[f] - 2.0*points[ index ].[f] + points[ [poff(index, 0, 0, 1, Nx, Ny, Nz)] ].[f] ) end
      value = rexpr value + b*( points[ [poff(index, 0, 0, -2, Nx, Ny, Nz)] ].[f] - 2.0*points[ index ].[f] + points[ [poff(index, 0, 0, 2, Nx, Ny, Nz)] ].[f] ) end
      value = rexpr value + c*( points[ [poff(index, 0, 0, -3, Nx, Ny, Nz)] ].[f] - 2.0*points[ index ].[f] + points[ [poff(index, 0, 0, 3, Nx, Ny, Nz)] ].[f] ) end
      value = rexpr onebydx*onebydx * (value) end
    end
  end
  return value
end

local function make_stencil_pattern_MND(points_c, points_e, f, index, a, b, c, Nx, Ny, Nz, onebydx, dir)
  local value
  if dir == 0 then      -- x direction stencil
      value = rexpr       - c*points_e[ [poff(index, -1, 0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value - b*points_c[ [poff(index, -1, 0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value - a*points_e[ [poff(index, -0, 0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + a*points_e[ [poff(index,  1, 0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + b*points_c[ [poff(index,  1, 0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + c*points_e[ [poff(index,  2, 0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr onebydx * ( value ) end
  elseif dir == 1 then  -- y direction stencil
      value = rexpr       - c*points_e[ [poff(index, 0, -1, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value - b*points_c[ [poff(index, 0, -1, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value - a*points_e[ [poff(index, 0, -0, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + a*points_e[ [poff(index, 0,  1, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + b*points_c[ [poff(index, 0,  1, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + c*points_e[ [poff(index, 0,  2, 0, Nx, Ny, Nz)] ].[f] end
      value = rexpr onebydx * ( value ) end
  elseif dir == 2 then  -- z direction stencil
      value = rexpr       - c*points_e[ [poff(index, 0, 0, -1, Nx, Ny, Nz)] ].[f] end
      value = rexpr value - b*points_c[ [poff(index, 0, 0, -1, Nx, Ny, Nz)] ].[f] end
      value = rexpr value - a*points_e[ [poff(index, 0, 0, -0, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + a*points_e[ [poff(index, 0, 0,  1, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + b*points_c[ [poff(index, 0, 0,  1, Nx, Ny, Nz)] ].[f] end
      value = rexpr value + c*points_e[ [poff(index, 0, 0,  2, Nx, Ny, Nz)] ].[f] end
      value = rexpr onebydx * ( value ) end
  end
  return value
end

local function make_stencil_x(r1, privileges_r1, f1, r2, privileges_r2, f2, Nx, Ny, Nz, onebydx, a, b, c, der)
  local rhs_x __demand(__inline) task rhs_x( [r1], [r2] )
  where
    [privileges_r1], [privileges_r2]
  do
    for i in r2 do
      r2[i].[f2] = [make_stencil_pattern(r1, f1, i, a, b, c, Nx, Ny, Nz, onebydx, 0, der)]
    end
  end
  return rhs_x
end

local function make_stencil_MND_x(r1, privileges_r1, f1, r2, privileges_r2, r3, privileges_r3, f3, Nx, Ny, Nz, onebydx, a, b, c, der)
  local rhs_x __demand(__inline) task rhs_x( [r1], [r2], [r3] )
  where
    [privileges_r1], [privileges_r2], [privileges_r3]
  do
    for i in r3 do
      r3[i].[f3] = [make_stencil_pattern_MND(r1, r2, f1, i, a, b, c, Nx, Ny, Nz, onebydx, 0, der)]
    end
  end
  return rhs_x
end

local function make_SolveXLU(r, privileges_r, f)
  local SolveXLU __demand(__inline) task SolveXLU( [r], LU : region(ispace(int3d), LU_struct) )
  where
    [privileges_r], reads(LU)
  do
    var bounds = [r].ispace.bounds
    var N = bounds.hi.x + 1
    var pr = LU.ispace.bounds.hi.y
    var pc = LU.ispace.bounds.hi.z
  
    for j = bounds.lo.y, bounds.hi.y+1 do
      for k = bounds.lo.z, bounds.hi.z+1 do
  
        -- Step 8
        [r][{1,j,k}].[f] = [r][{1,j,k}].[f] - LU[{1,pr,pc}].b*[r][{0,j,k}].[f]
        var sum1 : double = LU[{0,pr,pc}].k*[r][{0,j,k}].[f] + LU[{1,pr,pc}].k*[r][{1,j,k}].[f]
        var sum2 : double = LU[{0,pr,pc}].l*[r][{0,j,k}].[f] + LU[{1,pr,pc}].l*[r][{1,j,k}].[f]
  
        -- Step 9
        for i = 2,N-2 do
          [r][{i,j,k}].[f] = [r][{i,j,k}].[f] - LU[{i,pr,pc}].b*[r][{i-1,j,k}].[f] - LU[{i,pr,pc}].eg*[r][{i-2,j,k}].[f]
          sum1 += LU[{i,pr,pc}].k*[r][{i,j,k}].[f]
          sum2 += LU[{i,pr,pc}].l*[r][{i,j,k}].[f]
        end
  
        -- Step 10
        [r][{N-2,j,k}].[f] = [r][{N-2,j,k}].[f] - sum1
        [r][{N-1,j,k}].[f] = ( [r][{N-1,j,k}].[f] - sum2 - LU[{N-2,pr,pc}].l*[r][{N-2,j,k}].[f] )*LU[{N-1,pr,pc}].g
  
        -- Step 11
        [r][{N-2,j,k}].[f] = ( [r][{N-2,j,k}].[f] - LU[{N-2,pr,pc}].w*[r][{N-1,j,k}].[f] )*LU[{N-2,pr,pc}].g
        [r][{N-3,j,k}].[f] = ( [r][{N-3,j,k}].[f] - LU[{N-3,pr,pc}].v*[r][{N-2,j,k}].[f] - LU[{N-3,pr,pc}].w*[r][{N-1,j,k}].[f] )*LU[{N-3,pr,pc}].g
        [r][{N-4,j,k}].[f] = ( [r][{N-4,j,k}].[f] - LU[{N-4,pr,pc}].h*[r][{N-3,j,k}].[f] - LU[{N-4,pr,pc}].v*[r][{N-2,j,k}].[f] - LU[{N-4,pr,pc}].w*[r][{N-1,j,k}].[f] )*LU[{N-4,pr,pc}].g
        for i = N-5,-1,-1 do
          [r][{i,j,k}].[f] = ( [r][{i,j,k}].[f] - LU[{i,pr,pc}].h*[r][{i+1,j,k}].[f] - LU[{i,pr,pc}].ff*[r][{i+2,j,k}].[f] - LU[{i,pr,pc}].v*[r][{N-2,j,k}].[f] - LU[{i,pr,pc}].w*[r][{N-1,j,k}].[f] )*LU[{i,pr,pc}].g
        end
  
      end
    end
    return 1
  end
  return SolveXLU
end

function make_ddx(r_flux, f_flux, r_cnsr, f_cnsr, NX, NY, NZ, ONEBYDX, a, b, c)
  local privileges_r_flux   = regentlib.privilege(regentlib.reads,  r_flux, f_flux)
  local privileges_r_cnsr   = regentlib.privilege(regentlib.writes, r_cnsr, f_cnsr)

  local ComputeXRHS  = make_stencil_x(r_flux, privileges_r_flux, f_flux, r_cnsr, privileges_r_cnsr, f_cnsr, NX, NY, NZ, ONEBYDX, a, b, c, 1)
  local SolveXLU     = make_SolveXLU(r_cnsr, privileges_r_cnsr, f_cnsr)

  local task ddx( [r_flux],
                  [r_cnsr],
                  LU     : region(ispace(int3d), LU_struct) )
  where
    reads(LU), [privileges_r_flux], [privileges_r_cnsr]
  do
    [ComputeXRHS]([r_flux], [r_cnsr])
    var token = [SolveXLU]([r_cnsr],LU)
    token = r_cnsr[ [r_cnsr].ispace.bounds.lo ].[f_cnsr]
    return token
  end
  return ddx
end

function make_ddx_MND(r_flux, r_flux_e, f_flux, r_cnsr, f_cnsr, NX, NY, NZ, ONEBYDX, a, b, c)
  local privileges_r_flux   = regentlib.privilege(regentlib.reads,  r_flux,   f_flux)
  local privileges_r_flux_e = regentlib.privilege(regentlib.reads,  r_flux_e, f_flux)
  local privileges_r_cnsr   = regentlib.privilege(regentlib.writes, r_cnsr,   f_cnsr)

  local ComputeXRHS_MND  = make_stencil_MND_x(r_flux, privileges_r_flux, f_flux, r_flux_e, privileges_r_flux_e, r_cnsr, privileges_r_cnsr, f_cnsr, NX, NY, NZ, ONEBYDX, a, b, c, 1)
  local SolveXLU         = make_SolveXLU(r_cnsr, privileges_r_cnsr, f_cnsr)

  local task ddx_MND( [r_flux],
                      [r_flux_e],
                      [r_cnsr],
                      LU     : region(ispace(int3d), LU_struct) )
  where
    reads(LU), [privileges_r_flux], [privileges_r_flux_e], [privileges_r_cnsr]
  do
    [ComputeXRHS_MND]([r_flux], [r_flux_e], [r_cnsr])
    var token = [SolveXLU]([r_cnsr],LU)
    token = r_cnsr[ [r_cnsr].ispace.bounds.lo ].[f_cnsr]
    return token
  end
  return ddx_MND
end

-- task ddx_MND( r_flux   : region(ispace(int3d), conserved),
--               r_flux_e : region(ispace(int3d), conserved),
--               r_cnsr   : region(ispace(int3d), conserved),
--               LU       : region(ispace(int3d), LU_struct) )
-- where
--   reads(LU, r_flux.rho, r_flux_e.rho), reads writes(r_cnsr.rho)
-- do
--   ComputeXRHS_MND(r_flux, r_flux_e, r_cnsr)
--   var token = SolveXLU(r_cnsr,LU)
--   token = r_cnsr[r_cnsr.ispace.bounds.lo].rho
--   return token
-- end
