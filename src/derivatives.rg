import "regent"

require("fields")
local ip = require("interpolation")

-- Make the node and midpoint-node differencing tasks (Using pentadiagonal solver for this instead of tridiagonal solver)
alpha06d1 = 1.0/3.0
beta06d1  = 0.0
a06d1     = (14.0/ 9.0)/2.0
b06d1     = (1.0/ 9.0)/4.0
c06d1     = (0.0/100.0)/6.0

alpha06d1_1 = 3.0
p06d1_1     = 17.0/6.0
q06d1_1     = 3.0/2.0
r06d1_1     = 3.0/2.0
s06d1_1     = -1.0/6.0 

q06d1_2     = 3.0/4.0
alpha06d1_2 = 1.0/4.0

alpha06d1_2 = ((40.0*alpha06d1 - 1.0)*q06d1_1  + 7.0*(4.0*alpha06d1 -  1.0)*s06d1_1)/(16.0*(alpha06d1 + 2.0)*q06d1_1 + 8.0*(1.0 - 4.0*alpha06d1)*s06d1_1)
q06d1_2     = (1.0/3.0)*(alpha06d1_2 + 2.0)
r06d1_2     = (1.0/12.0)*(4.0*alpha06d1_2 - 1.0)
s06d1_2     =  0.0

w06d1_1 = (2.0*alpha06d1 + 1.0)/(2.0*(q06d1_1 + s06d1_1))
w06d1_2 = ((8.0*alpha06d1 + 7.0)*q06d1_1 - 6.0*(2.0*alpha06d1 + 1.0)*r06d1_1 + (8.0*alpha06d1 + 7.0)*s06d1_1)/(9.0*(q06d1_1 + s06d1_1))
w06d1_3 = (4.0*(alpha06d1 + 2.0)*q06d1_1 + 2*(1.0 - 4.0*alpha06d1)*s06d1_1) / (9.0*(q06d1_1 + s06d1_1))

alpha06d2 = 2.0/11.0
beta06d2  = 0.0
a06d2     = (12.0/11.0)/1.0
b06d2     = (3.0/11.0)/4.0
c06d2     = (0.0/100.0)/6.0

-- Compact MND finite difference
-- alpha06MND = -1.0/12.0
-- beta06MND  = 0.0
-- a06MND     = 16.0/9.0
-- b06MND     = (-17.0/18.0)/2.0
-- c06MND     = (0.0)/3.0

-- Compact staggered finite difference
alpha06MND = 9.0/62.0
beta06MND  = 0.0
a06MND     = 63.0/62.0
b06MND     = (0.0/18.0)/2.0
c06MND     = (17.0/62.0)/3.0

alpha06MND_LB = (40.0/31.0) * (9.0/80.0)
beta06MND_LB  = (40.0/31.0) * 0.0
a06MND_LB     = (40.0/31.0) * (1633.0/5376000.0)
b06MND_LB     = (40.0/31.0) * (9007.0/192000.0)
c06MND_LB     = (40.0/31.0) * (-29567.0/48000.0)
d06MND_LB     = (40.0/31.0) * (-65699.0/76800.0)
e06MND_LB     = (40.0/31.0) * (44033.0/24000.0)
f06MND_LB     = (40.0/31.0) * (-26353.0/38400.0)
g06MND_LB     = (40.0/31.0) * (104579.0/336000.0)
h06MND_LB     = (40.0/31.0) * (-27233.0/768000.0)

-- Explicit MND finite difference
-- alpha06MND = 0.0
-- beta06MND  = 0.0
-- a06MND     = 3.0/2.0
-- b06MND     = -3.0/10.0
-- c06MND     = 1.0/30.0

-- alpha06MND = 0.0
-- beta06MND  = 0.0
-- a06MND     = 1.0
-- b06MND     = 0.0
-- c06MND     = 0.0

-- alpha06MND_LB = 0.0
-- beta06MND_LB  = 0.0
-- a06MND_LB     = 0.0
-- b06MND_LB     = 0.0
-- c06MND_LB     = -1.0
-- d06MND_LB     = 0.0
-- e06MND_LB     = 1.0
-- f06MND_LB     = 0.0
-- g06MND_LB     = 0.0
-- h06MND_LB     = 0.0



-- Compact staggered finite difference flux difference form
alpha06MND_fd = (40.0/31.0) * (9.0/80.0)
beta06MND_fd  = (40.0/31.0) * 0.0
a06MND_fd     = (40.0/31.0) * (103.0/120.0)
b06MND_fd     = (40.0/31.0) * (0.0)
c06MND_fd     = (40.0/31.0) * (17.0/240.0)

alpha06MND_LB_fd = (80./71.0) * (9.0/80.0)
beta06MND_LB_fd  = (80./71.0) * 0.0
a06MND_LB_fd     = (80./71.0) * (-1633.0/5376000.0)
b06MND_LB_fd     = (80./71.0) * (-9007.0/192000.0)
c06MND_LB_fd     = (80./71.0) * (10989.0/16000.0)
d06MND_LB_fd     = (80./71.0) * (65699.0/76800.0)
e06MND_LB_fd     = (80./71.0) * (-7811.0/8000.0)
f06MND_LB_fd     = (80./71.0) * (26353.0/38400.0)
g06MND_LB_fd     = (80./71.0) * (-80779.0/336000.0)
h06MND_LB_fd     = (80./71.0) * (27233.0/768000.0)

-- alpha06MND_fd = 0.0
-- beta06MND_fd  = 0.0
-- a06MND_fd     = 1.0
-- b06MND_fd     = 0.0
-- c06MND_fd     = 0.0

-- alpha06MND_LB_fd = 0.0
-- beta06MND_LB_fd  = 0.0
-- a06MND_LB_fd     = 0.0
-- b06MND_LB_fd     = 0.0
-- c06MND_LB_fd     = 1.0
-- d06MND_LB_fd     = 0.0
-- e06MND_LB_fd     = 0.0
-- f06MND_LB_fd     = 0.0
-- g06MND_LB_fd     = 0.0
-- h06MND_LB_fd     = 0.0



function poff_wg( i, x, y, z, Nx, Ny, Nz )
  return rexpr int3d { x = (i.x - ip.n_ghosts + x + Nx)%Nx + ip.n_ghosts, y = (i.y - ip.n_ghosts + y + Ny)%Ny + ip.n_ghosts, z = (i.z - ip.n_ghosts + z + Nz)%Nz + ip.n_ghosts } end
end



task get_compact_matrix( mat      : region(ispace(int3d), LU_coeffs),
                         der      : int,
                         periodic : bool )
where
  writes( mat )
do
  var N = mat.ispace.bounds.hi.x + 1

  var e = beta06d1
  var a = alpha06d1
  var d = 1.0
  var c = alpha06d1
  var f = beta06d1

  if der == 2 then
    e = beta06d2
    a = alpha06d2
    d = 1.0
    c = alpha06d2
    f = beta06d2
  end

  for i in mat do
    mat[i].e = e
    mat[i].a = a
    mat[i].d = d
    mat[i].c = c
    mat[i].f = f

    -- Assuming tridiagonal :P
    if ((not periodic) and (i.x == 0)) then
      mat[i].e = 0.0
      mat[i].a = 0.0
    elseif ((not periodic) and (i.x == N-1)) then
      mat[i].c = 0.0
      mat[i].f = 0.0
    end
  end
end



task get_MND_matrix( mat      : region(ispace(int3d), LU_coeffs),
                     periodic : bool )
where
  writes( mat )
do
  var N = mat.ispace.bounds.hi.x + 1

  for i in mat do
    mat[i].e = beta06MND
    mat[i].a = alpha06MND
    mat[i].d = 1.0
    mat[i].c = alpha06MND
    mat[i].f = beta06MND

    if ((not periodic) and (i.x == 0)) then
      mat[i].e = 0.0
      mat[i].a = 0.0
      mat[i].c = alpha06MND_LB
      mat[i].f = beta06MND_LB
    elseif ((not periodic) and (i.x == N-1)) then
      mat[i].e = beta06MND_LB
      mat[i].a = alpha06MND_LB
      mat[i].c = 0.0
      mat[i].f = 0.0
    end
  end
end



task get_MND_matrix_fd( mat      : region(ispace(int3d), LU_coeffs),
                        periodic : bool )
where
  writes( mat )
do
  var N_e = mat.ispace.bounds.hi.x + 1

  for i in mat do
    mat[i].e = beta06MND_fd
    mat[i].a = alpha06MND_fd
    mat[i].d = 1.0
    mat[i].c = alpha06MND_fd
    mat[i].f = beta06MND_fd

    if ((not periodic) and (i.x == 0)) then
      mat[i].e = 0.0
      mat[i].a = 0.0
      mat[i].c = alpha06MND_LB_fd
      mat[i].f = beta06MND_LB_fd
    elseif ((not periodic) and (i.x == N_e-1)) then
      mat[i].e = beta06MND_LB_fd
      mat[i].a = alpha06MND_LB_fd
      mat[i].c = 0.0
      mat[i].f = 0.0
    end
  end
end



task get_LU_decomposition( LU  : region(ispace(int3d), LU_struct),
                           mat : region(ispace(int3d), LU_coeffs) )
where
  reads( mat ), reads writes( LU )
do
  var N : int64 = LU.ispace.bounds.hi.x + 1
  var pr = LU.ispace.bounds.hi.y
  var pc = LU.ispace.bounds.hi.z

  -- Step 1
  LU[{0,pr,pc}].g = mat[{0,pr,pc}].d
  LU[{1,pr,pc}].b = mat[{1,pr,pc}].a/LU[{0,pr,pc}].g
  LU[{0,pr,pc}].h = mat[{0,pr,pc}].c
  LU[{0,pr,pc}].k = mat[{N-2,pr,pc}].f/LU[{0,pr,pc}].g
  LU[{0,pr,pc}].w = mat[{0,pr,pc}].a
  LU[{0,pr,pc}].v = mat[{0,pr,pc}].e
  LU[{0,pr,pc}].l = mat[{N-1,pr,pc}].c/LU[{0,pr,pc}].g
  
  LU[{1,pr,pc}].g = mat[{1,pr,pc}].d - LU[{1,pr,pc}].b*LU[{0,pr,pc}].h
  LU[{1,pr,pc}].k = -LU[{0,pr,pc}].k*LU[{0,pr,pc}].h/LU[{1,pr,pc}].g
  LU[{1,pr,pc}].w = mat[{1,pr,pc}].e - LU[{1,pr,pc}].b*LU[{0,pr,pc}].w
  LU[{1,pr,pc}].v = -LU[{1,pr,pc}].b*LU[{0,pr,pc}].v
  LU[{1,pr,pc}].l = (mat[{N-1,pr,pc}].f - LU[{0,pr,pc}].l*LU[{0,pr,pc}].h) / LU[{1,pr,pc}].g
  LU[{1,pr,pc}].h = mat[{1,pr,pc}].c - LU[{1,pr,pc}].b*mat[{0,pr,pc}].f

  -- Step 2
  for i = 2,N-3 do
    LU[{i,pr,pc}].b = ( mat[{i,pr,pc}].a - ( mat[{i,pr,pc}].e/LU[{i-2,pr,pc}].g )*LU[{i-2,pr,pc}].h ) / LU[{i-1,pr,pc}].g
    LU[{i,pr,pc}].h = mat[{i,pr,pc}].c - LU[{i,pr,pc}].b*mat[{i-1,pr,pc}].f
    LU[{i,pr,pc}].g = mat[{i,pr,pc}].d - ( mat[{i,pr,pc}].e/LU[{i-2,pr,pc}].g )*mat[{i-2,pr,pc}].f - LU[{i,pr,pc}].b*LU[{i-1,pr,pc}].h
  end

  -- Step 3
  LU[{N-3,pr,pc}].b = ( mat[{N-3,pr,pc}].a - ( mat[{N-3,pr,pc}].e/LU[{N-5,pr,pc}].g )*LU[{N-5,pr,pc}].h ) / LU[{N-4,pr,pc}].g
  LU[{N-3,pr,pc}].g = mat[{N-3,pr,pc}].d - ( mat[{N-3,pr,pc}].e/LU[{N-5,pr,pc}].g )*mat[{N-5,pr,pc}].f - LU[{N-3,pr,pc}].b*LU[{N-4,pr,pc}].h

  -- Step 4
  for i = 2,N-4 do
    LU[{i,pr,pc}].k = -( LU[{i-2,pr,pc}].k*mat[{i-2,pr,pc}].f + LU[{i-1,pr,pc}].k*LU[{i-1,pr,pc}].h ) / LU[{i,pr,pc}].g
    LU[{i,pr,pc}].v = -( mat[{i,pr,pc}].e/LU[{i-2,pr,pc}].g )*LU[{i-2,pr,pc}].v - LU[{i,pr,pc}].b*LU[{i-1,pr,pc}].v
  end

  -- Step 5
  LU[{N-4,pr,pc}].k = ( mat[{N-2,pr,pc}].e - LU[{N-6,pr,pc}].k*mat[{N-6,pr,pc}].f - LU[{N-5,pr,pc}].k*LU[{N-5,pr,pc}].h ) / LU[{N-4,pr,pc}].g
  LU[{N-3,pr,pc}].k = ( mat[{N-2,pr,pc}].a - LU[{N-5,pr,pc}].k*mat[{N-5,pr,pc}].f - LU[{N-4,pr,pc}].k*LU[{N-4,pr,pc}].h ) / LU[{N-3,pr,pc}].g
  LU[{N-4,pr,pc}].v = mat[{N-4,pr,pc}].f  - ( mat[{N-4,pr,pc}].e/LU[{N-6,pr,pc}].g )*LU[{N-6,pr,pc}].v - LU[{N-4,pr,pc}].b*LU[{N-5,pr,pc}].v
  LU[{N-3,pr,pc}].v = mat[{N-3,pr,pc}].c - ( mat[{N-3,pr,pc}].e/LU[{N-5,pr,pc}].g )*LU[{N-5,pr,pc}].v - LU[{N-3,pr,pc}].b*LU[{N-4,pr,pc}].v
  LU[{N-2,pr,pc}].g = mat[{N-2,pr,pc}].d
  for i = 0,N-2 do
    LU[{N-2,pr,pc}].g -= LU[{i,pr,pc}].k*LU[{i,pr,pc}].v
  end

  -- Step 6
  for i = 2,N-3 do
    LU[{i,pr,pc}].w = -( mat[{i,pr,pc}].e/LU[{i-2,pr,pc}].g )*LU[{i-2,pr,pc}].w - LU[{i,pr,pc}].b*LU[{i-1,pr,pc}].w
    LU[{i,pr,pc}].l = -( LU[{i-2,pr,pc}].l*mat[{i-2,pr,pc}].f + LU[{i-1,pr,pc}].l*LU[{i-1,pr,pc}].h ) / LU[{i,pr,pc}].g
  end

  -- Step 7
  LU[{N-3,pr,pc}].w = mat[{N-3,pr,pc}].f - ( mat[{N-3,pr,pc}].e/LU[{N-5,pr,pc}].g )*LU[{N-5,pr,pc}].w - LU[{N-3,pr,pc}].b*LU[{N-4,pr,pc}].w
  LU[{N-2,pr,pc}].w = mat[{N-2,pr,pc}].c
  for i = 0,N-2 do
    LU[{N-2,pr,pc}].w -= LU[{i,pr,pc}].k*LU[{i,pr,pc}].w
  end
  LU[{N-3,pr,pc}].l = ( mat[{N-1,pr,pc}].e - LU[{N-5,pr,pc}].l*mat[{N-5,pr,pc}].f - LU[{N-4,pr,pc}].l*LU[{N-4,pr,pc}].h ) / LU[{N-3,pr,pc}].g
  LU[{N-2,pr,pc}].l = mat[{N-1,pr,pc}].a
  for i = 0,N-2 do
    LU[{N-2,pr,pc}].l -= LU[{i,pr,pc}].l*LU[{i,pr,pc}].v
  end
  LU[{N-2,pr,pc}].l = LU[{N-2,pr,pc}].l / LU[{N-2,pr,pc}].g
  LU[{N-1,pr,pc}].g = mat[{N-1,pr,pc}].d
  for i = 0,N-1 do
    LU[{N-1,pr,pc}].g -= LU[{i,pr,pc}].l*LU[{i,pr,pc}].w
  end

  -- Set eg = e/g
  for i = 2,N-2 do
    LU[{i,pr,pc}].eg = mat[{i,pr,pc}].e/LU[{i-2,pr,pc}].g
  end

  -- Set ff = f
  for i = 0,N-4 do
    LU[{i,pr,pc}].ff = mat[{i,pr,pc}].f
  end

  -- Set g = 1/g
  for i = 0,N do
    LU[{i,pr,pc}].g = 1.0/LU[{i,pr,pc}].g
  end

  return 1
end



local function make_stencil_pattern( points, f, index, nx, ny, nz, onebydx, dir, der )
  local value

  local a = a06d1
  local b = b06d1
  local c = c06d1

  if der == 1 then
    a = a06d1
    b = b06d1
    c = c06d1
  elseif der == 2 then
    a = a06d2
    b = b06d2
    c = c06d2
  end

  if dir == 0 then      -- x direction stencil
    if der == 1 then
      value = rexpr       - c*points[ [poff_wg(index, -3, 0, 0, nx, ny, nz)] ].[f] end
      value = rexpr value - b*points[ [poff_wg(index, -2, 0, 0, nx, ny, nz)] ].[f] end
      value = rexpr value - a*points[ [poff_wg(index, -1, 0, 0, nx, ny, nz)] ].[f] end
      value = rexpr value + a*points[ [poff_wg(index,  1, 0, 0, nx, ny, nz)] ].[f] end
      value = rexpr value + b*points[ [poff_wg(index,  2, 0, 0, nx, ny, nz)] ].[f] end
      value = rexpr value + c*points[ [poff_wg(index,  3, 0, 0, nx, ny, nz)] ].[f] end
      value = rexpr onebydx * ( value ) end
    elseif der == 2 then
      value = rexpr         a*( points[ [poff_wg(index, -1, 0, 0, nx, ny, nz)] ].[f] - 2.0*points[index].[f] + points[ [poff_wg(index, 1, 0, 0, nx, ny, nz)] ].[f] ) end
      value = rexpr value + b*( points[ [poff_wg(index, -2, 0, 0, nx, ny, nz)] ].[f] - 2.0*points[index].[f] + points[ [poff_wg(index, 2, 0, 0, nx, ny, nz)] ].[f] ) end
      value = rexpr value + c*( points[ [poff_wg(index, -3, 0, 0, nx, ny, nz)] ].[f] - 2.0*points[index].[f] + points[ [poff_wg(index, 3, 0, 0, nx, ny, nz)] ].[f] ) end
      value = rexpr onebydx*onebydx * (value) end
    end
  elseif dir == 1 then  -- y direction stencil
    if der == 1 then
      value = rexpr       - c*points[ [poff_wg(index, 0, -3, 0, nx, ny, nz)] ].[f] end
      value = rexpr value - b*points[ [poff_wg(index, 0, -2, 0, nx, ny, nz)] ].[f] end
      value = rexpr value - a*points[ [poff_wg(index, 0, -1, 0, nx, ny, nz)] ].[f] end
      value = rexpr value + a*points[ [poff_wg(index, 0,  1, 0, nx, ny, nz)] ].[f] end
      value = rexpr value + b*points[ [poff_wg(index, 0,  2, 0, nx, ny, nz)] ].[f] end
      value = rexpr value + c*points[ [poff_wg(index, 0,  3, 0, nx, ny, nz)] ].[f] end
      value = rexpr onebydx * ( value ) end
    elseif der == 2 then
      value = rexpr         a*( points[ [poff_wg(index, 0, -1, 0, nx, ny, nz)] ].[f] - 2.0*points[index].[f] + points[ [poff_wg(index, 0, 1, 0, nx, ny, nz)] ].[f] ) end
      value = rexpr value + b*( points[ [poff_wg(index, 0, -2, 0, nx, ny, nz)] ].[f] - 2.0*points[index].[f] + points[ [poff_wg(index, 0, 2, 0, nx, ny, nz)] ].[f] ) end
      value = rexpr value + c*( points[ [poff_wg(index, 0, -3, 0, nx, ny, nz)] ].[f] - 2.0*points[index].[f] + points[ [poff_wg(index, 0, 3, 0, nx, ny, nz)] ].[f] ) end
      value = rexpr onebydx*onebydx * (value) end
    end
  elseif dir == 2 then  -- z direction stencil
    if der == 1 then
      value = rexpr       - c*points[ [poff_wg(index, 0, 0, -3, nx, ny, nz)] ].[f] end
      value = rexpr value - b*points[ [poff_wg(index, 0, 0, -2, nx, ny, nz)] ].[f] end
      value = rexpr value - a*points[ [poff_wg(index, 0, 0, -1, nx, ny, nz)] ].[f] end
      value = rexpr value + a*points[ [poff_wg(index, 0, 0,  1, nx, ny, nz)] ].[f] end
      value = rexpr value + b*points[ [poff_wg(index, 0, 0,  2, nx, ny, nz)] ].[f] end
      value = rexpr value + c*points[ [poff_wg(index, 0, 0,  3, nx, ny, nz)] ].[f] end
      value = rexpr onebydx * ( value ) end
    elseif der == 2 then
      value = rexpr         a*( points[ [poff_wg(index, 0, 0, -1, nx, ny, nz)] ].[f] - 2.0*points[index].[f] + points[ [poff_wg(index, 0, 0, 1, nx, ny, nz)] ].[f] ) end
      value = rexpr value + b*( points[ [poff_wg(index, 0, 0, -2, nx, ny, nz)] ].[f] - 2.0*points[index].[f] + points[ [poff_wg(index, 0, 0, 2, nx, ny, nz)] ].[f] ) end
      value = rexpr value + c*( points[ [poff_wg(index, 0, 0, -3, nx, ny, nz)] ].[f] - 2.0*points[index].[f] + points[ [poff_wg(index, 0, 0, 3, nx, ny, nz)] ].[f] ) end
      value = rexpr onebydx*onebydx * (value) end
    end
  end
  return value
end



local function make_stencil_pattern_MND_interior( points_c, points_e, f, index, Nx, Ny, Nz, onebydx, dir, periodic )
  local value
  local Ne = 0
  local Nx_g = Nx + 2*ip.n_ghosts
  local Ny_g = Ny + 2*ip.n_ghosts
  local Nz_g = Nz + 2*ip.n_ghosts

  local a = a06MND
  local b = b06MND
  local c = c06MND

  if dir == 0 then      -- x direction stencil
      local index_e = rexpr int3d { x = index.x - ip.n_ghosts, y = index.y, z = index.z} end
      if periodic then
        Ne = Nx
      else
        Ne = Nx+1
      end

      value = rexpr       - b*points_c[ int3d {index.x - 1, index.y, index.z} ].[f] end
      value = rexpr value + b*points_c[ int3d {index.x + 1, index.y, index.z} ].[f] end

      value = rexpr value - c*points_e[ [poff(index_e, -1, 0, 0, Ne, Ny_g, Nz_g)] ].[f] end
      value = rexpr value - a*points_e[ [poff(index_e, -0, 0, 0, Ne, Ny_g, Nz_g)] ].[f] end
      value = rexpr value + a*points_e[ [poff(index_e,  1, 0, 0, Ne, Ny_g, Nz_g)] ].[f] end
      value = rexpr value + c*points_e[ [poff(index_e,  2, 0, 0, Ne, Ny_g, Nz_g)] ].[f] end
      value = rexpr onebydx * ( value ) end
  elseif dir == 1 then  -- y direction stencil
      local index_e = rexpr int3d { x = index.x, y = index.y - ip.n_ghosts, z = index.z} end
      if periodic then
        Ne = Ny
      else
        Ne = Ny+1
      end

      value = rexpr       - b*points_c[ int3d {index.x, index.y - 1, index.z} ].[f] end
      value = rexpr value + b*points_c[ int3d {index.x, index.y + 1, index.z} ].[f] end

      value = rexpr value - c*points_e[ [poff(index_e, 0, -1, 0, Nx_g, Ne, Nz_g)] ].[f] end
      value = rexpr value - a*points_e[ [poff(index_e, 0, -0, 0, Nx_g, Ne, Nz_g)] ].[f] end
      value = rexpr value + a*points_e[ [poff(index_e, 0,  1, 0, Nx_g, Ne, Nz_g)] ].[f] end
      value = rexpr value + c*points_e[ [poff(index_e, 0,  2, 0, Nx_g, Ne, Nz_g)] ].[f] end
      value = rexpr onebydx * ( value ) end
  elseif dir == 2 then  -- z direction stencil
      local index_e = rexpr int3d { x = index.x, y = index.y, z = index.z - ip.n_ghosts} end
      if periodic then
        Ne = Nz
      else
        Ne = Nz+1
      end

      value = rexpr       - b*points_c[ int3d {index.x, index.y, index.z - 1} ].[f] end
      value = rexpr value + b*points_c[ int3d {index.x, index.y, index.z + 1} ].[f] end

      value = rexpr value - c*points_e[ [poff(index_e, 0, 0, -1, Nx_g, Ny_g, Ne)] ].[f] end
      value = rexpr value - a*points_e[ [poff(index_e, 0, 0, -0, Nx_g, Ny_g, Ne)] ].[f] end
      value = rexpr value + a*points_e[ [poff(index_e, 0, 0,  1, Nx_g, Ny_g, Ne)] ].[f] end
      value = rexpr value + c*points_e[ [poff(index_e, 0, 0,  2, Nx_g, Ny_g, Ne)] ].[f] end
      value = rexpr onebydx * ( value ) end
  end
  return value
end



local function make_stencil_pattern_MND_LB( points_c, points_e, fs, index, onebydx, dir )
  local value

  local a = a06MND_LB
  local b = b06MND_LB
  local c = c06MND_LB
  local d = d06MND_LB
  local e = e06MND_LB
  local f = f06MND_LB
  local g = g06MND_LB
  local h = h06MND_LB

  if dir == 0 then      -- x direction stencil
      local index_e = rexpr int3d { x = index.x - ip.n_ghosts, y = index.y, z = index.z} end

      value = rexpr         a*points_c[ int3d {index.x   - 2, index.y,   index.z  } ].[fs] end
      value = rexpr value + b*points_c[ int3d {index.x   - 1, index.y,   index.z  } ].[fs] end
      value = rexpr value + c*points_e[ int3d {index_e.x + 0, index_e.y, index_e.z} ].[fs] end
      value = rexpr value + d*points_c[ int3d {index.x   - 0, index.y,   index.z  } ].[fs] end
      value = rexpr value + e*points_e[ int3d {index_e.x + 1, index_e.y, index_e.z} ].[fs] end
      value = rexpr value + f*points_c[ int3d {index.x   + 1, index.y,   index.z  } ].[fs] end
      value = rexpr value + g*points_e[ int3d {index_e.x + 2, index_e.y, index_e.z} ].[fs] end
      value = rexpr value + h*points_c[ int3d {index.x   + 2, index.y,   index.z  } ].[fs] end
      value = rexpr onebydx * ( value ) end
  elseif dir == 1 then  -- y direction stencil
      local index_e = rexpr int3d { x = index.x, y = index.y - ip.n_ghosts, z = index.z} end

      value = rexpr         a*points_c[ int3d {index.x,   index.y   - 2, index.z  } ].[fs] end
      value = rexpr value + b*points_c[ int3d {index.x,   index.y   - 1, index.z  } ].[fs] end
      value = rexpr value + c*points_e[ int3d {index_e.x, index_e.y + 0, index_e.z} ].[fs] end
      value = rexpr value + d*points_c[ int3d {index.x,   index.y   - 0, index.z  } ].[fs] end
      value = rexpr value + e*points_e[ int3d {index_e.x, index_e.y + 1, index_e.z} ].[fs] end
      value = rexpr value + f*points_c[ int3d {index.x,   index.y   + 1, index.z  } ].[fs] end
      value = rexpr value + g*points_e[ int3d {index_e.x, index_e.y + 2, index_e.z} ].[fs] end
      value = rexpr value + h*points_c[ int3d {index.x,   index.y   + 2, index.z  } ].[fs] end
      value = rexpr onebydx * ( value ) end
  elseif dir == 2 then  -- z direction stencil
      local index_e = rexpr int3d { x = index.x, y = index.y, z = index.z - ip.n_ghosts} end

      value = rexpr         a*points_c[ int3d {index.x,   index.y  , index.z   - 2} ].[fs] end
      value = rexpr value + b*points_c[ int3d {index.x,   index.y  , index.z   - 1} ].[fs] end
      value = rexpr value + c*points_e[ int3d {index_e.x, index_e.y, index_e.z + 0} ].[fs] end
      value = rexpr value + d*points_c[ int3d {index.x,   index.y  , index.z   - 0} ].[fs] end
      value = rexpr value + e*points_e[ int3d {index_e.x, index_e.y, index_e.z + 1} ].[fs] end
      value = rexpr value + f*points_c[ int3d {index.x,   index.y  , index.z   + 1} ].[fs] end
      value = rexpr value + g*points_e[ int3d {index_e.x, index_e.y, index_e.z + 2} ].[fs] end
      value = rexpr value + h*points_c[ int3d {index.x,   index.y  , index.z   + 2} ].[fs] end
      value = rexpr onebydx * ( value ) end
  end
  return value
end



local function make_stencil_pattern_MND_RB( points_c, points_e, fs, index, onebydx, dir )
  local value

  local a = a06MND_LB
  local b = b06MND_LB
  local c = c06MND_LB
  local d = d06MND_LB
  local e = e06MND_LB
  local f = f06MND_LB
  local g = g06MND_LB
  local h = h06MND_LB

  if dir == 0 then      -- x direction stencil
      local index_e = rexpr int3d { x = index.x - ip.n_ghosts, y = index.y, z = index.z} end

      value = rexpr       - a*points_c[ int3d {index.x   + 2, index.y,   index.z  } ].[fs] end
      value = rexpr value - b*points_c[ int3d {index.x   + 1, index.y,   index.z  } ].[fs] end
      value = rexpr value - c*points_e[ int3d {index_e.x + 1, index_e.y, index_e.z} ].[fs] end
      value = rexpr value - d*points_c[ int3d {index.x   + 0, index.y,   index.z  } ].[fs] end
      value = rexpr value - e*points_e[ int3d {index_e.x - 0, index_e.y, index_e.z} ].[fs] end
      value = rexpr value - f*points_c[ int3d {index.x   - 1, index.y,   index.z  } ].[fs] end
      value = rexpr value - g*points_e[ int3d {index_e.x - 1, index_e.y, index_e.z} ].[fs] end
      value = rexpr value - h*points_c[ int3d {index.x   - 2, index.y,   index.z  } ].[fs] end
      value = rexpr onebydx * ( value ) end
  elseif dir == 1 then  -- y direction stencil
      local index_e = rexpr int3d { x = index.x, y = index.y - ip.n_ghosts, z = index.z} end

      value = rexpr       - a*points_c[ int3d {index.x,   index.y   + 2, index.z  } ].[fs] end
      value = rexpr value - b*points_c[ int3d {index.x,   index.y   + 1, index.z  } ].[fs] end
      value = rexpr value - c*points_e[ int3d {index_e.x, index_e.y + 1, index_e.z} ].[fs] end
      value = rexpr value - d*points_c[ int3d {index.x,   index.y   + 0, index.z  } ].[fs] end
      value = rexpr value - e*points_e[ int3d {index_e.x, index_e.y - 0, index_e.z} ].[fs] end
      value = rexpr value - f*points_c[ int3d {index.x,   index.y   - 1, index.z  } ].[fs] end
      value = rexpr value - g*points_e[ int3d {index_e.x, index_e.y - 1, index_e.z} ].[fs] end
      value = rexpr value - h*points_c[ int3d {index.x,   index.y   - 2, index.z  } ].[fs] end
      value = rexpr onebydx * ( value ) end
  elseif dir == 2 then  -- z direction stencil
      local index_e = rexpr int3d { x = index.x, y = index.y, z = index.z - ip.n_ghosts} end

      value = rexpr       - a*points_c[ int3d {index.x,   index.y  , index.z   + 2} ].[fs] end
      value = rexpr value - b*points_c[ int3d {index.x,   index.y  , index.z   + 1} ].[fs] end
      value = rexpr value - c*points_e[ int3d {index_e.x, index_e.y, index_e.z + 1} ].[fs] end
      value = rexpr value - d*points_c[ int3d {index.x,   index.y  , index.z   + 0} ].[fs] end
      value = rexpr value - e*points_e[ int3d {index_e.x, index_e.y, index_e.z - 0} ].[fs] end
      value = rexpr value - f*points_c[ int3d {index.x,   index.y  , index.z   - 1} ].[fs] end
      value = rexpr value - g*points_e[ int3d {index_e.x, index_e.y, index_e.z - 1} ].[fs] end
      value = rexpr value - h*points_c[ int3d {index.x,   index.y  , index.z   - 2} ].[fs] end
      value = rexpr onebydx * ( value ) end
  end
  return value
end



local function make_stencil_pattern_MND_fd_interior( points_c, points_e, f, index, Nx, Ny, Nz, dir, periodic )
  local value
  local Ne = 0
  local Nx_g = Nx + 2*ip.n_ghosts
  local Ny_g = Ny + 2*ip.n_ghosts
  local Nz_g = Nz + 2*ip.n_ghosts

  local a = a06MND_fd
  local b = b06MND_fd
  local c = c06MND_fd

  if dir == 0 then      -- x direction stencil
      local index_c = rexpr int3d { x = index.x + ip.n_ghosts, y = index.y, z = index.z} end
      local index_e = rexpr index end

      if periodic then
        Ne = Nx
      else
        Ne = Nx + 1
      end

      value = rexpr         b*points_c[ int3d {index_c.x - 1, index_c.y, index_c.z} ].[f] end
      value = rexpr value + b*points_c[ int3d {index_c.x,     index_c.y, index_c.z} ].[f] end

      value = rexpr value + c*points_e[ [poff(index_e, -1, 0, 0, Ne, Ny_g, Nz_g)] ].[f] end
      value = rexpr value + a*points_e[ [poff(index_e,  0, 0, 0, Ne, Ny_g, Nz_g)] ].[f] end
      value = rexpr value + c*points_e[ [poff(index_e,  1, 0, 0, Ne, Ny_g, Nz_g)] ].[f] end
  elseif dir == 1 then  -- y direction stencil
      local index_c = rexpr int3d { x = index.x, y = index.y + ip.n_ghosts, z = index.z} end
      local index_e = rexpr index end

      if periodic then
        Ne = Ny
      else
        Ne = Ny + 1
      end

      value = rexpr         b*points_c[ int3d {index_c.x, index_c.y - 1, index_c.z} ].[f] end
      value = rexpr value + b*points_c[ int3d {index_c.x, index_c.y,     index_c.z} ].[f] end

      value = rexpr value + c*points_e[ [poff(index_e, 0, -1, 0, Nx_g, Ne, Nz_g)] ].[f] end
      value = rexpr value + a*points_e[ [poff(index_e, 0,  0, 0, Nx_g, Ne, Nz_g)] ].[f] end
      value = rexpr value + c*points_e[ [poff(index_e, 0,  1, 0, Nx_g, Ne, Nz_g)] ].[f] end
  elseif dir == 2 then  -- z direction stencil
      local index_c = rexpr int3d { x = index.x, y = index.y, z = index.z + ip.n_ghosts} end
      local index_e = rexpr index end

      if periodic then
        Ne = Nz
      else
        Ne = Nz + 1
      end

      value = rexpr         b*points_c[ int3d {index_c.x, index_c.y, index_c.z - 1} ].[f] end
      value = rexpr value + b*points_c[ int3d {index_c.x, index_c.y, index_c.z    } ].[f] end

      value = rexpr value + c*points_e[ [poff(index_e, 0, 0, -1, Nx_g, Ny_g, Ne)] ].[f] end
      value = rexpr value + a*points_e[ [poff(index_e, 0, 0,  0, Nx_g, Ny_g, Ne)] ].[f] end
      value = rexpr value + c*points_e[ [poff(index_e, 0, 0,  1, Nx_g, Ny_g, Ne)] ].[f] end
  end
  return value
end



local function make_stencil_pattern_MND_fd_LB( points_c, points_e, fs, index, dir )
  local value

  local a = a06MND_LB_fd
  local b = b06MND_LB_fd
  local c = c06MND_LB_fd
  local d = d06MND_LB_fd
  local e = e06MND_LB_fd
  local f = f06MND_LB_fd
  local g = g06MND_LB_fd
  local h = h06MND_LB_fd

  if dir == 0 then      -- x direction stencil
      local index_c = rexpr int3d { x = index.x + ip.n_ghosts, y = index.y, z = index.z} end
      local index_e = rexpr index end

      value = rexpr         a*points_c[ int3d {index_c.x - 2, index_c.y, index_c.z} ].[fs] end
      value = rexpr value + b*points_c[ int3d {index_c.x - 1, index_c.y, index_c.z} ].[fs] end
      value = rexpr value + c*points_e[ int3d {index_e.x + 0, index_e.y, index_e.z} ].[fs] end
      value = rexpr value + d*points_c[ int3d {index_c.x - 0, index_c.y, index_c.z} ].[fs] end
      value = rexpr value + e*points_e[ int3d {index_e.x + 1, index_e.y, index_e.z} ].[fs] end
      value = rexpr value + f*points_c[ int3d {index_c.x + 1, index_c.y, index_c.z} ].[fs] end
      value = rexpr value + g*points_e[ int3d {index_e.x + 2, index_e.y, index_e.z} ].[fs] end
      value = rexpr value + h*points_c[ int3d {index_c.x + 2, index_c.y, index_c.z} ].[fs] end
  elseif dir == 1 then  -- y direction stencil
      local index_c = rexpr int3d { x = index.x, y = index.y + ip.n_ghosts, z = index.z} end
      local index_e = rexpr index end

      value = rexpr         a*points_c[ int3d {index_c.x, index_c.y - 2, index_c.z} ].[fs] end
      value = rexpr value + b*points_c[ int3d {index_c.x, index_c.y - 1, index_c.z} ].[fs] end
      value = rexpr value + c*points_e[ int3d {index_e.x, index_e.y + 0, index_e.z} ].[fs] end
      value = rexpr value + d*points_c[ int3d {index_c.x, index_c.y - 0, index_c.z} ].[fs] end
      value = rexpr value + e*points_e[ int3d {index_e.x, index_e.y + 1, index_e.z} ].[fs] end
      value = rexpr value + f*points_c[ int3d {index_c.x, index_c.y + 1, index_c.z} ].[fs] end
      value = rexpr value + g*points_e[ int3d {index_e.x, index_e.y + 2, index_e.z} ].[fs] end
      value = rexpr value + h*points_c[ int3d {index_c.x, index_c.y + 2, index_c.z} ].[fs] end
  elseif dir == 2 then  -- z direction stencil
      local index_c = rexpr int3d { x = index.x, y = index.y, z = index.z + ip.n_ghosts} end
      local index_e = rexpr index end

      value = rexpr         a*points_c[ int3d {index_c.x, index_c.y, index_c.z - 2} ].[fs] end
      value = rexpr value + b*points_c[ int3d {index_c.x, index_c.y, index_c.z - 1} ].[fs] end
      value = rexpr value + c*points_e[ int3d {index_e.x, index_e.y, index_e.z + 0} ].[fs] end
      value = rexpr value + d*points_c[ int3d {index_c.x, index_c.y, index_c.z - 0} ].[fs] end
      value = rexpr value + e*points_e[ int3d {index_e.x, index_e.y, index_e.z + 1} ].[fs] end
      value = rexpr value + f*points_c[ int3d {index_c.x, index_c.y, index_c.z + 1} ].[fs] end
      value = rexpr value + g*points_e[ int3d {index_e.x, index_e.y, index_e.z + 2} ].[fs] end
      value = rexpr value + h*points_c[ int3d {index_c.x, index_c.y, index_c.z + 2} ].[fs] end
  end
  return value
end



local function make_stencil_pattern_MND_fd_RB( points_c, points_e, fs, index, dir )
  local value

  local a = a06MND_LB_fd
  local b = b06MND_LB_fd
  local c = c06MND_LB_fd
  local d = d06MND_LB_fd
  local e = e06MND_LB_fd
  local f = f06MND_LB_fd
  local g = g06MND_LB_fd
  local h = h06MND_LB_fd

  if dir == 0 then      -- x direction stencil
      local index_c = rexpr int3d { x = index.x + ip.n_ghosts - 1, y = index.y, z = index.z} end
      local index_e = rexpr index end

      value = rexpr         a*points_c[ int3d {index_c.x + 2, index_c.y, index_c.z} ].[fs] end
      value = rexpr value + b*points_c[ int3d {index_c.x + 1, index_c.y, index_c.z} ].[fs] end
      value = rexpr value + c*points_e[ int3d {index_e.x + 0, index_e.y, index_e.z} ].[fs] end
      value = rexpr value + d*points_c[ int3d {index_c.x + 0, index_c.y, index_c.z} ].[fs] end
      value = rexpr value + e*points_e[ int3d {index_e.x - 1, index_e.y, index_e.z} ].[fs] end
      value = rexpr value + f*points_c[ int3d {index_c.x - 1, index_c.y, index_c.z} ].[fs] end
      value = rexpr value + g*points_e[ int3d {index_e.x - 2, index_e.y, index_e.z} ].[fs] end
      value = rexpr value + h*points_c[ int3d {index_c.x - 2, index_c.y, index_c.z} ].[fs] end
  elseif dir == 1 then  -- y direction stencil
      local index_c = rexpr int3d { x = index.x, y = index.y + ip.n_ghosts - 1, z = index.z} end
      local index_e = rexpr index end

      value = rexpr         a*points_c[ int3d {index_c.x, index_c.y + 2, index_c.z} ].[fs] end
      value = rexpr value + b*points_c[ int3d {index_c.x, index_c.y + 1, index_c.z} ].[fs] end
      value = rexpr value + c*points_e[ int3d {index_e.x, index_e.y + 0, index_e.z} ].[fs] end
      value = rexpr value + d*points_c[ int3d {index_c.x, index_c.y + 0, index_c.z} ].[fs] end
      value = rexpr value + e*points_e[ int3d {index_e.x, index_e.y - 1, index_e.z} ].[fs] end
      value = rexpr value + f*points_c[ int3d {index_c.x, index_c.y - 1, index_c.z} ].[fs] end
      value = rexpr value + g*points_e[ int3d {index_e.x, index_e.y - 2, index_e.z} ].[fs] end
      value = rexpr value + h*points_c[ int3d {index_c.x, index_c.y - 2, index_c.z} ].[fs] end
  elseif dir == 2 then  -- z direction stencil
      local index_c = rexpr int3d { x = index.x, y = index.y, z = index.z + ip.n_ghosts - 1} end
      local index_e = rexpr index end

      value = rexpr         a*points_c[ int3d {index_c.x, index_c.y, index_c.z + 2} ].[fs] end
      value = rexpr value + b*points_c[ int3d {index_c.x, index_c.y, index_c.z + 1} ].[fs] end
      value = rexpr value + c*points_e[ int3d {index_e.x, index_e.y, index_e.z + 0} ].[fs] end
      value = rexpr value + d*points_c[ int3d {index_c.x, index_c.y, index_c.z + 0} ].[fs] end
      value = rexpr value + e*points_e[ int3d {index_e.x, index_e.y, index_e.z - 1} ].[fs] end
      value = rexpr value + f*points_c[ int3d {index_c.x, index_c.y, index_c.z - 1} ].[fs] end
      value = rexpr value + g*points_e[ int3d {index_e.x, index_e.y, index_e.z - 2} ].[fs] end
      value = rexpr value + h*points_c[ int3d {index_c.x, index_c.y, index_c.z - 2} ].[fs] end
  end
  return value
end



local function make_stencil( r1, privileges_r1, f1, r2, privileges_r2, f2, nx, ny, nz, onebydx, dir, der )
  local rhs __demand(__inline) task rhs( [r1], [r2] )
  where
    [privileges_r1], [privileges_r2]
  do
    for i in r2 do
      [r2][i].[f2] = [make_stencil_pattern(r1, f1, i, nx, ny, nz, onebydx, dir, der)]
    end
  end
  return rhs
end



local function make_stencil_MND( r1, privileges_r1, f1, r2, privileges_r2, r3, privileges_r3, f3, Nx, Ny, Nz, onebydx, dir, der, periodic )
  local rhs __demand(__inline) task rhs( [r1], [r2], [r3] )
  where
    [privileges_r1], [privileges_r2], [privileges_r3]
  do
    for i in r3.ispace do
      var idir = 0
      var N = 0
      if (dir == 0) then
        idir = i.x
        N = Nx
      elseif (dir == 1) then
        idir = i.y
        N = Ny
      elseif (dir == 2) then
        idir = i.z
        N = Nz
      end

      if ((not periodic) and (idir == ip.n_ghosts) ) then
        [r3][i].[f3] = [make_stencil_pattern_MND_LB(r1, r2, f1, i, onebydx, dir)]
      elseif ((not periodic) and (idir == N+ip.n_ghosts-1)) then
        [r3][i].[f3] = [make_stencil_pattern_MND_RB(r1, r2, f1, i, onebydx, dir)]
      else
        [r3][i].[f3] = [make_stencil_pattern_MND_interior(r1, r2, f1, i, Nx, Ny, Nz, onebydx, dir, periodic)]
      end

    end
  end
  return rhs
end



local function make_stencil_MND_fd( r1, privileges_r1, r2, privileges_r2, r3, privileges_r3, fs, Nx, Ny, Nz, dir, periodic )
  local rhs __demand(__inline) task rhs( [r1], [r2], [r3] )
  where
    [privileges_r1], [privileges_r2], [privileges_r3]
  do
    for i in r3.ispace do
      var idir = 0
      var N = 0
      if (dir == 0) then
        idir = i.x
        N = Nx + 1
      elseif (dir == 1) then
        idir = i.y
        N = Ny + 1
      elseif (dir == 2) then
        idir = i.z
        N = Nz + 1
      end

      if ((not periodic) and (idir == 0) ) then
        [r3][i].[fs] = [make_stencil_pattern_MND_fd_LB(r1, r2, fs, i, dir)]
      elseif ((not periodic) and (idir == N - 1)) then
        [r3][i].[fs] = [make_stencil_pattern_MND_fd_RB(r1, r2, fs, i, dir)]
      else
        [r3][i].[fs] = [make_stencil_pattern_MND_fd_interior(r1, r2, fs, i, Nx, Ny, Nz, dir, periodic)]
      end

    end
  end
  return rhs
end



local function make_SolveXLU( r, privileges_r, f, explicit )
  local SolveXLU __demand(__inline) task SolveXLU( [r],
                                                   LU : region(ispace(int3d), LU_struct),
                                                   n_ghosts : int64 )
  where
    [privileges_r], reads( LU )
  do
    if explicit == false then
      var bounds = [r].ispace.bounds
      var N = bounds.hi.x + 1 - n_ghosts
      var pr = LU.ispace.bounds.hi.y
      var pc = LU.ispace.bounds.hi.z
  
      for k = bounds.lo.z, bounds.hi.z+1 do
        for j = bounds.lo.y, bounds.hi.y+1 do
  
          -- Step 8
          [r][{1+n_ghosts,j,k}].[f] = [r][{1+n_ghosts,j,k}].[f] - LU[{1,pr,pc}].b*[r][{0+n_ghosts,j,k}].[f]
          var sum1 : double = LU[{0,pr,pc}].k*[r][{0+n_ghosts,j,k}].[f] + LU[{1,pr,pc}].k*[r][{1+n_ghosts,j,k}].[f]
          var sum2 : double = LU[{0,pr,pc}].l*[r][{0+n_ghosts,j,k}].[f] + LU[{1,pr,pc}].l*[r][{1+n_ghosts,j,k}].[f]
  
          -- Step 9
          for i = 2,N-2 do
            [r][{i+n_ghosts,j,k}].[f] = [r][{i+n_ghosts,j,k}].[f] - LU[{i,pr,pc}].b*[r][{i-1+n_ghosts,j,k}].[f] - LU[{i,pr,pc}].eg*[r][{i-2+n_ghosts,j,k}].[f]
            sum1 += LU[{i,pr,pc}].k*[r][{i+n_ghosts,j,k}].[f]
            sum2 += LU[{i,pr,pc}].l*[r][{i+n_ghosts,j,k}].[f]
          end
  
          -- Step 10
          [r][{N-2+n_ghosts,j,k}].[f] = [r][{N-2+n_ghosts,j,k}].[f] - sum1;
          [r][{N-1+n_ghosts,j,k}].[f] = ( [r][{N-1+n_ghosts,j,k}].[f] - sum2 - LU[{N-2,pr,pc}].l*[r][{N-2+n_ghosts,j,k}].[f] )*LU[{N-1,pr,pc}].g;
  
          -- Step 11
          [r][{N-2+n_ghosts,j,k}].[f] = ( [r][{N-2+n_ghosts,j,k}].[f] - LU[{N-2,pr,pc}].w*[r][{N-1+n_ghosts,j,k}].[f] )*LU[{N-2,pr,pc}].g;
          [r][{N-3+n_ghosts,j,k}].[f] = ( [r][{N-3+n_ghosts,j,k}].[f] - LU[{N-3,pr,pc}].v*[r][{N-2+n_ghosts,j,k}].[f] - LU[{N-3,pr,pc}].w*[r][{N-1+n_ghosts,j,k}].[f] )*LU[{N-3,pr,pc}].g;
          [r][{N-4+n_ghosts,j,k}].[f] = ( [r][{N-4+n_ghosts,j,k}].[f] - LU[{N-4,pr,pc}].h*[r][{N-3+n_ghosts,j,k}].[f] - LU[{N-4,pr,pc}].v*[r][{N-2+n_ghosts,j,k}].[f] - LU[{N-4,pr,pc}].w*[r][{N-1+n_ghosts,j,k}].[f] )*LU[{N-4,pr,pc}].g;
          for i = N-5,-1,-1 do
            [r][{i+n_ghosts,j,k}].[f] = ( [r][{i+n_ghosts,j,k}].[f] - LU[{i,pr,pc}].h*[r][{i+1+n_ghosts,j,k}].[f] - LU[{i,pr,pc}].ff*[r][{i+2+n_ghosts,j,k}].[f] - LU[{i,pr,pc}].v*[r][{N-2+n_ghosts,j,k}].[f] - LU[{i,pr,pc}].w*[r][{N-1+n_ghosts,j,k}].[f] )*LU[{i,pr,pc}].g;
          end
  
        end
      end
    end
    return 1
  end
  return SolveXLU
end



local function make_SolveYLU( r, privileges_r, f, explicit )
  local SolveYLU __demand(__inline) task SolveYLU( [r],
                                                   LU : region(ispace(int3d), LU_struct),
                                                   n_ghosts : int64 )
  where
    [privileges_r], reads( LU ) 
  do
    if explicit == false then
      var bounds = [r].ispace.bounds
      var N = bounds.hi.y + 1 - n_ghosts
      var pr = LU.ispace.bounds.hi.y
      var pc = LU.ispace.bounds.hi.z
  
      for k = bounds.lo.z, bounds.hi.z+1 do
        for i = bounds.lo.x, bounds.hi.x+1 do
  
          -- Step 8
          [r][{i,1+n_ghosts,k}].[f] = [r][{i,1+n_ghosts,k}].[f] - LU[{1,pr,pc}].b*[r][{i,0+n_ghosts,k}].[f]
          var sum1 : double = LU[{0,pr,pc}].k*[r][{i,0+n_ghosts,k}].[f] + LU[{1,pr,pc}].k*[r][{i,1+n_ghosts,k}].[f]
          var sum2 : double = LU[{0,pr,pc}].l*[r][{i,0+n_ghosts,k}].[f] + LU[{1,pr,pc}].l*[r][{i,1+n_ghosts,k}].[f]
  
          -- Step 9
          for j = 2,N-2 do
            [r][{i,j+n_ghosts,k}].[f] = [r][{i,j+n_ghosts,k}].[f] - LU[{j,pr,pc}].b*[r][{i,j-1+n_ghosts,k}].[f] - LU[{j,pr,pc}].eg*[r][{i,j-2+n_ghosts,k}].[f]
            sum1 += LU[{j,pr,pc}].k*[r][{i,j+n_ghosts,k}].[f]
            sum2 += LU[{j,pr,pc}].l*[r][{i,j+n_ghosts,k}].[f]
          end
  
          -- Step 10
          [r][{i,N-2+n_ghosts,k}].[f] = [r][{i,N-2+n_ghosts,k}].[f] - sum1;
          [r][{i,N-1+n_ghosts,k}].[f] = ( [r][{i,N-1+n_ghosts,k}].[f] - sum2 - LU[{N-2,pr,pc}].l*[r][{i,N-2+n_ghosts,k}].[f] )*LU[{N-1,pr,pc}].g;
  
          -- Step 11
          [r][{i,N-2+n_ghosts,k}].[f] = ( [r][{i,N-2+n_ghosts,k}].[f] - LU[{N-2,pr,pc}].w*[r][{i,N-1+n_ghosts,k}].[f] )*LU[{N-2,pr,pc}].g;
          [r][{i,N-3+n_ghosts,k}].[f] = ( [r][{i,N-3+n_ghosts,k}].[f] - LU[{N-3,pr,pc}].v*[r][{i,N-2+n_ghosts,k}].[f] - LU[{N-3,pr,pc}].w*[r][{i,N-1+n_ghosts,k}].[f] )*LU[{N-3,pr,pc}].g;
          [r][{i,N-4+n_ghosts,k}].[f] = ( [r][{i,N-4+n_ghosts,k}].[f] - LU[{N-4,pr,pc}].h*[r][{i,N-3+n_ghosts,k}].[f] - LU[{N-4,pr,pc}].v*[r][{i,N-2+n_ghosts,k}].[f] - LU[{N-4,pr,pc}].w*[r][{i,N-1+n_ghosts,k}].[f] )*LU[{N-4,pr,pc}].g
          for j = N-5,-1,-1 do
            [r][{i,j+n_ghosts,k}].[f] = ( [r][{i,j+n_ghosts,k}].[f] - LU[{j,pr,pc}].h*[r][{i,j+1+n_ghosts,k}].[f] - LU[{j,pr,pc}].ff*[r][{i,j+2+n_ghosts,k}].[f] - LU[{j,pr,pc}].v*[r][{i,N-2+n_ghosts,k}].[f] - LU[{j,pr,pc}].w*[r][{i,N-1+n_ghosts,k}].[f] )*LU[{j,pr,pc}].g
          end
  
        end
      end
    end
    return 1
  end
  return SolveYLU
end



local function make_SolveZLU( r, privileges_r, f, explicit )
  local SolveZLU __demand(__inline) task SolveZLU( [r],
                                                   LU : region(ispace(int3d), LU_struct),
                                                   n_ghosts : int64 )
  where
    [privileges_r], reads( LU )
  do
    if explicit == false then
      var bounds = [r].ispace.bounds
      var N = bounds.hi.z + 1 - n_ghosts
      var pr = LU.ispace.bounds.hi.y
      var pc = LU.ispace.bounds.hi.z
  
      for j = bounds.lo.y, bounds.hi.y+1 do
        for i = bounds.lo.x, bounds.hi.x+1 do
  
          -- Step 8
          [r][{i,j,n_ghosts+1}].[f] = [r][{i,j,n_ghosts+1}].[f] - LU[{1,pr,pc}].b*[r][{i,j,n_ghosts+0}].[f]
          var sum1 : double = LU[{0,pr,pc}].k*[r][{i,j,n_ghosts+0}].[f] + LU[{1,pr,pc}].k*[r][{i,j,n_ghosts+1}].[f]
          var sum2 : double = LU[{0,pr,pc}].l*[r][{i,j,n_ghosts+0}].[f] + LU[{1,pr,pc}].l*[r][{i,j,n_ghosts+1}].[f]
  
          -- Step 9
          for k = 2,N-2 do
            [r][{i,j,n_ghosts+k}].[f] = [r][{i,j,n_ghosts+k}].[f] - LU[{k,pr,pc}].b*[r][{i,j,n_ghosts+k-1}].[f] - LU[{k,pr,pc}].eg*[r][{i,j,n_ghosts+k-2}].[f]
            sum1 += LU[{k,pr,pc}].k*[r][{i,j,n_ghosts+k}].[f]
            sum2 += LU[{k,pr,pc}].l*[r][{i,j,n_ghosts+k}].[f]
          end
  
          -- Step 10
          [r][{i,j,n_ghosts+N-2}].[f] = [r][{i,j,n_ghosts+N-2}].[f] - sum1;
          [r][{i,j,n_ghosts+N-1}].[f] = ( [r][{i,j,n_ghosts+N-1}].[f] - sum2 - LU[{N-2,pr,pc}].l*[r][{i,j,n_ghosts+N-2}].[f] )*LU[{N-1,pr,pc}].g;
  
          -- Step 11
          [r][{i,j,n_ghosts+N-2}].[f] = ( [r][{i,j,n_ghosts+N-2}].[f] - LU[{N-2,pr,pc}].w*[r][{i,j,n_ghosts+N-1}].[f] )*LU[{N-2,pr,pc}].g;
          [r][{i,j,n_ghosts+N-3}].[f] = ( [r][{i,j,n_ghosts+N-3}].[f] - LU[{N-3,pr,pc}].v*[r][{i,j,n_ghosts+N-2}].[f] - LU[{N-3,pr,pc}].w*[r][{i,j,n_ghosts+N-1}].[f] )*LU[{N-3,pr,pc}].g;
          [r][{i,j,n_ghosts+N-4}].[f] = ( [r][{i,j,n_ghosts+N-4}].[f] - LU[{N-4,pr,pc}].h*[r][{i,j,n_ghosts+N-3}].[f] - LU[{N-4,pr,pc}].v*[r][{i,j,n_ghosts+N-2}].[f] - LU[{N-4,pr,pc}].w*[r][{i,j,n_ghosts+N-1}].[f] )*LU[{N-4,pr,pc}].g
          for k = N-5,-1,-1 do
            [r][{i,j,n_ghosts+k}].[f] = ( [r][{i,j,n_ghosts+k}].[f] - LU[{k,pr,pc}].h*[r][{i,j,n_ghosts+k+1}].[f] - LU[{k,pr,pc}].ff*[r][{i,j,n_ghosts+k+2}].[f] - LU[{k,pr,pc}].v*[r][{i,j,n_ghosts+N-2}].[f] - LU[{k,pr,pc}].w*[r][{i,j,n_ghosts+N-1}].[f] )*LU[{k,pr,pc}].g
          end
  
        end
      end
    end
    return 1
  end
  return SolveZLU
end



function make_ddx( r_func, f_func, r_der, f_der, NX, NY, NZ, ONEBYDX )
  local privileges_r_func = regentlib.privilege(regentlib.reads,  r_func, f_func)
  
  local reads_r_der      = regentlib.privilege(regentlib.reads,  r_der, f_der)
  local writes_r_der     = regentlib.privilege(regentlib.writes, r_der, f_der)
  local privileges_r_der = terralib.newlist({reads_r_der, writes_r_der})

  local ComputeXRHS = make_stencil(r_func, privileges_r_func, f_func, r_der, privileges_r_der, f_der, NX, NY, NZ, ONEBYDX, 0, 1)
  local SolveXLU    = make_SolveXLU(r_der, privileges_r_der, f_der, false)

  local ddx __demand(__inline) task ddx( [r_func],
                                         [r_der],
                                         LU     : region(ispace(int3d), LU_struct) )
  where                            
    reads( LU ), [privileges_r_func], [privileges_r_der]
  do
    [ComputeXRHS]([r_func], [r_der])
    var token = [SolveXLU]([r_der], LU, ip.n_ghosts)
    token = r_der[ [r_der].ispace.bounds.lo ].[f_der]
    return token
  end
  return ddx
end



function make_d2dx2( r_func, f_func, r_der, f_der, NX, NY, NZ, ONEBYDX )
  local privileges_r_func = regentlib.privilege(regentlib.reads,  r_func, f_func)
  
  local reads_r_der      = regentlib.privilege(regentlib.reads,  r_der, f_der)
  local writes_r_der     = regentlib.privilege(regentlib.writes, r_der, f_der)
  local privileges_r_der = terralib.newlist({reads_r_der, writes_r_der})

  local ComputeXRHS = make_stencil(r_func, privileges_r_func, f_func, r_der, privileges_r_der, f_der, NX, NY, NZ, ONEBYDX, 0, 2)
  local SolveXLU    = make_SolveXLU(r_der, privileges_r_der, f_der, false)

  local d2dx2 __demand(__inline) task d2dx2( [r_func],
                                             [r_der],
                                             LU     : region(ispace(int3d), LU_struct) )
  where                            
    reads( LU ), [privileges_r_func], [privileges_r_der]
  do
    [ComputeXRHS]([r_func], [r_der])
    var token = [SolveXLU]([r_der], LU, ip.n_ghosts)
    token = r_der[ [r_der].ispace.bounds.lo ].[f_der]
    return token
  end
  return d2dx2
end



function make_ddy( r_func, f_func, r_der, f_der, NX, NY, NZ, ONEBYDX )
  local privileges_r_func = regentlib.privilege(regentlib.reads,  r_func, f_func)
  
  local reads_r_der      = regentlib.privilege(regentlib.reads,  r_der, f_der)
  local writes_r_der     = regentlib.privilege(regentlib.writes, r_der, f_der)
  local privileges_r_der = terralib.newlist({reads_r_der, writes_r_der})

  local ComputeYRHS = make_stencil(r_func, privileges_r_func, f_func, r_der, privileges_r_der, f_der, NX, NY, NZ, ONEBYDX, 1, 1)
  local SolveYLU    = make_SolveYLU(r_der, privileges_r_der, f_der, false)

  local ddy __demand(__inline) task ddy( [r_func],
                                         [r_der],
                                         LU     : region(ispace(int3d), LU_struct) )
  where                            
    reads( LU ), [privileges_r_func], [privileges_r_der]
  do
    [ComputeYRHS]([r_func], [r_der])
    var token = [SolveYLU]([r_der], LU, ip.n_ghosts)
    token = r_der[ [r_der].ispace.bounds.lo ].[f_der]
    return token
  end
  return ddy
end



function make_d2dy2( r_func, f_func, r_der, f_der, NX, NY, NZ, ONEBYDX )
  local privileges_r_func = regentlib.privilege(regentlib.reads,  r_func, f_func)
  
  local reads_r_der      = regentlib.privilege(regentlib.reads,  r_der, f_der)
  local writes_r_der     = regentlib.privilege(regentlib.writes, r_der, f_der)
  local privileges_r_der = terralib.newlist({reads_r_der, writes_r_der})

  local ComputeYRHS = make_stencil(r_func, privileges_r_func, f_func, r_der, privileges_r_der, f_der, NX, NY, NZ, ONEBYDX, 1, 2)
  local SolveYLU    = make_SolveYLU(r_der, privileges_r_der, f_der, false)

  local d2dy2 __demand(__inline) task d2dy2( [r_func],
                                             [r_der],
                                             LU     : region(ispace(int3d), LU_struct) )
  where                            
    reads( LU ), [privileges_r_func], [privileges_r_der]
  do
    [ComputeYRHS]([r_func], [r_der])
    var token = [SolveYLU]([r_der], LU, ip.n_ghosts)
    token = r_der[ [r_der].ispace.bounds.lo ].[f_der]
    return token
  end
  return d2dy2
end



function make_ddz( r_func, f_func, r_der, f_der, NX, NY, NZ, ONEBYDX )
  local privileges_r_func = regentlib.privilege(regentlib.reads,  r_func, f_func)
  
  local reads_r_der      = regentlib.privilege(regentlib.reads,  r_der, f_der)
  local writes_r_der     = regentlib.privilege(regentlib.writes, r_der, f_der)
  local privileges_r_der = terralib.newlist({reads_r_der, writes_r_der})

  local ComputeZRHS = make_stencil(r_func, privileges_r_func, f_func, r_der, privileges_r_der, f_der, NX, NY, NZ, ONEBYDX, 2, 1)
  local SolveZLU    = make_SolveZLU(r_der, privileges_r_der, f_der, false)

  local ddz __demand(__inline) task ddz( [r_func],
                                         [r_der],
                                         LU     : region(ispace(int3d), LU_struct) )
  where                            
    reads( LU ), [privileges_r_func], [privileges_r_der]
  do
    [ComputeZRHS]([r_func], [r_der])
    var token = [SolveZLU]([r_der], LU, ip.n_ghosts)
    token = r_der[ [r_der].ispace.bounds.lo ].[f_der]
    return token
  end
  return ddz
end



function make_d2dz2( r_func, f_func, r_der, f_der, NX, NY, NZ, ONEBYDX )
  local privileges_r_func = regentlib.privilege(regentlib.reads,  r_func, f_func)
  
  local reads_r_der      = regentlib.privilege(regentlib.reads,  r_der, f_der)
  local writes_r_der     = regentlib.privilege(regentlib.writes, r_der, f_der)
  local privileges_r_der = terralib.newlist({reads_r_der, writes_r_der})

  local ComputeZRHS = make_stencil(r_func, privileges_r_func, f_func, r_der, privileges_r_der, f_der, NX, NY, NZ, ONEBYDX, 2, 2)
  local SolveZLU    = make_SolveZLU(r_der, privileges_r_der, f_der, false)

  local d2dz2 __demand(__inline) task d2dz2( [r_func],
                                             [r_der],
                                             LU     : region(ispace(int3d), LU_struct) )
  where                            
    reads( LU ), [privileges_r_func], [privileges_r_der]
  do
    [ComputeZRHS]([r_func], [r_der])
    var token = [SolveZLU]([r_der], LU, ip.n_ghosts)
    token = r_der[ [r_der].ispace.bounds.lo ].[f_der]
    return token
  end
  return d2dz2
end



function make_ddx_MND( r_func, r_func_e, f_func, r_der, f_der, NX, NY, NZ, ONEBYDX, periodic) 
  local privileges_r_func   = regentlib.privilege(regentlib.reads,  r_func,   f_func)
  local privileges_r_func_e = regentlib.privilege(regentlib.reads,  r_func_e, f_func)
  
  local reads_r_der      = regentlib.privilege(regentlib.reads,  r_der, f_der)
  local writes_r_der     = regentlib.privilege(regentlib.writes, r_der, f_der)
  local privileges_r_der = terralib.newlist({reads_r_der, writes_r_der})

  local ComputeXRHS_MND = make_stencil_MND(r_func, privileges_r_func, f_func, r_func_e, privileges_r_func_e, r_der, privileges_r_der, f_der, NX, NY, NZ, ONEBYDX, 0, 1, periodic)

  local explicit = false
  if alpha06MND == 0.0 and beta06MND == 0.0 and alpha06MND_LB == 0.0 and beta06MND_LB == 0.0 then
    explicit = true
    print("Using explicit MND scheme in X")
  end

  local SolveXLU = make_SolveXLU(r_der, privileges_r_der, f_der, explicit)

  local ddx_MND __demand(__inline) task ddx_MND( [r_func],
                                                 [r_func_e],
                                                 [r_der],
                                                 LU     : region(ispace(int3d), LU_struct) )
  where
    reads( LU ), [privileges_r_func], [privileges_r_func_e], [privileges_r_der]
  do
    [ComputeXRHS_MND]([r_func], [r_func_e], [r_der])
    var token = [SolveXLU]([r_der], LU, ip.n_ghosts)
    token = r_der[ [r_der].ispace.bounds.lo ].[f_der]
    return token
  end
  return ddx_MND
end



function make_ddy_MND( r_func, r_func_e, f_func, r_der, f_der, NX, NY, NZ, ONEBYDX, periodic )
  local privileges_r_func   = regentlib.privilege(regentlib.reads,  r_func,   f_func)
  local privileges_r_func_e = regentlib.privilege(regentlib.reads,  r_func_e, f_func)
  
  local reads_r_der      = regentlib.privilege(regentlib.reads,  r_der, f_der)
  local writes_r_der     = regentlib.privilege(regentlib.writes, r_der, f_der)
  local privileges_r_der = terralib.newlist({reads_r_der, writes_r_der})

  local ComputeYRHS_MND = make_stencil_MND(r_func, privileges_r_func, f_func, r_func_e, privileges_r_func_e, r_der, privileges_r_der, f_der, NX, NY, NZ, ONEBYDX, 1, 1, periodic)

  local explicit = false
  if alpha06MND == 0.0 and beta06MND == 0.0 and alpha06MND_LB == 0.0 and beta06MND_LB == 0.0 then
    explicit = true
    print("Using explicit MND scheme in Y")
  end

  local SolveYLU = make_SolveYLU(r_der, privileges_r_der, f_der, explicit)

  local ddy_MND __demand(__inline) task ddy_MND( [r_func],
                                                 [r_func_e],
                                                 [r_der],
                                                 LU     : region(ispace(int3d), LU_struct) )
  where
    reads( LU ), [privileges_r_func], [privileges_r_func_e], [privileges_r_der]
  do
    [ComputeYRHS_MND]([r_func], [r_func_e], [r_der])
    var token = [SolveYLU]([r_der], LU, ip.n_ghosts)
    token = r_der[ [r_der].ispace.bounds.lo ].[f_der]
    return token
  end
  return ddy_MND
end



function make_ddz_MND( r_func, r_func_e, f_func, r_der, f_der, NX, NY, NZ, ONEBYDX, periodic )
  local privileges_r_func   = regentlib.privilege(regentlib.reads,  r_func,   f_func)
  local privileges_r_func_e = regentlib.privilege(regentlib.reads,  r_func_e, f_func)
  
  local reads_r_der      = regentlib.privilege(regentlib.reads,  r_der, f_der)
  local writes_r_der     = regentlib.privilege(regentlib.writes, r_der, f_der)
  local privileges_r_der = terralib.newlist({reads_r_der, writes_r_der})

  local ComputeZRHS_MND = make_stencil_MND(r_func, privileges_r_func, f_func, r_func_e, privileges_r_func_e, r_der, privileges_r_der, f_der, NX, NY, NZ, ONEBYDX, 2, 1, periodic)

  local explicit = false
  if alpha06MND == 0.0 and beta06MND == 0.0 and alpha06MND_LB == 0.0 and beta06MND_LB == 0.0 then
    explicit = true
    print("Using explicit MND scheme in Z")
  end

  local SolveZLU = make_SolveZLU(r_der, privileges_r_der, f_der, explicit)

  local ddz_MND __demand(__inline) task ddz_MND( [r_func],
                                                 [r_func_e],
                                                 [r_der],
                                                 LU     : region(ispace(int3d), LU_struct) )
  where
    reads( LU ), [privileges_r_func], [privileges_r_func_e], [privileges_r_der]
  do
    [ComputeZRHS_MND]([r_func], [r_func_e], [r_der])
    var token = [SolveZLU]([r_der], LU, ip.n_ghosts)
    token = r_der[ [r_der].ispace.bounds.lo ].[f_der]
    return token
  end
  return ddz_MND
end



function make_ddx_MND_fd( r_func, r_func_e, r_func_ee, f_func, NX, NY, NZ, periodic )
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func,   f_func)
  local privileges_r_func_e = regentlib.privilege(regentlib.reads, r_func_e, f_func)
  
  local reads_r_func_ee      = regentlib.privilege(regentlib.reads,  r_func_ee, f_func)
  local writes_r_func_ee     = regentlib.privilege(regentlib.writes, r_func_ee, f_func)
  local privileges_r_func_ee = terralib.newlist({reads_r_func_ee, writes_r_func_ee})

  local ComputeXRHS_MND_fd = make_stencil_MND_fd(r_func, privileges_r_func, r_func_e, privileges_r_func_e, r_func_ee, privileges_r_func_ee, f_func, NX, NY, NZ, 0, periodic)

  local explicit = false
  if alpha06MND_fd == 0.0 and beta06MND_fd == 0.0 and alpha06MND_LB_fd == 0.0 and beta06MND_LB_fd == 0.0 then
    explicit = true
    print("Using explicit MND scheme (flux difference form) in X")
  end

  local SolveXLU = make_SolveXLU(r_func_ee, privileges_r_func_ee, f_func, explicit)

  local ddx_MND_fd __demand(__inline) task ddx_MND_fd( [r_func],
                                                       [r_func_e],
                                                       [r_func_ee],
                                                       LU : region(ispace(int3d), LU_struct) )
  where
    reads( LU ), [privileges_r_func], [privileges_r_func_e], [privileges_r_func_ee]
  do
    [ComputeXRHS_MND_fd]([r_func], [r_func_e], [r_func_ee])
    var token = [SolveXLU]([r_func_ee], LU, 0)
    token = r_func_ee[ [r_func_ee].ispace.bounds.lo ].[f_func]
    return token
  end
  return ddx_MND_fd
end



function make_ddy_MND_fd( r_func, r_func_e, r_func_ee, f_func, NX, NY, NZ, periodic )
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func,   f_func)
  local privileges_r_func_e = regentlib.privilege(regentlib.reads, r_func_e, f_func)

  local reads_r_func_ee      = regentlib.privilege(regentlib.reads,  r_func_ee, f_func)
  local writes_r_func_ee     = regentlib.privilege(regentlib.writes, r_func_ee, f_func)
  local privileges_r_func_ee = terralib.newlist({reads_r_func_ee, writes_r_func_ee})

  local ComputeYRHS_MND_fd = make_stencil_MND_fd(r_func, privileges_r_func, r_func_e, privileges_r_func_e, r_func_ee, privileges_r_func_ee, f_func, NX, NY, NZ, 1, periodic)

  local explicit = false
  if alpha06MND_fd == 0.0 and beta06MND_fd == 0.0 and alpha06MND_LB_fd == 0.0 and beta06MND_LB_fd == 0.0 then
    explicit = true
    print("Using explicit MND scheme (flux difference form) in Y")
  end

  local SolveYLU = make_SolveYLU(r_func_ee, privileges_r_func_ee, f_func, explicit)

  local ddy_MND_fd __demand(__inline) task ddy_MND_fd( [r_func],
                                                       [r_func_e],
                                                       [r_func_ee],
                                                       LU : region(ispace(int3d), LU_struct) )
  where
    reads( LU ), [privileges_r_func], [privileges_r_func_e], [privileges_r_func_ee]
  do
    [ComputeYRHS_MND_fd]([r_func], [r_func_e], [r_func_ee])
    var token = [SolveYLU]([r_func_ee], LU, 0)
    token = r_func_ee[ [r_func_ee].ispace.bounds.lo ].[f_func]
    return token
  end
  return ddy_MND_fd
end



function make_ddz_MND_fd( r_func, r_func_e, r_func_ee, f_func, NX, NY, NZ, periodic )
  local privileges_r_func   = regentlib.privilege(regentlib.reads, r_func,   f_func)
  local privileges_r_func_e = regentlib.privilege(regentlib.reads, r_func_e, f_func)

  local reads_r_func_ee      = regentlib.privilege(regentlib.reads,  r_func_ee, f_func)
  local writes_r_func_ee     = regentlib.privilege(regentlib.writes, r_func_ee, f_func)
  local privileges_r_func_ee = terralib.newlist({reads_r_func_ee, writes_r_func_ee})

  local ComputeZRHS_MND_fd = make_stencil_MND_fd(r_func, privileges_r_func, r_func_e, privileges_r_func_e, r_func_ee, privileges_r_func_ee, f_func, NX, NY, NZ, 2, periodic)

  local explicit = false
  if alpha06MND_fd == 0.0 and beta06MND_fd == 0.0 and alpha06MND_LB_fd == 0.0 and beta06MND_LB_fd == 0.0 then
    explicit = true
    print("Using explicit MND scheme (flux difference form) in Z")
  end
  local SolveZLU = make_SolveZLU(r_func_ee, privileges_r_func_ee, f_func, explicit)

  local ddz_MND_fd __demand(__inline) task ddz_MND_fd( [r_func],
                                                       [r_func_e],
                                                       [r_func_ee],
                                                       LU : region(ispace(int3d), LU_struct) )
  where
    reads( LU ), [privileges_r_func], [privileges_r_func_e], [privileges_r_func_ee]
  do
    [ComputeZRHS_MND_fd]([r_func], [r_func_e], [r_func_ee])
    var token = [SolveZLU]([r_func_ee], LU, 0)
    token = r_func_ee[ [r_func_ee].ispace.bounds.lo ].[f_func]
    return token
  end
  return ddz_MND_fd
end
