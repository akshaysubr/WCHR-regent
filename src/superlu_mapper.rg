-- Copyright 2017 Stanford University
--
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
--
--     http://www.apache.org/licenses/LICENSE-2.0
--
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

import "regent"

local csuperlu_mapper
do
  assert(os.getenv('LG_RT_DIR') ~= nil, "$LG_RT_DIR should be set!")
  -- local root_dir = arg[0]:match(".*/") or "./"
  -- local root_dir = "/home/akshays/Codes/WCHR-regent/src/"
  local root_dir = "/home/manlong/Codes/WCHR-regent/src/"
  local runtime_dir = os.getenv('LG_RT_DIR') .. "/"
  local runtime_dir = os.getenv('LG_RT_DIR') .. "/"
  local legion_dir = runtime_dir .. "legion/"
  local mapper_dir = runtime_dir .. "mappers/"
  local realm_dir = runtime_dir .. "realm/"
  local superlu_mapper_cc = root_dir .. "superlu_mapper.cc"
  local superlu_mapper_so
  if os.getenv('SAVEOBJ') == '1' then
    superlu_mapper_so = root_dir .. "libsuperlu_mapper.so"
  else
    superlu_mapper_so = os.tmpname() .. ".so" -- root_dir .. "mapper.so"
  end
  local cxx = os.getenv('CXX') or 'c++'

  local cxx_flags = "-O2 -Wall -Werror"
  if os.execute('test "$(uname)" = Darwin') == 0 then
    cxx_flags =
      (cxx_flags ..
         " -dynamiclib -single_module -undefined dynamic_lookup -fPIC")
  else
    cxx_flags = cxx_flags .. " -shared -fPIC"
  end

  local cmd = (cxx .. " " .. cxx_flags .. " -I " .. runtime_dir .. " " ..
                 " -I " .. mapper_dir .. " " .. " -I " .. legion_dir .. " " ..
                 " -I " .. realm_dir .. " " .. superlu_mapper_cc .. " -o " .. superlu_mapper_so)
  if os.execute(cmd) ~= 0 then
    print("Error: failed to compile " .. superlu_mapper_cc)
    assert(false)
  end
  terralib.linklibrary(superlu_mapper_so)
  csuperlu_mapper =
    terralib.includec("superlu_mapper.h", {"-I", root_dir, "-I", runtime_dir,
                                          "-I", mapper_dir, "-I", legion_dir,
                                          "-I", realm_dir})
end

return csuperlu_mapper

-- local c = regentlib.c
-- 
-- require("fields")
-- local superlu = require("superlu_util")
-- 
-- local terra get_base_pointer_3d(pr   : c.legion_physical_region_t[1],
--                                 fid  : c.legion_field_id_t[1],
--                                 rect : c.legion_rect_3d_t)
--   var subrect : c.legion_rect_3d_t
--   var offsets : c.legion_byte_offset_t[3]
--   var accessor = c.legion_physical_region_get_field_accessor_generic(pr[0], fid[0])
--   var base_pointer =
--     [&double](c.legion_accessor_generic_raw_rect_ptr_3d(
--       accessor, rect, &subrect, &(offsets[0])))
--   c.legion_accessor_generic_destroy(accessor)
--   return base_pointer
-- end
-- 
-- task initialize( r_prim_l_x : region(ispace(int3d), primitive),
--                  r_rhs_l_x  : region(ispace(int3d), primitive) )
-- where
--   reads writes(r_prim_l_x, r_rhs_l_x)
-- do
--   var bnds = r_prim_l_x.ispace.bounds
--   var nx = bnds.hi.x + 1
--   var ny = bnds.hi.y + 1
--   var nz = bnds.hi.z + 1
--   for i in r_prim_l_x do
--     r_prim_l_x[i].rho = 5*(i.x + nx*i.y + nx*ny*i.z)
--     r_prim_l_x[i].u   = 5*(i.x + nx*i.y + nx*ny*i.z) + 1
--     r_prim_l_x[i].v   = 5*(i.x + nx*i.y + nx*ny*i.z) + 2
--     r_prim_l_x[i].w   = 5*(i.x + nx*i.y + nx*ny*i.z) + 3
--     r_prim_l_x[i].p   = 5*(i.x + nx*i.y + nx*ny*i.z) + 4
-- 
--     r_rhs_l_x[i].rho = 5*(i.x + nx*i.y + nx*ny*i.z)
--     r_rhs_l_x[i].u   = 5*(i.x + nx*i.y + nx*ny*i.z) + 1
--     r_rhs_l_x[i].v   = 5*(i.x + nx*i.y + nx*ny*i.z) + 2
--     r_rhs_l_x[i].w   = 5*(i.x + nx*i.y + nx*ny*i.z) + 3
--     r_rhs_l_x[i].p   = 5*(i.x + nx*i.y + nx*ny*i.z) + 4
--   end
-- 
--   var p1 = get_base_pointer_3d(__physical(r_prim_l_x.rho), __fields(r_prim_l_x.rho), r_prim_l_x.bounds)
--   c.printf("=== initialize ===\n")
--   c.printf("%p\n", p1)
--   for i = 0, 5*nx*ny*nz do
--     c.printf("%.0f ", p1[i])
--   end
--   c.printf("\n==================\n")
-- end
-- 
-- task toplevel()
--   var Nx : int64 = 6
--   var Ny : int64 = 6
--   var Nz : int64 = 6
-- 
--   var Lx : double = 1
--   var Ly : double = 1
--   var Lz : double = 1
-- 
--   var dx : double = Lx / Nx
--   var dy : double = Ly / Ny
--   var dz : double = Lz / Nz
-- 
--   c.printf("================ Problem parameters ================\n")
--   c.printf("           Nx, Ny, Nz = %d, %d, %d\n", Nx, Ny, Nz)
--   c.printf("           Lx, Ly, Lz = %f, %f, %f\n", Lx, Ly, Lz)
--   c.printf("           dx, dy, dz = %f, %f, %f\n", dx, dy, dz)
--   c.printf("====================================================\n")
-- 
--   --------------------------------------------------------------------------------------------
--   --                       DATA STUCTURES
--   --------------------------------------------------------------------------------------------
--   var grid_e_x   = ispace(int3d, {x = Nx+1, y = Ny,   z = Nz  })  -- x cell edge index space
-- 
--   var r_prim_l_x = region(grid_e_x, primitive)  -- Primitive variables at left x cell edge
--   var r_rhs_l_x  = region(grid_e_x, primitive)  -- Store RHS for left interpolation in x
-- 
--   var pgrid_x    = ispace(int2d, {x = 1, y = 1}) -- Processor grid in x
--   var slu_x      = region(pgrid_x, superlu.c.superlu_vars_t) -- Super LU data structure for x interpolation
--   var matrix_l_x = region(pgrid_x, superlu.CSR_matrix) -- matrix data structure for x left interpolation
--   --------------------------------------------------------------------------------------------
--   --------------------------------------------------------------------------------------------
-- 
--   -- Initialize SuperLU stuff
--   matrix_l_x[{0,0}] = superlu.initialize_matrix_x(3./16., 5./8., 3./16., Nx, Ny, Nz)
--   
--   initialize(r_prim_l_x, r_rhs_l_x)
--   
--   superlu.initialize_superlu_vars( matrix_l_x[{0,0}], 5*(Nx+1)*Ny*Nz, r_rhs_l_x, r_prim_l_x, slu_x ) 
--   superlu.MatrixSolve( r_rhs_l_x, r_prim_l_x, matrix_l_x[{0,0}], Nx, Ny, Nz, slu_x ) 
-- end
-- 
-- if os.getenv('SAVEOBJ') == '1' then
--   local root_dir = arg[0]:match(".*/") or "./"
--   local link_flags = {"-L" .. root_dir, "-lsuperlu_mapper"}
--   regentlib.saveobj(toplevel, "superlu_mapper", "executable", csuperlu_mapper.register_mappers, link_flags)
-- else
--   regentlib.start(toplevel, csuperlu_mapper.register_mappers)
-- end
