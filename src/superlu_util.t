import "regent"

local superlu = {}
do
  local superlu_library = "-lsuperlu"
  local superlu_include_dir = "/opt/SuperLU_5.2.1"
  -- local root_dir = arg[0]:match(".*/") or "./"
  local root_dir = "/home/akshays/Codes/WCHR-regent/src/"
  local superlu_util_cc = root_dir .. "superlu_util.c"
  superlu_util_so = os.tmpname() .. ".so"
  local cc = os.getenv('CC') or 'cc'
  local cc_flags = "-g -O0 -Wall -Werror -std=c99"
  cc_flags = cc_flags .. " -I" .. superlu_include_dir
  local is_darwin = os.execute('test "$(uname)" = Darwin') == 0
  if is_darwin then
    cc_flags =
      (cc_flags ..
         " -dynamiclib -single_module -undefined dynamic_lookup -fPIC")
  else
    cc_flags = cc_flags .. " -shared -fPIC"
  end
  cc_flags = cc_flags .. " -lm -lblas " .. superlu_library 

  local cmd = (cc .. " " .. cc_flags .. " " .. superlu_util_cc .. " -o " .. superlu_util_so)
  -- print(cmd)

  if os.execute(cmd) ~= 0 then
    print("Error: failed to compile " .. superlu_util_cc)
    assert(false)
  end
  terralib.linklibrary(superlu_util_so)
  if is_darwin then
    terralib.linklibrary("libsuperlu.dylib")
  else
    terralib.linklibrary("libsuperlu.so")
  end
  superlu.c = terralib.includec("superlu_util.h", {"-I", root_dir, "-I", superlu_include_dir })
end

local c = regentlib.c

struct superlu.CSR_matrix {
  nzval  : &double,
  colind : &int,
  rowptr : &int,
  nnz    : int64,
}

terra superlu.initialize_matrix_x( alpha  : double,
                                   beta   : double,
                                   gamma  : double,
                                   nx     : int64,
                                   ny     : int64,
                                   nz     : int64 )
  var matrix : superlu.CSR_matrix

  var Nsize : int64 = (nx+1)*ny*nz
  matrix.nnz = (3*nx+2)*ny*nz
  matrix.rowptr = [&int] ( c.malloc ( (Nsize+1) * sizeof(int) ) )
  matrix.colind = [&int] ( c.malloc ( matrix.nnz * sizeof(int) ) )
  matrix.nzval  = [&double] ( c.malloc ( matrix.nnz * sizeof(double) ) )

  var xdim : int = nx+1

  var Avals : double[3]
  Avals[0] = alpha
  Avals[1] = beta
  Avals[2] = gamma

  var counter : int64 = 0
  matrix.rowptr[0] = counter

  for iz = 0, nz do
    for iy = 0, ny do
      for row = 0, nx do
        for j = 0, 3 do
          var col : int = row + j - 1
          var gcol : int64 = (col + nx)%nx + iy*xdim + iz*xdim*ny
          matrix.colind[counter] = gcol
          matrix.nzval [counter] = Avals[j]
          counter = counter + 1
        end
        var grow : int64 = row + iy*xdim + iz*xdim*ny
        matrix.rowptr[grow+1] = matrix.rowptr[grow] + 3
      end
      -- For the last point
      var gcol : int64 = nx + iy*xdim + iz*xdim*ny
      matrix.colind[counter] = gcol
      matrix.nzval [counter] = 1.0
      counter = counter + 1

      gcol = 0 + iy*xdim + iz*xdim*ny
      matrix.colind[counter] = gcol
      matrix.nzval [counter] = -1.0
      counter = counter + 1

      var grow : int64 = nx + iy*xdim + iz*xdim*ny
      matrix.rowptr[grow+1] = matrix.rowptr[grow] + 2
    end
  end

  return matrix
end

local terra get_base_pointer_2d(pr   : c.legion_physical_region_t,
                                fid  : c.legion_field_id_t,
                                rect : c.legion_rect_2d_t)
  var subrect : c.legion_rect_2d_t
  var offsets : c.legion_byte_offset_t[1]
  var accessor = c.legion_physical_region_get_field_accessor_generic(pr, fid)
  var base_pointer =
    [&superlu.c.superlu_vars_t](c.legion_accessor_generic_raw_rect_ptr_2d(
      accessor, rect, &subrect, &(offsets[0])))
  c.legion_accessor_generic_destroy(accessor)
  return base_pointer
end

local terra get_base_pointer_3d(pr   : c.legion_physical_region_t[1],
                                fid  : c.legion_field_id_t[1],
                                rect : c.legion_rect_3d_t)
  var subrect : c.legion_rect_3d_t
  var offsets : c.legion_byte_offset_t[3]
  var accessor = c.legion_physical_region_get_field_accessor_generic(pr[0], fid[0])
  var base_pointer =
    [&double](c.legion_accessor_generic_raw_rect_ptr_3d(
      accessor, rect, &subrect, &(offsets[0])))
  c.legion_accessor_generic_destroy(accessor)
  return base_pointer
end

terra superlu.initialize_superlu_vars( matrix : superlu.CSR_matrix,
                                       Nsize  : int64,
                                       pr1    : c.legion_physical_region_t[1],
                                       fid1   : c.legion_field_id_t[1],
                                       pr2    : c.legion_physical_region_t[1],
                                       fid2   : c.legion_field_id_t[1],
                                       rect   : c.legion_rect_3d_t,
                                       prv    : c.legion_physical_region_t,
                                       fidv   : c.legion_field_id_t,
                                       rectv  : c.legion_rect_2d_t )
  var b = get_base_pointer_3d(pr1, fid1, rect)
  var x = get_base_pointer_3d(pr2, fid2, rect)
  var vars = get_base_pointer_2d(prv, fidv, rectv)
  superlu.c.initialize_superlu_vars(matrix.nzval, matrix.colind, matrix.rowptr, Nsize, matrix.nnz, b, x, vars)
end

terra superlu.MatrixSolve( pr1    : c.legion_physical_region_t[1],
                           fid1   : c.legion_field_id_t[1],
                           pr2    : c.legion_physical_region_t[1],
                           fid2   : c.legion_field_id_t[1],
                           rect   : c.legion_rect_3d_t,
                           matrix : superlu.CSR_matrix,
                           nx     : int,
                           ny     : int,
                           nz     : int,
                           prv    : c.legion_physical_region_t,
                           fidv   : c.legion_field_id_t,
                           rectv  : c.legion_rect_2d_t )
  var b = get_base_pointer_3d(pr1, fid1, rect)
  var x = get_base_pointer_3d(pr2, fid2, rect)
  var vars = get_base_pointer_2d(prv, fidv, rectv)
  superlu.c.MatrixSolve(b, x, matrix.nzval, nx, ny, nz, vars)
end

return superlu
