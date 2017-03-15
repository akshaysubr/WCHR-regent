import "regent"

local superlu = {}
do
  local superlu_root = os.getenv('SUPERLU_PATH') or "/opt/SuperLU_5.2.1"
  --local superlu_library = "-lsuperlu"
  local superlu_library = "-L" .. superlu_root .. " -lsuperlu"
  local superlu_include_dir = superlu_root
  local root_dir = arg[0]:match(".*/") or "./"
  local superlu_util_cc = root_dir .. "superlu_util.c"
  superlu_util_so = os.tmpname() .. ".so"
  local cc = os.getenv('CC') or 'cc'
  local cc_flags = "-O3 -Wall -Werror -std=c99"
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
    terralib.linklibrary("libblas.so")
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

terra allocate_int( size: int64 )
  return [&int] ( c.malloc( size * sizeof(int) ) )
end

terra allocate_double( size: int64 )
  return [&double] ( c.malloc( size * sizeof(double) ) )
end

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

task superlu.initialize_matrix_char_x( matrix : region(ispace(int2d), superlu.CSR_matrix),
                                       alpha  : double,
                                       beta   : double,
                                       gamma  : double,
                                       nx     : int64,
                                       ny     : int64,
                                       nz     : int64 )
where
  reads writes (matrix)
do
  -- var matrix : superlu.CSR_matrix
  var pr = matrix.ispace.bounds.hi.x
  var pc = matrix.ispace.bounds.hi.y

  var Nsize : int64 = 5*(nx+1)*ny*nz
  matrix[{pr,pc}].nnz = (8*3*nx+10)*ny*nz
  matrix[{pr,pc}].rowptr = allocate_int( Nsize+1 ) -- [&int] ( c.malloc ( (Nsize+1) * sizeof(int) ) )
  matrix[{pr,pc}].colind = allocate_int( matrix[{pr,pc}].nnz ) -- [&int] ( c.malloc ( matrix[{pr,pc}].nnz * sizeof(int) ) )
  matrix[{pr,pc}].nzval  = allocate_double( matrix[{pr,pc}].nnz ) -- [&double] ( c.malloc ( matrix[{pr,pc}].nnz * sizeof(double) ) )

  var dim : int = nx+1

  var Avals : double[3]
  Avals[0] = alpha
  Avals[1] = beta
  Avals[2] = gamma

  -- c.printf("In initialize_matrix_char_x\n")

  matrix[{pr,pc}].rowptr[0] = 0

  for iz = 0, nz do
    for iy = 0, ny do
      for brow = 0, nx do
        var grow : int64 = 5*brow + iy*5*dim + iz*5*dim*ny
        var bcounter : int64 = 8*3*brow + iy*(8*3*nx+10) + iz*(8*3*nx+10)*ny

        -- rho
        for j = 0, 3 do
          var bcol : int64 = brow + j - 1
          var gcol : int64 = ((bcol + nx)%nx)*5 + iy*5*dim + iz*5*dim*ny -- Top of the block
          var counter : int64 = bcounter + 2*j

          matrix[{pr,pc}].colind[counter] = gcol+1
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (-0.5)

          matrix[{pr,pc}].colind[counter+1] = gcol+4
          matrix[{pr,pc}].nzval [counter+1] = Avals[j] * (0.5)
        end
        bcounter = bcounter + 3*2
        matrix[{pr,pc}].rowptr[grow+1] = bcounter

        -- u
        for j = 0, 3 do
          var bcol : int64 = brow + j - 1
          var gcol : int64 = ((bcol + nx)%nx)*5 + iy*5*dim + iz*5*dim*ny -- Top of the block
          var counter : int64 = bcounter + 2*j

          matrix[{pr,pc}].colind[counter] = gcol
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (1.0)

          matrix[{pr,pc}].colind[counter+1] = gcol+4
          matrix[{pr,pc}].nzval [counter+1] = Avals[j] * (-1.0)
        end
        bcounter = bcounter + 3*2
        matrix[{pr,pc}].rowptr[grow+2] = bcounter

        -- v
        for j = 0, 3 do
          var bcol : int64 = brow + j - 1
          var gcol : int64 = ((bcol + nx)%nx)*5 + iy*5*dim + iz*5*dim*ny -- Top of the block
          var counter : int64 = bcounter + j

          matrix[{pr,pc}].colind[counter] = gcol+2
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (1.0)
        end
        bcounter = bcounter + 3*1
        matrix[{pr,pc}].rowptr[grow+3] = bcounter

        -- w
        for j = 0, 3 do
          var bcol : int64 = brow + j - 1
          var gcol : int64 = ((bcol + nx)%nx)*5 + iy*5*dim + iz*5*dim*ny -- Top of the block
          var counter : int64 = bcounter + j

          matrix[{pr,pc}].colind[counter] = gcol+3
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (1.0)
        end
        bcounter = bcounter + 3*1
        matrix[{pr,pc}].rowptr[grow+4] = bcounter

        -- p
        for j = 0, 3 do
          var bcol : int64 = brow + j - 1
          var gcol : int64 = ((bcol + nx)%nx)*5 + iy*5*dim + iz*5*dim*ny -- Top of the block
          var counter : int64 = bcounter + 2*j

          matrix[{pr,pc}].colind[counter] = gcol+1
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (0.5)

          matrix[{pr,pc}].colind[counter+1] = gcol+4
          matrix[{pr,pc}].nzval [counter+1] = Avals[j] * (0.5)
        end
        bcounter = bcounter + 3*2
        matrix[{pr,pc}].rowptr[grow+5] = bcounter
      end
      -- For the last point
      
      var bcounter : int64 = 8*3*nx + iy*(8*3*nx+10) + iz*(8*3*nx+10)*ny
      for pvar = 0,5 do
        var gcol : int64 = 5*nx+pvar + iy*5*dim + iz*5*dim*ny
        var counter : int64 = bcounter + 2*pvar
        matrix[{pr,pc}].colind[counter] = gcol
        matrix[{pr,pc}].nzval [counter] = 1.0

        gcol = 5*0+pvar + iy*5*dim + iz*5*dim*ny
        matrix[{pr,pc}].colind[counter+1] = gcol
        matrix[{pr,pc}].nzval [counter+1] = -1.0

        var grow : int64 = 5*nx+pvar + iy*5*dim + iz*5*dim*ny
        matrix[{pr,pc}].rowptr[grow+1] = counter + 2
      end
    end
  end

  -- c.printf("nzval = numpy.array([")
  -- for i = 0,matrix[{pr,pc}].nnz do
  --   c.printf("%g, ", matrix[{pr,pc}].nzval[i])
  -- end
  -- c.printf("])\n\n")
  -- c.printf("colind = numpy.array([")
  -- for i = 0,matrix[{pr,pc}].nnz do
  --   c.printf("%d, ", matrix[{pr,pc}].colind[i])
  -- end
  -- c.printf("])\n\n")
  -- c.printf("rowptr = numpy.array([")
  -- for i = 0,Nsize+1 do
  --   c.printf("%d, ", matrix[{pr,pc}].rowptr[i])
  -- end
  -- c.printf("])\n\n")
end

task superlu.initialize_matrix_char_y( matrix : region(ispace(int2d), superlu.CSR_matrix),
                                       alpha  : double,
                                       beta   : double,
                                       gamma  : double,
                                       nx     : int64,
                                       ny     : int64,
                                       nz     : int64 )
where
  reads writes (matrix)
do
  -- var matrix : superlu.CSR_matrix
  var pr = matrix.ispace.bounds.hi.x
  var pc = matrix.ispace.bounds.hi.y

  var Nsize : int64 = 5*(ny+1)*nx*nz
  matrix[{pr,pc}].nnz = (8*3*ny+10)*nx*nz
  matrix[{pr,pc}].rowptr = allocate_int( Nsize+1 ) -- [&int] ( c.malloc ( (Nsize+1) * sizeof(int) ) )
  matrix[{pr,pc}].colind = allocate_int( matrix[{pr,pc}].nnz ) -- [&int] ( c.malloc ( matrix[{pr,pc}].nnz * sizeof(int) ) )
  matrix[{pr,pc}].nzval  = allocate_double( matrix[{pr,pc}].nnz ) -- [&double] ( c.malloc ( matrix[{pr,pc}].nnz * sizeof(double) ) )

  var dim : int = ny+1

  var Avals : double[3]
  Avals[0] = alpha
  Avals[1] = beta
  Avals[2] = gamma

  -- c.printf("In initialize_matrix_char_y\n")

  matrix[{pr,pc}].rowptr[0] = 0

  for iz = 0, nz do
    for iy = 0, ny do
      for brow = 0, nx do
        var grow : int64 = 5*brow + iy*5*nx + iz*5*dim*nx
        var bcounter : int64 = 8*3*brow + iy*(8*3)*nx + iz*(8*3*ny+10)*nx

        -- rho
        for j = 0, 3 do
          var bcol : int64 = iy + j - 1
          var gcol : int64 = brow*5 + ((bcol+ny)%ny)*5*nx + iz*5*dim*nx -- Top of the block
          var counter : int64 = bcounter + 2*j

          matrix[{pr,pc}].colind[counter] = gcol+2
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (-0.5)

          matrix[{pr,pc}].colind[counter+1] = gcol+4
          matrix[{pr,pc}].nzval [counter+1] = Avals[j] * (0.5)
        end
        bcounter = bcounter + 3*2
        matrix[{pr,pc}].rowptr[grow+1] = bcounter

        -- u
        for j = 0, 3 do
          var bcol : int64 = iy + j - 1
          var gcol : int64 = brow*5 + ((bcol+ny)%ny)*5*nx + iz*5*dim*nx -- Top of the block
          var counter : int64 = bcounter + j

          matrix[{pr,pc}].colind[counter] = gcol+1
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (1.0)
        end
        bcounter = bcounter + 3*1
        matrix[{pr,pc}].rowptr[grow+2] = bcounter

        -- v
        for j = 0, 3 do
          var bcol : int64 = iy + j - 1
          var gcol : int64 = brow*5 + ((bcol+ny)%ny)*5*nx + iz*5*dim*nx -- Top of the block
          var counter : int64 = bcounter + 2*j

          matrix[{pr,pc}].colind[counter] = gcol
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (1.0)

          matrix[{pr,pc}].colind[counter+1] = gcol+4
          matrix[{pr,pc}].nzval [counter+1] = Avals[j] * (-1.0)
        end
        bcounter = bcounter + 3*2
        matrix[{pr,pc}].rowptr[grow+3] = bcounter

        -- w
        for j = 0, 3 do
          var bcol : int64 = iy + j - 1
          var gcol : int64 = brow*5 + ((bcol+ny)%ny)*5*nx + iz*5*dim*nx -- Top of the block
          var counter : int64 = bcounter + j

          matrix[{pr,pc}].colind[counter] = gcol+3
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (1.0)
        end
        bcounter = bcounter + 3*1
        matrix[{pr,pc}].rowptr[grow+4] = bcounter

        -- p
        for j = 0, 3 do
          var bcol : int64 = iy + j - 1
          var gcol : int64 = brow*5 + ((bcol+ny)%ny)*5*nx + iz*5*dim*nx -- Top of the block
          var counter : int64 = bcounter + 2*j

          matrix[{pr,pc}].colind[counter] = gcol+2
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (0.5)

          matrix[{pr,pc}].colind[counter+1] = gcol+4
          matrix[{pr,pc}].nzval [counter+1] = Avals[j] * (0.5)
        end
        bcounter = bcounter + 3*2
        matrix[{pr,pc}].rowptr[grow+5] = bcounter
      end

    end
      
    for brow = 0, nx do
      -- For the last point
      var bcounter : int64 = 10*brow + ny*(8*3)*nx + iz*(8*3*ny+10)*nx
      for pvar = 0,5 do
        var gcol : int64 = pvar + brow*5 + ny*5*nx + iz*5*dim*nx
        var counter : int64 = bcounter + 2*pvar
        matrix[{pr,pc}].colind[counter] = gcol
        matrix[{pr,pc}].nzval [counter] = 1.0

        gcol = pvar + brow*5 + 0*5*nx + iz*5*dim*nx
        matrix[{pr,pc}].colind[counter+1] = gcol
        matrix[{pr,pc}].nzval [counter+1] = -1.0

        var grow : int64 = pvar + 5*brow + ny*5*nx + iz*5*dim*nx
        matrix[{pr,pc}].rowptr[grow+1] = counter + 2
      end
    end
  end

  -- c.printf("nzval = numpy.array([")
  -- for i = 0,matrix[{pr,pc}].nnz do
  --   c.printf("%g, ", matrix[{pr,pc}].nzval[i])
  -- end
  -- c.printf("])\n\n")
  -- c.printf("colind = numpy.array([")
  -- for i = 0,matrix[{pr,pc}].nnz do
  --   c.printf("%d, ", matrix[{pr,pc}].colind[i])
  -- end
  -- c.printf("])\n\n")
  -- c.printf("rowptr = numpy.array([")
  -- for i = 0,Nsize+1 do
  --   c.printf("%d, ", matrix[{pr,pc}].rowptr[i])
  -- end
  -- c.printf("])\n\n")
end

task superlu.initialize_matrix_char_z( matrix : region(ispace(int2d), superlu.CSR_matrix),
                                       alpha  : double,
                                       beta   : double,
                                       gamma  : double,
                                       nx     : int64,
                                       ny     : int64,
                                       nz     : int64 )
where
  reads writes (matrix)
do
  -- var matrix : superlu.CSR_matrix
  var pr = matrix.ispace.bounds.hi.x
  var pc = matrix.ispace.bounds.hi.y

  var Nsize : int64 = 5*(nz+1)*nx*ny
  matrix[{pr,pc}].nnz = (8*3*nz+10)*nx*ny
  matrix[{pr,pc}].rowptr = allocate_int( Nsize+1 ) -- [&int] ( c.malloc ( (Nsize+1) * sizeof(int) ) )
  matrix[{pr,pc}].colind = allocate_int( matrix[{pr,pc}].nnz ) -- [&int] ( c.malloc ( matrix[{pr,pc}].nnz * sizeof(int) ) )
  matrix[{pr,pc}].nzval  = allocate_double( matrix[{pr,pc}].nnz ) -- [&double] ( c.malloc ( matrix[{pr,pc}].nnz * sizeof(double) ) )

  var dim : int = nz+1

  var Avals : double[3]
  Avals[0] = alpha
  Avals[1] = beta
  Avals[2] = gamma

  -- c.printf("In initialize_matrix_char_z\n")

  matrix[{pr,pc}].rowptr[0] = 0

  for iz = 0, nz do
    for iy = 0, ny do
      for brow = 0, nx do
        var grow : int64 = 5*brow + iy*5*nx + iz*5*nx*ny
        var bcounter : int64 = 8*3*brow + iy*(8*3)*nx + iz*(8*3)*ny*nx

        -- rho
        for j = 0, 3 do
          var bcol : int64 = iz + j - 1
          var gcol : int64 = brow*5 + iy*5*nx + ((bcol+nz)%nz)*5*nx*ny -- Top of the block
          var counter : int64 = bcounter + 2*j

          matrix[{pr,pc}].colind[counter] = gcol+3
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (-0.5)

          matrix[{pr,pc}].colind[counter+1] = gcol+4
          matrix[{pr,pc}].nzval [counter+1] = Avals[j] * (0.5)
        end
        bcounter = bcounter + 3*2
        matrix[{pr,pc}].rowptr[grow+1] = bcounter

        -- u
        for j = 0, 3 do
          var bcol : int64 = iz + j - 1
          var gcol : int64 = brow*5 + iy*5*nx + ((bcol+nz)%nz)*5*nx*ny -- Top of the block
          var counter : int64 = bcounter + j

          matrix[{pr,pc}].colind[counter] = gcol+1
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (1.0)
        end
        bcounter = bcounter + 3*1
        matrix[{pr,pc}].rowptr[grow+2] = bcounter

        -- v
        for j = 0, 3 do
          var bcol : int64 = iz + j - 1
          var gcol : int64 = brow*5 + iy*5*nx + ((bcol+nz)%nz)*5*nx*ny -- Top of the block
          var counter : int64 = bcounter + j

          matrix[{pr,pc}].colind[counter] = gcol+2
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (1.0)
        end
        bcounter = bcounter + 3*1
        matrix[{pr,pc}].rowptr[grow+3] = bcounter

        -- w
        for j = 0, 3 do
          var bcol : int64 = iz + j - 1
          var gcol : int64 = brow*5 + iy*5*nx + ((bcol+nz)%nz)*5*nx*ny -- Top of the block
          var counter : int64 = bcounter + 2*j

          matrix[{pr,pc}].colind[counter] = gcol
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (1.0)

          matrix[{pr,pc}].colind[counter+1] = gcol+4
          matrix[{pr,pc}].nzval [counter+1] = Avals[j] * (-1.0)
        end
        bcounter = bcounter + 3*2
        matrix[{pr,pc}].rowptr[grow+4] = bcounter

        -- p
        for j = 0, 3 do
          var bcol : int64 = iz + j - 1
          var gcol : int64 = brow*5 + iy*5*nx + ((bcol+nz)%nz)*5*nx*ny -- Top of the block
          var counter : int64 = bcounter + 2*j

          matrix[{pr,pc}].colind[counter] = gcol+3
          matrix[{pr,pc}].nzval [counter] = Avals[j] * (0.5)

          matrix[{pr,pc}].colind[counter+1] = gcol+4
          matrix[{pr,pc}].nzval [counter+1] = Avals[j] * (0.5)
        end
        bcounter = bcounter + 3*2
        matrix[{pr,pc}].rowptr[grow+5] = bcounter

      end
    end
  end

  for iy = 0, ny do
    for brow = 0, nx do
      -- For the last point
      var bcounter : int64 = 10*brow + 10*iy*nx + nz*(8*3)*ny*nx
      for pvar = 0,5 do
        var gcol : int64 = pvar + brow*5 + iy*5*nx + nz*5*nx*ny
        var counter : int64 = bcounter + 2*pvar
        matrix[{pr,pc}].colind[counter] = gcol
        matrix[{pr,pc}].nzval [counter] = 1.0

        gcol = pvar + brow*5 + iy*5*nx + 0*5*nx*ny
        matrix[{pr,pc}].colind[counter+1] = gcol
        matrix[{pr,pc}].nzval [counter+1] = -1.0

        var grow : int64 = pvar + brow*5 + iy*5*nx + nz*5*nx*ny
        matrix[{pr,pc}].rowptr[grow+1] = counter + 2
      end
    end
  end

  -- c.printf("nzval = numpy.array([")
  -- for i = 0,matrix[{pr,pc}].nnz do
  --   c.printf("%g, ", matrix[{pr,pc}].nzval[i])
  -- end
  -- c.printf("])\n\n")
  -- c.printf("colind = numpy.array([")
  -- for i = 0,matrix[{pr,pc}].nnz do
  --   c.printf("%d, ", matrix[{pr,pc}].colind[i])
  -- end
  -- c.printf("])\n\n")
  -- c.printf("rowptr = numpy.array([")
  -- for i = 0,Nsize+1 do
  --   c.printf("%d, ", matrix[{pr,pc}].rowptr[i])
  -- end
  -- c.printf("])\n\n")
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

__demand(__external)
task superlu.initialize_superlu_vars( matrix : superlu.CSR_matrix,
                                      Nsize  : int64,
                                      r_rhs  : region(ispace(int3d), primitive),
                                      r_sol  : region(ispace(int3d), primitive),
                                      slu    : region(ispace(int2d), superlu.c.superlu_vars_t) )
where
  reads(r_rhs), writes(r_sol), reads writes(slu)
do
  var b = get_base_pointer_3d(__physical(r_rhs.rho), __fields(r_rhs.rho), r_rhs.bounds)
  var x = get_base_pointer_3d(__physical(r_sol.rho), __fields(r_sol.rho), r_sol.bounds)
  var vars = get_base_pointer_2d(__physical(slu)[0], __fields(slu)[0], slu.bounds)
  superlu.c.initialize_superlu_vars(matrix.nzval, matrix.colind, matrix.rowptr, Nsize, matrix.nnz, b, x, vars)
end

task superlu.init_superlu_vars( matrix : region(ispace(int2d), superlu.CSR_matrix),
                                Nsize  : int64,
                                r_rhs  : region(ispace(int3d), primitive),
                                r_sol  : region(ispace(int3d), primitive),
                                slu    : region(ispace(int2d), superlu.c.superlu_vars_t) )
where
  reads(matrix, r_rhs), writes(r_sol), reads writes(slu)
do
  var pr = matrix.ispace.bounds.hi.x
  var pc = matrix.ispace.bounds.hi.y

  -- c.printf("In init_superlu_vars\n")

  superlu.initialize_superlu_vars( matrix[{pr,pc}], Nsize, r_rhs, r_sol, slu )
end

__demand(__external)
task superlu.MatrixSolve( r_rhs  : region(ispace(int3d), primitive),
                          r_sol  : region(ispace(int3d), primitive),
                          matrix : superlu.CSR_matrix,
                          nx     : int,
                          ny     : int,
                          nz     : int,
                          slu    : region(ispace(int2d), superlu.c.superlu_vars_t) )
where
  reads(r_rhs), writes(r_sol), reads writes(slu)
do
  var b = get_base_pointer_3d(__physical(r_rhs.rho), __fields(r_rhs.rho), r_rhs.bounds)
  var x = get_base_pointer_3d(__physical(r_sol.rho), __fields(r_sol.rho), r_sol.bounds)
  var vars = get_base_pointer_2d(__physical(slu)[0], __fields(slu)[0], slu.bounds)
  superlu.c.MatrixSolve(b, x, matrix.nzval, nx, ny, nz, vars)
end

return superlu
