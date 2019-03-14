import "regent"

local use_io = os.getenv('USE_IO') == '1' or false
if not use_io then
  task read_primitive( r_prim     : region(ispace(int3d), primitive),
                       fileprefix : rawstring,
                       vizcount   : int32,
                       n_ghosts   : int64,
                       pencil     : int2d )
  where
    writes( r_prim )
  do
    return 1
  end

  task read_coords( r_coords   : region(ispace(int3d), coordinates),
                    fileprefix : rawstring,
                    n_ghosts   : int64,
                    pencil     : int2d )
  where
    writes( r_coords )
  do
    return 1
  end

  task write_primitive( r_prim     : region(ispace(int3d), primitive),
                        fileprefix : rawstring,
                        vizcount   : int32,
                        n_ghosts   : int64,
                        pencil     : int2d )
  where
    reads( r_prim )
  do
    return 1
  end

  task write_coords( r_coords   : region(ispace(int3d), coordinates),
                     fileprefix : rawstring,
                     n_ghosts   : int64,
                     pencil     : int2d )
  where
    reads( r_coords )
  do
    return 1
  end

  return false
end

require("fields")

local c  = regentlib.c
local st = terralib.includec("string.h")

local hdf5_path = os.getenv('HDF_ROOT')
local hdf5_include = hdf5_path .. "/include"

hdf5 = terralib.includec("hdf5.h", {"-I", hdf5_include})

-- there's some funny business in hdf5.h that prevents terra from being able to
--  see some of the #define's, so we fix it here, and hope the HDF5 folks don't
--  change the internals very often...
hdf5.H5F_ACC_TRUNC = 2
hdf5.H5T_STD_I32LE = hdf5.H5T_STD_I32LE_g
hdf5.H5T_STD_I64LE = hdf5.H5T_STD_I64LE_g
hdf5.H5T_IEEE_F64LE = hdf5.H5T_IEEE_F64LE_g
hdf5.H5P_DEFAULT = 0



terra generate_hdf5_file(filename : rawstring, NX : int64, NY : int64, NZ : int64)
  var fid = hdf5.H5Fcreate(filename, hdf5.H5F_ACC_TRUNC, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)

  var dims : hdf5.hsize_t[3]
  dims[0] = NZ
  dims[1] = NY
  dims[2] = NX
  var did = hdf5.H5Screate_simple(3, dims, [&uint64](0))

  var ds1id = hdf5.H5Dcreate2(fid, "rho", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds1id)

  var ds2id = hdf5.H5Dcreate2(fid, "u", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds2id)

  var ds3id = hdf5.H5Dcreate2(fid, "v", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds3id)

  var ds4id = hdf5.H5Dcreate2(fid, "w", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds4id)

  var ds5id = hdf5.H5Dcreate2(fid, "p", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds5id)

  hdf5.H5Sclose(did)
  hdf5.H5Fclose(fid)
end



terra file_exists( filename : rawstring )
  var file_handle = c.fopen(filename, 'r')
  if ( file_handle == nil ) then
    return false
  else
    return true
  end
end



local terra read_grid( nx : &int64, ny : &int64, nz : &int64, f : &c.FILE )
  c.fscanf(f, "%d", nx)
  c.fscanf(f, "%d", ny)
  c.fscanf(f, "%d", nz)
end



local terra read_coord_data( x : &double, y : &double, z : &double,
                             f : &c.FILE )
  c.fscanf(f, "%lf", x)
  c.fscanf(f, "%lf", y)
  c.fscanf(f, "%lf", z)
end



local terra read_primitive_data( rho : &double, u : &double, v : &double, w : &double, p : &double,
                                 f : &c.FILE )
  c.fscanf(f, "%lf", rho)
  c.fscanf(f, "%lf", u)
  c.fscanf(f, "%lf", v)
  c.fscanf(f, "%lf", w)
  c.fscanf(f, "%lf", p)
end



task write_primitive( r_prim     : region(ispace(int3d), primitive),
                      fileprefix : rawstring,
                      vizcount   : int32,
                      n_ghosts   : int64,
                      pencil     : int2d )
where
  reads( r_prim )
do
  var Nx = r_prim.ispace.bounds.hi.x - r_prim.ispace.bounds.lo.x + 1
  var Ny = r_prim.ispace.bounds.hi.y - r_prim.ispace.bounds.lo.y + 1
  var Nz = r_prim.ispace.bounds.hi.z - r_prim.ispace.bounds.lo.z + 1

  -- var vizstr : &int8 = [&int8] ( c.malloc(21*8) )
  -- c.sprintf(vizstr, "%04d_px%04d_pz%04d.h5", vizcount, pencil.x, pencil.y)
  var vizstr : &int8 = [&int8] ( c.malloc(22*8) )
  c.sprintf(vizstr, "%04d_px%04d_pz%04d.dat", vizcount, pencil.x - 1, pencil.y - 1)
  var filename : &int8 = [&int8] ( c.malloc(256*8) )
  st.strcpy(filename, fileprefix)
  st.strcat(filename, vizstr)
  c.free(vizstr)

  c.printf("\nWriting visualization file %s\n", filename)

  var file_handle = c.fopen(filename, 'w')

  var bounds_c = r_prim.ispace.bounds
  c.fprintf(file_handle, "%d %d %d\n", bounds_c.lo.x - n_ghosts, bounds_c.lo.y - n_ghosts, bounds_c.lo.z - n_ghosts )
  c.fprintf(file_handle, "%d %d %d\n", bounds_c.hi.x - n_ghosts, bounds_c.hi.y - n_ghosts, bounds_c.hi.z - n_ghosts )

  for k = bounds_c.lo.z, bounds_c.hi.z+1 do
    for j = bounds_c.lo.y, bounds_c.hi.y+1 do
      for i = bounds_c.lo.x, bounds_c.hi.x+1 do
        c.fprintf(file_handle, "%26.16e %26.16e %26.16e %26.16e %26.16e\n", r_prim[{i,j,k}].rho, r_prim[{i,j,k}].u, r_prim[{i,j,k}].v, r_prim[{i,j,k}].w, r_prim[{i,j,k}].p)
      end
    end
  end

  c.fclose(file_handle) 
  c.free(filename)

  -- generate_hdf5_file(filename, Nx, Ny, Nz)

  -- var grid = ispace(int3d, {x = Nx, y = Ny, z = Nz})
  -- var tmp_prim = region(grid, primitive)

  -- attach(hdf5, tmp_prim.{rho, u, v, w, p}, filename, regentlib.file_read_write)
  -- copy(r_prim.{rho, u, v, w, p}, tmp_prim.{rho, u, v, w, p})
  -- detach(hdf5, tmp_prim.{rho, u, v, w, p})

  -- c.free(filename)
  -- __delete(tmp_prim)
  return 1
end



task read_primitive( r_prim     : region(ispace(int3d), primitive),
                     fileprefix : rawstring,
                     vizcount   : int32,
                     n_ghosts   : int64,
                     pencil     : int2d )
where
  writes( r_prim )
do
  var Nx = r_prim.ispace.bounds.hi.x - r_prim.ispace.bounds.lo.x + 1
  var Ny = r_prim.ispace.bounds.hi.y - r_prim.ispace.bounds.lo.y + 1
  var Nz = r_prim.ispace.bounds.hi.z - r_prim.ispace.bounds.lo.z + 1

  -- var vizcount : int = 0
  -- var found_restart = false
  -- while not found_restart do
  --   var vizstr : &int8 = [&int8] ( c.malloc(22*8) )
  --   var filename : &int8 = [&int8] ( c.malloc(256*8) )
  --   c.sprintf(vizstr, "%04d_px%04d_pz%04d.dat", vizcount, pencil.x, pencil.y)
  --   st.strcpy(filename, fileprefix)
  --   st.strcat(filename, vizstr)
  --   c.printf("Trying to read %s\n", filename)
  --   found_restart = (not file_exists(filename))
  --   vizcount += 1
  --   c.free(vizstr)
  --   c.free(filename)
  -- end
  -- vizcount -= 2
  -- if vizcount < 0 then
  --   c.printf("Error: Could not find any restart file\n")
  --   c.abort()
  -- end

  var vizstr : &int8 = [&int8] ( c.malloc(22*8) )
  var filename : &int8 = [&int8] ( c.malloc(256*8) )
  c.sprintf(vizstr, "%04d_px%04d_pz%04d.dat", vizcount, pencil.x - 1, pencil.y - 1)
  st.strcpy(filename, fileprefix)
  st.strcat(filename, vizstr)

  var file_handle : &c.FILE
  if file_exists( filename ) then
    file_handle = c.fopen(filename, 'r')
  else
    c.printf("Error: Failed to open \"%s\"\n", filename)
    c.abort()
  end
  c.printf("\nReading visualization file %s\n", filename)
  c.free(vizstr)

  var ix : int64[1]
  var iy : int64[1]
  var iz : int64[1]
  ix[0], iy[0], iz[0] = 0, 0, 0

  var bounds_c = r_prim.ispace.bounds
  -- c.fprintf(file_handle, "%d %d %d\n", bounds_c.lo.x, bounds_c.lo.y, bounds_c.lo.z )
  read_grid(ix, iy, iz, file_handle)
  regentlib.assert( (ix[0] == bounds_c.lo.x - n_ghosts), "Grid size in restart file does not match partition size in x" )
  regentlib.assert( (iy[0] == bounds_c.lo.y - n_ghosts), "Grid size in restart file does not match partition size in y" )
  regentlib.assert( (iz[0] == bounds_c.lo.z - n_ghosts), "Grid size in restart file does not match partition size in z" )

  -- c.fprintf(file_handle, "%d %d %d\n", bounds_c.hi.x, bounds_c.hi.y, bounds_c.hi.z )
  read_grid(ix, iy, iz, file_handle)
  regentlib.assert( (ix[0] == bounds_c.hi.x - n_ghosts), "Grid size in restart file does not match partition size in x" )
  regentlib.assert( (iy[0] == bounds_c.hi.y - n_ghosts), "Grid size in restart file does not match partition size in y" )
  regentlib.assert( (iz[0] == bounds_c.hi.z - n_ghosts), "Grid size in restart file does not match partition size in z" )

  var rho_dat : double[1]
  var u_dat   : double[1]
  var v_dat   : double[1]
  var w_dat   : double[1]
  var p_dat   : double[1]
  rho_dat[0], u_dat[0], v_dat[0], w_dat[0], p_dat[0] = 0., 0., 0., 0., 0.

  for k = bounds_c.lo.z, bounds_c.hi.z+1 do
    for j = bounds_c.lo.y, bounds_c.hi.y+1 do
      for i = bounds_c.lo.x, bounds_c.hi.x+1 do
        -- c.fprintf(file_handle, "%26.16e %26.16e %26.16e %26.16e %26.16e\n", r_prim[{i,j,k}].rho, r_prim[{i,j,k}].u, r_prim[{i,j,k}].v, r_prim[{i,j,k}].w, r_prim[{i,j,k}].p)
        read_primitive_data( rho_dat, u_dat, v_dat, w_dat, p_dat, file_handle)
        r_prim[{i,j,k}].rho = rho_dat[0]
        r_prim[{i,j,k}].u   = u_dat[0]
        r_prim[{i,j,k}].v   = v_dat[0]
        r_prim[{i,j,k}].w   = w_dat[0]
        r_prim[{i,j,k}].p   = p_dat[0]
      end
    end
  end

  c.fclose(file_handle) 
  c.free(filename)

  return vizcount
end



terra generate_hdf5_file_coords(filename : rawstring, NX : int64, NY : int64, NZ : int64)
  var fid = hdf5.H5Fcreate(filename, hdf5.H5F_ACC_TRUNC, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)

  var dims : hdf5.hsize_t[3]
  dims[0] = NZ
  dims[1] = NY
  dims[2] = NX
  var did = hdf5.H5Screate_simple(3, dims, [&uint64](0))

  var ds1id = hdf5.H5Dcreate2(fid, "x_c", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds1id)

  var ds2id = hdf5.H5Dcreate2(fid, "y_c", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds2id)

  var ds3id = hdf5.H5Dcreate2(fid, "z_c", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  hdf5.H5Dclose(ds3id)

  hdf5.H5Sclose(did)
  hdf5.H5Fclose(fid)
end



task write_coords( r_coords   : region(ispace(int3d), coordinates),
                   fileprefix : rawstring,
                   n_ghosts   : int64,
                   pencil     : int2d )
where
  reads( r_coords )
do
  var Nx = r_coords.ispace.bounds.hi.x - r_coords.ispace.bounds.lo.x + 1
  var Ny = r_coords.ispace.bounds.hi.y - r_coords.ispace.bounds.lo.y + 1
  var Nz = r_coords.ispace.bounds.hi.z - r_coords.ispace.bounds.lo.z + 1

  -- var vizstr : &int8 = [&int8] ( c.malloc(23*8) )
  -- c.sprintf(vizstr, "coords_px%04d_pz%04d.h5", pencil.x, pencil.y)
  var vizstr : &int8 = [&int8] ( c.malloc(24*8) )
  c.sprintf(vizstr, "coords_px%04d_pz%04d.dat", pencil.x - 1, pencil.y - 1)
  var filename : &int8 = [&int8] ( c.malloc(256*8) )
  st.strcpy(filename, fileprefix)
  st.strcat(filename, vizstr)
  c.free(vizstr)

  var file_handle = c.fopen(filename, 'w')

  var bounds_c = r_coords.ispace.bounds
  c.fprintf(file_handle, "%d %d %d\n", bounds_c.lo.x - n_ghosts, bounds_c.lo.y - n_ghosts, bounds_c.lo.z - n_ghosts )
  c.fprintf(file_handle, "%d %d %d\n", bounds_c.hi.x - n_ghosts, bounds_c.hi.y - n_ghosts, bounds_c.hi.z - n_ghosts )

  for k = bounds_c.lo.z, bounds_c.hi.z+1 do
    for j = bounds_c.lo.y, bounds_c.hi.y+1 do
      for i = bounds_c.lo.x, bounds_c.hi.x+1 do
        c.fprintf(file_handle, "%26.16e %26.16e %26.16e\n", r_coords[{i,j,k}].x_c, r_coords[{i,j,k}].y_c, r_coords[{i,j,k}].z_c )
      end
    end
  end

  c.fclose(file_handle) 
  c.free(filename)

  -- generate_hdf5_file_coords(filename, Nx, Ny, Nz)

  -- var grid = ispace(int3d, {x = Nx, y = Ny, z = Nz})
  -- var tmp_coords = region(grid, coordinates)

  -- attach(hdf5, tmp_coords.{x_c, y_c, z_c}, filename, regentlib.file_read_write)
  -- copy(r_coords.{x_c, y_c, z_c}, tmp_coords.{x_c, y_c, z_c})
  -- detach(hdf5, tmp_coords.{x_c, y_c, z_c})

  -- c.free(filename)
  -- __delete(tmp_coords)

  return 1
end



task read_coords( r_coords   : region(ispace(int3d), coordinates),
                  fileprefix : rawstring,
                  n_ghosts   : int64,
                  pencil     : int2d )
where
  writes( r_coords )
do
  var Nx = r_coords.ispace.bounds.hi.x - r_coords.ispace.bounds.lo.x + 1
  var Ny = r_coords.ispace.bounds.hi.y - r_coords.ispace.bounds.lo.y + 1
  var Nz = r_coords.ispace.bounds.hi.z - r_coords.ispace.bounds.lo.z + 1

  var vizstr : &int8 = [&int8] ( c.malloc(24*8) )
  c.sprintf(vizstr, "coords_px%04d_pz%04d.dat", pencil.x - 1, pencil.y - 1)
  var filename : &int8 = [&int8] ( c.malloc(256*8) )
  st.strcpy(filename, fileprefix)
  st.strcat(filename, vizstr)
  c.free(vizstr)

  var file_handle : &c.FILE
  if file_exists( filename ) then
    file_handle = c.fopen(filename, 'r')
  else
    c.printf("Error: Failed to open \"%s\"\n", filename)
    c.abort()
  end

  var bounds_c = r_coords.ispace.bounds

  var ix : int64[1]
  var iy : int64[1]
  var iz : int64[1]
  ix[0], iy[0], iz[0] = 0, 0, 0

  -- c.fprintf(file_handle, "%d %d %d\n", bounds_c.lo.x, bounds_c.lo.y, bounds_c.lo.z )
  read_grid(ix, iy, iz, file_handle)
  regentlib.assert( (ix[0] == bounds_c.lo.x - n_ghosts), "Grid size in restart file does not match partition size in x" )
  regentlib.assert( (iy[0] == bounds_c.lo.y - n_ghosts), "Grid size in restart file does not match partition size in y" )
  regentlib.assert( (iz[0] == bounds_c.lo.z - n_ghosts), "Grid size in restart file does not match partition size in z" )

  -- c.fprintf(file_handle, "%d %d %d\n", bounds_c.hi.x, bounds_c.hi.y, bounds_c.hi.z )
  read_grid(ix, iy, iz, file_handle)
  regentlib.assert( (ix[0] == bounds_c.hi.x - n_ghosts), "Grid size in restart file does not match partition size in x" )
  regentlib.assert( (iy[0] == bounds_c.hi.y - n_ghosts), "Grid size in restart file does not match partition size in y" )
  regentlib.assert( (iz[0] == bounds_c.hi.z - n_ghosts), "Grid size in restart file does not match partition size in z" )

  var x_dat : double[1]
  var y_dat : double[1]
  var z_dat : double[1]
  x_dat[0], y_dat[0], z_dat[0] = 0., 0., 0.

  for k = bounds_c.lo.z, bounds_c.hi.z+1 do
    for j = bounds_c.lo.y, bounds_c.hi.y+1 do
      for i = bounds_c.lo.x, bounds_c.hi.x+1 do
        -- c.fprintf(file_handle, "%26.16e %26.16e %26.16e\n", r_coords[{i,j,k}].x_c, r_coords[{i,j,k}].y_c, r_coords[{i,j,k}].z_c )
        read_coord_data(x_dat, y_dat, z_dat, file_handle)
        r_coords[{i,j,k}].x_c = x_dat[0]
        r_coords[{i,j,k}].y_c = y_dat[0]
        r_coords[{i,j,k}].z_c = z_dat[0]
      end
    end
  end

  c.fclose(file_handle) 
  c.free(filename)

  return 1
end

return true
