import "regent"

local use_io = os.getenv('USE_IO') == '1' or false
if not use_io then
  task write_primitive( r_prim     : region(ispace(int3d), primitive),
                        fileprefix : rawstring,
                        vizcount   : int32,
                        pencil     : int2d )
  where
    reads( r_prim )
  do
    return 1
  end

  task write_coords( r_coords   : region(ispace(int3d), coordinates),
                     fileprefix : rawstring,
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

task write_primitive( r_prim     : region(ispace(int3d), primitive),
                      fileprefix : rawstring,
                      vizcount   : int32,
                      pencil     : int2d )
where
  reads( r_prim )
do
  var Nx = r_prim.ispace.bounds.hi.x - r_prim.ispace.bounds.lo.x + 1
  var Ny = r_prim.ispace.bounds.hi.y - r_prim.ispace.bounds.lo.y + 1
  var Nz = r_prim.ispace.bounds.hi.z - r_prim.ispace.bounds.lo.z + 1

  var vizstr : &int8 = [&int8] ( c.malloc(21*8) )
  c.sprintf(vizstr, "%04d_px%04d_pz%04d.h5", vizcount, pencil.x, pencil.y)
  var filename : &int8 = [&int8] ( c.malloc(256*8) )
  st.strcpy(filename, fileprefix)
  st.strcat(filename, vizstr)
  c.free(vizstr)

  c.printf("\nWriting visualization file %s\n", filename)

  generate_hdf5_file(filename, Nx, Ny, Nz)

  var grid = ispace(int3d, {x = Nx, y = Ny, z = Nz})
  var tmp_prim = region(grid, primitive)

  attach(hdf5, tmp_prim.{rho, u, v, w, p}, filename, regentlib.file_read_write)
  copy(r_prim.{rho, u, v, w, p}, tmp_prim.{rho, u, v, w, p})
  detach(hdf5, tmp_prim.{rho, u, v, w, p})

  c.free(filename)
  __delete(tmp_prim)
  return 1
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
                   pencil     : int2d )
where
  reads( r_coords )
do
  var Nx = r_coords.ispace.bounds.hi.x - r_coords.ispace.bounds.lo.x + 1
  var Ny = r_coords.ispace.bounds.hi.y - r_coords.ispace.bounds.lo.y + 1
  var Nz = r_coords.ispace.bounds.hi.z - r_coords.ispace.bounds.lo.z + 1

  var vizstr : &int8 = [&int8] ( c.malloc(23*8) )
  c.sprintf(vizstr, "coords_px%04d_pz%04d.h5", pencil.x, pencil.y)
  var filename : &int8 = [&int8] ( c.malloc(256*8) )
  st.strcpy(filename, fileprefix)
  st.strcat(filename, vizstr)
  c.free(vizstr)

  generate_hdf5_file_coords(filename, Nx, Ny, Nz)

  var grid = ispace(int3d, {x = Nx, y = Ny, z = Nz})
  var tmp_coords = region(grid, coordinates)

  attach(hdf5, tmp_coords.{x_c, y_c, z_c}, filename, regentlib.file_read_write)
  copy(r_coords.{x_c, y_c, z_c}, tmp_coords.{x_c, y_c, z_c})
  detach(hdf5, tmp_coords.{x_c, y_c, z_c})

  c.free(filename)
  __delete(tmp_coords)

  return 1
end

return true
