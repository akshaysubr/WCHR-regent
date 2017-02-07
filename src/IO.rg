import "regent"

require("fields")

local c = terralib.includec("assert.h")

hdf5 = terralib.includec("hdf5.h", {"-I", "/opt/hdf5-GNU-openmpi/include", "-I", "/usr/lib/openmpi/include" })
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
  --c.assert(fid > 0)

  var dims : hdf5.hsize_t[3]
  dims[0] = NZ
  dims[1] = NY
  dims[2] = NX
  var did = hdf5.H5Screate_simple(3, dims, [&uint64](0))
  --c.assert(did > 0)

  var ds1id = hdf5.H5Dcreate2(fid, "rho", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  --c.assert(ds1id > 0)
  hdf5.H5Dclose(ds1id)

  var ds2id = hdf5.H5Dcreate2(fid, "u", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  --c.assert(ds2id > 0)
  hdf5.H5Dclose(ds2id)

  var ds3id = hdf5.H5Dcreate2(fid, "v", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  --c.assert(ds3id > 0)
  hdf5.H5Dclose(ds3id)

  var ds4id = hdf5.H5Dcreate2(fid, "w", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  --c.assert(ds4id > 0)
  hdf5.H5Dclose(ds4id)

  var ds5id = hdf5.H5Dcreate2(fid, "p", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  --c.assert(ds5id > 0)
  hdf5.H5Dclose(ds5id)

  hdf5.H5Sclose(did)
  hdf5.H5Fclose(fid)
end

terra generate_hdf5_file_coords(filename : rawstring, NX : int64, NY : int64, NZ : int64)
  var fid = hdf5.H5Fcreate(filename, hdf5.H5F_ACC_TRUNC, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  --c.assert(fid > 0)

  var dims : hdf5.hsize_t[3]
  dims[0] = NZ
  dims[1] = NY
  dims[2] = NX
  var did = hdf5.H5Screate_simple(3, dims, [&uint64](0))
  --c.assert(did > 0)

  var ds1id = hdf5.H5Dcreate2(fid, "x_c", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  --c.assert(ds1id > 0)
  hdf5.H5Dclose(ds1id)

  var ds2id = hdf5.H5Dcreate2(fid, "y_c", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  --c.assert(ds2id > 0)
  hdf5.H5Dclose(ds2id)

  var ds3id = hdf5.H5Dcreate2(fid, "z_c", hdf5.H5T_IEEE_F64LE, did,
                              hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT, hdf5.H5P_DEFAULT)
  --c.assert(ds3id > 0)
  hdf5.H5Dclose(ds3id)

  hdf5.H5Sclose(did)
  hdf5.H5Fclose(fid)
end
