import "regent"

local c = regentlib.c

require("fields")
require("IO")
local superlu = require("superlu_util")

__demand(__inline)
task get_char_values( r_prim_c : region(ispace(int3d), primitive),
                      idx      : int3d,
                      Nx       : int64,
                      Ny       : int64,
                      Nz       : int64)
where
  reads(r_prim_c)
do
  var char_values : double[6][5]

  for i = -3, 3 do
    char_values[0][i+3] = r_prim_c[ [poff(idx, i, 0, 0, Nx, Ny, Nz)] ].rho
    char_values[1][i+3] = r_prim_c[ [poff(idx, i, 0, 0, Nx, Ny, Nz)] ].u
    char_values[2][i+3] = r_prim_c[ [poff(idx, i, 0, 0, Nx, Ny, Nz)] ].v
    char_values[3][i+3] = r_prim_c[ [poff(idx, i, 0, 0, Nx, Ny, Nz)] ].w
    char_values[4][i+3] = r_prim_c[ [poff(idx, i, 0, 0, Nx, Ny, Nz)] ].p
  end

  return char_values
end

-- __demand(__inline)
-- task get_beta( values : double[5][6],
--                eq     : int32 )
terra get_beta( values : double[6][5],
                eq     : int32 )
  var beta : double[4]

  beta[0] = 1.0/3.0*(values[eq][0]*(4.0*values[eq][0] - 19.0*values[eq][1] + 11.0*values[eq][2])
            + values[eq][1]*(25.0*values[eq][1] - 31.0*values[eq][2]) + 10.0*values[eq][2]*values[eq][2])
  
  beta[1] = 1.0/3.0*(values[eq][1]*(4.0*values[eq][1] - 13.0*values[eq][2] + 5.0*values[eq][3])
            + 13.0*values[eq][2]*(values[eq][2] - values[eq][3]) + 4.0*values[eq][3]*values[eq][3])
  
  beta[2] = 1.0/3.0*(values[eq][2]*(10.0*values[eq][2] - 31.0*values[eq][3] + 11.0*values[eq][4])
            + values[eq][3]*(25.0*values[eq][3] - 19.0*values[eq][4]) + 4.0*values[eq][4]*values[eq][4])
  
  beta[3] = (values[eq][0]*((525910327.0/232243200.0)*values[eq][0] - (4562164630.0/232243200.0)*values[eq][1]
            + (7799501420.0/232243200.0)*values[eq][2] - (6610694540.0/232243200.0)*values[eq][3] + (2794296070.0/232243200.0)*values[eq][4]
            - (472758974.0/232243200.0)*values[eq][5]) + 5.0*values[eq][1]*((2146987907.0/232243200.0)*values[eq][1] 
            - (7722406988.0/232243200.0)*values[eq][2] + (6763559276.0/232243200.0)*values[eq][3] - (2926461814.0/232243200.0)*values[eq][4] 
            + (503766638.0/232243200.0)*values[eq][5]) + 20.0*values[eq][2]*((1833221603.0/232243200.0)*values[eq][2] 
            - (3358664662.0/232243200.0)*values[eq][3] + (1495974539.0/232243200.0)*values[eq][4] - (263126407.0/232243200.0)*values[eq][5]) 
            + 20.0*values[eq][3]*((1607794163.0/232243200.0)*values[eq][3] - (1486026707.0/232243200.0)*values[eq][4]
            + (268747951.0/232243200.0)*values[eq][5]) + 5.0*values[eq][4]*((1432381427.0/232243200.0)*values[eq][4] 
            - (536951582.0/232243200.0)*values[eq][5]) + (263126407.0/232243200.0)*values[eq][5]*values[eq][5])

  return beta
end

-- __demand(__inline)
-- task get_nonlinear_weights( values : double[5][6] )
terra get_nonlinear_weights( values : double[6][5] )
  var nlweights : double[4][5]

  for eq = 0, 5 do 
    var beta = get_beta(values, eq)

    var epsilon : double = 1.0e-10
    
    -- Compute the nonlinear weights
    var d = array(1.0/4.0, 1.0/2.0, 1.0/4.0, 0.0)
    var d_central = array(1.0/8.0, 3.0/8.0, 3.0/8.0, 1.0/8.0)
    var sum : double = 0.0
    for i = 0, 4 do
      -- nlweights[eq][i] = d[i] / ( (beta[i] + epsilon)*(beta[i] + epsilon) )  -- Hard coded q = 2 here to avoid using power
      nlweights[eq][i] = d_central[i]
      sum = sum + nlweights[eq][i]
    end

    c.printf("get_nonlinear_weights: \n")
    for i = 0, 4 do
      nlweights[eq][i] = nlweights[eq][i] / sum
      c.printf("%g ", nlweights[eq][i])
    end
    c.printf("\n")
  end

  return nlweights
end

-- __demand(__inline)
-- task get_coefficients( nlweights : double[5][4] )
terra get_coefficients( nlweights : double[4][5] )

  var lcoeff0 = array(3.0/4.0, 1.0/4.0, 0.0/4.0, 1.0/4.0, 3.0/4.0, 0.0/4.0, 0.0/4.0)
  var lcoeff1 = array(1.0/4.0, 3.0/4.0, 0.0/4.0, 0.0/4.0, 3.0/4.0, 1.0/4.0, 0.0/4.0)
  var lcoeff2 = array(0.0/4.0, 3.0/4.0, 1.0/4.0, 0.0/4.0, 1.0/4.0, 3.0/4.0, 0.0/4.0)
  var lcoeff3 = array(0.0/4.0, 1.0/4.0, 3.0/4.0, 0.0/4.0, 0.0/4.0, 3.0/4.0, 1.0/4.0)

  var coeffs : double[7][5]

  for eq = 0, 5 do
    c.printf("get_coefficients: \n")
    for i = 0, 4 do
      c.printf("%g ", nlweights[eq][i])
    end
    c.printf("\n")
  end
  c.printf("\n")

  for eq = 0, 5 do
    c.printf("get_coefficients: \n")
    for i = 0, 7 do
      coeffs[eq][i] = lcoeff0[i]*nlweights[eq][0] + lcoeff1[i]*nlweights[eq][1] + lcoeff2[i]*nlweights[eq][2] + lcoeff3[i]*nlweights[eq][3]
      c.printf("%g ", coeffs[eq][i])
    end
    c.printf("\n")
  end

  return coeffs
end

task WCHR_interpolation_x( r_prim_c   : region(ispace(int3d), primitive),
                           r_prim_l_x : region(ispace(int3d), primitive),
                           r_prim_r_x : region(ispace(int3d), primitive),
                           matrix_l_x : region(ispace(int2d), superlu.CSR_matrix),
                           matrix_r_x : region(ispace(int2d), superlu.CSR_matrix),
                           slu_x      : region(ispace(int2d), superlu.c.superlu_vars_t) )
where
  reads(r_prim_c, matrix_l_x, matrix_r_x, slu_x), writes(r_prim_l_x, r_prim_r_x, matrix_l_x, matrix_r_x, slu_x)
do

  var nx = r_prim_c.ispace.bounds.hi.x - r_prim_c.ispace.bounds.lo.x + 1
  var ny = r_prim_c.ispace.bounds.hi.y - r_prim_c.ispace.bounds.lo.y + 1
  var nz = r_prim_c.ispace.bounds.hi.z - r_prim_c.ispace.bounds.lo.z + 1

  var r_rhs = region(ispace(int3d, {x = nx+1, y = ny, z = nz}), primitive)  

  var xdim : int64 = nx+1

  var counter : int64 = 0

  var pr = matrix_l_x.ispace.bounds.hi.x
  var pc = matrix_l_x.ispace.bounds.hi.y
  matrix_l_x[{pr,pc}].rowptr[0] = counter

  var bounds_c = r_prim_c.ispace.bounds
  var bounds_x = r_prim_l_x.ispace.bounds

  regentlib.assert(bounds_c.lo.x == 0, "Can only perform X interpolation in the X pencil")
  regentlib.assert(bounds_x.lo.x == 0, "Can only perform X interpolation in the X pencil")

  for iz = 0, nz do
    for iy = 0, ny do
      for row = 0, nx do
        var eq : int32 = 0  -- Interpolation only for rho right now

        var char_values : double[6][5] = get_char_values(r_prim_c, int3d({x = row, y = iy + bounds_c.lo.y, z = iz + bounds_c.lo.z}), nx, ny, nz)
        var nlweights = get_nonlinear_weights(char_values)
        var coeffs = get_coefficients(nlweights)

        c.printf("WCHR_interpolation_x: \n")
        for ci = 0, 7 do
          c.printf("%g ", coeffs[eq][ci])
        end
        c.printf("\n\n")

        -- coeffs[eq][0], coeffs[eq][1], coeffs[eq][2] = 3.0/16.0, 5.0/8.0, 3.0/16.0
        for j = 0, 3 do
          var col : int = row + j - 1
          var gcol : int64 = (col + nx)%nx + iy*xdim + iz*xdim*ny
          matrix_l_x[{pr,pc}].colind[counter] = gcol
          matrix_l_x[{pr,pc}].nzval [counter] = coeffs[eq][j]
          counter = counter + 1
        end
        var grow : int64 = row + iy*xdim + iz*xdim*ny
        matrix_l_x[{pr,pc}].rowptr[grow+1] = matrix_l_x[{pr,pc}].rowptr[grow] + 3

        r_rhs[{row,iy,iz}].rho = coeffs[eq][3] * char_values[eq][1]
                               + coeffs[eq][4] * char_values[eq][2]
                               + coeffs[eq][5] * char_values[eq][3]
                               + coeffs[eq][6] * char_values[eq][4]
        -- r_rhs[{row,iy,iz}].rho =  1.0/32.0 * char_values[eq][1]
        --                        + 15.0/32.0 * char_values[eq][2]
        --                        + 15.0/32.0 * char_values[eq][3]
        --                        +  1.0/32.0 * char_values[eq][4]
      end
      -- For the last point (enforcing periodicity)
      var gcol : int64 = nx + iy*xdim + iz*xdim*ny
      matrix_l_x[{pr,pc}].colind[counter] = gcol
      matrix_l_x[{pr,pc}].nzval [counter] = 1.0
      counter = counter + 1

      gcol = 0 + iy*xdim + iz*xdim*ny
      matrix_l_x[{pr,pc}].colind[counter] = gcol
      matrix_l_x[{pr,pc}].nzval [counter] = -1.0
      counter = counter + 1

      var grow : int64 = nx + iy*xdim + iz*xdim*ny
      matrix_l_x[{pr,pc}].rowptr[grow+1] = matrix_l_x[{pr,pc}].rowptr[grow] + 2

      r_rhs[{nx,iy,iz}].rho = 0.0
    end
  end

  -- write_primitive(r_rhs, "rhs_l_x", 0)

  superlu.MatrixSolve(__physical(r_rhs.rho), __fields(r_rhs.rho),
                      __physical(r_prim_l_x.rho), __fields(r_prim_l_x.rho), r_prim_l_x.bounds,
                      matrix_l_x[{pr,pc}], nx, ny, nz,
                      __physical(slu_x)[0], __fields(slu_x)[0], slu_x.bounds )

  __delete(r_rhs)
end
