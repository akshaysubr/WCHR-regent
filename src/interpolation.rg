import "regent"

require("fields")
local superlu = require("superlu_util")

__demand(__inline)
task get_char_values( r_prim_c : region(ispace(int3d), primitive),
                      idx : int3d,
                      Nx : int64,
                      Ny : int64,
                      Nz : int64)
where
  reads(r_prim_c)
do
  var char_values : double[5][6]

  for i = -3,2 do
    char_values[0][i+3] = r_prim_c[ [poff(idx, i, 0, 0, Nx, Ny, Nz)] ].rho
    char_values[1][i+3] = r_prim_c[ [poff(idx, i, 0, 0, Nx, Ny, Nz)] ].u
    char_values[2][i+3] = r_prim_c[ [poff(idx, i, 0, 0, Nx, Ny, Nz)] ].v
    char_values[3][i+3] = r_prim_c[ [poff(idx, i, 0, 0, Nx, Ny, Nz)] ].w
    char_values[4][i+3] = r_prim_c[ [poff(idx, i, 0, 0, Nx, Ny, Nz)] ].p
  end
end

__demand(__inline)
task get_beta( values : double[5][6],
               eq     : int32 )
do
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
end

__demand(__inline)
task get_nonlinear_weights( values : double[5][6], eq : int32 )
do
  var nlweights : double[5][4] = 0.0

  for eq = 0, 4 do 
    var beta = get_beta(values, eq)

    var epsilon : double = 1.0e-10
    
    -- Compute the nonlinear weights
    var d = array([1./4., 1./2., 1./4., 0.])
    var sum : double = 0.0
    for i = 0, 3 do
      nlweights[eq][i] = d[i] / ( (beta[i] + epsilon)*(beta[i] + epsilon) )
      sum += nlweights[eq][i]
    end

    for i = 0, 3 do
      nlweights[eq][i] = nlweights[eq][i] / sum
    end
  end
  return nlweights
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

  __demand(__vectorize)
  for idx in r_prim_l_x do
    var char_values : double[5][6] = get_char_values(r_prim_c, idx)

    var nlweights = get_nonlinear_weights(char_values[eq])
    var coeffs = get_coefficients(nlweights)

    -- Create matrix
    -- Create RHS
  end
end
