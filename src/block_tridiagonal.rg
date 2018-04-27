import "regent"

local c       = regentlib.c
local cstdlib = terralib.includec("stdlib.h")

require("fields")
require("SOE")
local problem = require("problem")

local periodic_x = problem.periodic_x
local periodic_y = problem.periodic_y
local periodic_z = problem.periodic_z

-- All matrices here are in Fortran order

-- terra get_Rinv( rho : double, sos : double, Rinv : &double )
--   Rinv[0 + 3*0] = 0.; Rinv[0 + 3*1] = -0.5*rho*sos; Rinv[0 + 3*2] = 0.5;
--   Rinv[1 + 3*0] = 1.; Rinv[1 + 3*1] = 0.;           Rinv[1 + 3*2] = -1./(sos*sos);
--   Rinv[2 + 3*0] = 0.; Rinv[2 + 3*1] =  0.5*rho*sos; Rinv[2 + 3*2] = 0.5;
-- end

local terra multiply_diagonal_l( matrix : &double, d0 : double, d1 : double, d2 : double )
  matrix[0 + 3*0] = d0*matrix[0 + 3*0]; matrix[0 + 3*1] = d0*matrix[0 + 3*1]; matrix[0 + 3*2] = d0*matrix[0 + 3*2];
  matrix[1 + 3*0] = d1*matrix[1 + 3*0]; matrix[1 + 3*1] = d1*matrix[1 + 3*1]; matrix[1 + 3*2] = d1*matrix[1 + 3*2];
  matrix[2 + 3*0] = d2*matrix[2 + 3*0]; matrix[2 + 3*1] = d2*matrix[2 + 3*1]; matrix[2 + 3*2] = d2*matrix[2 + 3*2];
end
multiply_diagonal_l:setinlined(true)


local __demand(__inline) task multiply_diagonal_l_r( r : region(ispace(int3d), double[9]), i : int3d(double[9], r), d0 : double, d1 : double, d2 : double )
where
  reads writes (r)
do
  (@i)[0 + 3*0] = d0*((@i)[0 + 3*0]); (@i)[0 + 3*1] = d0*((@i)[0 + 3*1]); (@i)[0 + 3*2] = d0*((@i)[0 + 3*2]);
  (@i)[1 + 3*0] = d1*((@i)[1 + 3*0]); (@i)[1 + 3*1] = d1*((@i)[1 + 3*1]); (@i)[1 + 3*2] = d1*((@i)[1 + 3*2]);
  (@i)[2 + 3*0] = d2*((@i)[2 + 3*0]); (@i)[2 + 3*1] = d2*((@i)[2 + 3*1]); (@i)[2 + 3*2] = d2*((@i)[2 + 3*2]);
end

local terra multiply_diagonal_r( matrix : &double, d0 : double, d1 : double, d2 : double )
  matrix[0 + 3*0] = d0*matrix[0 + 3*0]; matrix[0 + 3*1] = d1*matrix[0 + 3*1]; matrix[0 + 3*2] = d2*matrix[0 + 3*2];
  matrix[1 + 3*0] = d0*matrix[1 + 3*0]; matrix[1 + 3*1] = d1*matrix[1 + 3*1]; matrix[1 + 3*2] = d2*matrix[1 + 3*2];
  matrix[2 + 3*0] = d0*matrix[2 + 3*0]; matrix[2 + 3*1] = d1*matrix[2 + 3*1]; matrix[2 + 3*2] = d2*matrix[2 + 3*2];
end
multiply_diagonal_r:setinlined(true)

local terra mult_matrix_vector( matrix : &double, vector : &double )
  var output : double[3]

  output[0] = matrix[0 + 3*0]*vector[0] + matrix[0 + 3*1]*vector[1] + matrix[0 + 3*2]*vector[2]
  output[1] = matrix[1 + 3*0]*vector[0] + matrix[1 + 3*1]*vector[1] + matrix[1 + 3*2]*vector[2]
  output[2] = matrix[2 + 3*0]*vector[0] + matrix[2 + 3*1]*vector[1] + matrix[2 + 3*2]*vector[2]

  return output
end
mult_matrix_vector:setinlined(true)

local terra mult_matrix_matrix( matrix1 : &double, matrix2 : &double, output : &double )
  -- var output : double[9]

  for i = 0,3 do
    for j = 0,3 do
      output[i+3*j] = 0.
      for k = 0,3 do
        output[i+3*j] = output[i+3*j] + matrix1[i+3*k] * matrix2[k+3*j]
      end
    end
  end 

  -- return output
end
mult_matrix_matrix:setinlined(true)


local __demand(__inline) task mult_matrix_matrix_r( matrix1 : &double, matrix2 : &double, r : region(ispace(int3d), double[9]), ii : int3d(double[9], r) )
where
  reads writes (r)
do
  -- var output : double[9]

  for i = 0,3 do
    for j = 0,3 do
      (r[ii])[i+3*j] = 0.
      for k = 0,3 do
        (r[ii])[i+3*j] = (r[ii])[i+3*j] + matrix1[i+3*k] * matrix2[k+3*j]
      end
    end
  end 

  -- return output
end

local terra print_matrix( matrix : &double )

  for i = 0,3 do
    for j = 0,3 do
	  c.printf(" %11.8f ", i, j, matrix[i + 3*j])
    end
    c.printf("\n")
  end
  c.printf("\n")
end
print_matrix:setinlined(true)

local terra invert_matrix( matrix : &double )
  -- var ipiv  : int[3]
  -- var lwork : int[1]
  -- var work  : double[9]
  -- var info  : int[1]
  -- var N     : int[1]

  -- N[0] = 3
  -- lwork[0] = 9

  -- lapack.dgetrf_(N,N,matrix,N,ipiv,info);
  -- if info[0] ~= 0 then
  --   print_matrix( matrix )
  -- end
  -- regentlib.assert(info[0] == 0, "DGETRF did not work as expected. Check for errors.")

  -- lapack.dgetri_(N,matrix,N,ipiv,work,lwork,info);
  -- regentlib.assert(info[0] == 0, "DGETRI did not work as expected. Check for errors.")

  var inv_det = matrix[0 + 3*0]*(matrix[1 + 3*1]*matrix[2 + 3*2]-matrix[1 + 3*2]*matrix[2 + 3*1]) 
              - matrix[0 + 3*1]*(matrix[1 + 3*0]*matrix[2 + 3*2]-matrix[2 + 3*0]*matrix[1 + 3*2]) 
              + matrix[0 + 3*2]*(matrix[1 + 3*0]*matrix[2 + 3*1]-matrix[2 + 3*0]*matrix[1 + 3*1])
  inv_det = 1. / inv_det

  var inv : double[9]

  inv[0 + 3*0] = (matrix[1 + 3*1]*matrix[2 + 3*2]-matrix[1 + 3*2]*matrix[2 + 3*1])
  inv[1 + 3*0] = (matrix[1 + 3*2]*matrix[2 + 3*0]-matrix[1 + 3*0]*matrix[2 + 3*2])
  inv[2 + 3*0] = (matrix[1 + 3*0]*matrix[2 + 3*1]-matrix[1 + 3*1]*matrix[2 + 3*0])

  inv[0 + 3*1] = (matrix[0 + 3*2]*matrix[2 + 3*1]-matrix[0 + 3*1]*matrix[2 + 3*2])
  inv[1 + 3*1] = (matrix[0 + 3*0]*matrix[2 + 3*2]-matrix[0 + 3*2]*matrix[2 + 3*0])
  inv[2 + 3*1] = (matrix[0 + 3*1]*matrix[2 + 3*0]-matrix[0 + 3*0]*matrix[2 + 3*1])

  inv[0 + 3*2] = (matrix[0 + 3*1]*matrix[1 + 3*2]-matrix[0 + 3*2]*matrix[1 + 3*1])
  inv[1 + 3*2] = (matrix[0 + 3*2]*matrix[1 + 3*0]-matrix[0 + 3*0]*matrix[1 + 3*2])
  inv[2 + 3*2] = (matrix[0 + 3*0]*matrix[1 + 3*1]-matrix[0 + 3*1]*matrix[1 + 3*0])

  for ii = 0,9 do
    matrix[ii] = inv[ii] * inv_det
  end

end
invert_matrix:setinlined(true)


local __demand(__inline) task invert_matrix_r( r : region(ispace(int3d), double[9]), i : int3d(double[9], r) )
where
  reads writes(r)
do

  var inv_det = ((@i)[0 + 3*0]) * ( ((@i)[1 + 3*1]) * ((@i)[2 + 3*2]) - ((@i)[1 + 3*2]) * ((@i)[2 + 3*1]) ) 
              - ((@i)[0 + 3*1]) * ( ((@i)[1 + 3*0]) * ((@i)[2 + 3*2]) - ((@i)[2 + 3*0]) * ((@i)[1 + 3*2]) ) 
              + ((@i)[0 + 3*2]) * ( ((@i)[1 + 3*0]) * ((@i)[2 + 3*1]) - ((@i)[2 + 3*0]) * ((@i)[1 + 3*1]) )
  inv_det = 1. / inv_det

  var inv : double[9]

  inv[0 + 3*0] = ((@i)[1 + 3*1]*(@i)[2 + 3*2]-(@i)[1 + 3*2]*(@i)[2 + 3*1])
  inv[1 + 3*0] = ((@i)[1 + 3*2]*(@i)[2 + 3*0]-(@i)[1 + 3*0]*(@i)[2 + 3*2])
  inv[2 + 3*0] = ((@i)[1 + 3*0]*(@i)[2 + 3*1]-(@i)[1 + 3*1]*(@i)[2 + 3*0])

  inv[0 + 3*1] = ((@i)[0 + 3*2]*(@i)[2 + 3*1]-(@i)[0 + 3*1]*(@i)[2 + 3*2])
  inv[1 + 3*1] = ((@i)[0 + 3*0]*(@i)[2 + 3*2]-(@i)[0 + 3*2]*(@i)[2 + 3*0])
  inv[2 + 3*1] = ((@i)[0 + 3*1]*(@i)[2 + 3*0]-(@i)[0 + 3*0]*(@i)[2 + 3*1])

  inv[0 + 3*2] = ((@i)[0 + 3*1]*(@i)[1 + 3*2]-(@i)[0 + 3*2]*(@i)[1 + 3*1])
  inv[1 + 3*2] = ((@i)[0 + 3*2]*(@i)[1 + 3*0]-(@i)[0 + 3*0]*(@i)[1 + 3*2])
  inv[2 + 3*2] = ((@i)[0 + 3*0]*(@i)[1 + 3*1]-(@i)[0 + 3*1]*(@i)[1 + 3*0])

  for ii = 0,9 do
    (@i)[ii] = inv[ii] * inv_det
  end

end

local terra axpby( x : &double, y : &double, a : double, b : double, N : int )
  for i = 0,N do
    x[i] = a*x[i] + b*y[i]
  end
end
axpby:setinlined(true)

local __demand(__inline) task axpby_r( r : region(ispace(int3d), double[9]), ii : int3d(double[9], r), y : &double, a : double, b : double, N : int )
where
  reads writes (r)
do
  for i = 0,N do
    (@ii)[i] = a*(@ii)[i] + b*y[i]
  end
end

local terra print_vector( vector : &double )

  for i = 0,3 do
    c.printf(" %11.8f ", vector[i])
    c.printf("\n")
  end
end

local terra random_number()
  return ( [double](cstdlib.rand()) / [double](cstdlib.RAND_MAX + 1.) )
end

__demand(__inline)
task solve_block_tridiagonal_x( alpha   : region( ispace(int3d), coeffs ),
                                beta    : region( ispace(int3d), coeffs ),
                                gamma   : region( ispace(int3d), coeffs ),
                                rho_avg : region( ispace(int3d), double ),
                                sos_avg : region( ispace(int3d), double ),
                                sol     : region( ispace(int3d), primitive ),
                                d       : region( ispace(int3d), double[9] ),
                                Uinv    : region( ispace(int3d), double[9] ) )
                                -- d       : region( ispace(int3d), &double ),
                                -- Uinv    : region( ispace(int3d), &double ) )
where
  reads(alpha.{_0,_1,_4}, beta.{_0,_1,_4}, gamma.{_0,_1,_4}, rho_avg, sos_avg, sol.{rho,u,p}, d, Uinv),
  writes(sol.{rho,u,p}, d, Uinv, beta.{_0,_1,_4})
do

  var bounds = sol.ispace.bounds
  var N : int = bounds.hi.x + 1
  if periodic_x then
    N = bounds.hi.x -- Don't solve for last edge if periodic
  end

  for k = bounds.lo.z, bounds.hi.z+1 do
    for j = bounds.lo.y, bounds.hi.y+1 do

      if periodic_x then
        -- If periodic, make correction for Sherman-Morrison
        beta[{0,j,k}]._0 = beta[{0,j,k}]._0 + alpha[{0,j,k}]._0
        beta[{0,j,k}]._1 = beta[{0,j,k}]._1 + alpha[{0,j,k}]._1
        beta[{0,j,k}]._4 = beta[{0,j,k}]._4 + alpha[{0,j,k}]._4

        beta[{N-1,j,k}]._0 = beta[{N-1,j,k}]._0 + gamma[{N-1,j,k}]._0
        beta[{N-1,j,k}]._1 = beta[{N-1,j,k}]._1 + gamma[{N-1,j,k}]._1
        beta[{N-1,j,k}]._4 = beta[{N-1,j,k}]._4 + gamma[{N-1,j,k}]._4
      end

      -- Forward elimination
      get_Rinv_r( rho_avg[{0,j,k}], sos_avg[{0,j,k}], d, unsafe_cast(int3d(double[9], d), int3d {0,j,k}) )
      multiply_diagonal_l_r( d, unsafe_cast(int3d(double[9], d), int3d {0,j,k}), beta[{0,j,k}]._0, beta[{0,j,k}]._1, beta[{0,j,k}]._4 )
      invert_matrix_r( d, unsafe_cast(int3d(double[9], d), int3d {0,j,k}) )

      sol[{0,j,k}].rho = -sol[{0,j,k}].rho
      sol[{0,j,k}].u   = -sol[{0,j,k}].u
      sol[{0,j,k}].p   = -sol[{0,j,k}].p

      if periodic_x then
        get_Rinv_r( rho_avg[{0,j,k}], sos_avg[{0,j,k}], Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {0,j,k}) )
        multiply_diagonal_l_r( Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {0,j,k}), alpha[{0,j,k}]._0, alpha[{0,j,k}]._1, alpha[{0,j,k}]._4 )
      end

      for i = 1,N do
        var Rinv_i : double[9]
        get_Rinv( rho_avg[{i,j,k}], sos_avg[{i,j,k}], Rinv_i ) -- Get Rinv at i

        var mat : double[9]
        mult_matrix_matrix( Rinv_i, d[{i-1,j,k}], mat )
        multiply_diagonal_l( mat, alpha[{i,j,k}]._0, alpha[{i,j,k}]._1, alpha[{i,j,k}]._4 )

        var gammaRinv_im1 : double[9]
        get_Rinv( rho_avg[{i-1,j,k}], sos_avg[{i-1,j,k}], gammaRinv_im1 ) -- Get Rinv at i-1
        multiply_diagonal_l( gammaRinv_im1, gamma[{i-1,j,k}]._0, gamma[{i-1,j,k}]._1, gamma[{i-1,j,k}]._4 )

        multiply_diagonal_l( Rinv_i, beta[{i,j,k}]._0, beta[{i,j,k}]._1, beta[{i,j,k}]._4 )

        mult_matrix_matrix_r( mat, gammaRinv_im1, d, unsafe_cast(int3d(double[9], d), int3d {i,j,k}) )
        axpby_r( d, unsafe_cast(int3d(double[9], d), int3d {i,j,k}), Rinv_i, -1., 1., 9 ) -- Delta_i = beta_i - alpha_i Delta_i gamma_i-1
        invert_matrix_r( d, unsafe_cast(int3d(double[9], d), int3d {i,j,k}) )

        var prim : double[3] = array( sol[{i-1,j,k}].rho, sol[{i-1,j,k}].u, sol[{i-1,j,k}].p )
        var cprime = mult_matrix_vector( mat, prim )
        sol[{i,j,k}].rho = -sol[{i,j,k}].rho - cprime[0]
        sol[{i,j,k}].u   = -sol[{i,j,k}].u   - cprime[1]
        sol[{i,j,k}].p   = -sol[{i,j,k}].p   - cprime[2]

        if periodic_x then
          mult_matrix_matrix_r( mat, Uinv[{i-1,j,k}], Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,j,k}) )
          multiply_diagonal_l_r( Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,j,k}), -1., -1., -1. )
        end
      end

      if periodic_x then
        var Rinv : double[9]
        get_Rinv( rho_avg[{N-1,j,k}], sos_avg[{N-1,j,k}], Rinv ) -- Get Rinv
        multiply_diagonal_l( Rinv, gamma[{N-1,j,k}]._0, gamma[{N-1,j,k}]._1, gamma[{N-1,j,k}]._4 )
        axpby_r( Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {N-1,j,k}), Rinv, 1., -1., 9 )
      end

      -- Back substitution
      var prim : double[3] = array( sol[{N-1,j,k}].rho, sol[{N-1,j,k}].u, sol[{N-1,j,k}].p )
      var cprime = mult_matrix_vector( d[{N-1,j,k}], prim )
      sol[{N-1,j,k}].rho = - cprime[0]
      sol[{N-1,j,k}].u   = - cprime[1]
      sol[{N-1,j,k}].p   = - cprime[2]

      if periodic_x then
        var tmp : double[9]
        mult_matrix_matrix( d[{N-1,j,k}], Uinv[{N-1,j,k}], tmp )
        for ii = 0,9 do
          (Uinv[{N-1,j,k}])[ii] = -tmp[ii]
        end
      end

      for i = N-1,0,-1 do
        var gammaRinv_im1 : double[9]
        get_Rinv( rho_avg[{i-1,j,k}], sos_avg[{i-1,j,k}], gammaRinv_im1 ) -- Get Rinv at i-1
        multiply_diagonal_l( gammaRinv_im1, gamma[{i-1,j,k}]._0, gamma[{i-1,j,k}]._1, gamma[{i-1,j,k}]._4 )

        prim = array( sol[{i,j,k}].rho, sol[{i,j,k}].u, sol[{i,j,k}].p )
        cprime = mult_matrix_vector( gammaRinv_im1, prim )    

        prim = array( sol[{i-1,j,k}].rho, sol[{i-1,j,k}].u, sol[{i-1,j,k}].p )
        axpby( prim, cprime, 1., 1., 3 )

        cprime = mult_matrix_vector( d[{i-1,j,k}], prim )
        sol[{i-1,j,k}].rho = - cprime[0]
        sol[{i-1,j,k}].u   = - cprime[1]
        sol[{i-1,j,k}].p   = - cprime[2]

        if periodic_x then
          var tmp : double[9]
          mult_matrix_matrix( gammaRinv_im1, Uinv[{i,j,k}], tmp )
          axpby( tmp, Uinv[{i-1,j,k}], -1., -1., 9 )
          mult_matrix_matrix_r( d[{i-1,j,k}], tmp, Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i-1,j,k}) )
        end -- periodic
      end

      -- Sherman-Morrison correction
      if periodic_x then
        var M : double[9] = array(1., 0., 0., 0., 1., 0., 0., 0., 1.) -- Identity matrix
        axpby( M, Uinv[{0,j,k}],   1., 1., 9 )
        axpby( M, Uinv[{N-1,j,k}], 1.,-1., 9 )
        invert_matrix( M )
        prim = array( sol[{0,j,k}].rho - sol[{N-1,j,k}].rho, sol[{0,j,k}].u - sol[{N-1,j,k}].u, sol[{0,j,k}].p - sol[{N-1,j,k}].p )
        var corrfact : double[3] = mult_matrix_vector( M, prim )

        for i = 0,N do
          cprime = mult_matrix_vector( Uinv[{i,j,k}], corrfact )
          sol[{i,j,k}].rho = sol[{i,j,k}].rho - cprime[0]
          sol[{i,j,k}].u   = sol[{i,j,k}].u   - cprime[1]
          sol[{i,j,k}].p   = sol[{i,j,k}].p   - cprime[2]
        end

        -- Copy last edge if periodic
        sol[{N,j,k}].rho = sol[{0,j,k}].rho
        sol[{N,j,k}].u   = sol[{0,j,k}].u
        sol[{N,j,k}].p   = sol[{0,j,k}].p
      end

    end
  end

end


function make_solve_tridiagonal_x( fi, fn )
  -- fi: field index in coeff region (_1, _2, _3 ...)
  -- fn: field name in primitive region (u, v, w ...)

  local solve_tridiagonal_x __demand(__inline) task solve_tridiagonal_x( alpha   : region( ispace(int3d), coeffs ),
                                                                         beta    : region( ispace(int3d), coeffs ),
                                                                         gamma   : region( ispace(int3d), coeffs ),
                                                                         sol     : region( ispace(int3d), primitive ),
                                                                         Uinv    : region( ispace(int3d), double) )
  where
    reads (alpha.[fi], gamma.[fi]), reads writes ( beta.[fi], sol.[fn], Uinv )
  do
  
    var bounds = sol.ispace.bounds
    var N : int = bounds.hi.x + 1
    if periodic_x then
      N = bounds.hi.x -- Don't solve for last edge if periodic
    end
  
    for k = bounds.lo.z, bounds.hi.z+1 do
      for j = bounds.lo.y, bounds.hi.y+1 do
  
        var beta0 : double = beta[{0,j,k}].[fi]
        if periodic_x then
          -- Sherman-Morrison reduced matrix
          Uinv[{0,j,k}] = -beta[{0,j,k}].[fi]
          beta[{N-1,j,k}].[fi] = beta[{N-1,j,k}].[fi] + gamma[{N-1,j,k}].[fi] * alpha[{0,j,k}].[fi] / beta[{0,j,k}].[fi]
          beta[{0,j,k}].[fi] = 2.*beta[{0,j,k}].[fi]
        end

        -- Forward substitution
        beta[{0,j,k}].[fi] = 1./beta[{0,j,k}].[fi]
        sol[{0,j,k}].[fn] = beta[{0,j,k}].[fi] * sol[{0,j,k}].[fn]
        if periodic_x then
          Uinv[{0,j,k}] =  beta[{0,j,k}].[fi] * Uinv[{0,j,k}]
        end

        for i = 1,N do
          beta[{i,j,k}].[fi] = 1./( beta[{i,j,k}].[fi] - alpha[{i,j,k}].[fi] * beta[{i-1,j,k}].[fi] * gamma[{i-1,j,k}].[fi] )
          sol[{i,j,k}].[fn] = beta[{i,j,k}].[fi] * (sol[{i,j,k}].[fn] - alpha[{i,j,k}].[fi] * sol[{i-1,j,k}].[fn])
          if periodic_x then
            Uinv[{i,j,k}] = beta[{i,j,k}].[fi] * ( - alpha[{i,j,k}].[fi] * Uinv[{i-1,j,k}])
          end

        end

        if periodic_x then
          Uinv[{N-1,j,k}] = Uinv[{N-1,j,k}] + beta[{N-1,j,k}].[fi] * ( gamma[{N-1,j,k}].[fi] )
        end
        
        -- Backward elimination
        for i = N-1,0,-1 do
          sol[{i-1,j,k}].[fn] = sol[{i-1,j,k}].[fn] - beta[{i-1,j,k}].[fi] * gamma[{i-1,j,k}].[fi] * sol[{i,j,k}].[fn]
          if periodic_x then
            Uinv[{i-1,j,k}] = Uinv[{i-1,j,k}] - beta[{i-1,j,k}].[fi] * gamma[{i-1,j,k}].[fi] * Uinv[{i,j,k}]
          end
        end
  
        -- Sherman-Morrison correction
        if periodic_x then
          var corr_factor : double = ( sol[{0,j,k}].[fn] - (alpha[{0,j,k}].[fi] / beta0) * sol[{N-1,j,k}].[fn] )
          corr_factor = corr_factor / ( 1. + Uinv[{0,j,k}] - (alpha[{0,j,k}].[fi] / beta0) * Uinv[{N-1,j,k}] )

          for i = 0,N do
            sol[{i,j,k}].[fn] = sol[{i,j,k}].[fn] - corr_factor * Uinv[{i,j,k}]
          end

          -- Copy last edge if periodic
          sol[{N,j,k}].[fn] = sol[{0,j,k}].[fn]
        end

      end
    end
  
  end
  return solve_tridiagonal_x
end

__demand(__inline)
task solve_block_tridiagonal_y( alpha   : region( ispace(int3d), coeffs ),
                                beta    : region( ispace(int3d), coeffs ),
                                gamma   : region( ispace(int3d), coeffs ),
                                rho_avg : region( ispace(int3d), double ),
                                sos_avg : region( ispace(int3d), double ),
                                sol     : region( ispace(int3d), primitive ),
                                d       : region( ispace(int3d), double[9] ),
                                Uinv    : region( ispace(int3d), double[9] ) )
                                -- d       : region( ispace(int3d), &double ),
                                -- Uinv    : region( ispace(int3d), &double ) )
where
  reads(alpha.{_0,_2,_4}, beta.{_0,_2,_4}, gamma.{_0,_2,_4}, rho_avg, sos_avg, sol.{rho,v,p}, d, Uinv),
  writes(sol.{rho,v,p}, d, Uinv, beta.{_0,_2,_4})
do

  var bounds = sol.ispace.bounds
  var N : int = bounds.hi.y + 1
  if periodic_y then
    N = bounds.hi.y -- Don't solve for last edge if periodic
  end

  for k = bounds.lo.z, bounds.hi.z+1 do

    if periodic_y then
      for i = bounds.lo.x, bounds.hi.x+1 do
        -- If periodic, make correction for Sherman-Morrison
        beta[{i,0,k}]._0 = beta[{i,0,k}]._0 + alpha[{i,0,k}]._0
        beta[{i,0,k}]._2 = beta[{i,0,k}]._2 + alpha[{i,0,k}]._2
        beta[{i,0,k}]._4 = beta[{i,0,k}]._4 + alpha[{i,0,k}]._4

        beta[{i,N-1,k}]._0 = beta[{i,N-1,k}]._0 + gamma[{i,N-1,k}]._0
        beta[{i,N-1,k}]._2 = beta[{i,N-1,k}]._2 + gamma[{i,N-1,k}]._2
        beta[{i,N-1,k}]._4 = beta[{i,N-1,k}]._4 + gamma[{i,N-1,k}]._4
      end
    end

    -- Forward elimination
    for i = bounds.lo.x, bounds.hi.x+1 do
      get_Rinv_r( rho_avg[{i,0,k}], sos_avg[{i,0,k}], d, unsafe_cast(int3d(double[9], d), int3d {i,0,k}) )
      multiply_diagonal_l_r( d, unsafe_cast(int3d(double[9], d), int3d {i,0,k}), beta[{i,0,k}]._0, beta[{i,0,k}]._2, beta[{i,0,k}]._4 )
      invert_matrix_r( d, unsafe_cast(int3d(double[9], d), int3d {i,0,k}) )
    end

    for i = bounds.lo.x, bounds.hi.x+1 do
      sol[{i,0,k}].rho = -sol[{i,0,k}].rho
      sol[{i,0,k}].v   = -sol[{i,0,k}].v
      sol[{i,0,k}].p   = -sol[{i,0,k}].p
    end

    if periodic_y then
      for i = bounds.lo.x, bounds.hi.x+1 do
        get_Rinv_r( rho_avg[{i,0,k}], sos_avg[{i,0,k}], Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,0,k}) )
        multiply_diagonal_l_r( Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,0,k}), alpha[{i,0,k}]._0, alpha[{i,0,k}]._2, alpha[{i,0,k}]._4 )
      end
    end

    for j = 1,N do
      for i = bounds.lo.x, bounds.hi.x+1 do
        var Rinv_i : double[9]
        get_Rinv( rho_avg[{i,j,k}], sos_avg[{i,j,k}], Rinv_i ) -- Get Rinv at i

        var mat : double[9]
        mult_matrix_matrix( Rinv_i, d[{i,j-1,k}], mat )
        multiply_diagonal_l( mat, alpha[{i,j,k}]._0, alpha[{i,j,k}]._2, alpha[{i,j,k}]._4 )

        var gammaRinv_im1 : double[9]
        get_Rinv( rho_avg[{i,j-1,k}], sos_avg[{i,j-1,k}], gammaRinv_im1 ) -- Get Rinv at i-1
        multiply_diagonal_l( gammaRinv_im1, gamma[{i,j-1,k}]._0, gamma[{i,j-1,k}]._2, gamma[{i,j-1,k}]._4 )

        multiply_diagonal_l( Rinv_i, beta[{i,j,k}]._0, beta[{i,j,k}]._2, beta[{i,j,k}]._4 )

        mult_matrix_matrix_r( mat, gammaRinv_im1, d, unsafe_cast(int3d(double[9], d), int3d {i,j,k}) )
        axpby_r( d, unsafe_cast(int3d(double[9], d), int3d {i,j,k}), Rinv_i, -1., 1., 9 ) -- Delta_i = beta_i - alpha_i Delta_i gamma_i-1
        invert_matrix_r( d, unsafe_cast(int3d(double[9], d), int3d {i,j,k}) )

        var prim : double[3] = array( sol[{i,j-1,k}].rho, sol[{i,j-1,k}].v, sol[{i,j-1,k}].p )
        var cprime = mult_matrix_vector( mat, prim )
        sol[{i,j,k}].rho = -sol[{i,j,k}].rho - cprime[0]
        sol[{i,j,k}].v   = -sol[{i,j,k}].v   - cprime[1]
        sol[{i,j,k}].p   = -sol[{i,j,k}].p   - cprime[2]

        if periodic_y then
          mult_matrix_matrix_r( mat, Uinv[{i,j-1,k}], Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,j,k}) )
          multiply_diagonal_l_r( Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,j,k}), -1., -1., -1. )
        end
      end
    end

    if periodic_y then
      for i = bounds.lo.x, bounds.hi.x+1 do
        var Rinv : double[9]
        get_Rinv( rho_avg[{i,N-1,k}], sos_avg[{i,N-1,k}], Rinv ) -- Get Rinv
        multiply_diagonal_l( Rinv, gamma[{i,N-1,k}]._0, gamma[{i,N-1,k}]._2, gamma[{i,N-1,k}]._4 )
        axpby_r( Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,N-1,k}), Rinv, 1., -1., 9 )
      end
    end

    for i = bounds.lo.x, bounds.hi.x+1 do
      -- Back substitution
      var prim : double[3] = array( sol[{i,N-1,k}].rho, sol[{i,N-1,k}].v, sol[{i,N-1,k}].p )
      var cprime = mult_matrix_vector( d[{i,N-1,k}], prim )
      sol[{i,N-1,k}].rho = - cprime[0]
      sol[{i,N-1,k}].v   = - cprime[1]
      sol[{i,N-1,k}].p   = - cprime[2]
    end

    if periodic_y then
      for i = bounds.lo.x, bounds.hi.x+1 do
        var tmp : double[9]
        mult_matrix_matrix( d[{i,N-1,k}], Uinv[{i,N-1,k}], tmp )
        for ii = 0,9 do
          (Uinv[{i,N-1,k}])[ii] = -tmp[ii]
        end
      end
    end

    for j = N-1,0,-1 do
      for i = bounds.lo.x, bounds.hi.x+1 do
        var gammaRinv_im1 : double[9]
        get_Rinv( rho_avg[{i,j-1,k}], sos_avg[{i,j-1,k}], gammaRinv_im1 ) -- Get Rinv at i-1
        multiply_diagonal_l( gammaRinv_im1, gamma[{i,j-1,k}]._0, gamma[{i,j-1,k}]._2, gamma[{i,j-1,k}]._4 )

        var prim : double[3] = array( sol[{i,j,k}].rho, sol[{i,j,k}].v, sol[{i,j,k}].p )
        var cprime = mult_matrix_vector( gammaRinv_im1, prim )    

        prim = array( sol[{i,j-1,k}].rho, sol[{i,j-1,k}].v, sol[{i,j-1,k}].p )
        axpby( prim, cprime, 1., 1., 3 )

        cprime = mult_matrix_vector( d[{i,j-1,k}], prim )
        sol[{i,j-1,k}].rho = - cprime[0]
        sol[{i,j-1,k}].v   = - cprime[1]
        sol[{i,j-1,k}].p   = - cprime[2]

        if periodic_y then
          var tmp : double[9]
          mult_matrix_matrix( gammaRinv_im1, Uinv[{i,j,k}], tmp )
          axpby( tmp, Uinv[{i,j-1,k}], -1., -1., 9 )
          mult_matrix_matrix_r( d[{i,j-1,k}], tmp, Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,j-1,k}) )
        end -- periodic
      end
    end

    -- Sherman-Morrison correction
    if periodic_y then
      for i = bounds.lo.x, bounds.hi.x+1 do
        var M : double[9] = array(1., 0., 0., 0., 1., 0., 0., 0., 1.) -- Identity matrix
        axpby( M, Uinv[{i,0,k}],   1., 1., 9 )
        axpby( M, Uinv[{i,N-1,k}], 1.,-1., 9 )
        invert_matrix( M )
        var prim : double[3] = array( sol[{i,0,k}].rho - sol[{i,N-1,k}].rho, sol[{i,0,k}].v - sol[{i,N-1,k}].v, sol[{i,0,k}].p - sol[{i,N-1,k}].p )
        var corrfact : double[3] = mult_matrix_vector( M, prim )

        for j = 0,N do
          var cprime = mult_matrix_vector( Uinv[{i,j,k}], corrfact )
          sol[{i,j,k}].rho = sol[{i,j,k}].rho - cprime[0]
          sol[{i,j,k}].v   = sol[{i,j,k}].v   - cprime[1]
          sol[{i,j,k}].p   = sol[{i,j,k}].p   - cprime[2]
        end

        -- Copy last edge if periodic
        sol[{i,N,k}].rho = sol[{i,0,k}].rho
        sol[{i,N,k}].v   = sol[{i,0,k}].v
        sol[{i,N,k}].p   = sol[{i,0,k}].p
      end
    end

  end

end

function make_solve_tridiagonal_y( fi, fn )
  -- fi: field index in coeff region (_1, _2, _3 ...)
  -- fn: field name in primitive region (u, v, w ...)

  local solve_tridiagonal_y __demand(__inline) task solve_tridiagonal_y( alpha   : region( ispace(int3d), coeffs ),
                                                                         beta    : region( ispace(int3d), coeffs ),
                                                                         gamma   : region( ispace(int3d), coeffs ),
                                                                         sol     : region( ispace(int3d), primitive ),
                                                                         Uinv    : region( ispace(int3d), double) )
  where
    reads (alpha.[fi], gamma.[fi]), reads writes ( beta.[fi], sol.[fn], Uinv )
  do
  
    var bounds = sol.ispace.bounds
    var N : int = bounds.hi.y + 1
    if periodic_y then
      N = bounds.hi.y -- Don't solve for last edge if periodic
    end
  
    for k = bounds.lo.z, bounds.hi.z+1 do
  
      for i = bounds.lo.x, bounds.hi.x+1 do
        var beta0 : double = beta[{i,0,k}].[fi]
        if periodic_y then
            -- Sherman-Morrison reduced matrix
            Uinv[{i,0,k}] = -beta[{i,0,k}].[fi]
            beta[{i,N-1,k}].[fi] = beta[{i,N-1,k}].[fi] + gamma[{i,N-1,k}].[fi] * alpha[{i,0,k}].[fi] / beta[{i,0,k}].[fi]
            beta[{i,0,k}].[fi] = 2.*beta[{i,0,k}].[fi]
        end

        -- Forward substitution
        beta[{i,0,k}].[fi] = 1./beta[{i,0,k}].[fi]
        sol[{i,0,k}].[fn] = beta[{i,0,k}].[fi] * sol[{i,0,k}].[fn]
        if periodic_y then
          Uinv[{i,0,k}] =  beta[{i,0,k}].[fi] * Uinv[{i,0,k}]
        end
      end

      for j = 1,N do
        for i = bounds.lo.x, bounds.hi.x+1 do
          beta[{i,j,k}].[fi] = 1./( beta[{i,j,k}].[fi] - alpha[{i,j,k}].[fi] * beta[{i,j-1,k}].[fi] * gamma[{i,j-1,k}].[fi] )
          sol[{i,j,k}].[fn] = beta[{i,j,k}].[fi] * (sol[{i,j,k}].[fn] - alpha[{i,j,k}].[fi] * sol[{i,j-1,k}].[fn])
          if periodic_y then
            Uinv[{i,j,k}] = beta[{i,j,k}].[fi] * ( - alpha[{i,j,k}].[fi] * Uinv[{i,j-1,k}])
          end
        end
      end

      if periodic_y then
        for i = bounds.lo.x, bounds.hi.x+1 do
          Uinv[{i,N-1,k}] = Uinv[{i,N-1,k}] + beta[{i,N-1,k}].[fi] * ( gamma[{i,N-1,k}].[fi] )
        end
      end
      
      -- Backward elimination
      for j = N-1,0,-1 do
        for i = bounds.lo.x, bounds.hi.x+1 do
          sol[{i,j-1,k}].[fn] = sol[{i,j-1,k}].[fn] - beta[{i,j-1,k}].[fi] * gamma[{i,j-1,k}].[fi] * sol[{i,j,k}].[fn]
          if periodic_y then
            Uinv[{i,j-1,k}] = Uinv[{i,j-1,k}] - beta[{i,j-1,k}].[fi] * gamma[{i,j-1,k}].[fi] * Uinv[{i,j,k}]
          end
        end
      end
  
      -- Sherman-Morrison correction
      if periodic_y then
        for i = bounds.lo.x, bounds.hi.x+1 do
          var beta0 = (1./beta[{i,0,k}].[fi])/2.
          var corr_factor : double = ( sol[{i,0,k}].[fn] - (alpha[{i,0,k}].[fi] / beta0) * sol[{i,N-1,k}].[fn] )
          corr_factor = corr_factor / ( 1. + Uinv[{i,0,k}] - (alpha[{i,0,k}].[fi] / beta0) * Uinv[{i,N-1,k}] )

          for j = 0,N do
            sol[{i,j,k}].[fn] = sol[{i,j,k}].[fn] - corr_factor * Uinv[{i,j,k}]
          end

          -- Copy last edge if periodic
          sol[{i,N,k}].[fn] = sol[{i,0,k}].[fn]
        end
      end

    end
  
  end
  return solve_tridiagonal_y
end


__demand(__inline)
task solve_block_tridiagonal_z( alpha   : region( ispace(int3d), coeffs ),
                                beta    : region( ispace(int3d), coeffs ),
                                gamma   : region( ispace(int3d), coeffs ),
                                rho_avg : region( ispace(int3d), double ),
                                sos_avg : region( ispace(int3d), double ),
                                sol     : region( ispace(int3d), primitive ),
                                d       : region( ispace(int3d), double[9] ),
                                Uinv    : region( ispace(int3d), double[9] ) )
                                -- d       : region( ispace(int3d), &double ),
                                -- Uinv    : region( ispace(int3d), &double ) )
where
  reads(alpha.{_0,_3,_4}, beta.{_0,_3,_4}, gamma.{_0,_3,_4}, rho_avg, sos_avg, sol.{rho,w,p}, d, Uinv),
  writes(sol.{rho,w,p}, d, Uinv, beta.{_0,_3,_4})
do

  var bounds = sol.ispace.bounds
  var N : int = bounds.hi.z + 1
  if periodic_z then
    N = bounds.hi.z -- Don't solve for last edge if periodic
  end


  if periodic_z then
    for j = bounds.lo.y, bounds.hi.y+1 do
      for i = bounds.lo.x, bounds.hi.x+1 do
        -- If periodic, make correction for Sherman-Morrison
        beta[{i,j,0}]._0 = beta[{i,j,0}]._0 + alpha[{i,j,0}]._0
        beta[{i,j,0}]._3 = beta[{i,j,0}]._3 + alpha[{i,j,0}]._3
        beta[{i,j,0}]._4 = beta[{i,j,0}]._4 + alpha[{i,j,0}]._4

        beta[{i,j,N-1}]._0 = beta[{i,j,N-1}]._0 + gamma[{i,j,N-1}]._0
        beta[{i,j,N-1}]._3 = beta[{i,j,N-1}]._3 + gamma[{i,j,N-1}]._3
        beta[{i,j,N-1}]._4 = beta[{i,j,N-1}]._4 + gamma[{i,j,N-1}]._4
      end
    end
  end

  -- Forward elimination
  for j = bounds.lo.y, bounds.hi.y+1 do
    for i = bounds.lo.x, bounds.hi.x+1 do
      get_Rinv_r( rho_avg[{i,j,0}], sos_avg[{i,j,0}], d, unsafe_cast(int3d(double[9], d), int3d {i,j,0}) )
      multiply_diagonal_l_r( d, unsafe_cast(int3d(double[9], d), int3d {i,j,0}), beta[{i,j,0}]._0, beta[{i,j,0}]._3, beta[{i,j,0}]._4 )
      invert_matrix_r( d, unsafe_cast(int3d(double[9], d), int3d {i,j,0}) )
    end
  end

  for j = bounds.lo.y, bounds.hi.y+1 do
    for i = bounds.lo.x, bounds.hi.x+1 do
      sol[{i,j,0}].rho = -sol[{i,j,0}].rho
      sol[{i,j,0}].w   = -sol[{i,j,0}].w
      sol[{i,j,0}].p   = -sol[{i,j,0}].p
    end
  end

  if periodic_z then
    for j = bounds.lo.y, bounds.hi.y+1 do
      for i = bounds.lo.x, bounds.hi.x+1 do
        get_Rinv_r( rho_avg[{i,j,0}], sos_avg[{i,j,0}], Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,j,0}) )
        multiply_diagonal_l_r( Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,j,0}), alpha[{i,j,0}]._0, alpha[{i,j,0}]._3, alpha[{i,j,0}]._4 )
      end
    end
  end

  for k = 1,N do
    for j = bounds.lo.y, bounds.hi.y+1 do
      for i = bounds.lo.x, bounds.hi.x+1 do
        var Rinv_i : double[9]
        get_Rinv( rho_avg[{i,j,k}], sos_avg[{i,j,k}], Rinv_i ) -- Get Rinv at i

        var mat : double[9]
        mult_matrix_matrix( Rinv_i, d[{i,j,k-1}], mat )
        multiply_diagonal_l( mat, alpha[{i,j,k}]._0, alpha[{i,j,k}]._3, alpha[{i,j,k}]._4 )

        var gammaRinv_im1 : double[9]
        get_Rinv( rho_avg[{i,j,k-1}], sos_avg[{i,j,k-1}], gammaRinv_im1 ) -- Get Rinv at i-1
        multiply_diagonal_l( gammaRinv_im1, gamma[{i,j,k-1}]._0, gamma[{i,j,k-1}]._3, gamma[{i,j,k-1}]._4 )

        multiply_diagonal_l( Rinv_i, beta[{i,j,k}]._0, beta[{i,j,k}]._3, beta[{i,j,k}]._4 )

        mult_matrix_matrix_r( mat, gammaRinv_im1, d, unsafe_cast(int3d(double[9], d), int3d {i,j,k}) )
        axpby_r( d, unsafe_cast(int3d(double[9], d), int3d {i,j,k}), Rinv_i, -1., 1., 9 ) -- Delta_i = beta_i - alpha_i Delta_i gamma_i-1
        invert_matrix_r( d, unsafe_cast(int3d(double[9], d), int3d {i,j,k}) )

        var prim : double[3] = array( sol[{i,j,k-1}].rho, sol[{i,j,k-1}].w, sol[{i,j,k-1}].p )
        var cprime = mult_matrix_vector( mat, prim )
        sol[{i,j,k}].rho = -sol[{i,j,k}].rho - cprime[0]
        sol[{i,j,k}].w   = -sol[{i,j,k}].w   - cprime[1]
        sol[{i,j,k}].p   = -sol[{i,j,k}].p   - cprime[2]

        if periodic_z then
          mult_matrix_matrix_r( mat, Uinv[{i,j,k-1}], Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,j,k}) )
          multiply_diagonal_l_r( Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,j,k}), -1., -1., -1. )
        end

      end
    end
  end

  if periodic_z then
    for j = bounds.lo.y, bounds.hi.y+1 do
      for i = bounds.lo.x, bounds.hi.x+1 do
        var Rinv : double[9]
        get_Rinv( rho_avg[{i,j,N-1}], sos_avg[{i,j,N-1}], Rinv ) -- Get Rinv
        multiply_diagonal_l( Rinv, gamma[{i,j,N-1}]._0, gamma[{i,j,N-1}]._3, gamma[{i,j,N-1}]._4 )
        axpby_r( Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,j,N-1}), Rinv, 1., -1., 9 )
      end
    end
  end

  -- Back substitution
  for j = bounds.lo.y, bounds.hi.y+1 do
    for i = bounds.lo.x, bounds.hi.x+1 do
      var prim : double[3] = array( sol[{i,j,N-1}].rho, sol[{i,j,N-1}].w, sol[{i,j,N-1}].p )
      var cprime = mult_matrix_vector( d[{i,j,N-1}], prim )
      sol[{i,j,N-1}].rho = - cprime[0]
      sol[{i,j,N-1}].w   = - cprime[1]
      sol[{i,j,N-1}].p   = - cprime[2]
    end
  end

  if periodic_z then
    for j = bounds.lo.y, bounds.hi.y+1 do
      for i = bounds.lo.x, bounds.hi.x+1 do
        var tmp : double[9]
        mult_matrix_matrix( d[{i,j,N-1}], Uinv[{i,j,N-1}], tmp )
        for ii = 0,9 do
          (Uinv[{i,j,N-1}])[ii] = -tmp[ii]
        end
      end
    end
  end

  for k = N-1,0,-1 do
    for j = bounds.lo.y, bounds.hi.y+1 do
      for i = bounds.lo.x, bounds.hi.x+1 do
        var gammaRinv_im1 : double[9]
        get_Rinv( rho_avg[{i,j,k-1}], sos_avg[{i,j,k-1}], gammaRinv_im1 ) -- Get Rinv at i-1
        multiply_diagonal_l( gammaRinv_im1, gamma[{i,j,k-1}]._0, gamma[{i,j,k-1}]._3, gamma[{i,j,k-1}]._4 )

        var prim : double[3] = array( sol[{i,j,k}].rho, sol[{i,j,k}].w, sol[{i,j,k}].p )
        var cprime = mult_matrix_vector( gammaRinv_im1, prim )    

        prim = array( sol[{i,j,k-1}].rho, sol[{i,j,k-1}].w, sol[{i,j,k-1}].p )
        axpby( prim, cprime, 1., 1., 3 )

        cprime = mult_matrix_vector( d[{i,j,k-1}], prim )
        sol[{i,j,k-1}].rho = - cprime[0]
        sol[{i,j,k-1}].w   = - cprime[1]
        sol[{i,j,k-1}].p   = - cprime[2]

        if periodic_z then
          var tmp : double[9]
          mult_matrix_matrix( gammaRinv_im1, Uinv[{i,j,k}], tmp )
          axpby( tmp, Uinv[{i,j,k-1}], -1., -1., 9 )
          mult_matrix_matrix_r( d[{i,j,k-1}], tmp, Uinv, unsafe_cast(int3d(double[9], Uinv), int3d {i,j,k-1}) )
        end -- periodic
      end
    end
  end

  -- Sherman-Morrison correction
  if periodic_z then
    for j = bounds.lo.y, bounds.hi.y+1 do
      for i = bounds.lo.x, bounds.hi.x+1 do
        var M : double[9] = array(1., 0., 0., 0., 1., 0., 0., 0., 1.) -- Identity matrix
        axpby( M, Uinv[{i,j,0}],   1., 1., 9 )
        axpby( M, Uinv[{i,j,N-1}], 1.,-1., 9 )
        invert_matrix( M )
        var prim : double[3] = array( sol[{i,j,0}].rho - sol[{i,j,N-1}].rho, sol[{i,j,0}].w - sol[{i,j,N-1}].w, sol[{i,j,0}].p - sol[{i,j,N-1}].p )
        var corrfact : double[3] = mult_matrix_vector( M, prim )

        for k = 0,N do
          var cprime = mult_matrix_vector( Uinv[{i,j,k}], corrfact )
          sol[{i,j,k}].rho = sol[{i,j,k}].rho - cprime[0]
          sol[{i,j,k}].w   = sol[{i,j,k}].w   - cprime[1]
          sol[{i,j,k}].p   = sol[{i,j,k}].p   - cprime[2]
        end

        -- Copy last edge if periodic
        sol[{i,j,N}].rho = sol[{i,j,0}].rho
        sol[{i,j,N}].w   = sol[{i,j,0}].w
        sol[{i,j,N}].p   = sol[{i,j,0}].p
      end
    end
  end

end

function make_solve_tridiagonal_z( fi, fn )
  -- fi: field index in coeff region (_1, _2, _3 ...)
  -- fn: field name in primitive region (u, v, w ...)

  local solve_tridiagonal_z __demand(__inline) task solve_tridiagonal_z( alpha   : region( ispace(int3d), coeffs ),
                                                                         beta    : region( ispace(int3d), coeffs ),
                                                                         gamma   : region( ispace(int3d), coeffs ),
                                                                         sol     : region( ispace(int3d), primitive ),
                                                                         Uinv    : region( ispace(int3d), double) )
  where
    reads (alpha.[fi], gamma.[fi]), reads writes ( beta.[fi], sol.[fn], Uinv )
  do
  
    var bounds = sol.ispace.bounds
    var N : int = bounds.hi.z + 1
    if periodic_z then
      N = bounds.hi.z -- Don't solve for last edge if periodic
    end
  
    for j = bounds.lo.y, bounds.hi.y+1 do
      for i = bounds.lo.x, bounds.hi.x+1 do
  
        var beta0 : double = beta[{i,j,0}].[fi]
        if periodic_z then
          -- Sherman-Morrison reduced matrix
          Uinv[{i,j,0}] = -beta[{i,j,0}].[fi]
          beta[{i,j,N-1}].[fi] = beta[{i,j,N-1}].[fi] + gamma[{i,j,N-1}].[fi] * alpha[{i,j,0}].[fi] / beta[{i,j,0}].[fi]
          beta[{i,j,0}].[fi] = 2.*beta[{i,j,0}].[fi]
        end

        -- Forward substitution
        beta[{i,j,0}].[fi] = 1./beta[{i,j,0}].[fi]
        sol[{i,j,0}].[fn] = beta[{i,j,0}].[fi] * sol[{i,j,0}].[fn]
        if periodic_z then
          Uinv[{i,j,0}] =  beta[{i,j,0}].[fi] * Uinv[{i,j,0}]
        end

      end
    end

    for k = 1,N do
      for j = bounds.lo.y, bounds.hi.y+1 do
        for i = bounds.lo.x, bounds.hi.x+1 do
          beta[{i,j,k}].[fi] = 1./( beta[{i,j,k}].[fi] - alpha[{i,j,k}].[fi] * beta[{i,j,k-1}].[fi] * gamma[{i,j,k-1}].[fi] )
          sol[{i,j,k}].[fn] = beta[{i,j,k}].[fi] * (sol[{i,j,k}].[fn] - alpha[{i,j,k}].[fi] * sol[{i,j,k-1}].[fn])
          if periodic_z then
            Uinv[{i,j,k}] = beta[{i,j,k}].[fi] * ( - alpha[{i,j,k}].[fi] * Uinv[{i,j,k-1}])
          end
        end
      end
    end

    if periodic_z then
      for j = bounds.lo.y, bounds.hi.y+1 do
        for i = bounds.lo.x, bounds.hi.x+1 do
          Uinv[{i,j,N-1}] = Uinv[{i,j,N-1}] + beta[{i,j,N-1}].[fi] * ( gamma[{i,j,N-1}].[fi] )
        end
      end
    end
    
    -- Backward elimination
    for k = N-1,0,-1 do
      for j = bounds.lo.y, bounds.hi.y+1 do
        for i = bounds.lo.x, bounds.hi.x+1 do
          sol[{i,j,k-1}].[fn] = sol[{i,j,k-1}].[fn] - beta[{i,j,k-1}].[fi] * gamma[{i,j,k-1}].[fi] * sol[{i,j,k}].[fn]
          if periodic_z then
            Uinv[{i,j,k-1}] = Uinv[{i,j,k-1}] - beta[{i,j,k-1}].[fi] * gamma[{i,j,k-1}].[fi] * Uinv[{i,j,k}]
          end
        end
      end
    end
  
    -- Sherman-Morrison correction
    if periodic_z then
      for j = bounds.lo.y, bounds.hi.y+1 do
        for i = bounds.lo.x, bounds.hi.x+1 do
          var beta0 = (1./beta[{i,j,0}].[fi])/2.
          var corr_factor : double = ( sol[{i,j,0}].[fn] - (alpha[{i,j,0}].[fi] / beta0) * sol[{i,j,N-1}].[fn] )
          corr_factor = corr_factor / ( 1. + Uinv[{i,j,0}] - (alpha[{i,j,0}].[fi] / beta0) * Uinv[{i,j,N-1}] )

          for k = 0,N do
            sol[{i,j,k}].[fn] = sol[{i,j,k}].[fn] - corr_factor * Uinv[{i,j,k}]
          end

          -- Copy last edge if periodic
          sol[{i,j,N}].[fn] = sol[{i,j,0}].[fn]
        end
      end
    end

  
  end
  return solve_tridiagonal_z
end

