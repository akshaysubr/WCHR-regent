import "regent"

local lapack = terralib.includecstring [[
extern void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

extern void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
]]

if os.execute("bash -c \"[ `uname` == 'Darwin' ]\"") == 0 then
  terralib.linklibrary("libblas.dylib")
  terralib.linklibrary("liblapack.dylib")
else
  terralib.linklibrary("libblas.so")
  terralib.linklibrary("liblapack.so")
end

local c       = regentlib.c
local cmath   = terralib.includec("math.h")
local cstdlib = terralib.includec("stdlib.h")

fspace coeffs {
  _0 : double,
  _1 : double,
  _2 : double,
}

fspace primitive {
  rho : double,
  u   : double,
  p   : double,
}

-- All matrices here are in Fortran order

terra get_Rinv( rho : double, sos : double, Rinv : &double )
  Rinv[0 + 3*0] = 0.; Rinv[0 + 3*1] = -0.5*rho*sos; Rinv[0 + 3*2] = 0.5;
  Rinv[1 + 3*0] = 1.; Rinv[1 + 3*1] = 0.;           Rinv[1 + 3*2] = -1./(sos*sos);
  Rinv[2 + 3*0] = 0.; Rinv[2 + 3*1] =  0.5*rho*sos; Rinv[2 + 3*2] = 0.5;
end

local terra multiply_diagonal_l( matrix : &double, d0 : double, d1 : double, d2 : double )
  matrix[0 + 3*0] = d0*matrix[0 + 3*0]; matrix[0 + 3*1] = d0*matrix[0 + 3*1]; matrix[0 + 3*2] = d0*matrix[0 + 3*2];
  matrix[1 + 3*0] = d1*matrix[1 + 3*0]; matrix[1 + 3*1] = d1*matrix[1 + 3*1]; matrix[1 + 3*2] = d1*matrix[1 + 3*2];
  matrix[2 + 3*0] = d2*matrix[2 + 3*0]; matrix[2 + 3*1] = d2*matrix[2 + 3*1]; matrix[2 + 3*2] = d2*matrix[2 + 3*2];
end

local terra multiply_diagonal_r( matrix : &double, d0 : double, d1 : double, d2 : double )
  matrix[0 + 3*0] = d0*matrix[0 + 3*0]; matrix[0 + 3*1] = d1*matrix[0 + 3*1]; matrix[0 + 3*2] = d2*matrix[0 + 3*2];
  matrix[1 + 3*0] = d0*matrix[1 + 3*0]; matrix[1 + 3*1] = d1*matrix[1 + 3*1]; matrix[1 + 3*2] = d2*matrix[1 + 3*2];
  matrix[2 + 3*0] = d0*matrix[2 + 3*0]; matrix[2 + 3*1] = d1*matrix[2 + 3*1]; matrix[2 + 3*2] = d2*matrix[2 + 3*2];
end

local terra mult_matrix_vector( matrix : &double, vector : &double )
  var output : double[3]

  output[0] = matrix[0 + 3*0]*vector[0] + matrix[0 + 3*1]*vector[1] + matrix[0 + 3*2]*vector[2]
  output[1] = matrix[1 + 3*0]*vector[0] + matrix[1 + 3*1]*vector[1] + matrix[1 + 3*2]*vector[2]
  output[2] = matrix[2 + 3*0]*vector[0] + matrix[2 + 3*1]*vector[1] + matrix[2 + 3*2]*vector[2]

  return output
end

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

local terra invert_matrix( matrix : &double )
  var ipiv  : int[3]
  var lwork : int[1]
  var work  : double[9]
  var info  : int[1]
  var N     : int[1]

  N[0] = 3
  lwork[0] = 9

  lapack.dgetrf_(N,N,matrix,N,ipiv,info);
  regentlib.assert(info[0] == 0, "DGETRF did not work as expected. Check for errors.")

  lapack.dgetri_(N,matrix,N,ipiv,work,lwork,info);
  regentlib.assert(info[0] == 0, "DGETRI did not work as expected. Check for errors.")
end

local terra print_vector( vector : &double )

  for i = 0,3 do
    c.printf(" %11.8f ", vector[i])
    c.printf("\n")
  end
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

local terra random_number()
  return ( [double](cstdlib.rand()) / [double](cstdlib.RAND_MAX + 1.) )
end

local terra axpby( x : &double, y : &double, a : double, b : double, N : int )
  for i = 0,N do
    x[i] = a*x[i] + b*y[i]
  end
end

task solve_block_tridiagonal( alpha   : region( ispace(int3d), coeffs ),
                              beta    : region( ispace(int3d), coeffs ),
                              gamma   : region( ispace(int3d), coeffs ),
                              rho_avg : region( ispace(int3d), double ),
                              sos_avg : region( ispace(int3d), double ),
                              sol     : region( ispace(int3d), primitive ),
                              d       : region( ispace(int3d), &double ) )
where
  reads(alpha, beta, gamma, rho_avg, sos_avg, sol, d), writes(sol, d)
do
  var bounds = sol.ispace.bounds
  var N : int = bounds.hi.x + 1

  for k = bounds.lo.z, bounds.hi.z+1 do
    for j = bounds.lo.y, bounds.hi.y+1 do

      -- Forward substitution
      get_Rinv( rho_avg[{0,j,k}], sos_avg[{0,j,k}], d[{0,j,k}] )
      multiply_diagonal_l( d[{0,j,k}], beta[{0,j,k}]._0, beta[{0,j,k}]._1, beta[{0,j,k}]._2 )
      invert_matrix( d[{0,j,k}] )

      sol[{0,j,k}].rho = -sol[{0,j,k}].rho
      sol[{0,j,k}].u   = -sol[{0,j,k}].u
      sol[{0,j,k}].p   = -sol[{0,j,k}].p

      for i = 1,N do
        var Rinv_i : double[9]
        get_Rinv( rho_avg[{i,j,k}], sos_avg[{i,j,k}], Rinv_i ) -- Get Rinv at i

        var mat : double[9]
        mult_matrix_matrix( Rinv_i, d[{i-1,j,k}], mat )
        multiply_diagonal_l( mat, alpha[{i,j,k}]._0, alpha[{i,j,k}]._1, alpha[{i,j,k}]._2 )

        var gammaRinv_im1 : double[9]
        get_Rinv( rho_avg[{i-1,j,k}], sos_avg[{i-1,j,k}], gammaRinv_im1 ) -- Get Rinv at i-1
        multiply_diagonal_l( gammaRinv_im1, gamma[{i-1,j,k}]._0, gamma[{i-1,j,k}]._1, gamma[{i-1,j,k}]._2 )

        multiply_diagonal_l( Rinv_i, beta[{i,j,k}]._0, beta[{i,j,k}]._1, beta[{i,j,k}]._2 )

        mult_matrix_matrix( mat, gammaRinv_im1, d[{i,j,k}] )
        axpby( d[{i,j,k}], Rinv_i, -1., 1., 9 ) -- Delta_i = beta_i - alpha_i Delta_i gamma_i-1
        invert_matrix( d[{i,j,k}] )

        var prim : double[3] = array( sol[{i-1,j,k}].rho, sol[{i-1,j,k}].u, sol[{i-1,j,k}].p )
        var cprime = mult_matrix_vector( mat, prim )
        sol[{i,j,k}].rho = -sol[{i,j,k}].rho - cprime[0]
        sol[{i,j,k}].u   = -sol[{i,j,k}].u   - cprime[1]
        sol[{i,j,k}].p   = -sol[{i,j,k}].p   - cprime[2]

      end

      -- Back substitution
      var prim : double[3] = array( sol[{N-1,j,k}].rho, sol[{N-1,j,k}].u, sol[{N-1,j,k}].p )
      var cprime = mult_matrix_vector( d[{N-1,j,k}], prim )
      sol[{N-1,j,k}].rho = - cprime[0]
      sol[{N-1,j,k}].u   = - cprime[1]
      sol[{N-1,j,k}].p   = - cprime[2]

      for i = N-1,0,-1 do
        var gammaRinv_im1 : double[9]
        get_Rinv( rho_avg[{i-1,j,k}], sos_avg[{i-1,j,k}], gammaRinv_im1 ) -- Get Rinv at i-1
        multiply_diagonal_l( gammaRinv_im1, gamma[{i-1,j,k}]._0, gamma[{i-1,j,k}]._1, gamma[{i-1,j,k}]._2 )

        prim = array( sol[{i,j,k}].rho, sol[{i,j,k}].u, sol[{i,j,k}].p )
        cprime = mult_matrix_vector( gammaRinv_im1, prim )    

        prim = array( sol[{i-1,j,k}].rho, sol[{i-1,j,k}].u, sol[{i-1,j,k}].p )
        axpby( prim, cprime, 1., 1., 3 )

        cprime = mult_matrix_vector( d[{i-1,j,k}], prim )
        sol[{i-1,j,k}].rho = - cprime[0]
        sol[{i-1,j,k}].u   = - cprime[1]
        sol[{i-1,j,k}].p   = - cprime[2]
      end

    end
  end

end

local task print_sol( sol : region( ispace(int3d), primitive ) )
where
  reads(sol)
do
  var N : int = [int](sol.ispace.bounds.hi.x) + 1

  c.printf("rho = numpy.array([ ")
  for i = 0,N do
    c.printf(" %20.16e, ", sol[{i,3,6}].rho)
  end
  c.printf("\b\b])\n\n")

  c.printf("u = numpy.array([ ")
  for i = 0,N do
    c.printf(" %20.16e, ", sol[{i,3,6}].u)
  end
  c.printf("\b\b])\n\n")

  c.printf("p = numpy.array([ ")
  for i = 0,N do
    c.printf(" %20.16e, ", sol[{i,3,6}].p)
  end
  c.printf("\b\b])\n\n")
end

terra allocate_double( size : int64 )
  return [&double] ( c.malloc( size * sizeof(double) ) )
end

task main()

  var nx : int = 8
  var ny : int = 8
  var nz : int = 8

  var alpha = region( ispace(int3d, {x = nx, y = ny, z = nz}), coeffs )
  var beta  = region( ispace(int3d, {x = nx, y = ny, z = nz}), coeffs )
  var gamma = region( ispace(int3d, {x = nx, y = ny, z = nz}), coeffs )

  fill(alpha.{_0,_1,_2}, 3./16.)
  fill(beta.{_0,_1,_2}, 5./8.)
  fill(gamma.{_0,_1,_2}, 3./16.)

  var sol = region( ispace(int3d, {x = nx, y = ny, z = nz}), primitive )
  for i in sol do
    sol[i].rho = random_number()
    sol[i].u   = random_number()
    sol[i].p   = random_number()
  end
  print_sol( sol )

  var rho_avg = region( ispace(int3d, {x = nx, y = ny, z = nz}), double )
  var sos_avg = region( ispace(int3d, {x = nx, y = ny, z = nz}), double )
  for i in rho_avg do
    rho_avg[i] = 0.5 + 0.5*random_number()
    sos_avg[i] = 0.5 + 0.5*random_number()
  end

  c.printf("rho_avg = numpy.array([ ")
  for i = 0,nx do
    c.printf(" %20.16e, ", rho_avg[{i,3,6}])
  end
  c.printf("\b\b])\n\n")

  c.printf("sos_avg = numpy.array([ ")
  for i = 0,nx do
    c.printf(" %20.16e, ", sos_avg[{i,3,6}])
  end
  c.printf("\b\b])\n\n")

  var d = region( ispace(int3d, {x = nx, y = ny, z = nz}), &double )
  for i in d do
    d[i] = allocate_double(9)
  end

  solve_block_tridiagonal( alpha, beta, gamma, rho_avg, sos_avg, sol, d )
  print_sol( sol )
end

regentlib.start(main)
