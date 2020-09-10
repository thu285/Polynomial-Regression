# Final Project CPS 338

## Aiden & Thu

## Description:

This project is to implement linear regression and polynomial regression using Householder matrix and backsubstitution. The functions are written in C++ and can be used in R using the library Rcpp. The file contains the following functions:
 
 ### linfit338: do linear regression, output a list containing coefficients of the line and two vectors FX and FY that aids in the production of smoothed regression lines.

### polyfit338: do polynomial regression, output a list containing coefficients of the line and
two vectors FX and FY that aids in the production of smoothed regression lines.

### house: compute the householder vectors v, returns a list containing the vector v and theta.

### qr338: does QR decomposition, return the matrix R with 0s substituted by values of v.

### hls: computes the coefficients that guarantee least-squares using matrix R and back substitution.

### some other helpful helper functions as follows:
  + meanvec: taking average of a vector
  + matmult: multiply two matrices
  + inprod: calculate inner product of vectors
  + outprod: calculate outer product of vectors
  + tranpose: return the transpose of a matrix
  + makeid: returns an identity matrix with given dimension.
