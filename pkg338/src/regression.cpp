/* Author: Aiden & Thu
 *
 * Description:
 *
 * The file contains the following functions:
 *
 * - linfit338: do linear regression, output a list containing coefficients of the line and
 * two vectors FX and FY that aids in the production of smoothed regression lines.
 *
 * - polyfit338: do polynomial regression, output a list containing coefficients of the line and
 * two vectors FX and FY that aids in the production of smoothed regression lines.
 *
 * - house: compute the householder vectors v, returns a list containing the vector v and theta.
 *
 * - qr338: does QR decomposition, return the matrix R with 0s substituted by values of v.
 *
 * - hls: computes the coefficients that guarantee least-squares using matrix R and back
 * substution.
 *
 * - some other helpful helper functions as follows:
 *  + meanvec: taking average of a vector
 *  + matmult: multiply two matrices
 *  + inprod: calculate inner product of vectors
 *  + outprod: calculate outer product of vectors
 *  + tranpose: return the transpose of a matrix
 *  + makeid: returns an identity matrix with given dimension.
 */

#include <Rcpp.h>
using namespace Rcpp;

/* taking average of a vector
 * param: vector x
 * output: scalar value represent average of the elements in the vector
 */
double meanvec(NumericVector x) {
  int n = x.size();
  double tot = 0;
  for (int i = 0; i < n; i++) {
    tot += x[i];
  }
  tot /= n;
  return tot;
}



/* multiplying two matrices
 * param: two matrices a and b
 * output: a single matrix c that is product of the two
 */
NumericMatrix matmult(NumericMatrix a, NumericMatrix b) {
  int m = a.nrow();
  int n = a.ncol();
  int p = b.ncol();
  
  NumericMatrix c(m,p);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < p; j++) {
      c(i,j) = 0;
      for (int k = 0; k < n; k++) {
        c(i,j) += a(i,k) * b(k,j);
      }
    }
  }
  return c;
}



/* inner product
 * param: two vectors u and v
 * output: scalar value that represents inner product of the two
 */
double inprod(NumericVector u, NumericVector v) {
  int m = u.size();
  double val = 0;
  for (int i = 0; i < m; i++) {
    val += u[i] * v[i];
  }
  return val;
}



/* transpose matrix
 * param: matrix a
 * output: matrix b which is the transpose matrix of a
 */
NumericMatrix transpose(NumericMatrix a){
  int m = a.nrow();
  int n = a.ncol();
  NumericMatrix b(n,m);
  for (int i = 0; i < m; i++){
    for (int j = 0; j < n; j++){
      b(j,i) = a(i,j);
    }
  }
  return b;
}



/* identity matrix
 * param: a number n
 * output: the identity matrix of dimension NxN
 */
NumericMatrix makeid(int n) {
  NumericMatrix x(n,n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) { 
      x(i,j) = 0;
    }
  } for (int i = 0; i < n; i++) {
    x(i,i) = 1;
  }
  return x;
}



/* outer product
 * param: two vectors u and v
 * output: a matrix that stores the outer product of the two
 */
NumericMatrix outprod(NumericVector u, NumericVector v) {
  int m = u.size();
  int n = v.size();
  NumericMatrix mat(m,n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++){
      mat(i,j) = u[i] * v[j];
    }
  }
  return mat;
}



/* household vector
 * param: vector x where x represents values of predictor variable
 * output: a List containing "v" which is the householder vector and the value "theta"
 */
List house(NumericVector x) {
  // size of vector
  int n = x.size();
  // sigma
  NumericVector esx(n-1);
  for (int i = 0; i < n-1; i++) {
    esx[i] = x[i+1];
  }
  double sig = inprod(esx,esx);
  // v
  NumericVector v(n);
  v[0] = 1;
  for (int i = 1; i < n; i++) {
    v[i] = x[i];
  }
  // making theta initially
  double theta = 0;
  
  if ((sig == 0) && (x[0] >= 0)) {
    theta = 0;
  } else if ((sig == 0) && (x[0] < 0)) {
    theta = -2;
  } else {
    double mu = sqrt(pow(x[0], 2.0) + sig);
    if (x[0] <= 0) {
      v[0] = x[0] - mu;
    } else {
      v[0] = (-1 * sig) / (x[0] + mu);
    }
    theta = 2 * pow(v[0], 2.0) / (sig + pow(v[0], 2.0));
    for (int i = 1; i < n; i++) {
      v[i] /= v[0];
    }
    v[0] = 1;
  }
  List vt;
  vt["v"] = v;
  vt["theta"] = theta;
  return vt;
}



/* qr decomposition
 * param: matrix H which is the basis matrix for poly regression
 * output: matrix R which replaces 0s with values of the householder vector
 */
NumericMatrix qr338(NumericMatrix H) {
  int n = H.nrow();
  int m = H.ncol();
  if (m > n) stop("This matrix does not have a QR decomposition.");
  // making matrix R
  NumericMatrix R(n,m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      R(i,j) = H(i,j);
    }
  }
  // updating R for each column
  for (int j = 0; j < m; j++) {
    // subsetting correct column vector
    NumericVector x(n-j);
    for (int i = 0; i < n-j; i++) {
      x[i] = R(i+j,j);
    }
    // householder vector and constant
    List par = house(x);
    NumericVector v = as<NumericVector>(par["v"]);
    double theta = as<double>(par["theta"]);
    // creating householder matrix
    NumericMatrix I = makeid(n-j);
    NumericMatrix O = outprod(v,v);
    NumericMatrix T(n-j,n-j);
    for (int k = 0; k < n-j; k++) {
      for (int l = 0; l < n-j; l++) {
        T(k,l) = I(k,l) - (theta * O(k,l));
      }
    }
    // creating submatrix of R
    NumericMatrix A(n-j,m-j);
    for (int k = 0; k < n-j; k++) {
      for (int l = 0; l < m-j; l++) {
        A(k,l) = R(k+j,l+j);
      }
    }
    // updating R
    NumericMatrix B = matmult(T,A);
    for (int k = 0; k < n-j; k++) {
      for (int l = 0; l < m-j; l++) {
        R(k+j,l+j) = B(k,l);
      }
    }
    // adding values of householder vector to H
    NumericVector subV = v[Range(1,n-j-1)];
    if (j < n-1) {
      for (int k = 0; k < n-j-1; k++) {
        R(k+j+1,j) = subV[k];
      }
    }
  }
  return R;
}



/* back substitution
 * param: matrix mat
 * output: vector x that stores the result of the equation
 */
NumericVector backSub(NumericMatrix mat){ 
  int N = mat.nrow();
  NumericVector x(N);  // array to store solution
  
  // start calculating from last equation up to the first
  for (int i = N-1; i >= 0; i--) 
  { 
    // start with the RHS of the equation
    x[i] = mat(i,N);
    // Initialize j to i+1 since matrix is upper triangular
    for (int j=i+1; j<N; j++) 
    { 
      /* subtract all the lhs values 
       * except the coefficient of the variable 
       * whose value is being calculated */
      x[i] -= mat(i,j)*x[j]; 
    } 
    
    /* divide the RHS by the coefficient of the 
     unknown being calculated */
    x[i] = x[i]/mat(i,i); 
  } 
  return x;
} 



/* creating a basis matrix for polynomials
 * param: vector data which stores predictor variable values, degree of the desired polynomial
 * output: basis matrix H
 */
NumericMatrix createMatrixH(NumericVector data, int degree){
  int nrows = data.size();
  NumericMatrix H(nrows,degree+1);
  NumericVector c(degree+1);
  for (int i = 0; i < degree+1; i++){
    for (int j = 0; j < nrows; j++){
      H(j,i) = 1;
      if (i > 0){
        H(j,i) = pow(data[j],i);
      }
    }
  }
  return H;
}



/* least squares for matrix
 * param: basis matrix H and vector y which stores response variable values.
 */
NumericVector hls(NumericMatrix H, NumericVector y) {
  NumericMatrix R = qr338(H);
  int N = H.nrow();
  int M = H.ncol();
  for (int j = 0; j < M; j++){
    NumericVector Rcol =  R(_,j);
    NumericVector subRcol = Rcol[Range(j+1,N-1)]; // N-j-1 elements
    int nrowV = N-j; // (N-1) - (j+1) + 1 + 1 = N - j elements
    NumericMatrix v(nrowV,1);
    v(0,0) = 1;
    for (int i = 0; i < nrowV-1; i++){
      v(i+1,0) = subRcol[i];
    }
    NumericMatrix vtrans = transpose(v);
    NumericMatrix vprod = matmult(vtrans,v);
    double theta = 2/vprod(0,0);
    NumericMatrix subY(N-j,1);
    for (int i = 0; i < N-j; i++){
      subY(i,0) = y[i+j];
    }
    NumericMatrix prod = matmult(vtrans,subY);
    double scalar = theta*prod(0,0);
    for (int i = 0; i < N-j; i++){
      y[i+j] = subY(i,0) - scalar*v(i,0);
    }
  }
  // create the matrix for system equation
  NumericMatrix mat(M,M+1);
  for (int i = 0; i < M; i++){
    for (int j = 0; j < M; j++){
      mat(i,j) = R(i,j);
    }
    mat(i,M) = y(i);
  }
  NumericVector beta = backSub(mat); // solve for the coefficients using back substitution
  
  return beta;
}

/* fitting polynomial regression
 * param:   vector data which stores predictor variable values
 *          vector response which stores response variable values
 *          degree which represents the desired degree of the polynomial
 * output:  a List containing "coef" which represents coefficients of the polynomial,
 *          vector FX and FY that aid in producing the smoothed predicted curve
 */
// [[Rcpp::export]]
List polyfit338(NumericVector data, NumericVector response, int degree){
  NumericMatrix basis = createMatrixH(data,degree);
  NumericMatrix basisCopy = clone(basis);
  NumericVector yCopy = clone(response);
  NumericVector beta = hls(basisCopy,yCopy);
  /*
   find the range of the predictor variable to output the fitted line stretching
   from the smallest data point to the largest data point
   */
  double min = data[0];
  double max = data[0];
  for (int i = 1; i < data.size(); i++){
    if (data[i] < min){
      min = data[i];
    }
    if (data[i] > max){
      max = data[i];
    }
  }
  double step = 0.1;
  int size = int((max-min)*10)+1;
  // these two FX and FY vectors are for outputting smoothed predicted curve
  NumericVector FX(size);
  NumericVector FY(size);
  for (int i = 0; i < size; i++){
    FX(i) = min + step*i;
    FY(i) = 0;
    for (int j = 0; j <= degree; j++){
      FY(i) += beta[j]*pow(FX(i),j);
    }
  }
  // return the coefficients and the vectors for smoothed regression line in a list
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("coef") = beta, Rcpp::Named("X")=FX,Rcpp::Named("Y")=FY);
  return res;
}

/* fitting linear regression
 * param:   vector x which stores predictor variable values
 *          vector y which stores response variable values
 * output:  a List containing "coef" which represents coefficients of the line,
 *          vector FX and FY that aid in producing the smoothed predicted line
 */
// [[Rcpp::export]]
List linfit338(NumericVector x, NumericVector y) {
  int n = x.size();
  double xm = meanvec(x);
  double ym = meanvec(y);
  NumericVector xa(n);
  for (int i = 0; i < n; i++) {
    xa[i] = x[i] - xm;
  }
  NumericVector ya(n);
  for (int i = 0; i < n; i++) {
    ya[i] = y[i] - ym;
  }
  double b1 = inprod(xa, ya) / inprod(xa, xa);
  double b0 = ym - (b1 * xm);
  NumericVector beta = NumericVector::create(b0, b1);
  /*
   find the range of the predictor variable to output the fitted line stretching
   from the smallest data point to the largest data point
   */
  double min = x[0];
  double max = x[0];
  for (int i = 1; i < x.size(); i++){
    if (x[i] < min){
      min = x[i];
    }
    if (x[i] > max){
      max = x[i];
    }
  }
  double step = 0.1;
  int size = int((max-min)*10)+1;
  // these two FX and FY vectors are for outputting smoothed predicted curve
  NumericVector FX(size);
  NumericVector FY(size);
  for (int i = 0; i < size; i++){
    FX(i) = min + step*i;
    FY(i) = beta[0] + beta[1]*FX(i);
  }
  // return the coefficients and the vectors for smoothed regression line in a list
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("coef") = beta, Rcpp::Named("X")=FX,Rcpp::Named("Y")=FY);
  return res;
}
