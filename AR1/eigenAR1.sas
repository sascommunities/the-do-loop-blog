/* SAS program to accompany the article 
   "An explicit formula for eigenvalues of the AR(1) correlation matrix"
   by Rick Wicklin, published 03MAR2025 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2025/03/03/formula-eigenvalues-ar1-corr.html 

   This program shows how to compute the k_th eigenvalue and eigenvector
   of a (n x n) AR1(rho) correlation matrix by using the formulas in 
   P. J. Sherman, "On the Eigenstructure of the AR(1) Covariance," 
      2023 IEEE Statistical Signal Processing Workshop (SSP), Hanoi, Vietnam, 2023, pp. 6-10.
      https://ieeexplore.ieee.org/abstract/document/10208005
*/


/* First, visualize the rational trig function P_n(x; rho) on the interval [0, pi].
   The roots of these functions are used in the formulas for the eigenvalues
   and eigenvectors */
proc iml;
/* if you need to explicitly form an AR1 matrix, use this function
   https://blogs.sas.com/content/iml/2022/12/14/heterogeneous-covariance-matrices.html
*/
start AR1Corr(dim, rho);
   u = cuprod(j(1,dim-1,rho)); /* cumulative product */
   return( toeplitz(1 || u) );
finish;

/* set parameters for example */
n = 4;
rho = 0.5;
A = AR1Corr(n, rho);
print A;

/* scalar formula for rational trig function */
/*
start SineRational(w, rho, n);
   if w=0 then 
      return (n+1) -2*rho*n + rho##2*(n-1); 
   p = sin( (n+1)*w ) -2*rho*sin(n*w) + rho##2*sin((n-1)*w);
   return( p / sin(w) );
finish;
*/

/* Eqn 1: Create a vectorize version that avoids the undefined value at w=0 */
start SineRational(w, rho, n);
   y = j(nrow(w), ncol(w), .);
   idx = loc(w = 0);
   if ncol(idx)>0 then 
      /* limit as w->0 of sin(a*w)/sin(w) = a */
      y[idx] = (n+1) -2*rho*n + rho##2*(n-1); 
   idx = loc(w ^= 0);
   if ncol(idx)>0 then do;
      z = w[idx];
      p = sin( (n+1)*z ) -2*rho*sin(n*z) + rho##2*sin((n-1)*z);
      y[idx] = p / sin(z);
   end;
   return y;
finish;

/* visualize the function, its roots, and the intervals on which the roots are found */
pi = constant('pi');
w = do(0, pi, 0.01);
y = SineRational(w, rho, n);

endpts = (0:n)*pi / (n+1);
print endpts;
title "Graph of Trigonometric Function for n=4 and rho=0.5";
title2 "One Root in Each Interval";
call series(w, y) grid={x,y} other="refline 0; 
     refline 0 0.6283185 1.2566371 1.8849556 2.5132741 /axis=x lineattrs=(color=DarkRed);
     xaxis grid values=(0 to 3 by 0.25) valueshint";
QUIT;


/****************************************/
/* MAIN COMPUTATION: SHERMAN'S FORMULAS */
/****************************************/


proc iml;
/* Eqn 1: Define the rational trig function that avoids the undefined value at w=0 */
start SineRational(w, rho, n);
   if w=0 then 
      return (n+1) -2*rho*n + rho##2*(n-1); 
   p = sin( (n+1)*w ) -2*rho*sin(n*w) + rho##2*sin((n-1)*w);
   return( p / sin(w) );
finish;
/* A one-variable version of the rational trig function; put parameters into GLOBAL clause */
start SineRationalFun(w) global( g_rho, g_dim );
   return(  SineRational(w, g_rho, g_dim) ); 
finish;

/* Get the k_th roots of a rational trig function. 
   By default, get all eigenvalues when seq=1:n. The last arg 
   determines a subset of roots. 
   For example, 1:3 will get the first 3 roots. */
start GetSineRationalRoots(rho, n, seq=1:n) global( g_rho, g_dim );
   g_rho = rho; g_dim = n;   /* copy parameters to global vars */
   pi = constant('pi');
   endpts = (0:n)*pi / (n+1);
   intervals = endpts[1:n] || endpts[2:(n+1)];  /* get all intervals */
   intervals = intervals[seq,];      /* get only the specified roots */
   roots = froot("SineRationalFun", intervals);
   return roots;
finish;

/* Get eigenvalues of an (n x n) AR1(rho) matrix.
   By default, get all eigenvalues 1:n. Specify the last arg 
   to get only a subset. For example, 1:3 will get only the top 3. */
start EigenvalsAR1(rho, n, seq=1:n);
    w = GetSineRationalRoots(rho, n, seq);
    lambda = (1-rho##2) / (1 - 2*rho*cos(w) + rho##2);
    return lambda;
finish;

/* Example: find roots and eigenvalues for n=4 and rho=0.5 */
n = 4;
rho = 0.5;
w = GetSineRationalRoots(rho, n);
print w[L="w (Roots)" r=('w1':'w4')];

/* find the eigenvalues */
eval = EigenvalsAR1(rho, n);
print eval[L="Eigenvalues" r=('lambda1':'lambda4')];

/* verify by calling EIGVAL */
start AR1Corr(dim, rho);
   u = cuprod(j(1,dim-1,rho)); /* cumulative product */
   return( toeplitz(1 || u) );
finish;

A = AR1Corr(n, rho);
ev2 = eigval(A);
diff = max(abs(eval-ev2));
print diff[L='max diff'];



/* Example 2: find only top two eigenvals */
seq = 1:2;
/*
w = GetSineRationalRoots(rho, n, seq);
print w[L="w (Roots)" r=seq];
*/
eval2 = EigenvalsAR1(rho, n, seq);
print eval2[L="Eigenvalues" r=seq];

/* Example 3: find first and last eigenvalues */
seq = {1 4};
eval3 = EigenvalsAR1(rho, n, seq);
print eval3[L="Eigenvalues" r=seq];
   

/* ----- EIGENVECTORS ----- */

/* helper functions to return constants for the eigenvector formulas */
start GetPhi(w, rho);
   return atan( rho*sin(w) / (1 - rho*cos(w)) );
finish;
start GetC(w, rho);
   return sqrt(1 - rho*cos(w) + rho##2);
finish;

/* Return the eigenvectors of an (n x n) AR(1) correlation matrix.
   By default, get all n eigenvectors. Specify the last arg 
   to get only a subset. For example, 1:3 will get only the top 3. */
start EigenvectorsAR1(rho, n, seq=1:n);
   w = GetSineRationalRoots(rho, n, seq);
   w = rowvec(w);   /* change to row vector */
   phi = GetPhi(w, rho);
   c = GetC(w, rho);
   /* use the outer product to form matrix where columns are eigenvectors */
   theta = T(1:n) * w + phi;
   return( c # sin( theta ) );
finish;

/* Example: compute eigenvalues for n=4 */
evec = EigenvectorsAR1(rho, n);
print evec[c=('v1':'v4')];

/* Check: are they eigenvectors? */
eval = EigenvalsAR1(rho, n);
diff = A*evec - eval`#evec;
print diff[L='Zero matrix'];

/* Check: Do the formulas agree with the EIGVEC call in IML? */
/* by convention, eigenvectors are standardardized to have unit length */
evec_std = evec / sqrt(evec[##,]);
print evec_std;
evec2 = eigvec(A);
print evec2;

/* The vectors might point in different directions. If so, change sign */
do k = 1 to n;
   if sign(evec_std[1,k]) ^= sign(evec2[1,k]) then 
      evec_std[,k] = -evec_std[,k];
end;

/* now they should be the same */
diff = max(abs(evec_std-evec2));
print diff[L='max diff'];


/* Example 2: find only top two eigenvectors */
seq = 1:2;
ev2 = EigenvectorsAR1(rho, n, seq);
print ev2[c=seq];

/* Example 3: find first and last eigenvectors */
seq = {1 4};
ev3 = EigenvectorsAR1(rho, n, seq);
print ev3[c=seq];

QUIT;
