
/* ---------------------------------------------------------- */
/* ------------------ Test Cases for cdftvn ----------------- */
/* ---------------------------------------------------------- */
/* The tests in this file are as follows:
   1. Uncorrelated Case: Sigma=I(3) should yield the product of univariate CDFs.
   2. Left-tail Orthant Probability: For b=0, the formula involves arcsin of correlations.
   3. Rank-1 Correlation Matrix: R = D + v*v' has a known integral formula.
   4. Genz and Bretz (2009) Numerical Example on p 4-5, Eqn 1.5: A specific case with a known solution.
   5. Partial Independence: If X1 is independent of (X2, X3), then CDF factorizes.
   6. Highly Correlated Case: Test stability with rho close to 1 using rank-1 as benchmark.
   7. Bug case: Need to manually switch the sign of the integral from CALL QUAD when integrating on [c,d] and c > d.
   8. Iterate all possible signs of a 3x3 correlation matrix.
   */
proc iml;
load module=_all_;

EPSILON = 1e-9; /* Tolerance for numerical comparisons */

/* Test 1: If Sigma=I(3), then prob is the product of the uncorrelated univariate CDFs */
b = {-1  0  2,
      0  0  0,
      1  2  1};
Sigma = I(3);
prob = cdftvn(b, Sigma);
correct = cdf("Normal", b[,1]) # cdf("Normal", b[,2]) # cdf("Normal", b[,3]);
maxDiff = max(abs(prob-correct));
if maxDiff > EPSILON then 
   print "--- ERROR in Test 1 ---", maxDiff prob correct;
else
   print "--- Test 1 passes ---";

/* Test 2: When the domain of integration is {x1 < 0, x2 < 0, and x3 < 0} then 
   there is a formula: 1/8 + 1/(4*pi) *(arsin(rho_12)+arsin(rho_13)+arsin(rho_23)) */
b = {0  0  0};
Sigma = {1    0.6 0.2,
         0.6  1   0.4,
         0.2 0.4  1};
prob = cdftvn(b, Sigma);
rho = Sigma[{2 3 6}];
pi = constant('pi');
correct = 1/8+ 1/(4*pi) *sum(arsin(rho));
maxDiff = max(abs(prob-correct));
if maxDiff > EPSILON then 
   print "--- ERROR in Test 2 ---", maxDiff prob correct;
else
   print "--- Test 2 passes ---";

/* Test 3: Genz and Bretz (2009, p. 16) show an exact solution for 
   R = D + V*V` where V is kxm matrix. When m=1, then R = "diagonal + rank-1".
   For m=1, Eqn 2.17 shows that the probability reduces to the integral on [0,1] of
   a product of functions that involve the 1-D CDF. Eqn 2.17 is for rectangular regions.
   When a=-Infinity, the second term vanishes.   
*/
/* Helper Function for Rank-1 Integration (Tests 3 and 6)
   Implements the 1-D integral for R = D + v*v' [Genz and Bretz 2009, Eqn 2.17]
   where D[i,i] = 1 - v[i]**2 and |v[i]| < 1.
*/
start Rank1Integrand(t) GLOBAL(g_lambda, g_b);
   if t < constant('maceps') then return 1;
   if t > 1 - 2*constant('maceps') then return 0;
   z = quantile("Normal", t);
   /* Compute Phi((b - lambda*z) / sqrt(1 - lambda^2)) */
   w = (g_b - g_lambda*z) / sqrt(1 - g_lambda##2);
   return( prod(cdf("Normal", w)) );
finish;


/* Test 3: Rank-1 Correlation Matrix
   A rank-1 correlation matrix has the form R = D + v*v'.
   The probability can be computed via a 1-D integral over [0,1].
*/
v = {0.8, 0.7, 0.6};
Sigma = v*v` + diag(1 - v##2);
b = {1.0 0.5 0.0};
/* Compute "correct" value via direct numerical integration of the rank-1 formula */
g_lambda = v; g_b = colvec(b);
call quad(correct, "Rank1Integrand", {0 1});
prob = cdftvn(b, Sigma);
maxDiff = max(abs(prob-correct));
if maxDiff > EPSILON then 
   print "--- ERROR in Test 3 ---", maxDiff prob correct;
else 
   print "--- Test 3 passes ---";

/* Test 4: Genz and Bretz (2009, p. 4-5) Numerical Example
   This is a specific test case provided in the literature with a known solution.
   Correct probability: 0.827984897456834.
*/
b = {1 4 2};
/* 
Sigma = {1.0          0.6           0.3333333333,
         0.6          1.0           0.7333333333,
         0.3333333333 0.7333333333  1.0         };
*/
Sigma = I(3);
Sigma[1,2] = 3/5;      Sigma[2,1] = Sigma[1,2];
Sigma[1,3] = 1/3;      Sigma[3,1] = Sigma[1,3];
Sigma[2,3] = 11/15;    Sigma[3,2] = Sigma[2,3];

prob = cdftvn(b, Sigma);
correct = 0.827984897456834;
maxDiff = max(abs(prob-correct));
if maxDiff > EPSILON then 
   print "--- ERROR in Test 4 ---", maxDiff prob correct;
else 
   print "--- Test 4 passes ---";

/* Test 5: Partial Independence
   If X1 is independent of (X2, X3), then the 3-D CDF is the 
   product of the univariate CDF and the bivariate CDF.
*/
b = {0.5 1.0 1.5};
Sigma = {1.0  0.0  0.0,
         0.0  1.0  0.5,
         0.0  0.5  1.0};
prob = cdftvn(b, Sigma);
/* Correct value via independence property */
correct = cdf("Normal", 0.5) * probbnrm(1.0, 1.5, 0.5);
maxDiff = max(abs(prob-correct));
if maxDiff > EPSILON then 
   print "--- ERROR in Test 5 ---", maxDiff prob correct;
else 
   print "--- Test 5 passes ---";

/* Test 6: Highly Correlated Case (Equicorrelation)
   Test stability with high correlation (rho=0.9) using the Rank-1 formula as benchmark.
*/
rho_val = 0.9;
Sigma = j(3,3,rho_val); 
Sigma[{1 5 9}] = 1;
b = {1.0 1.0 1.0};
/* Benchmark using rank-1 integral where lambda = sqrt(0.9) */
g_lambda = j(3,1,sqrt(rho_val)); 
g_b = colvec(b);
call quad(correct, "Rank1Integrand", {0 1});
prob = cdftvn(b, Sigma);
/* Quadrature near rho=1 is more challenging; check for reasonable precision */
maxDiff = max(abs(prob-correct));
if maxDiff > EPSILON then 
   print "--- ERROR in Test 6 ---", maxDiff prob correct;
else 
   print "--- Test 6 passes ---";

/* Test 7: 
   A. Sigma has negative correlations
   B. Largest magnitude is R[1,3]
   C. Use a noncentral MVN distribution */
mu ={ 0.8  0    1.65};
Sigma={14.3 -1.2 -4.4, 
       -1.2  5.2 -1.4, 
       -4.4 -1.4  9.1};
b  ={1.75  0.1 -1};
prob = cdftvn(b, Sigma, mu);

start MonteCarloEstimate(N, b, Sigma, mu={0 0 0});
   X = randnormal(N, mu, Sigma);
   inRegion = (X[,1] < b[1] & X[,2] < b[2] & X[,3] < b[3]);
   MC_Est = mean(inRegion);
   return MC_Est;
finish;

call randseed(1);
correct = MonteCarloEstimate(1E6, b, Sigma, mu={0 0 0}); /* Monte Carlo estimate = 0.022556 */
maxDiff = max(abs(prob-correct));
if maxDiff > 1E-3 then 
   print "--- ERROR in Test 7 ---", maxDiff prob correct;
else 
   print "--- Test 7 passes ---";

/* Test 8: Iterate over all combinations of signs for a correlation matrix.
   Largest magnitude is R[1,3]
*/
R0 = { 1    0.2  0.4, 
       0.2  1    0.1, 
       0.4  0.1  1   };
b  = {-1  0.05 -2};
/* there are 8 possible combinations for the signs of the 3 corr coefficients */
signs = { 1  1  1    1  1  1   1  1 1, 
          1 -1  1   -1  1  1   1  1 1, 
          1  1 -1    1  1  1  -1  1 1, 
          1  1  1    1  1 -1   1 -1 1, 
          1 -1 -1   -1  1  1  -1  1 1, 
          1 -1  1   -1  1 -1   1 -1 1, 
          1  1 -1    1  1 -1  -1 -1 1, 
          1 -1 -1   -1  1 -1  -1 -1 1 };
do i = 1 to nrow(signs);
   S = shape(signs[i,], 3, 3);
   R = R0 # S;
   prob = cdftvn(b, R);
   correct = MonteCarloEstimate(5E5, b, R);
   maxDiff = max(abs(prob-correct));
   if maxDiff > 1E-3 then 
      print "--- ERROR in Test 8 ---", maxDiff prob correct;
   else do;
      msg = cat("--- Test 8.",char(i,1)," passes ---");
      print (msg);
end;

print "--- DONE ---";