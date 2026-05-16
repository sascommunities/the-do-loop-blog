options ps=32000 nodate nonumber;

proc iml;
load module=_all_;

/* --- TEST SUITE --- */ 

/* Helper module to format test results */
start check_test(test_name, prob, correct, tol=0.001);
   maxDiff = max(abs(prob-correct));
   if maxDiff > tol then do;
      msg = cat("--- ",test_name, " FAILS ---");
      print msg[L=""], maxDiff prob correct;
   end;
   else do; 
      msg = cat("--- ",test_name, " passes ---");
      print msg[L=""];
   end;
finish;

/* Use Monte Carlo simulation to estimate the probability that a 
   MVN random variable is less than a specified value in each coordinate.
   Sigma is a kxk covariance matrix; mu is an option row vector with k elements.
   The row vector b specifies the upper limits of integration. 
   The function returns an estimate of 
   P(X1<b[1] & X2<b[2] & ... & Xk<b[k] | X~MVN(mu, Sigma))
   by simulating N random variates from MVN(mu, Sigma) and returning the proportion
   that are in the specified region.
*/
start MC_CDFMVN(N, b, Sigma, mu=j(1,ncol(Sigma),0));
   X = randnormal(N, mu, Sigma);
   inRegion = j(N,1,1);
   do i = 1 to ncol(b);
      v = (X[,i] < b[i]);
      inRegion = inRegion & v;
   end;
   return mean(inRegion);
finish;

/* Monte Carlo simulation of P( L < X < U | X~MVN(mu,Sigma) )
   where L and U are row vectors and a missing value 
   represents -Infinity in L and represents +Infinity in U */
start MC_PROBMVN(N, Lcov, Ucov, Sigma, mu=j(1,ncol(Sigma),0));
   /* standardize the parameters to the correlation scale */
   L = Xform_Limits_Cov2Corr(Lcov, Sigma, mu);
   U = Xform_Limits_Cov2Corr(Ucov, Sigma, mu);
   R = cov2corr(Sigma);
   X = randnormal(N, j(1,ncol(R),0), R);
   inRegion = j(N,1,1);
   /* Note: The '<' operator works correctly with a missing value
      on the left. The expression (a<y & y<b) is correct if a=.;
      However, if b=., you need to use (a<y). */
   do i = 1 to ncol(L);
      if U[i]=. then 
         v = (L[i] < X[,i]);
      else
         v = ((L[i] < X[,i]) & (X[,i] < U[i]));
      inRegion = inRegion & v;
   end;
   return mean(inRegion);
finish;

/* return a list with the MC est and a 95% CL. The list looks like
   [prob, lower95, upper95] */
start MC_PROBMVN_CL(N, L, U, Sigma, mu=j(1,ncol(Sigma),0));
   prob_MC = MC_PROBMVN(N, L, U, Sigma);
   SE_MC = sqrt( prob_MC * (1-prob_MC)/N );
   Lower95 = prob_MC - 1.96*SE_MC;
   Upper95 = prob_MC + 1.96*SE_MC;
   return( [prob_MC, Lower95, Upper95] );
finish;


   /* basic validation test */
call randseed(12345);
testName = "Test 0: 5-D Identity Matrix; [a,b]=[-2,2] in all coordinates";
R = i(5);
n = ncol(R);
lower = j(1,n,-2);
upper = j(1,n, 2);
prob = probmvn_mod(lower, upper, R);
correct = prod( probuvn_mod(colvec(lower), colvec(upper)));
run check_test(testName, prob, correct);


/***********************************************************/
/* FIRST, test with lower limits set to -Infinity (missing).
   We have many tests for the CDF. Reuse them. */
/***********************************************************/

/* 1. 4-D Identity Matrix */
testName = "Test 1: 4-D Identity Matrix; CDF";
R = I(4);
n = ncol(R);
upper = {0 -1 -2 3};
lower = j(1, n, .); 
prob = probmvn_mod(lower, upper, R);
correct = prod(cdf("Normal", upper));
run check_test(testName, prob, correct);


/* 2. 5-D Identity Matrix */
testName = "Test 2: 5-D Identity Matrix; CDF";
R = I(5);
upper = {1 1 1 1 1};
n = ncol(R);
lower = j(1, n, .); 
prob = probmvn_mod(lower, upper, R);
correct = prod(cdf("Normal", upper));
run check_test(testName, prob, correct);

/* 3. 5-D Rank-1 Update (Equicorrelated) 
   Let R = 0.5*I + 0.5*11'. This is equivalent to rho=0.5 
   For rho=0.5, b=0, the prob is 1/(dim+1) = 1/6 = 0.16666...
*/
testName = "Test 3: 5-D Equicorrelated (rho=0.5); CDF";
v = j(5,1, sqrt(0.5));
R = 0.5*I(5) + v*v`;
upper = {0 0 0 0 0};
n = ncol(R);
lower = j(1, n, .); 
prob = probmvn_mod(lower, upper, R);
correct = 1/6;
run check_test(testName, prob, correct);

/* 4. 5-D Min Matrix */
testName = "Test 4: 5-D Min Matrix; CDF";
Sigma = {1 1 1 1 1, 
         1 2 2 2 2, 
         1 2 3 3 3, 
         1 2 3 4 4, 
         1 2 3 4 5};
upper = 0:4;
n = ncol(Sigma);
lower = j(1, n, .); 
prob = probmvn_mod(lower, upper, Sigma);
correct = 0.4597946;
run check_test(testName, prob, correct);

/* 4a. 5-D Min Matrix G&B p. 5, Eqn 1.8 */
testName = "Test 4a: 5-D Min Matrix G & B p. 5, Eqn 1.8, Rectangular domain";
L = -1:-5;
U = 2:6;
prob = probmvn_mod(L, U, Sigma);
correct = 0.7615;
run check_test(testName, prob, correct);

/* 5. 8-D Min Matrix */
testName = "Test 5: 8-D Min Matrix; CDF";
Sigma = {1 1 1 1 1 1 1 1, 
         1 2 2 2 2 2 2 2, 
         1 2 3 3 3 3 3 3, 
         1 2 3 4 4 4 4 4, 
         1 2 3 4 5 5 5 5, 
         1 2 3 4 5 6 6 6, 
         1 2 3 4 5 6 7 7, 
         1 2 3 4 5 6 7 8};
upper = 0:7;
n = ncol(Sigma);
lower = j(1, n, .); 
prob = probmvn_mod(lower, upper, Sigma);
correct = 0.4590496;
run check_test(testName, prob, correct);

/* 5a: 8-D Min Matrix G&B p. 5, Eqn 1.9 */
testName = "Test 5a: 8-D Min Matrix G & B p. 5, Eqn 1.9, Rectangular domain";
L = -1:-8;
U = 2:9;
prob = probmvn_mod(L, U, Sigma);
correct = 0.7595;
run check_test(testName, prob, correct);

/* 6. A kxk equicorrelated matrix with rho=0.5 and b=0.
      Theoretical prob is 1/(k+1) 
*/
testName = "Test 6: Several equicorrelated matrices (rho=0.5); CDF";
kk = {9, 10, 15, 20};
kk = {9, 10, 15};
do i=1 to nrow(kk); 
   n = kk[i];
   print "--- Dimension =" n[L=""] "---";
   v = j(n,1,sqrt(0.5));
   R = 0.5*I(n) + v*v`;
   upper = j(1,n,0);
   lower = j(1, n, .); 
   prob = probmvn_mod(lower, upper, R);
   correct = 1/(n+1);
   run check_test(testName, prob, correct);
end;

/* 7. Rank-1 singular correlation matrix.
   If R is a matrix of all 1s, then X1=X2=...=Xn.
   P(X1 < 1, X2 < 2, ..., X8 < 8) = P(X1 < min(b)) = P(X1 < 1).
*/
testName = "Test 7: 8-D Singular (All 1s); CDF";
n = 8;
R = j(n,n,1);
upper = 1:n;
lower = j(1, n, .); 
prob = probmvn_mod(lower, upper, R);
correct = cdf("Normal", min(b));
run check_test(testName, prob, correct);

/* 8. Block-diagonal correlation matrix. If the blocks are 1x1 and 2x2,
      then the MVN probability is the product of univariate and bivariate probs.
*/
testName = "Test 8: Block-Diagonal Matrix; CDF";
/* Test 9: 5-D Block Diagonal Decomposition */
/* Construct the Block Diagonal Correlation Matrix */
R = { 1.0  0.5  0.0  0.0  0.0,
       0.5  1.0  0.0  0.0  0.0,
       0.0  0.0  1.0  0.0  0.0,
       0.0  0.0  0.0  1.0 -0.3,
       0.0  0.0  0.0 -0.3  1.0 };
n = ncol(R);
upper = {0.5  0.8  0.4  1.0 -0.2};
lower = j(1, n, .); 
prob = probmvn_mod(lower, upper, R);
/* Calculate the correct value using the product of components */
p1 = probbnrm(0.5, 0.8, 0.5);   /* Bivariate Block 1 */
p2 = cdf("Normal", 0.4);       /* Univariate Block 2 */
p3 = probbnrm(1.0, -0.2, -0.3); /* Bivariate Block 3 */
correct = p1 * p2 * p3;
run check_test(testName, prob, correct);

/* 9. 3-D Negative correlation matrix and orthant probability */
testName = "Test 9: 3-D Negative Equicorrelated; CDF";
R = {1.0 -0.6 -0.4, 
    -0.6  1.0 -0.2, 
    -0.4 -0.2 1.0};
upper = {0 0 0};
n = ncol(R);
lower = j(1, n, .); 
prob = probmvn_mod(lower, upper, R);
/* Exact formula for q=3, b=0:
   P = 1/8 + 1/(4pi) * sum(arsin(rho_ij))
*/
rho_ij = R[{2 3 6}];
correct = 1/8 + 1/(4*constant("pi")) * sum(arsin( rho_ij ));
run check_test(testName, prob, correct);


/***********************************************************/
/* SECOND, test with structured problems for which the exact answer is known. 
   10. Independent ractangles: If R=I(n), then the MVN probability is the 
      product of univariate probabilities.
   11. Singular correlation: If R is entirely 1s, it is a rank-1 matrix, so
      X_1 = X_2 = ... = X_n = Z \sim N(0,1). The probability that all 
      X_i fall within their respective (L_i, U_i) bounds is the 
      probability that Z falls in the tightest intersecting interval.
   12. Positive Orthant (Equicorrelated): For an equicorrelated matrix with 
      rho=0.5, we know the left-tailed orthant probability 
      P(X \le 0) = 1/(d+1). 
      By symmetry, -X \sim MVN(0, R). Therefore, the positive orthant 
      probability P(X > 0) is identical to P(X < 0). Test this by 
      setting L = 0 and U = missing.
   13. Block-Diagonal Matrix: The probability of 1-D and 2-D rectangular blocks
      can be computed exactly by using the probuvn_std function for 1-D blocks
      and using probbvn_std for 2-D blocks. The overall probability is the 
      product of the block probabilities.
   14. 3-D Mixed Orthant via Reflection: For a general 3-D correlation matrix, 
      compute P(X_1 < 0, X_2 > 0, X_3 < 0). This requires 
      L = {-Infty, 0, -Infty} and U = {0, Infty, 0}. 
      If you substitute Y_2 = -X_2, you can reframe this as a standard lower orthant problem: 
      P(X_1 < 0, Y_2 < 0, X_3 < 0). The distribution of this modified vector is MVN(0, R2), 
      where the signs of row 2 and column 2 in the correlation matrix have been flipped. 
      You can then use the known analytical 3-D arcsin formula on R2.
*/
/***********************************************************/
/* Test 10: 4-D Identity Matrix, Rectangular Limits */
test_name = "Test 10: 4-D Identity Matrix, Rectangular Limits";
R = I(4);
L = {-1.0  0.0 -0.5 -2.0};
U = { 1.0  2.0  0.5  0.0};
correct = prod( probuvn_mod(colvec(L), colvec(U)) );
prob = probmvn_mod(L, U, R);
run check_test(test_name, prob, correct);

/* Test 11: 8-D Singular (All 1s) Rectangular Limits */
test_name = "Test 11: 8-D Singular (All 1s), Rectangular Limits";
k = 8;
R = j(k,k,1);
L = {-1.0 -2.0 -0.5 -3.0 -4.0 -1.5 -2.5 -1.0};
U = { 2.0  1.5  3.0  2.5  1.0  4.0  3.5  1.0};
/* The probability reduces to P( max(L) < Z < min(U) ) */
correct = prod( probuvn_mod(max(L), min(U)) );
prob = probmvn_mod(L, U, R);
run check_test(test_name, prob, correct);

/* Test 12: 5-D Equicorrelated Positive Orthant (rho=0.5) */
test_name = "Test 12: 5-D Equicorrelated Positive Orthant";
k = 5;
v = j(k,1,sqrt(0.5));
R = 0.5*I(k) + v*v`;
L = j(1,k,0);
U = j(1,k,.); /* infinity for upper limits */
correct = 1/(k+1);
prob = probmvn_mod(L, U, R);
run check_test(test_name, prob, correct);

/* Test 13: 5-D Block-Diagonal, Rectangular Limits */
test_name = "Test 13: 5-D Block-Diagonal, Rectangular Limits";
R = { 1.0  0.5  0.0  0.0  0.0,
      0.5  1.0  0.0  0.0  0.0,
      0.0  0.0  1.0  0.0  0.0,
      0.0  0.0  0.0  1.0 -0.3,
      0.0  0.0  0.0 -0.3  1.0 };
L = {-1.0  0.0 -0.5 -1.0 -0.5};
U = { 1.0  1.0  0.5  2.0  1.5};
p1 = probbvn_std(L[,1:2], U[,1:2], 0.5);   * Bivariate Block 1 (Cols 1,2);
p2 = probuvn_mod(L[,3], U[,3]);            * Univariate Block 2 (Col 3);
p3 = probbvn_std(L[,4:5], U[,4:5], -0.3);  * Bivariate Block 3 (Cols 4,5);
correct = p1 * p2 * p3;
prob = probmvn_mod(L, U, R);
run check_test(test_name, prob, correct);

/* Test 14: 3-D Mixed Orthant with General Correlation */
test_name = "Test 14: 3-D Mixed Orthant";
/* Original general correlation matrix */
R0 = { 1.0  0.6  0.4,
       0.6  1.0  0.2,
       0.4  0.2  1.0 };
/* Compute P(X1 < 0, X2 > 0, X3 < 0) */
L = {.  0 .};
U = {0  . 0};
prob = probmvn_mod(L, U, R0);
/* By flipping the sign of X2, this becomes P(X1 < 0, Y2 < 0, X3 < 0).
   The new correlation matrix R_star flips signs on row 2 and col 2. */
rhos = R0[{2,3,6}];
rhos[1] = -rhos[1]; /* rho_12 */
rhos[3] = -rhos[3]; /* rho_23 */
correct = 1/8 + 1/(4*constant("pi")) * sum(arsin(rhos));
run check_test(test_name, prob, correct);

/* Test 15: Repeat previous test for all combinations of signs of the off-diagonal correlations. */
/* there are 8 possible combinations for the signs of the 3 corr coefficients */
/* Compute P(X1 < 0, X2 > 0, X3 < 0) */
L = {.  0 .};
U = {0  . 0};
signs = { 1  1  1    1  1  1   1  1 1, /* we already did first row */
          1 -1  1   -1  1  1   1  1 1, 
          1  1 -1    1  1  1  -1  1 1, 
          1  1  1    1  1 -1   1 -1 1, 
          1 -1 -1   -1  1  1  -1  1 1, 
          1 -1  1   -1  1 -1   1 -1 1, 
          1  1 -1    1  1 -1  -1 -1 1, 
          1 -1 -1   -1  1 -1  -1 -1 1 };
do i = 2 to nrow(signs);
   S = shape(signs[i,], 3, 3);
   R = R0 # S;
   prob = probmvn_mod(L, U, R);
   rhos = R[{2,3,6}];
   rhos[1] = -rhos[1]; /* rho_12 */
   rhos[3] = -rhos[3]; /* rho_23 */
   correct = 1/8 + 1/(4*constant("pi")) * sum(arsin(rhos));
   test_name = "Test 15: 3-D Mixed Orthant, Version " + char(i,2);
   run check_test(test_name, prob, correct);
end;

/* Test 16: Orthant probability in 4-D. See
   https://blogs.sas.com/content/iml/2026/05/18/4d-orthant-probability-mvn.html
*/
/* The base case: P2. For a scalar correlation, rho, returning the orthant probability, P2.
   See https://blogs.sas.com/content/iml/2026/05/11/mvn-orthant-probability.html */
start P2(rho);
  return( 0.25 + arsin(rho) / (2*constant("pi")) );
finish;
 
/* The sum of the I4 integrands. Instead of evaluating the three integrals separately,
   evaluate the sum of the integrals and perform one numerical integration.   
   Evaluates the Sun (1988) residual partial correlation matrix.
   Whereas Genz/Bretz Eqn 2.8-2.9 uses a general covariance matrix, this function assumes 
   R is a correlation matrix, so R[i,i]=1  */
start ProbInt4(x) global(R_global);
  R = R_global;
  pi = constant("pi");
  x2 = x # x;
  sum = 0;
 
  /* The sum of three terms: i=2, 3, and 4 */
  do i = 2 to 4;
    /* Determine the indices (j, k) of the two remaining variables */
    if i=2 then do;      j=3; k=4; end;
    else if i=3 then do; j=2; k=4; end;
    else do;             j=2; k=3; end;
 
    /* R_ii term: The quadratic denominator */
    w = 1 - R[1,i]##2 * x2;          /* R[i,i]=1 for all i */
    if w< 1E-12 then w = 1E-12;      /* Numerical safeguard */
 
    /* Compute the conditional variances and covariance. The are a_11, a_12, etc, in Sun and Asano */
    a_jj = 1      -R[1,j]##2 * x2   -((R[i,j] - R[1,i]*R[1,j]*x2)##2) / w;
    a_kk = 1      -R[1,k]##2 * x2   -((R[i,k] - R[1,i]*R[1,k]*x2)##2) / w;
    a_jk = R[j,k] -R[1,j]*R[1,k]*x2 -(R[i,j] - R[1,i]*R[1,j]*x2)*(R[i,k] - R[1,i]*R[1,k]*x2) / w;
 
    /* rho is the off-diagonal element of the conditional correlation matrix */
    rho = a_jk / sqrt(a_jj * a_kk);
 
    /* ensure rho is in [-1,1]. See
       https://blogs.sas.com/content/iml/2026/02/04/clip-values.html */
    rho = ( -1 <> (rho >< 1) ); /* clamp to [-1,1] */
 
    /* Call P2 function to get the I2 term */
    I2_val = 2 * pi * ( P2(rho) - 0.25 );     /* this actually simplifies to ARSIN(rho) :-) */
 
    /* Add the i-th term to the sum */
    sum = sum + (R[1,i] / sqrt(w)) * I2_val;
  end;
  return(sum);
finish;
 
/* The main P4 function, which evaluates the Childs-Sun formula for 4 dimensions */
start P4_Childs(R) global(R_global);
  R_global = R; /* Set the global correlation matrix for QUAD */
 
  /* Integrate ProbInt4 over the interval [0, 1] */
  call quad(I4, "ProbInt4", {0 1});
 
  /* Apply the Childs (1967) formula */
  pi = constant("pi");
  sum_arsin = sum( arsin(R[{2 3 4 7 8 12}]) );  /* sum of arcsine for all 6 off-diagonal elements */
  prob = 1/16 + sum_arsin / (8*pi) + I4 / (4 * pi**2);  /* fixes the typo in Genz and Bretz, Eqn 2.8 */
  return(prob);
finish;
 
/* --- Test all 16 orthant probabilities. We compare Childs formul in each octant to
   a call to probmvn_mod. See
   https://blogs.sas.com/content/iml/2026/05/11/mvn-orthant-probability.html
   https://blogs.sas.com/content/iml/2026/05/18/4d-orthant-probability-mvn.html
*/
R0 = {1.0  0.5  0.3  0.2,
      0.5  1.0  0.4  0.3,
      0.3  0.4  1.0  0.5,
      0.2  0.3  0.5  1.0};
L0 = {. . . .};
U0 = {0 0 0 0};
/* Generate all 16 combinations of sign vectors for 4-D orthants. See
   https://blogs.sas.com/content/iml/2026/06/01/generate-all-combinations-of-signs.html */
Signs = expandgrid({1 -1}, {1 -1}, {1 -1}, {1 -1});
do i = 1 to nrow(Signs);
   s = Signs[i,];
   R = s # R0 # s`;
   correct = P4_Childs(R);
   /* now reverse the limits of integration every location where there is a -1 sign */
   L = L0;
   U = U0;
   flip_idx = loc( s = -1 );
   if ncol(flip_idx) > 0 then do;
      L[flip_idx] = U0[ flip_idx ];
      U[flip_idx] = L0[ flip_idx ];
   end;
   prob = probmvn_mod(L, U, R0);
   test_name = "Test 16: 4x4 Orthant Probability (Childs, 1967): Orthant=" + putn(i,"f2.");
   run check_test(test_name, prob, correct);
end;

/* Test 17: A series of rectangular regions for 5-D.
   Include all possible ranges: 
   (-Infty, Infty), (-Infty, b), (a, Infty), (a, b) 
   The final test has all variables in (-Infty, Infty) so answer should be 1! */
test_name = "--- Test 17: 36 Different 5-D Rectangular Regions:";
print test_name[L=""] "only failures will be printed ---";
R = {1    -0.25  0.15 -0.35 -0.15 ,
    -0.25  1    -0.4   0.55  0.35 ,
     0.15 -0.4   1     0.05 -0.55 ,
    -0.35  0.55  0.05  1     0.1  ,
    -0.15  0.35 -0.55  0.1   1    };

L_Block = {.M .M .M .M .M,
           .M -1 -1 -1 .M,
           -2 -1 .M -1 -1,
           -3 .M -1 -2 -1,
           -1 -3 -2 -1 -2,
           -2 -1 -3 -2 -1  };
U_Block = {.I  1  0  1 .I,
            2  1 .I  1  0,
            3  2  1 .I  1,
            0  1  3  2  1,
            1  3  2  1  0,
           .I .I .I .I .I  };
N = 1E6;
do i = 1 to nrow(L_Block);
   L = L_Block[i,];
   do j = 1 to nrow(U_Block);
      U = U_Block[j,];
      prob = probmvn_mod(L, U, R);
      /* Compute the MC estimate; see if the QMC value is in the 95% CI */
      correct_list = MC_PROBMVN_CL(N, L, U, R);
      correct = correct_list$1;
      lower95 = correct_list$2;
      upper95 = correct_list$3;
      if prob < lower95 | prob > upper95 then do;
         run check_test(test_name, prob, correct);
         print L, U;
      end;
   end;
end;
print test_name[L=""] "DONE ---";
