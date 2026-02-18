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

call randseed(12345);

print "--- Starting Test Suite for CDFMVN ---";

/* 1. 4-D Identity Matrix */
test_name = "Test 1: 4-D Identity Matrix";
R = I(4);
b = {0 -1 -2 3};
correct = cdf("Normal", 0) * cdf("Normal", -1) * cdf("Normal", -2) * cdf("Normal", 3);
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);

/* 2. 5-D Identity Matrix */
test_name = "Test 2: 5-D Identity Matrix";
R = I(5);
b = {1 1 1 1 1};
correct = cdf("Normal", 1)**5;
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);

/* 3. 5-D Rank-1 Update (Equicorrelated) 
   Let R = 0.5*I + 0.5*11'. This is equivalent to rho=0.5 
   For rho=0.5, b=0, the prob is 1/(dim+1) = 1/6 = 0.16666...
*/
test_name = "Test 3: 5-D Equicorrelated (rho=0.5)";
v = j(5,1, sqrt(0.5));
R = 0.5*I(5) + v*v`;
b = {0 0 0 0 0};
correct = 1/6;
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);

/* 4. 5-D Min Matrix */
test_name = "Test 4: 5-D Min Matrix";
Sigma = {1 1 1 1 1, 
         1 2 2 2 2, 
         1 2 3 3 3, 
         1 2 3 4 4, 
         1 2 3 4 5};
b = 0:4;
correct = 0.4597946;
prob = cdfmvn(b, Sigma);
run check_test(test_name, prob, correct);

/* 5. 8-D Min Matrix */
test_name = "Test 5: 8-D Min Matrix";
Sigma = {1 1 1 1 1 1 1 1, 
         1 2 2 2 2 2 2 2, 
         1 2 3 3 3 3 3 3, 
         1 2 3 4 4 4 4 4, 
         1 2 3 4 5 5 5 5, 
         1 2 3 4 5 6 6 6, 
         1 2 3 4 5 6 7 7, 
         1 2 3 4 5 6 7 8};
b = 0:7;
correct = 0.4590496;
prob = cdfmvn(b, Sigma);
run check_test(test_name, prob, correct);

/* 6. A kxk equicorrelated matrix with rho=0.5 and b=0.
      Theoretical prob is 1/(k+1) 
*/
kk = {9, 12, 15};
do i=1 to nrow(kk); 
   k = kk[i];
   v = j(k,1,sqrt(0.5));
   R = 0.5*I(k) + v*v`;
   b = j(1,k,0);
   correct = 1/(k+1);
   prob = cdfmvn(b, R);
   /* for large matrices, reduce the desired precision */
   test_name = cat("Test 6: Equicorrelated (rho=0.5), dim=",strip(char(k,2)));
   run check_test(test_name, prob, correct, 0.01);
end;


/* 7. Rank-1 singular correlation matrix.
   If R is a matrix of all 1s, then X1=X2=...=Xn.
   P(X1 < 1, X2 < 2, ..., X8 < 8) = P(X1 < min(b)) = P(X1 < 1).
*/
test_name = "Test 7: 8-D Singular (All 1s)";
k = 8;
R = j(k,k,1);
b = 1:k;
correct = cdf("Normal", min(b));
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);


/* 8. Block-diagonal correlation matrix. If the blocks are 1x1 and 2x2,
      then the MVN probability is the product of univariate and bivariate probs.
*/
test_name = "Test 8: Block-Diagonal Matrix";
/* Unit Test: 5-D Block Diagonal Correlation Matrix */
R = { 1.0  0.5  0.0  0.0  0.0,
       0.5  1.0  0.0  0.0  0.0,
       0.0  0.0  1.0  0.0  0.0,
       0.0  0.0  0.0  1.0 -0.3,
       0.0  0.0  0.0 -0.3  1.0 };
b = {0.5  0.8  0.4  1.0 -0.2};
/* Calculate the correct value using the product of components */
p1 = probbnrm(0.5, 0.8, 0.5);   /* Bivariate Block 1 */
p2 = cdf("Normal", 0.4);       /* Univariate Block 2 */
p3 = probbnrm(1.0, -0.2, -0.3); /* Bivariate Block 3 */
correct = p1 * p2 * p3;
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);

/* 9. 3-D Negative correlation matrix and orthant probability */
test_name = "Test 9: 3-D Negative Equicorrelated";
R = {1.0 -0.6 -0.4, 
    -0.6  1.0 -0.2, 
    -0.4 -0.2 1.0};
b = {0 0 0};
/* Exact formula for q=3, b=0:
   P = 1/8 + 1/(4pi) * sum(arsin(rho_ij))
*/
rho_ij = R[{2 3 6}];
correct = 1/8 + 1/(4*constant("pi")) * sum(arsin( rho_ij ));
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);

/* 10. 4-D example from Bretz */
test_name = "Test 10: 4-D Example from Bretz";
R={1          0.7071068    0          0,
   0.7071068  1            0.5        0,
   0          0.5          1          0.3333333,
   0          0            0.3333333  1};
b={1 1 1 1};
correct = 0.583122;
prob = cdfmvn(b, R);
run check_test(test_name, prob, correct);


test_name = "Test 11: 4-D Diagonal Covariance+Mean";
Sigma = diag({1, 4, 9, 16});
mu    = {10 20 30 40};
b     = mu + rowvec(sqrt(vecdiag(Sigma)));  /* Evaluates at z=1 for all */
correct = cdf("Normal", 1)**4;
prob    = cdfmvn(b, Sigma, mu);
run check_test(test_name, prob, correct);


test_name = "Test 12: 4-D Block-Diagonal Covariance+Mean";
Sigma = { 2  1   0   0,
          1  2   0   0,
          0  0   3 -1.5,
          0  0 -1.5  3 };
mu    = {-5  5  0 10};
b     = mu + rowvec(sqrt(vecdiag(Sigma)));
/* Block 1 rho: 1/2 = 0.5. Block 2 rho: -1.5/3 = -0.5 */
correct = probbnrm(1, 1, 0.5) * probbnrm(1, 1, -0.5);
prob    = cdfmvn(b, Sigma, mu);
run check_test(test_name, prob, correct);


test_name = "Test 13: 5-D Block-Diagonal Covariance+Mean";
Sigma = { 4  2  0  0  0,
          2  4  0  0  0,
          0  0  9  0  0,
          0  0  0  5  4,
          0  0  0  4  5 };
mu    = {1  2  3  4  5};
b     = mu + rowvec(sqrt(vecdiag(Sigma)));
/* Block 1 rho: 0.5. Block 2: Univariate. Block 3 rho: 4/5 = 0.8 */
correct = probbnrm(1, 1, 0.5) * cdf("Normal", 1) * probbnrm(1, 1, 0.8);
prob    = cdfmvn(b, Sigma, mu);
run check_test(test_name, prob, correct);


test_name = "Test 14: 6-D Block-Diagonal Covariance+Mean";
Sigma = { 2  1.6  0  0  0   0,
          1.6 2   0  0  0   0,
          0  0    3  0  0   0,
          0  0    0  4  0   0,
          0  0    0  0  5 -2.5,
          0  0    0  0 -2.5  5 };
mu    = {10  20  30  40  50  60};
b     = mu + rowvec(sqrt(vecdiag(Sigma)));
/* Block 1 rho: 0.8. Blocks 2/3: Univariate. Block 4 rho: -0.5 */
correct = probbnrm(1, 1, 0.8) * cdf("Normal", 1)**2 * probbnrm(1, 1, -0.5);
prob    = cdfmvn(b, Sigma, mu);
call check_test(test_name, prob, correct);
 

/* PERFORMANCE */
   
/* TEST PERFORMANCE for kxk correlation matrix */
/*
kk = T(6:16);
elapsed_time = j(nrow(kk),1,.);
do i=1 to nrow(kk); 
   k = kk[i];
   v = j(k,1,sqrt(0.5));
   R = 0.5*I(k) + v*v`;
   b = j(1,k,0);
   correct = 1/(k+1);
   nRep = 3;
   t0 = time();
   do rep=1 to nRep;
      prob = cdfmvn(b, R);
   end;
   elapsed_time[i] = (time() - t0) / nRep;
end;
title "Performance of Orthant Probability in k Dimensions";
call series(kk, elapsed_time) grid={x y} label={"Dimension" "Elapsed Time"};
*/
/* evaluate the performance of 4-D qmc_eval function as (n,p) increase */
start TimePerfQMC(_=0);
q = 4;
b = {0 0 0 0}; 
C = {1 0 0 0, 
     1 1 0 0, 
     1 1 1 0, 
     1 1 1 1};
vec = 0:q-2;

/* The sequence of primes and the specific h-values for q=4 (row 3 of the mat matrix) */
p_vector = {157 313 619 1249 2503 5003 10007 20011};
h_vector = { 46  93 178  136  652 1476  2325  2894};

/* Define the grid of n and p values to test. 
   Using indices 4, 6, and 8 from p_vector (1249, 5003, and 20011) 
   and n values of 10, 30, and 50 to ensure measurable times. 
*/
n_values = {10 30 50}; 
test_indices = {4 6 8}; 

num_rows = ncol(n_values) * ncol(test_indices);
results = j(num_rows, 5, .);
row = 1;
do i = 1 to ncol(n_values);
    n = n_values[i];
    do j = 1 to ncol(test_indices);
        index = test_indices[j];
        p = p_vector[index];
        h = h_vector[index];
        z = mod(j(1, q-1, h)##vec, p); /* Compute the generator vector z */
        t0 = time(); 
           call qmc_eval(intval, varsum, n, p, z, q, b, C);
        elapsed = time() - t0;
        results[row, ] = n || p || (n*p) || elapsed || intval;
        row = row + 1;
    end;
end;

/* print the results */
cnames = {"n" "p" "Total_Evals" "Time (s)"};
print results[colname=cnames format=BEST10.];
finish;

*run TimePerfQMC();
/*
results = {
10 1249 12490 0.11199999 0.27341795 
10 5003 50030 0.44099998 0.27343176 
10 20011 200110 2.23399997 0.27343443 
30 1249 37470 0.37000012 0.27344875 
30 5003 150090 1.7329998 0.27343717 
30 20011 600330 6.21500015 0.27343613 
50 1249 62450 0.61999989 0.27343094 
50 5003 250150 2.43400002 0.27343614 
50 20011 1000550 10.694 0.27343765 
};
results = shape(results, 0, 5);
cnames = {"n" "p" "Total_Evals" "Time (s)"};
title "QMC Performance";
call scatter(results[,3], results[,4]) grid={x y} label={'Total Evaluations' 'Time (s)'};
*/

QUIT;