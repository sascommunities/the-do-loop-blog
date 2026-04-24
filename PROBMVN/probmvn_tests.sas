proc iml;
load module=_all_;


   /*********************/
   /* call main program */
   /*********************/
   /*
   Read-Only Constants: 
   eps, sqtwpi, nl, plim, klim, minsmp, p, c, psqt, mxdim, mxhsum, bb.

   Read/Write State Variables: 
   hisum, olds, nn, covars, done, eone, infi, a, b, y.
   You must reset all state variables in mvndnt so that you can call mvn_dist multiple times.
   */
   hisum = .;
   olds = 0;
   mxdim = 80;
   mxhsum = 50;
   bb = 2;
   psqt={1.414213562373 1.732050807569 2.236067977500 2.645751311065 3.316624790355 3.605551275464 4.123105625618 4.358898943541 4.795831523313 5.385164807135 5.567764362830 6.082762530298 6.403124237433 6.557438524302 6.855654600401 7.280109889281 7.681145747869 7.810249675907 8.185352771872 8.426149773176 8.544003745318 8.888194417316 9.110433579144 9.433981132057 9.848857801796 10.04987562112 10.14889156509 10.34408043279 10.44030650891 10.63014581273 11.26942766958 11.44552314226 11.70469991072 11.78982612255 12.20655561573 12.28820572744 12.52996408614 12.76714533480 12.92284798332 13.15294643797 13.37908816026 13.45362404707 13.82027496109 13.89244398945 14.03566884762 14.10673597967 14.52583904633 14.93318452307 15.06651917332 15.13274595042 15.26433752247 15.45962483374 15.52417469626 15.84297951775 16.03121954188 16.21727474023 16.40121946686 16.46207763315 16.64331697709 16.76305461424 16.82260384126 17.11724276862 17.52141546794 17.63519208855 17.69180601295 17.80449381476 18.19340539866 18.35755975069 18.62793601020 18.68154169227 18.78829422806 18.94729532150 19.15724406067 19.31320791583 19.46792233393 19.57038579078 19.72308292332 19.92485884517 20.02498439450 20.22374841616};
   eps = 1e-10;
   sqtwpi = 2.506628274631000502415765284811045253;
   plim = 25;
   klim = 20;
   minsmp = 8;
   p = { 31 47 73 113 173 263 397 593 907 1361 2053 3079 4621 6947 10427 15641 23473 35221 52837 79259 118891 178349 267523 401287 601942};
   c = { 12 9 9 13 12 12 12 12 12 12 12 12 3 3 3 12 7 7 12, 13 11 17 10 15 15 15 15 15 15 22 15 15 6 6 6 15 15 9 , 27 28 10 11 11 20 11 11 28 13 13 28 13 13 13 14 14 14 14 , 35 27 27 36 22 29 29 20 45 5 5 5 21 21 21 21 21 21 21 , 64 66 28 28 44 44 55 67 10 10 10 10 10 10 38 38 10 10 10 , 111 42 54 118 20 31 31 72 17 94 14 14 11 14 14 14 94 10 10 , 163 154 83 43 82 92 150 59 76 76 47 11 11 100 131 116 116 116 116 , 246 189 242 102 250 250 102 250 280 118 196 118 191 215 121 121 49 49 49 , 347 402 322 418 215 220 339 339 339 337 218 315 315 315 315 167 167 167 167 , 505 220 601 644 612 160 206 206 206 422 134 518 134 134 518 652 382 206 158 , 794 325 960 528 247 247 338 366 847 753 753 236 334 334 461 711 652 381 381 , 1189 888 259 1082 725 811 636 965 497 497 1490 1490 392 1291 508 508 1291 1291 508 , 1763 1018 1500 432 1332 2203 126 2240 1719 1284 878 1983 266 266 266 266 747 747 127 , 2872 3233 1534 2941 2910 393 1796 919 446 919 919 1117 103 103 103 103 103 103 103 , 4309 3758 4034 1963 730 642 1502 2246 3834 1511 1102 1102 1522 1522 3427 3427 3928 915 915 , 6610 6977 1686 3819 2314 5647 3953 3614 5115 423 423 5408 7426 423 423 487 6227 2660 6227 , 9861 3647 4073 2535 3430 9865 2830 9328 4320 5913 10365 8272 3706 6186 7806 7806 7806 8610 2563 , 10327 7582 7124 8214 9600 10271 10193 10800 9086 2365 4409 13812 5661 9344 9344 10362 9344 9344 8585 , 19540 19926 11582 11113 24585 8726 17218 419 4918 4918 4918 15701 17710 4037 4037 15808 11401 19398 25950 , 34566 9579 12654 26856 37873 38806 29501 17271 3663 10763 18955 1298 26560 17132 17132 4753 4753 8713 18624 , 31929 49367 10982 3527 27066 13226 56010 18911 40574 20767 20767 9686 47603 47603 11736 11736 41601 12888 32948 , 40701 69087 77576 64590 39397 33179 10858 38935 43129 35468 35468 2196 61518 61518 27945 70975 70975 86478 86478 , 103650 125480 59978 46875 77172 83021 126904 14541 56299 43636 11655 52680 88549 29804 101894 113675 48040 113675 34987 , 165843 90647 59925 189541 67647 74795 68365 167485 143918 74912 167289 75517 8148 172106 126159 35867 35867 35867 121694 , 130365 236711 110235 125699 56483 93735 234469 60549 1291 93937 245291 196061 258647 162489 176631 204895 73353 172319 28881};

   /* Global constants for the Richtmyer generators and lattice rules */
   nl = 100;

   /* --- TEST SUITE --- */
/* TEMP: helper module to call mvn_dist with fewer parameters. Eventually this will 
   be the main module that users call.  */
start probmvn_mod(L, U, R);
   /* eventually we won't to pass these parameters because they are 
      derivable from other parameters. We can create them inside other functions,
      if needed. */
   n = ncol(R);
   infin = GetInfinityFlag(L, U);
   maxpts = 2000*n*n*n; 
   /* Eventually, we will make abstol a constant inside mvn_dist. 
      We will also releps because we only want to use abs errors */
   abseps = .0001; 
   releps = 0; 
   run mvn_dist(n, L, U, infin, R, maxpts, abseps, releps,
                error, value, nevals, inform );
   return(value);
finish;


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
prob = probmvn_mod(L, U, R0);
/* Compute P(X1 < 0, X2 > 0, X3 < 0) */
L = {.  0 .};
U = {0  . 0};
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
   test_name = "Test 15: 3-D Mixed Orthant, Version " + char(i);
   run check_test(test_name, prob, correct);
end;
