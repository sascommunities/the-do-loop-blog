/* Analytical Tests for special cases of bivariate probabilities */
proc iml;
start cdfmvn(b,R);
   return probbnrm(b[1], b[2], R[1,2]);
finish;

/* Analytical Test 1: Orthant probabilities */
pi = constant('pi');
R = {1 0.5, 0.5 1};
b = {0 0};
prob = cdfmvn(b, R);
correct = 1/4 + arsin(R[1,2]) / (2*pi);
if abs(prob - correct) > 1E-14 then 
   print "ERROR: Orthant probability", rho prob correct;

/* Analytical Test 2: Diagonal matrices ==> Identity correlation */
b = {-1.2 0.7};
R = I(2);
prob = cdfmvn(b, R);
correct = cdf('normal',b[1]) # cdf('normal',b[2]);
if abs(prob - correct) > 1E-14 then 
   print "ERROR: Diagonal matrix", prob correct;

/* Comparison Test with PROBBNRM */
b = { 1.2 -0.7};
R = {1 -0.34, -0.34 1};
prob = cdfmvn(b, R);
correct = probbnrm(b[1], b[2], R[1,2]);
if abs(prob - correct) > 1E-14 then 
   print "ERROR: Diagonal matrix", prob correct;

print "--- DONE ---";

rho_vec = do(-0.9, 0.9, 0.1);
R = I(2);
/* Analytical Test 1: Orthant probabilities */
b = {0 0};
do i = 1 to ncol(rho_vec);
   rho = rho_vec[i];
   R[1,2] = rho;
   R[2,1] = R[1,2];
   prob = cdfmvn(b, R);
   correct = 1/4 + arsin(rho) / (2*pi);
   if abs(prob - correct) > 1E-14 then 
      print "ERROR: Orthant probability", rho prob correct;
end;

/* Analytical Test 2: Diagonal matrices ==> Identity correlation */
xy = expandgrid( -2:2, -2:2 );
R = I(2);
do j = 1 to nrow(xy);
   b = xy[j,];
   prob = cdfmvn(b, R);
   correct = cdf('normal',b[1]) # cdf('normal',b[2]);
   if abs(prob - correct) > 1E-14 then 
      print "ERROR: Diagonal matrix", prob correct;
end;

/* Comparison Test with PROBBNRM: Correlation {1 rho, rho 1} */
rho_vec = do(-0.9, 0.9, 0.1);
xy = expandgrid( -2:2, -2:2 );
R = I(2);
do j = 1 to nrow(xy);
   b = xy[j,];
   do i = 1 to ncol(rho_vec);
      rho = rho_vec[i];
      R[1,2] = rho;
      R[2,1] = R[1,2];
      prob = cdfmvn(b, R);
      correct = probbnrm(b[1],b[2], R[1,2]);
      if abs(prob - correct) > 1E-14 then 
         print "ERROR: Orthant probability" , rho pro correct;
   end;
end;

print "--- DONE ---";
