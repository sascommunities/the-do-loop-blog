/* Unit tests for the common validation functions defined in mvn_validate.sas.
   Covers:
   1. IsSym   -- matrix symmetry check
   2. IsSPD   -- symmetric positive definite check
   3. IsCorr  -- correlation matrix check

   Test convention: a failing test prints a message beginning with "FAIL:".
   A passing test is silent. After each section, "End <section>" is printed.
   If the log contains no "FAIL:" lines, all tests passed.
*/
proc iml;
load module=_all_;

***************************************************;
print "Start unit tests for IsSym";
***************************************************;

* Non-square input -- expect 0 ;
S = 1:5;
t = IsSym(S);
if t^=0 then print "FAIL: IsSym should return 0 for a non-square matrix", S;

* Square but non-symmetric -- expect 0 ;
S = {4 2, 1 3};
t = IsSym(S);
if t^=0 then print "FAIL: IsSym should return 0 for a non-symmetric matrix", S;

* Symmetric matrix -- expect 1 ;
S = {4 2, 2 3};
t = IsSym(S);
if t^=1 then print "FAIL: IsSym should return 1 for a symmetric matrix", S;

print "End unit tests for IsSym";

***************************************************;
print "Start unit tests for IsSPD";
***************************************************;

* Symmetric but not positive definite -- expect 0 ;
R = {1.0 0.6 0.9,
     0.6 1.0 0.9,
     0.9 0.9 1.0};
t = IsSPD(R);
if t^=0 then print "FAIL: IsSPD should return 0 for a non-PD matrix", R;

* Valid SPD matrix -- expect 1 ;
R = {1.0 0.6 0.4,
     0.6 1.0 0.3,
     0.4 0.3 1.0};
t = IsSPD(R);
if t^=1 then print "FAIL: IsSPD should return 1 for a valid SPD matrix", R;

print "End unit tests for IsSPD";

***************************************************;
print "Start unit tests for IsCorr";
***************************************************;

* SPD but diagonal not all 1 -- expect 0 ;
S = {4 2, 2 3};
t = IsCorr(S);
if t^=0 then print "FAIL: IsCorr should return 0 when diagonal is not all 1", S;

* Valid correlation matrix -- expect 1 ;
S = {1.0 0.6 0.4,
     0.6 1.0 0.3,
     0.4 0.3 1.0};
t = IsCorr(S);
if t^=1 then print "FAIL: IsCorr should return 1 for a valid correlation matrix", S;

print "End unit tests for IsCorr";

***************************************************;
print "Start unit tests for Xform_Limits_Cov2Corr";
***************************************************;

Sigma = {1 1 1 1 1,
         1 2 2 2 2,
         1 2 3 3 3,
         1 2 3 4 4,
         1 2 3 4 5 };
mu = {1 2 3 4 5};
b = {-1 2 2 3 7};
/* standardize the parameters to the correlation scale */
c = Xform_Limits_Cov2Corr(b, Sigma, mu);
R = cov2corr(Sigma);

c_correct = {-2 0 -0.57735 -0.5 0.8944272 };
if max(abs(c - c_correct)) > 1e-6 then print "FAIL: Xform_Limits_Cov2Corr did not produce expected standardized limits", c;
R_correct = {1 0.7071068 0.5773503 0.5 0.4472136 ,
0.7071068 1 0.8164966 0.7071068 0.6324555 ,
0.5773503 0.8164966 1 0.8660254 0.7745967 ,
0.5 0.7071068 0.8660254 1 0.8944272 ,
0.4472136 0.6324555 0.7745967 0.8944272 1 
};
if max(abs(R - R_correct)) > 1e-6 then print "FAIL: Xform_Limits_Cov2Corr did not produce expected correlation matrix", R;
print "End unit tests for Xform_Limits_Cov2Corr";

QUIT;
