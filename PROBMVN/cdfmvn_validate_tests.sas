/* This file contins te following tests for validation of parameters used for computing the CDF:
   1. Test Sigma for symmetry
   2. Test Sigma for positive definiteness
   3. Test whether Sigma is a correlation matrix (versus a covariance)
   4. Test dimensions of arguments to cdftvn
   5. If the problem is specified on a covariance scale, test that the problem is correctly 
      converted to the correlation scale.
*/
proc iml;
load module=_all_;

***************************************************;
print "Start unit tests for validation of matrices";
***************************************************;
*************;
* IsSym tests;
*************;
S = 1:5;
t = IsSym(S);
if t^=0 then print "IsSym test failed", S;

S = {4 2, 1 3};
t = IsSym(S);  
if t^=0 then print "IsSym test failed", S;

S = {4 2, 2 3};
t = IsSym(S);  
if t^=1 then print "IsSym test failed", S;

*************;
* IsSPD tests;
*************;
R = {1.0 0.6 0.9,
     0.6 1.0 0.9,
     0.9 0.9 1.0};
t = IsSpd(R);  
if t^=0 then print "SPD test failed", R;

R = {1.0 0.6 0.4,
     0.6 1.0 0.3,
     0.4 0.3 1.0};
t = IsSpd(R);  
if t^=1 then print "SPD test failed", R;

**************;
* IsCorr tests;
**************;
S = {4 2, 2 3};
t = IsCorr(S);  
if t^=0 then print "IsCorr test failed", S;

S = {1.0 0.6 0.4,
     0.6 1.0 0.3,
     0.4 0.3 1.0};
t = IsCorr(S);  
if t^=1 then print "IsCorr test failed", S;

print "End unit tests";

***************************************************;
* unit tests for Programming Pattern Validation;
***************************************************;
print "Start unit tests for Parameter Validation";

* 1. Test Dimensional Compatibility ;
b = {1 1 1}; Sigma = j(2,2,1); mu = {0 0 0};
t = IsValidParmsTVN(b, Sigma, mu); * Mismatch b and Sigma ;
if t^=0 then print "Test Failed: Dimension mismatch b/Sigma not caught";

* 2. Test Mean (mu) Compatibility ;
b = {1 1 1}; Sigma = I(3); mu = {0 0};
t = IsValidParmsTVN(b, Sigma, mu); * Mismatch mu and Sigma ;
if t^=0 then print "Test Failed: Dimension mismatch mu/Sigma not caught";

* 3. Test Missing Values in b ;
b = {1 . 1}; Sigma = I(3); mu = {0 0 0};
t = IsValidParmsTVN(b, Sigma, mu);
if t^=0 then print "Test Failed: Missing value in b not caught";

* 4. Test Non-3D input for TVN specific check ;
b = {1 1}; Sigma = I(2); mu = {0 0};
t = IsValidParmsTVN(b, Sigma, mu);
if t^=0 then print "Test Failed: Non-3D input allowed in TVN check";

* 5. Test Mathematical failure propagation (Symmetry) ;
b = {1 1 1}; Sigma = {1 0.5 0.2, 0.1 1 0.3, 0.2 0.3 1}; mu = {0 0 0};
t = IsValidParmsTVN(b, Sigma, mu);
if t^=0 then print "Test Failed: Non-symmetric Sigma allowed";

* 6. Test Mathematical failure propagation (Non-SPD) ;
b = {1 1 1}; Sigma = {1 0.9 0.9, 0.9 1 0.1, 0.9 0.1 1}; mu = {0 0 0};
t = IsValidParmsTVN(b, Sigma, mu);
if t^=0 then print "Test Failed: Non-SPD Sigma allowed";

* 7. Test Valid Case ;
b = {0 0 0}; Sigma = {1 0.5 0.5, 0.5 1 0.5, 0.5 0.5 1}; mu = {0 0 0};
t = IsValidParmsTVN(b, Sigma, mu);
if t^=1 then print "Test Failed: Valid parameters rejected";
print "End unit tests";




* unit tests for Covariance to Correlation conversion;
print "Start unit tests for Convert_Cov_Limits";

* Test Case 1: Variances = {4, 9, 1}, Covariances such that rho = 0.5 ;
Sigma = {4 3 1,
         3 9 1.5,
         1 1.5 1};
b_orig = {2 3 1};   * This should scale to {1 1 1} since sigma = {2 3 1} ;
b = b_orig;
R = Sigma;
run Convert_Cov_Limits(b, R);
correct = {1 0.5 0.5, 0.5 1 0.5, 0.5 0.5 1}; * Expected R should have 0.5 in off-diagonals ;
if max(abs(R - correct)) > 1e-10 then 
   print "Test Failed: Correlation matrix R not correct", R;
correct = {1 1 1};    * Expected b should be {1 1 1} ;
if max(abs(b - correct)) > 1e-10 then 
   print "Test Failed: Scaled limits b not correct", b;

* Test Case 2: b is a matrix of row vectors ;
b_orig = {2 3 1, 4 6 2}; 
b = b_orig;
R = Sigma;
run Convert_Cov_Limits(b, R);
correct = {1 1 1, 2 2 2};
if max(abs(b - correct)) > 1e-10 then 
   print "Test Failed: Scaled limits b not correct", b;
print "End unit tests";
QUIT;

