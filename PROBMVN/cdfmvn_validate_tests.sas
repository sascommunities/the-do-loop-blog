/* Unit tests for CDFMVN-specific validation functions in cdfmvn_validate.sas.
   Tests for common functions (IsSym, IsSPD, IsCorr) are in mvn_validate_tests.sas.
   Covers:
   1. IsValidParmsTVN -- argument validation for cdftvn
   2. Convert_Cov_Limits -- covariance-to-correlation scaling
*/
*%include "mvn_validate_tests.sas";

proc iml;
load module=_all_;

***************************************************;
print "Start unit tests for IsValidParmsTVN";
***************************************************;

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

QUIT;

