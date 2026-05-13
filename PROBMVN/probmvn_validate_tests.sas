/* Unit tests for PROBMVN-specific validation functions in probmvn_validate.sas.
   Tests for common functions (IsSym, IsSPD, IsCorr) are in mvn_validate_tests.sas.
   Covers:
   1. IsValidRectLimits    -- lower/upper limit vector consistency;
                              missing values (representing +/-Infinity) must be accepted
   2. IsValidParmsPROBMVN  -- full argument validation for PROBMVN_MOD

   Test convention: a failing test prints a message beginning with "FAIL:".
   A passing test is silent. After each section, "End <section>" is printed.
   If the log contains no "FAIL:" lines, all tests passed.
*/
proc iml;
load module=_all_;

***************************************************;
print "Start unit tests for IsValidRectLimits";
***************************************************;

* L and U with different lengths -- expect 0 ;
L = {-1 -1};
U = {1 1 1};
t = IsValidRectLimits(L, U);
if t^=0 then print "FAIL: IsValidRectLimits should return 0 when L and U have different lengths";

* L[i] > U[i] for one non-missing pair -- expect 0 ;
L = {-1 2 -1};
U = { 1 1  1};
t = IsValidRectLimits(L, U);
if t^=0 then print "FAIL: IsValidRectLimits should return 0 when L[i] > U[i]";

* L[i] = U[i] (point mass) is allowed -- expect 1 ;
L = {0 0 0};
U = {0 1 2};
t = IsValidRectLimits(L, U);
if t^=1 then print "FAIL: IsValidRectLimits should return 1 when L[i] = U[i] (point mass)";

* Missing values in L are allowed (represent -Infinity) -- expect 1 ;
L = {. -1 .};
U = {1  1 1};
t = IsValidRectLimits(L, U);
if t^=1 then print "FAIL: IsValidRectLimits should return 1 when L has missing values";

* Missing values in U are allowed (represent +Infinity) -- expect 1 ;
L = {-1 -1 -1};
U = { 1  . .};
t = IsValidRectLimits(L, U);
if t^=1 then print "FAIL: IsValidRectLimits should return 1 when U has missing values";

* Both L and U missing for a dimension -- expect 1 (entire real line) ;
L = {. -1};
U = {.  1};
t = IsValidRectLimits(L, U);
if t^=1 then print "FAIL: IsValidRectLimits should return 1 when both L[i] and U[i] are missing";

* All valid, no missing values -- expect 1 ;
L = {-2 -2 -2};
U = { 2  2  2};
t = IsValidRectLimits(L, U);
if t^=1 then print "FAIL: IsValidRectLimits should return 1 for a fully valid L/U pair";

print "End unit tests for IsValidRectLimits";

***************************************************;
print "Start unit tests for IsValidParmsPROBMVN";
***************************************************;

* Reference valid inputs used throughout this section ;
n = 3;
L_ok  = {-2 -2 -2};
U_ok  = { 2  2  2};
mu_ok = { 0  0  0};
S_ok  = {1.0 0.5 0.2,
         0.5 1.0 0.3,
         0.2 0.3 1.0};

* 1. Non-symmetric Sigma -- expect 0 ;
S_bad = {1.0 0.5 0.2,
         0.1 1.0 0.3,
         0.2 0.3 1.0};
t = IsValidParmsPROBMVN(L_ok, U_ok, S_bad, mu_ok);
if t^=0 then print "FAIL 1: Non-symmetric Sigma should be rejected";

* 2. Missing value in Sigma -- expect 0 ;
S_bad = S_ok;
S_bad[1,2] = .;
t = IsValidParmsPROBMVN(L_ok, U_ok, S_bad, mu_ok);
if t^=0 then print "FAIL 2: Missing value in Sigma should be rejected";

* 3. Missing value in mu -- expect 0 ;
mu_bad = {0 . 0};
t = IsValidParmsPROBMVN(L_ok, U_ok, S_ok, mu_bad);
if t^=0 then print "FAIL 3: Missing value in mu should be rejected";

* 4. L dimension mismatch with Sigma -- expect 0 ;
L_bad = {-2 -2};
t = IsValidParmsPROBMVN(L_bad, U_ok, S_ok, mu_ok);
if t^=0 then print "FAIL 4: L dimension mismatch with Sigma should be rejected";

* 5. U dimension mismatch with Sigma -- expect 0 ;
U_bad = {2 2};
t = IsValidParmsPROBMVN(L_ok, U_bad, S_ok, mu_ok);
if t^=0 then print "FAIL 5: U dimension mismatch with Sigma should be rejected";

* 6. mu dimension mismatch with Sigma -- expect 0 ;
mu_bad = {0 0};
t = IsValidParmsPROBMVN(L_ok, U_ok, S_ok, mu_bad);
if t^=0 then print "FAIL 6: mu dimension mismatch with Sigma should be rejected";

* 7. L[i] > U[i] for a non-missing pair -- expect 0 ;
L_bad = {-2 3 -2};
t = IsValidParmsPROBMVN(L_bad, U_ok, S_ok, mu_ok);
if t^=0 then print "FAIL 7: L[i] > U[i] should be rejected";

* 8. Non-positive-definite Sigma (singular) -- expect 0 ;
S_bad = {1.0 0.6 0.95,
         0.6 1.0 0.9,
         0.95 0.9 1.0};
t = IsValidParmsPROBMVN(L_ok, U_ok, S_bad, mu_ok);
if t^=0 then print "FAIL 8: Non-positive-definite Sigma should be rejected";

* 9. Dimension n > 100 -- expect 0 ;
n_big = 101;
L_big  = j(1, n_big, -2);
U_big  = j(1, n_big,  2);
S_big  = I(n_big);
mu_big = j(1, n_big,  0);
t = IsValidParmsPROBMVN(L_big, U_big, S_big, mu_big);
if t^=0 then print "FAIL 9: Dimension > 100 should be rejected";

* 10. Valid case: all non-missing, 3-D SPD Sigma -- expect 1 ;
t = IsValidParmsPROBMVN(L_ok, U_ok, S_ok, mu_ok);
if t^=1 then print "FAIL 10: Valid parameters should be accepted";

* 11. Valid case: missing values in L (represent -Infinity) -- expect 1 ;
L_miss = {. -2 .};
t = IsValidParmsPROBMVN(L_miss, U_ok, S_ok, mu_ok);
if t^=1 then print "FAIL 11: Missing values in L should be accepted";

* 12. Valid case: missing values in U (represent +Infinity) -- expect 1 ;
U_miss = { 2 . 2};
t = IsValidParmsPROBMVN(L_ok, U_miss, S_ok, mu_ok);
if t^=1 then print "FAIL 12: Missing values in U should be accepted";

* 13. Valid case: non-zero mean vector -- expect 1 ;
mu_nonzero = {1 -1 0.5};
t = IsValidParmsPROBMVN(L_ok, U_ok, S_ok, mu_nonzero);
if t^=1 then print "FAIL 13: Non-zero mu should be accepted with valid Sigma";

* 14. Valid case: 5-D identity covariance -- expect 1 ;
n5 = 5;
t = IsValidParmsPROBMVN(j(1,n5,-2), j(1,n5,2), I(n5), j(1,n5,0));
if t^=1 then print "FAIL 14: 5-D identity covariance should be accepted";

print "End unit tests for IsValidParmsPROBMVN";
QUIT;
