/* SAS program to accompany the article 
   "On the computation of the Wilcoxon signed rank statistic"
   by Rick Wicklin, published 19JUL2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/07/19/wilcoxon-signed-rank.html

   This article discusses the Wilcoxon signed rank test and explains
   why there are multiple test statistics that you can use to 
   implement the test.
*/

/* Wilcoxon Signed Rank Test 
   https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test

   Often used for matched data to test whether 
   pre/post data has significant difference in location,
   You take the difference (pre - post) and ask whether 
   the distribution of the difference has mean 0.
   This is the same test as the paired t test, but the
   signed rank test does not assume anything about the 
   distribution of the data (or the difference).
*/

/* data from the PROC TTEST documentation:
   https://go.documentation.sas.com/doc/en/pgmsascdc/v_021/statug/statug_ttest_examples03.htm
   Systolic blood pressures of 12 men before and after a treatment.
*/
data BP;
input SBP_Before SBP_After @@;
Diff = SBP_Before - SBP_After;
datalines;
120 128   124 131   130 131   118 127
140 132   128 125   140 141   135 137
126 118   130 132   126 129   127 135
;
 
proc ttest data=BP plots=none;
   paired SBP_Before*SBP_After;
run;
/* t test is not significant. No evidence that 
   the treatment significantly changes the systolic blood pressure.
*/

proc univariate data=BP mu0=0;
   var Diff;
   ods select TestsForLocation;
run;


/* NIST
https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/signrank.htm
"The signed rank test weakens the assumption of normality 
of the paired t-test to an assumption of symmetry. 
The signed rank test is more powerful than a sign test 
(it takes the magnitude of the differences into account as well as the sign), 
but it has stronger assumptions than the sign test. 
So if your data is at approximately symmetric, then the signed rank test 
is preferred to the sign test. However, if the symmetry assumption is not 
reasonable, the sign test is preferred."

The symmetry is only used for computing p-values b/c 
we assume each datum/difference has P=1/2 of being pos or negative.

This also explains why sometimes the test is called 
a test for a difference of means whereas most often it is called a test
for a difference of medians. For a symmetric distribution,
the mean equals the median.
*/

proc iml;
use BP; read all var "Diff" into y; close;
mu0 = 0;

/* test statistic for Wilcoxon Signed Rank Test
   according to formula from PROC UNIVARIATE documentation */
start SignedRankTest_SAS(y, mu0=0);
   x = y - mu0;
   x = x[loc(x ^= 0)];            /* 1. Exclude if x=mu0 */
   nT = nrow(x);                  /* count remaining obs */
   R = ranktie(abs(x));           /* 2. Rank abs(x) with ties=MEAN */
   idx = loc(x > 0);
   Tplus = sum(R[idx]);           /* 3. Sum ranks for x > mu0 */
   print x R (R#(R>=0))[L='R#Ind'] 
         (cusum(R#(x>=0)))[L='cusum(R#(x>0))']
         (cusum(R#(R>=0)))[L='cusum(R#Ind)']
         , Tplus (-nT*(nT+1)/4)[L='-nt.../4'];
   S = Tplus - nT*(nT+1)/4;       /* 4. Test statistic for PROC UNIVARIATE */
   return S;
finish;

/* compute test statistic used in PROC UNIVARIATE */
S = SignedRankTest_SAS(y, mu0);
print "Signed Rank Test: S = " S;

/* The signed rank formula in Wikipedia is a different statistic, 
   T = sum of the signed ranks */
start SignedRankTest(y, mu0=0);
   x = y - mu0;
   SR = j(nrow(x), 1, 0);         /* allocate vector; assign 0 if x=mu0 */
   idx = loc(x ^= 0);             /* 1. Exclude if x=mu0 */
   SR[idx] = sign(x[idx])#ranktie(abs(x[idx])); /* 2. Signed ranks of abs(x) */
   T = sum(SR);                   /* 3. test statistic is sum of signed ranks */
   print x SR (cusum(SR));
   return T;
finish;

T = SignedRankTest(y, mu0);
print "Signed Rank Test: T = " T;


/* Wiki goes on to say that T and Tplus are related:
   Tplus = T/2 + nT(nT+1)/4
   And because S = Tplus - nT*(nT+1)/4, we obtain
   S = T/2  or T = 2S
*/

/* Since the T and Tplus are related by an linear transformation,
   the distribution of one is merely a scaled and translated
  copy of the other. You can use either statistic to compute 
  things like p-values and confidence intervals.

  T is useful for some theoretical calculations because 
  E[T] = 0 under the null hypothesis, whereas
  E[Tplus] = n(n+1)/4. This is why the test statistic
  is S = Tplus - n(n+1)/4 and not Tplus itself.
*/
QUIT;
