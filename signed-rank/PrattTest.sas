/* SAS program to accompany the article 
   "Modifications of the Wilcoxon signed rank test and exact p-values"
   by Rick Wicklin, published 24JUL2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/07/24/exact-signed-rank-pratt.html

   This program shows how to download and call the %SIGNEDRANK macro.
   The macro is by 
   Leuchs AK, Neuhäuser M. (2010) "A SAS/IML algorithm for exact nonparametric paired tests."
   GMS Med Inform Biom Epidemiol. 6(1)
   Freely available from: http://www.egms.de/en/journals/mibe/2010-6/mibe000104.shtml

   This version on the GitHub site is slightly modified and reformatted by 
   Rick Wicklin 19JUL2023.
   See 
   https://blogs.sas.com/content/iml/2023/07/24/exact-signed-rank-pratt.html

SYNTAX:
%signedrank(data, var_diff, test, alternative, round=4);

*data: data set 
*var_diff: name of variable (difference) to evaluate 
*test (optional): Specification of the test 
    - 'Signed' for Wilcoxon's signed rank test 
    - 'Pratt' for modified test according to Pratt 
    - 'Original' for test based on original data 
    DEFAULT: TEST='Signed' 
*alternative (optional): Specification of the alternative 
    - 'less' one-sided test with H1: mu<0 
    - 'greater' one-sided test with H1: mu>0 
    - 'two' two-sided test with H1: mu ^= 0
    - 'all' all three tests 
    DEFAULT: ALTERNATIVE='two' 
*round (optional): specifies the number of decimal places to round on 
    (rounding only for the test based on original data) 
    DEFAULT: ROUND=4 
*/

/* define GitHub repo (source) and the local directory (destination) */
options dlcreatedir;                       /* create RepoPath directory if it doesn't exist */
%let gitURL = https://github.com/sascommunities/the-do-loop-blog/;  /* Git repo to copy */
%let RepoPath = %sysfunc(getoption(WORK))/BlogRepo;                 /* location to put copy */
 
/* Clone the GitHub repo into RepoPath; if exists, skip download */
data _null_;
if fileexist("&RepoPath.") then do;
   put 'Repository already exists; skipping the clone operation';
end;
else do;
   put "Cloning repository from &gitURL";
   /* NOTE: use GITFN_CLONE for 9.4M5; use GIT_CLONE for 9.4M6 and for Viya */
   rc = gitfn_clone("&gitURL", "&RepoPath." ); 
end;
run;
 
/* read the SAS program that defines the %SIGNEDRANK macro */
%include "&RepoPath/signed-rank/SignedRankMacro.sas";

options pagesize=32000;
/* data (Table 1) in Leuchs and Neuhäuser (2010) */
data Leucocytes;
input Leucocytes_Before Leucocytes_After @@;
Diff = Leucocytes_Before - Leucocytes_After;
label Leucocytes_Before="Number of Leucocytes/h (baseline)"
      Leucocytes_After ="Number of Leucocytes/h (after treatment)"
      Diff = "Difference (Baseline - After)";
datalines;
1.1 0.3   3.4 0.4   2.5 0.2   4.5 0.2   5.1 0.3
7.3 2.8   4.9 4.9   3.3 0.5   2.2 4.2   6.0 6.0
;

title "Signed Rank Test";
title2 "One-sided test: mu > 0";
%signedrank(Leucocytes, Diff, 'Signed', 'Greater');
title2 "Two-sided test: mu ^= 0";
%signedrank(Leucocytes, Diff); /* defaults='Signed' and 'two' */

title "Pratt's Modification of the Signed Rank Test";
title2 "One-sided test: mu > 0";
%signedrank(Leucocytes, Diff, 'Pratt', 'Greater');
title2 "Two-sided test: mu ^= 0";
%signedrank(Leucocytes, Diff, 'Pratt'); /* default='two' */


title "Signed Rank Test";
title2 "Two-sided test: mu ^= 0";
proc univariate data=Leucocytes mu0=0;
   var Diff;
   ods select TestsForLocation;
run;

/* a larger example: N=30, and nt=26 nonzero differences */
data Test;
input Before After @@;
Diff = Before - After; 
datalines;
12 13   12 13   13 13   11 13  15 13
14 13   13 13   14 14   15 13  16 18
13 15   14 16   17 18   12 14  15 14
14 16   15 17   16 18   14 15  14 17
17 17   12 15   14 16   13 18  18 17
13 15   13 12   14 13   16 12  15 13
;

proc univariate data=Test mu0=0;
   var Diff;
   ods select TestsForLocation;
run;

title "Signed Rank Test";
title2 "Two-sided test: mu ^= 0";
%signedrank(Test, Diff); /* defaults='Signed' and 'two' */

title "Pratt's Modification of the Signed Rank Test";
title2 "Two-sided test: mu ^= 0";
%signedrank(Test, Diff, 'Pratt'); /* default='two' */


