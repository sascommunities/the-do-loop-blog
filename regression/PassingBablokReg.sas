/* SAS program to accompany the article 
   "Passing-Bablok regression in SAS"
   by Rick Wicklin, published 14FEB2022 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2022/02/14/passing-bablok-regression-sas.html

   This program shows how to perform Passing-Bablok regression in SAS.
   The program uses the SAS/IML language to implement the method in

   Passing and Bablok (1983) "A New Biometrical Procedure for Testing the 
   Equality of Measurements from Two Different Analytical Methods"
   J. Clin. Chem. Clin. Biochem.
   Vol. 21, 1983, pp. 709-720
   https://edoc.hu-berlin.de/bitstream/handle/18452/11511/d.pdf?sequence=1

   For an overview of PB reg, see the documentation in the NCSS software:
   https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Passing-Bablok_Regression_for_Method_Comparison.pdf
*/

/* store SAS/IML helper functions */
%let INFTY = 1E12;         /* define a big number to use as "infinity" for vertical slopes */
proc iml;
/* Extract only the values X[i,j] where i < j.
   https://blogs.sas.com/content/iml/2012/08/16/extract-the-lower-triangular-elements-of-a-matrix.html */
   start StrictLowerTriangular(X);
   v = cusum( 1 || (ncol(X):2) );
   return( remove(vech(X), v) );
finish;
/* Compute pairwise differences:
   https://blogs.sas.com/content/iml/2011/02/09/computing-pairwise-differences-efficiently-of-course.html */ 
start PairwiseDiff(x);   /* x is a column vector */
   n = nrow(x);   
   m = shape(x, n, n);   /* duplicate data */
   diff = m` - m;        /* pairwise differences */
   return( StrictLowerTriangular(diff) );
finish;

/* Implement the Passing-Bablok (1983) regression method */
start PassingBablok(_x, _y, alpha=0.05, DEBUG=0);
   keepIdx = loc(_x^=. & _y^=.);      /* 1. remove missing values */
   x = _x[keepIdx]; y = _y[keepIdx];
   nObs = nrow(x);                    /* sample size of nonmissing values */

   /* 1. Compute all the pairwise slopes. For vertical, use large positive or negative values. */
   DX = T(PairwiseDiff(x));           /* x[i] - x[j] */
   DY = T(PairwiseDiff(y));           /* y[i] - y[j] */

   /* "Identical pairs of measurements with x_i=x_j and y_i=y_j
      do not contribute to the estimation of beta." (PB, p 711) */
   idx = loc(DX^=0 | DY^=0);          /* exclude 0/0 slopes (DX=0 and DY=0) */
   G = DX[idx] || DY[idx] || j(ncol(idx), 1, .); /* last col is DY/DX */
   idx = loc( G[,1] ^= 0 );           /* find where DX ^= 0 */
   G[idx,3] = G[idx,2] / G[idx,1];    /* form Sij = DY / DX for DX^=0 */

   /* Special cases:
      1) x_i=x_j, but y_i^=y_j   ==> Sij = +/- infinity
      2) x_i^=x_j, but y_i=y_j   ==> Sij = 0
      "The occurrence of any of these special cases has a probability 
      of zero (experimental data should exhibit these cases very rarely)."

      Passing and Bablock do not specify how to handle, but one way is
      to keep these values within the set of valid slopes. Use large 
      positive or negative numbers for +infty and -infty, respectively.
   */
   posIdx = loc(G[,1]=0 & G[,2]>0);    /* DX=0 & DY>0 */
   if ncol(posIdx) then G[posIdx,3] =  &INFTY;
   negIdx = loc(G[,1]=0 & G[,2]<0);    /* DX=0 & DY<0 */
   if ncol(negIdx) then G[negIdx,3] = -&INFTY;
   if DEBUG then print (T(1:nrow(G)) || G)[c={"N" "DX" "DY" "S"}];

   /* 2. Exclude slopes that are exactly -1 */
   /* We want to be able to interchange the X and Y variables. 
      "For reasons of symmetry, any Sij with a value of -1 is also discarded." */
   idx = loc(G[,3] ^= -1); 
   S = G[idx,3];
   N = nrow(S);                      /* number of valid slopes */
   call sort(S);                     /* sort slopes ==> pctl computations easier */
   if DEBUG then print (S`)[L='sort(S)'];

   /* 3. Estimate the slope as the shifted median of the pairwise slopes. */
   K = sum(S < -1);                  /* shift median by K values */
   b = choose(mod(N,2),              /* estimate of slope */
              S[(N+1)/2+K ],         /* if N odd, shifted median */
              mean(S[N/2+K+{0,1}])); /* if N even, average of adjacent values */
   /* 4. Estimate the intercept as the median of the values {yi – bxi}. */
   a = median(y - b*x);              /* estimate of the intercept, a */

   /* 5. Estimate confidence intervals. */
   C = quantile("Normal", 1-alpha/2) * sqrt(nObs*(nObs-1)*(2*nObs+5)/18);
   M1 = round((N - C)/2);
   M2 = N - M1 + 1;
   if DEBUG then print K C M1 M2 (M1+K)[L="M1+K"] (M2+k)[L="M1+K"];
   CI_b = {. .};                     /* initialize CI to missing */
   if 1 <= M1+K & M1+K <= N then CI_b[1] = S[M1+K];
   if 1 <= M2+K & M2+K <= N then CI_b[2] = S[M2+K];
   CI_a = median(y - CI_b[,2:1]@x);  /* aL=med(y-bU*x); aU=med(y-bL*x); */

   /* return a list that contains the results */
   outList = [#'a'=a, #'b'=b, #'CI_a'=CI_a, #'CI_b'=CI_b]; /* pack into list */
   return (outList);
finish;

store module=(StrictLowerTriangular PairwiseDiff PassingBablok);
quit;


/* data posted by @soulmate at
   https://communities.sas.com/t5/SAS-Communities-Library/An-example-of-Passing-Bablok-Regression-in-SAS/tac-p/792402
   which has MedCalc estimates 
Var       Estimate  Lower95 Upper95
Intercept 0.02860   0.02203 0.03832
Slope     0.9117    0.8971  0.9310
*/
data PBData;
label x="Method 1" y="Method 2";
input x y @@;
datalines;
0.07 0.07  0.04 0.10  0.07 0.08  0.05 0.10  0.06 0.08
0.06 0.07  0.03 0.05  0.04 0.05  0.05 0.08  0.05 0.10
0.12 0.21  0.10 0.17  0.03 0.11  0.18 0.24  0.29 0.33
0.05 0.05  0.04 0.05  0.11 0.21  0.38 0.36  0.04 0.05
0.03 0.06  0.06 0.07  0.46 0.27  0.08 0.05  0.02 0.09
0.12 0.23  0.16 0.25  0.07 0.10  0.07 0.07  0.04 0.06
0.11 0.23  0.23 0.21  0.15 0.20  0.05 0.07  0.04 0.14
0.05 0.11  0.17 0.26  8.43 8.12  1.86 1.47  0.17 0.27
8.29 7.71  5.98 5.26  0.57 0.43  0.07 0.10  0.06 0.11
0.21 0.24  8.29 7.84  0.28 0.23  8.25 7.79  0.18 0.25
0.17 0.22  0.15 0.21  8.33 7.67  5.58 5.05  0.09 0.14
1.00 3.72  0.09 0.13  3.88 3.54  0.51 0.42  4.09 3.64
0.89 0.71  5.61 5.18  4.52 4.15  0.09 0.12  0.05 0.05
0.05 0.05  0.04 0.07  0.05 0.07  0.09 0.11  0.03 0.05
2.27 2.08  1.50 1.21  5.05 4.46  0.22 0.25  2.13 1.93
0.05 0.08  4.09 3.61  1.46 1.13  1.20 0.97  0.02 0.05
1.00 1.00  3.39 3.11  1.00 0.05  2.07 1.83  6.68 6.06
3.00 2.97  0.06 0.09  7.17 6.55  1.00 1.00  2.00 0.05
2.91 2.45  3.92 3.36  1.00 1.00  7.20 6.88  6.42 5.83
2.38 2.04  1.97 1.76  4.72 4.37  1.64 1.25  5.48 4.77
5.54 5.16  0.02 0.06
;

%let DSName = PBData;    /* name of data set */
%let XVar = x;           /* name of X variable (Method 1) */
%let YVar = y;           /* name of Y variable (Method 2) */

ods graphics / width=480px height=480px;
title "Measurements from Two Methods";
proc sgplot data=&DSName;
   scatter x=&XVar y=&YVar;
   lineparm x=0 y=0 slope=1 / clip;    /* identity line: Intercept=0, Slope=1 */
   xaxis grid;
   yaxis grid;
run;

proc iml;
load module=(StrictLowerTriangular PairwiseDiff PassingBablok);
use &DSName; 
   read all var {"&XVar"} into x;
   read all var {"&YVar"} into y;
close;
R = PassingBablok(x, y);
ParameterEstimates = (R$'a' || R$'CI_a') // 
                     (R$'b' || R$'CI_b');
print ParameterEstimates[c={'Estimate' 'LCL' 'UCL'} r={'Intercept' 'Slope'}];
QUIT;

/* --------------------------------------------- */
/* For convenience, make a little macro          */
/* --------------------------------------------- */
%macro PBReg(DSName, XVar, YVar);
   proc iml;
   load module=(StrictLowerTriangular PairwiseDiff PassingBablok);
   use &DSName; 
      read all var {"&XVar"} into x;
      read all var {"&YVar"} into y;
   close;
   R = PassingBablok(x, y);
   ParameterEstimates = (R$'a' || R$'CI_a') // 
                        (R$'b' || R$'CI_b');
   print ParameterEstimates[c={'Estimate' 'LCL' 'UCL'} r={'Intercept' 'Slope'}];
   QUIT;
%mend;

/* --------------------------------------------- */
/* Show that PB regression is robust to outliers */
/* --------------------------------------------- */

%let DSName = PBOutlier;    /* name of data set */
data Outlier;
x = 24; y = 35;
run;
data &DSName;
set PBData Outlier;
run;

title "Add Extreme Outlier";
proc sgplot data=PBOutlier;
   scatter x=x y=y;
   lineparm x=0 y=0 slope=1 / clip;    /* identity line: Intercept=0, Slope=1 */
   xaxis grid;
   yaxis grid;
run;

%PBReg(PBOutlier, x, y);


/* --------------------------------------------- */
/* Show an example for which PB determines that the methods are equivalent */
/* --------------------------------------------- */

%let DSName = PBData2;    /* name of data set */
%let XVar = m1;           /* name of X variable (Method 1) */
%let YVar = m2;           /* name of Y variable (Method 2) */

data &DSName;
label x="Method 1" y="Method 2";
input m1 m2 @@;
datalines;
10.1 10.1  10.5 10.6  12.8 13.1   8.7  8.7  10.8 11.5 
10.6 10.4   9.6  9.9  11.3 11.3   7.7  7.5  11.5 11.9 
 7.9  7.8   9.9  9.7   9.3  9.9  16.3 16.2   9.4  9.8 
11.0 11.8   9.1  8.7  12.2 11.9  13.4 12.6   9.2  9.5 
 8.8  9.1   9.3  9.2   8.5  8.7   9.6  9.7  13.5 13.1 
 9.4  9.1   9.5  9.6  10.8 11.3  14.6 14.8   9.7  9.3 
16.4 16.5   8.1  8.6   8.3  8.1   9.5  9.6  20.3 20.4 
 8.6  8.5   9.1  8.8   8.8  8.7   9.2  9.9   8.1  8.1 
 9.2  9.0  17.0 17.0  10.2 10.0  10.0  9.8   6.5  6.6 
 7.9  7.6  11.3 10.7  14.2 14.1  11.9 12.7   9.9  9.4 
;

title "Equivalent Measurements";
proc sgplot data=&DSName;
   scatter x=&XVar y=&YVar;
   lineparm x=0 y=0 slope=1 / clip;
   xaxis grid;
   yaxis grid;
run;

%PBReg( &DSName, &XVar, &YVar );


/* -------------------------------------------------------- */
/* STOP HERE UNLESS YOU HAVE THE mcr PACKAGE in R INSTALLED */
/* -------------------------------------------------------- */


/* --------------------------------------------- */
/* Compare output with mcr package in R */
/* --------------------------------------------- */
title "Compare with mcr Package in R";

proc iml;
load module=(StrictLowerTriangular PairwiseDiff PassingBablok);

use &DSName; 
   read all var {"&XVar"} into x;
   read all var {"&YVar"} into y;
close;

R = PassingBablok(x, y);
ParameterEstimates = (R$'a' || R$'CI_a') // 
                     (R$'b' || R$'CI_b');
print ParameterEstimates[c={'Estimate' 'LCL' 'UCL'} r={'Intercept' 'Slope'}];

/* compare with R */
call ExportMatrixToR(x, "x"); call ExportMatrixToR(y, "y");
submit/ R;
#install.packages(mcr, dependencies = TRUE)
library(mcr)
m <- mcreg(x, y, method.reg="PaBa", method.ci="analytical", 
            slope.measure="radian", 
            mref.name="X",mtest.name="Y",na.rm=TRUE)
printSummary(m)
coef <- as.data.frame( getCoefficients(m) )
endsubmit;

call ImportMatrixFromR(RCoef, 'coef');
RResult = RCoef[, {1 3 4}];
Diff = ParameterEstimates - RResult;
print Diff[c={'Estimate' 'LCL' 'UCL'} r={'Intercept' 'Slope'}];

QUIT;

/* --------------------------------------------- */
/* Compare output with creatin data from mcr package in R */
/* --------------------------------------------- */

/* Data from the mcr package in R. From the doc:
"This data set gives the blood and serum preoperative creatinine measurements in 110 heart surgery
patients. ... Serum and plasma creatinin measurements in mg/dL
for each sample."
*/
data Creatin;
label x="Serum Creatin" y="Plasma Creatin";
input x y @@;
datalines;
0.82 0.79  1.83 1.62  1.39 1.36  0.81 1.3   1.72 1.88 
3.23 3.35  1.23 1.06  1.37 1.34  1.72 1.56  0.96 0.94 
0.76 0.69  1.15 1.16  0.95 0.71  1 0.83     1.52 1.44 
1.56 1.26  2.45 2.36  1.85 1.9   0.89 0.95  0.82 0.77 
1.01 0.82  0.96 0.88  3.38 3.42  1.31 1.15  1.1 0.96 
1.26 1.12  1.15 1.02  1.2 1.29   1.01 0.87  1.33 1.13 
0.87 0.82  1.66 1.33  1.41 1.14  1.08 0.93  1.17 1.09 
0.82 .     1.77 1.89  1.33 1.57  1.23 1.29  1.21 1.13 
1.2 1.11   1.49 1.61  0.89 0.79  1.09 1.19  1.03 0.96 
1.07 1.25  0.89 0.87  1.39 1.36  0.96 0.86  1.06 1.03 
1.17 0.86  0.9 0.8    1.39 1.22  0.91 0.86  2.28 2.25 
0.9 0.87   0.83 .     1.34 1.3   1.32 1.04  1.02 0.92 
1.04 1.01  1.27 1.21  1.08 1.11  1.1 1.21   1.39 1.54 
1.07 0.99  0.79 0.79  1.02 1.06  0.94 1.02  1.99 2.15 
0.7 0.56   0.84 1.01  0.66 0.59  1.1 1.17   0.9 0.82 
2.06 2.09  1.31 1.44  1.03 1.22  1.61 1.77  1.52 1.6 
1.11 1.13  0.85 1.06  1.2 1.38   1.09 1.24  1.27 1.29 
2.06 2.08  1.21 1.31  0.77 0.82  1.06 1.28  1.02 1.15 
1 1.04     1.47 1.51  0.97 1.13  1.7 1.84   1.78 1.94 
0.93 1.32  0.92 1.36  1.59 1.88  0.68 0.73  0.7 0.78 
1.08 1.01  0.83 0.9   1.17 1.3   0.98 0.96  0.77 0.61 
0.91 1.27  1.1 1.05   0.8 1.12   0.85 0.72  0.82 0.87 
;

title "Creatin Data";
proc sgplot data=Creatin noautolegend;
   scatter x=x y=y;
   lineparm x=0 y=0 slope=1 / clip;    /* identity line: Intercept=0, Slope=1 */
   xaxis grid;
   yaxis grid;
run;

proc iml;
submit / R;
library("mcr")
data(creatinine,package="mcr")
x <- creatinine$serum.crea
y <- creatinine$plasma.crea

# Passing-Bablok regression fit.
# The confidence intervals are calculated by using the analytical method.
model <- mcreg(x,y,method.reg="PaBa", method.ci="analytical", 
               slope.measure="radian", 
               mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE)
printSummary(model)
endsubmit;

/*
PASSING-BABLOK REGRESSION FIT:                                  
                                                                
                 EST SE        LCI       UCI                    
Intercept -0.1171729 NA -0.2001149 -0.020000                    
Slope      1.0880089 NA  1.0000000  1.173005                    
*/

call ImportMatrixFromR(x, 'x');
call ImportMatrixFromR(y, 'y');

load module=(StrictLowerTriangular PairwiseDiff PassingBablok);
R = PassingBablok(x, y);
ParameterEstimates = (R$'a' || R$'CI_a') // 
                     (R$'b' || R$'CI_b');
print ParameterEstimates[c={'Estimate' 'LCL' 'UCL'} r={'Intercept' 'Slope'}];
QUIT;


