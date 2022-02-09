/* SAS program to accompany the article 
   "Deming regression for comparing different measurement methods"
   by Rick Wicklin, published 7JAN2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/01/07/deming-regression-sas.html

   This program shows how to perform Deming regression in SAS.

   A complete statistical presentation is available from 
   K. Linnet (1993), "Evaluation of Regression Procedures for Methods 
      Comparison Studies", CLIN.CHEM. 39(3),424-432 
      Appendix in http://clinchem.aaccjnls.org/content/clinchem/39/3/424.full.pdf

   -  or  -
   "Deming Regression" chapter in _NCSS Statistical Software_
   https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Deming_Regression.pdf
   -  or  -
   Maximum likelihood formula in Anders Christian Jensen  MethComp package
   link to 
   https://r-forge.r-project.org
   at end of Wikipedia article on Deming Regression:
   https://en.wikipedia.org/wiki/Deming_regression
*/

/* Demming regression */
data BloodTest;
label x="micrograms per deciliter" y="kiloOhms";
input x y @@;
subjid = _N_;
datalines;
169.0 45.5 130.8 33.4 109.0 23.8 94.1 19.8 86.3 20.4 78.4 18.7 
 76.1 16.1  72.2 16.7  70.0 11.9 69.8 14.6 69.5 10.6 68.7 12.7 67.3 16.9 
174.7 57.8 137.9 39.0 114.6 30.4 99.8 21.1 90.1 21.7 85.1 25.2
 80.7 20.6 78.1 19.3  77.8 20.9  76.0 18.2 77.8 18.3 74.2 15.7 73.1 13.9 
182.5 55.5 144.0 38.7 123.8 35.1 107.6 30.6 96.9 25.7 92.8 19.2 
 87.2 22.4  86.3 18.4  84.4 20.7  83.7 20.6 83.3 20.0 83.9 18.8 82.7 21.8 
160.8 49.9 122.7 32.2 102.6 19.2 86.6 14.7 76.1 16.6 69.6 18.8 
 66.7  7.4  64.4  8.2  63.0 15.5 61.7 13.7 61.2 9.2 62.4 12.0 58.4 15.2 
171.3 48.7 136.3 36.1 111.9 28.6 96.5 21.8 90.3 25.6 82.9 16.8 
 78.1 14.1  76.5 14.2  73.5 11.9 74.4 17.7 73.9 17.6 71.9 10.2 72.0 15.6 
;


title "Deming Regression";
title2 "Gold Standard (X) vs New Method (Y)";
proc sgplot data=BloodTest noautolegend;
   scatter x=x y=y;
   lineparm x=0 y=-10.56415 slope=0.354463 / clip; /* Deming regression estimates */
   xaxis grid label="Lab Test (micrograms per deciliter)";
   yaxis grid label="New Device (kiloohms)";
run;

/* Deming Regression */
/* See also 
   Njoya and Hemyari (2017) "Application of Deming Regression in 
            Molecular Diagnostics using a SAS Macro"
   https://www.lexjansen.com/pharmasug/2017/TT/PharmaSUG-2017-TT12.pdf

   Estimate linear model 
       y = b0 + b1*x
   when x and y both contain measurement error.
   Deal, Pate, and El Rouby handle the case where 
   lambda = var(y) / var(x) = 1, but you can estimate lambda from the data.

   https://en.wikipedia.org/wiki/Deming_regression
   https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Deming_Regression.pdf
*/
proc iml;
start Deming(XY, lambda=);
   /* Equations from https://en.wikipedia.org/wiki/Deming_regression */
   m = mean(XY);
   xMean = m[1]; yMean = m[2];
   S = cov(XY);
   Sxx = S[1,1]; Sxy = S[1,2]; Syy = S[2,2];
   /* if lambda is specified (eg, lambda=1), use it. Otherwise, estimate. */
   if IsEmpty(lambda) then
      delta = Sxx / Syy;        /* estimate of ratio of variance */
   else delta = lambda;
   c = Syy - delta*Sxx;
   b1 = (c + sqrt(c**2 + 4*delta*Sxy**2)) / (2*Sxy);
   b0 = yMean - b1*xMean;
   return (b0 || b1);
finish;

/* Test the program on the blood test data */
use BloodTest; read all var {x y} into XY; close;
b = Deming(XY);
print b[c={'Intercept' 'Slope'} L="Deming Regression"];
/*
   S = cov(XY);
   Sxx = S[1,1]; Sxy = S[1,2]; Syy = S[2,2];
   lambda = Sxx / Syy;        * estimate of ratio of variance;
print lambda;
b = Deming(XY, lambda);
print b[c={'Intercept' 'Slope'}];

b = Deming(XY, 1);
print b[c={'Intercept' 'Slope'}];
*/

b0 = b[1]; b1 = b[2];
x=XY[,1]; y = XY[,2];
lineStmt = "lineparm x=0 y=" + char(b0) + " slope=" + char(b1) + "/clip;";
call scatter(x, y) other=lineStmt grid={x y};

/* Helper modules for jackknife estimates of standard error and CI for parameters */
/* See https://blogs.sas.com/content/iml/2017/06/21/jackknife-estimate-standard-error-sas.html
*/
/* return the vector {1,2,...,i-1, i+1,...,n}, which excludes the scalar value i */ 
start SeqExclude(n,i);
   if i=1 then return 2:n;
   if i=n then return 1:n-1;
   return (1:i-1) || (i+1:n);
finish;
 
/* return the i_th jackknife sample for (n x p) matrix X */
start JackSamp(X,i);
   return X[ SeqExclude(nrow(X), i), ];  /* return data without i_th row */
finish;

/* 1. Compute T = statistic on original data */
T = b;

/* 2. Compute statistic on each leave-one-out jackknife sample */
n = nrow(XY);
T_LOO = j(n,2,.);             /* LOO = "Leave One Out" */
do i = 1 to n;
   J = JackSamp(XY,i);
   T_LOO[i,] = Deming(J); 
end;
 
/* 3. compute mean of the LOO statistics */
T_Avg = mean( T_LOO );  
 
/* 4. Compute jackknife estimates of standard error and CI */
stdErrJack = sqrt( (n-1)/n * (T_LOO - T_Avg)[##,] );
alpha = 0.05;
tinv = quantile("T", 1-alpha/2, n-2); /* use df=n-2 b/c both x and y are estimated */
Lower = T - tinv#stdErrJack;
Upper = T + tinv#stdErrJack;
result = T` || T_Avg` || stdErrJack` || Lower` || Upper`;
print result[F=7.3 c={"Estimate" "Mean Jackknife Estimate" "Std Error" 
             "Lower 95% CL" "Upper 95% CL"} r={'Intercept' 'Slope'}];


/* you might want to know the predicted (x,y) values for each x or y values */
delta = var(x) / var(y);        /* estimate of ratio of variance */
xPred = x + (b1/(b1**2 + delta))*(y - b0 - b1*x);
yPred = b0 + b1*xPred;
create DemingOut var {x y xPred yPred};
append;
close;

QUIT;

/* When the variance of X and Y are equal, then delta=1
   and Deming regression is euivalent to orthogonal regression
   which minimizes the sum of squared perpendicular distances from the 
   data points to the regression line. When the variances are not equal,
   the predicted X and Y values are not projections onto the line. */

proc sgplot data=DemingOut;
   label x="micrograms per deciliter" y="kiloOhms";
   scatter x=x y=y;
   lineparm x=0 y=-10.56415 slope=0.354463 / clip;
   *vector x=xPred y=yPred / xorigin=x yorigin=y;
   xaxis grid;
   yaxis grid;
run;

