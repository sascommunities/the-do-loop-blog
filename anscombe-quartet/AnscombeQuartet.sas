/* SAS program to accompany the article 
   "Create your own version of Anscombe's quartet: 
          Dissimilar data that have similar statistics"
   by Rick Wicklin, published 17APR2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/04/17/create-version-of-anscombes-quartet.html

   This program shows how to create data that are similar to 
   Anscombe's quartet
   https://en.wikipedia.org/wiki/Anscombe%27s_quartet

   Create 
   1. Linear data with normal errors
   2. Quadratic data, no errors
   The remaining Anscombe data sets are
   3. Linear with one extreme outlier
   4. Constant X, varying Y, except for one outlier
*/

/* THE (original) Anscombe's quartet data:
   For all pairs (x, y1), (x, y2), (x, y3), and (x4, y4)
   the correlation between is 0.816	
   and the linear regression line is y = 3 + 0.5x
*/
data Anscombe;
input x y1 y2 y3 x4 y4;
datalines;
10.0    8.04    9.14    7.46    8.0     6.58
8.0     6.95    8.14    6.77    8.0     5.76
13.0    7.58    8.74    12.74   8.0     7.71
9.0     8.81    8.77    7.11    8.0     8.84
11.0    8.33    9.26    7.81    8.0     8.47
14.0    9.96    8.10    8.84    8.0     7.04
6.0     7.24    6.13    6.08    8.0     5.25
4.0     4.26    3.10    5.39    19.0    12.50
12.0    10.84   9.13    8.15    8.0     5.56
7.0     4.82    7.26    6.42    8.0     7.91
5.0     5.68    4.74    5.73    8.0     6.89
;

ods graphics / width=300px height=200px;
ods layout gridded columns=2 advance=table;
title "Anscombe's Quartet: Linear Data Set";
proc sgplot data=Anscombe noautolegend;
   scatter x=x y=y1;
   lineparm x=4 y=5 slope=0.5;
run;
title "Anscombe's Quartet: Quadratic Data Set";
proc sgplot data=Anscombe noautolegend;
   scatter x=x y=y2;
   lineparm x=4 y=5 slope=0.5;
run;
ods layout end;


proc iml;
/* 1. Create first data set (randomly) */
call randseed(12345);
x = T( do(4, 14, 0.2) );                              /* evenly spaced X */
eps = round( randfun(nrow(x), "Normal", 0, 1), 0.01); /* normal error */
y = 3 + 0.5*x + eps;                                  /* linear Y + error */
 
/* Helper function. Return paremater estimates for linear regression. Args are col vectors */
start LinearReg(Y, tX);
   X = j(nrow(tX), 1, 1) || tX;
   b = solve(X`*X, X`*Y);       /* solve normal equation */
   return b;
finish;
 
/* 2. compute target statistics. Other data sets have to match these */
targetB = LinearReg(y, x);          /* compute regression estimates */
targetCorr = corr(y||x)[2];         /* compute sample correlation */
print (targetB`||targetCorr)[c={'b0' 'b1' 'corr'} F=5.3 L="Target"];


/* Define system of simultaneous equations:
   https://blogs.sas.com/content/iml/2018/02/28/solve-system-nonlinear-equations-sas.html */
/* This function returns linear regression estimates (b0, b1) and correlation for a choice of beta */
start LinearFitCorr(beta) global(x);
   y2 = beta[1] + beta[2]*x + beta[3]*x##2;    /* construct quadratic Y */
   b = LinearReg(y2, x);      /* linear fit */
   corr = corr(y2||x)[2];     /* sample corr */
   return ( b` || corr);      /* return row vector */
finish;
 
/* This function returns the vector quantity (beta - target). 
   Find value that minimizes Sum | F_i(beta)-Target+i |^2 */
start Func(beta) global(targetB, targetCorr);
   target = rowvec(targetB) || targetCorr;
   G = LinearFitCorr(beta) - target;
   return( G );              /* return row vector */
finish;
 
/* 3. Let's solve for quadratic parameters so that same 
   linear fit and correlation occur */
beta0 = {-5 1 -0.1};         /* initial guess */
con = {.  .  .,              /* constraint matrix */
       0  .  0};             /* quadratic term is negative */
optn = ncol(beta0) || 0;     /* LS with 3 components || amount of printing */
/* minimize sum( beta[i] - target[i])**2 */
call nlphqn(rc, beta, "Func", beta0, optn) blc=con;  /* LS solution */
print beta[L="Optimal beta"];
 
/* How nearly does the solution solve the problem? Did we match the target values? */
Y2Stats = LinearFitCorr(beta);
print Y2Stats[c={'b0' 'b1' 'corr'} F=5.3];


/* 4. Visualize the new versions of the first two data set
      in Anscombe's quartet */
y2 = beta[1] + beta[2]*x + beta[3]*x##2;
create Anscombe2 var {x y y2};
append;
close;
QUIT;

ods layout gridded columns=2 advance=table;
proc sgplot data=Anscombe2 noautolegend;
   scatter x=X y=y;
   lineparm x=0 y=3.561 slope=0.447 / clip;
run;
proc sgplot data=Anscombe2 noautolegend;
   scatter x=x y=y2;
   lineparm x=0 y=3.561 slope=0.447 / clip;
run;
ods layout end;
