/* SAS program to accompany the article by Rick Wicklin on The DO Loop blog:
   "The QR algorithm for least-squares regression"
   published 12JUL2021 
   https://blogs.sas.com/content/iml/2021/07/12/qr-least-squares.html
   and
   "Compare computational methods for least squares regression"
   published 14JUL2021 
   https://blogs.sas.com/content/iml/2021/07/14/performance-ls-regression.html
*/

/*
   This program shows how to use various linear algebra routines to form and 
   solve the normal equations: (X`*X)*b = X`*y. The solution, b, is an estimate
   of the regression coefficients for least squares regression.

   This article solves the LS estimates by using the SOLVE, INV, and two kinds of QR calls.
*/

options ps=32000;

/* Create dummies as GLMSELECT design matrix */
proc glmselect data=Sashelp.Class outdesign=DesignMat;
   class Sex;
   model Weight = Height Sex Height*Sex/ selection=none;
   output out=GLMOut Pred=Pred;
run;

proc reg data=DesignMat plots=none;
   model Weight = Height Sex_F Height_Sex_F;
run;

proc iml;
use DesignMat;
   read all var {'Intercept' 'Height' 'Sex_F' 'Height_Sex_F'} into X;
   read all var {'Weight'} into Y;
close;

/* form normal equations */
XpX = X`*X;
Xpy = X`*y;

/* an efficient numerical solution */
b = solve(XpX, Xpy);
print b[L="SOLVE" F=D10.4];

/* an inefficient numerical solution */
Ainv = inv(XpX);  /* explicitly form the inverse matrix */
b = Ainv*Xpy;
print b[L="INV" F=D10.4];

/* similarly, it is inefficient to factor A=QR and solve
   the general equation R*b = Q`*c */
call QR(Q, R, piv, lindep, XpX);
c = Q`*Xpy;
b = trisolv(1, R, c, piv); /* equivalent to b = inv(R)*Q`*c */
print b[L="General QR" F=D10.4];

/* more efficient: solve the specific system with RHS c */
call QR(c, R, piv, lindep, XpX, ,Xpy);  /* c = Q`*Xpy */
b = trisolv(1, R, c, piv);
print b[L="Specific QR" F=D10.4];
print R, piv;

/* interestingly, you don't even need to form the normal equations.
   The QR algorithm can work with non-square and nonsymmetric matrices.
   So, you can work directly with the design matrix and the observed responses.
   However, this is very inefficient because the matrices are so much bigger
   than the normal equations. */
m = ncol(X);
call QR(Qty, R, piv, lindep, X, , y);
c = QTy[1:m];    /* we only need to first m rows of Q`*y */
b = trisolv(1, R, c, piv);
print b[L="Direct QR" F=D10.4];
QUIT;


/* SAS program to accompany the article by Rick Wicklin on The DO Loop blog:
   "Compare computational methods for least squares regression"
   published 14JUL2021 
   https://blogs.sas.com/content/iml/2021/07/14/performance-ls-regression.html

   This article compares the performance of five different algorithms.
   Four are based on the normal equations (SOLVE, INV, and two QR calls)
   and one uses QR to directly find a LS solution to the overdetermined system X`*b=y.
*/

proc iml;
/********************************/
/* simulation for large problem */
/********************************/
n = 1e5;                              /* number of observations */
nVars = {10 25 50 75 100 250 500};
call randseed(12345);

   /* TRICK: to eliminate the time required to load the DLLs, 
      solve the equations once before timing them */
   /*****************************************/
      numVars = nVars[nrow(nvars)];             /* number of variables */
      beta = randfun(numVars, "Uniform",-1,1);  /* regression coefficients */
      X = randfun(n||numVars, "Normal");
      X[,1] = 1;                                /* overwrite intercept column */
      eps = randfun(n, "Normal");               /* random errors */
      y = X*beta + eps;                         /* observed responses */
      XpX = X`*X;
      Xpy = X`*y;
      b = solve(XpX, Xpy);
      Ainv = inv(XpX);  
      b = Ainv*Xpy;
      call QR(Q, R, piv, lindep, XpX);
   /*****************************************/


NReps = 2;  /* repeat each computation this many times and record the average */
method=BlankStr(15); m=.; Time=.;
create Timing var {'Method' 'numVars' 'Time'};
do i = 1 to ncol(nVars);
   /* simulate regression data */
   numVars = nVars[i];                       /* number of variables */
   beta = randfun(numVars, "Uniform",-1,1);  /* regression coefficients */
   X = randfun(n||numVars, "Normal");
   X[,1] = 1;                                /* overwrite intercept column */
   eps = randfun(n, "Normal");               /* random errors */
   y = X*beta + eps;                         /* observed responses */

   /* time each method */
   method = "Normal Eqns";
   t0 = time();
   do k=1 to nReps;
      XpX = X`*X;
      Xpy = X`*y;
   end;
   Time = (time() - t0)/NReps;
   append;
   tNormal = Time;  /* store how long it takes to form normal eqns */

   method = "SOLVE";
   t0 = time();
   do k=1 to nReps;
      b = solve(XpX, Xpy);
   end;
   Time = (time() - t0)/NReps + tNormal;
   append;

   method = "INV";
   t0 = time();
   do k=1 to nReps;
      Ainv = inv(XpX);  
      b = Ainv*Xpy;
   end;
   Time = (time() - t0)/NReps + tNormal;
   append;

   method = "General QR";
   t0 = time();
   do k=1 to nReps;
      call QR(Q, R, piv, lindep, XpX);
      c = Q`*Xpy;
      b = trisolv(1, R, c, piv);
   end;
   Time = (time() - t0)/NReps + tNormal;
   append;

   method = "Specific QR";
   t0 = time();
   do k=1 to nReps;
      call QR(c, R, piv, lindep, XpX, ,Xpy);  /* c = Q`*Xpy */
      b = trisolv(1, R, c, piv);
   end;
   Time = (time() - t0)/NReps + tNormal;
   append;

   if numVars <= 100 then do;
   method = "Direct QR";
   t0 = time();
   do k=1 to nReps;
      m = ncol(X);
      call QR(Qty, R, piv, lindep, X, , y);
      b = trisolv(1, R, QTy[1:m], piv);
   end;
   Time = (time() - t0)/NReps;
   append;
   end;
end;
close;
QUIT;

title "Solving Matrix Equations";
title2 "Least-Squares Regression Estimates (n=100,000)";
proc sgplot data=Timing;
   series x=numVars y=Time / group=Method lineattrs=(thickness=2);
   keylegend / position=E sortorder=reverseauto;
   xaxis grid;
   yaxis grid max=10;
run;


proc print data=Timing;
run;
/*
proc print data=Timing noobs;
   where numVars=500;
run;
*/
proc print data=Timing noobs;
   where method="Direct QR";
run;

data Diff;
retain T0;
set Timing(where=(numVars=500));
if method="Normal Eqns" then do;
   T0 = Time;
   SolveTime=.;
end;
else SolveTime = Time - T0;
prop = SolveTime/Time;
run;
proc print data=Diff(drop=prop t0) noobs; 
run;


/**********************************/
/* For fun, look at PROC REG on a 100,000 x 500 data set.
   Of course, PROC REG does more than compute the estimates. It also 
   provides standard errors, p-values, and other estimates such as 
   sums-of-squares, MSE, F-values, R-squared stastistics, and more.
*/
/**********************************/
proc iml;
n = 1e5;                              /* number of observations */
numVars = 500;
call randseed(12345);
   beta = randfun(numVars, "Uniform",-1,1);  /* regression coefficients */
   X = randfun(n||numVars, "Normal");
   X[,1] = 1;                                /* overwrite intercept column */
   eps = randfun(n, "Normal");               /* random errors */
   y = X*beta + eps;                         /* observed responses */

   create In from X Y; append from X y; close;
QUIT;

proc reg data=In noprint plots=none outest=est;
model col501 = col2-col500;
quit;


