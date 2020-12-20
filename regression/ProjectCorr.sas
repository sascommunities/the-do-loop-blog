/* SAS programs to accompany the articles
   "Find a vector that has a specified correlation with another vector"
   by Rick Wicklin, published 17DEC2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/12/17/generate-correlated-vector.html

   and
 
   "Create a response variable that has a specified R-square value"
   by Rick Wicklin, published 21DEC2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/12/20/response-specified-r-square.html

   Given a vector, x, and a correlation, rho, find y such that
   corr(x,y) = rho.

   Recall that the correlation is the inner product between standardized 
   versions of the vectors. It is also equal to cos(theta), where theta
   is the angle between the vectors:
   https://blogs.sas.com/content/iml/2017/09/05/7-ways-view-correlation.html
*/

proc iml;

/* HELPER FUNCTIONS */
/* center a column vector by subtracting its mean */
start Center(v);
   return ( v - mean(v) );
finish;
/* create a unit vector in the direction of a column vector */
start UnitVec(v);
   return ( v / norm(v) );
finish;

/* Given x and rho, find a vector y such that corr(x,y) = rho.
   The initial guess can be almost any vector that 
   is not in span(x), orthog to span(x), and not in span(1) */
start CorrVec1(x, rho, guess);
   /* 1. Center the x and z vectors. Scale them to unit length. */
   u = UnitVec( Center(x)     );
   z = UnitVec( Center(guess) );

   /* 2. Project z onto the span(u) and the orthog complement of span(u) */
   w = (z`*u) * u;
   wPerp = z - w;       

   /* 3. The requirement that cos(theta)=rho results in a right triangle 
         where y (the hypotenuse) has unit length and the legs
         have lengths rho and sqrt(1-rho^2), respectively */
   v1 = rho * UnitVec(w);
   v2 = sqrt(1 - rho**2) * UnitVec(wPerp);
   y = v1 + v2;

   /* 4. Check the sign of y`*u. Flip the sign of y, if necessary */
   if sign(y`*u) ^= sign(rho) then 
      y = -y;
   return ( y );
finish;



/***************************************************/
/* First article and program                       */
/***************************************************/

/* Example: Call the CorrVec1 function */
x = {1,2,3};
rho = 0.543;
guess = {0, 1, -1};  
y = CorrVec1(x, rho, guess);
corr = corr(x||y);
print corr;

/* call for a wide variety of positive and negative correlations */
do i = 1 to 20;
  r = randfun(1, "Uniform", -1, 1);
  guess = randfun(nrow(x), "Normal");
  y = CorrVec1(x, r, guess);
  c = corr(x||y)[2];
  if abs(r-c)>1e-6 then print r c;
end;

/* because corr is a relationship between standardized vectors, 
   you can translate and scale Y any way you want */
y2 = 100 + 23*y;
corr = corr(x||y2);
print corr;

use sashelp.class;
   read all var {"Height"} into X;
close;
rho = 0.678;
call randseed(123);
guess = randfun(nrow(x), "Normal");
y = CorrVec1(x, rho, guess);

mean = 100;
std = 23*sqrt(nrow(x)-1);
v = mean + std*y;
title "Correlation = 0.678";
title2 "Random Normal Vector";
call scatter(X, v) grid={x y};

v = mean + std*y;
print (mean(v)) (std(v)) (var(v));
print (norm(std*y));

guess = randfun(nrow(x), "Expon");
y = CorrVec1(x, rho, guess);
title2 "Random Exponential Vector";
call scatter(X, mean + std*y);

mean = 100;
std = 23*sqrt(nrow(x)-1);
call randseed(123, 1);
r = {-0.75 -0.25 0.25 0.75};
create CorrOut var {'rho' 'x' 'y'};
do i = 1 to ncol(r);
   guess = randfun(nrow(x), "Normal");
   y = mean + std*CorrVec1(x, r[i], guess);
   rho = j(nrow(x), 1, r[i]);
   append;
end;
close;

submit;
   title "Random Vectors with Specified Correlations";
   proc sgpanel data=CorrOut;
      panelby rho / columns=2;
      scatter x=x y=y;
      rowaxis grid;
      colaxis grid;
   run;
endsubmit;


/***************************************************/
/* Second article and program                      */
/***************************************************/

/* You can compute a vector that has a certain angle with 
   a higher-dimensional subspace. Use OLS to do the projection of the 
   guess onto the span(X1, X2, ..., Xk). Call the projection w.
   Then find Y such that corr(Y, w) = rho

   1. Make almost any guess, z
   2. Project z onto span(1, X1, X2, X_k)
   3. Use previous function to find a vector correlated with z_hat
*/

/* Define or load the modules from 
   https://blogs.sas.com/content/iml/2020/12/17/generate-correlated-vector.html
*/

/* read some data X1, X2, ... into columns of a matrix, X */
use sashelp.class;
   read all var {"Height" "Weight" "Age"} into X;  /* read data into (X1,X2,X3) */
close;

/* Least-squares fit = Project Y onto span(1,X1,X2,...,Xk) */
start OLSPred(y, _x);
   X = j(nrow(_x), 1, 1) || _x;
   b = solve(X`*X, X`*y);
   yhat = X*b;
   return yhat;
finish;

/* specify the desired correlation between Y and \hat{Y}. Equiv: R-square = rho^2 */
rho = 0.543;        

call randseed(123, 1);
guess = randfun(nrow(X), "Normal");   /* 1. make random guess */
w = OLSPred(guess, X);                /* 2. w is in Span(1,X1,X2,...) */
Y = CorrVec1(w, rho, guess);          /* 3. Find Y such that corr(Y,w) = rho   */
/* optional: you can scale Y anyway you want ... */
/* in regression, R-square is squared correlation between Y and YHat */
corr = corr(Y||w)[2];
R2 = rho**2;                          
PRINT rho corr R2;


/* Write to a data set, then call PROC REG */
Z = Y || X;
create SimCorr from Z[c={Y X1 X2 X3}];
append from Z;
close;
QUIT;


proc reg data=SimCorr plots=none;
   model Y = X1 X2 X3;
   ods select FitStatistics ParameterEstimates;
quit;

