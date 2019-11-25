/* SAS program to accompany the articles
   "Evaluate a quadratic polynomial in SAS" 
   by Rick Wicklin, published 25NOV2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/11/25/quadratic-polynomial-sas.html
   -and-
   "Evaluate a function on a linear subspace"
   by Rick Wicklin, published 27NOV2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/11/27/function-on-linear-subspace.html

   This program shows 
   1. How to use matrix computations to efficiently compute a quadratic polynomial.
   2. How to use ideas in vector calculus to restrict a function to a linear subspace.
*/


/* 1: Evaluate a quadratic polynomial by using matrix computations in SAS 
   "Evaluate a quadratic polynomial in SAS" 
   https://blogs.sas.com/content/iml/2019/11
*/
proc iml;
/* Evaluate  f(x) = 0.5 * x` * Q * x + L`*x + const, where
   Q is p x p symmetric matrix, 
   L and x are col vector with p elements.
   This version evaluates ONE vector x and returns a scalar value. */
start EvalQuad(x, Q, L, const=0);
   return 0.5 * x`*Q*x + L`*x + const;
finish;

/* compute Q and L for f(x,y)= 9*x##2 + x#y + 4*y##2) - 12*x - 4*y + 6 */
Q = {18 1,             /* matrix of second derivatives */
      1 8};
L = { -12, -4};        /* use column vectors for L and x */
const = 6;
x0 = {0, -1};          /* evaluate polynomial at this point */
f = EvalQuad(x0, Q, L, const);
print f[L="f(0,-1)"];
/* find function value at global minimum */
xx = {0.6433566, 0.4195804};
f = EvalQuad(xx, Q, L, const);
print f[L="f(opt)"];

/* Evaluate the quadratic function at each column of X and return a row vector. */
start EvalQuadVec(X, Q, L, const=0);
   f = j(1, ncol(X), .);
   do i = 1 to ncol(X);
      v = X[,i];
      f[i] = 0.5 * v`*Q*v + L`*v + const;
   end;
   return f;
finish;

/*    X1  X2  X3 X4  X5  */
vx = {-1 -0.5 0  0.5 1 ,
      -3 -2  -1  0   1 };
f = EvalQuadVec(vx, Q, L, const=0);
print (vx // f)[r={'x' 'y' 'f(x,y)'} c=('X1':'X5')];


x = do(-1, 1, 0.1); 
y = do(-3, 1.5, 0.1);
xy = expandgrid(x, y);             /* 966 x 2 matrix */
f = EvalQuadVec(xy`, Q, L, const); /* evaluate polynomial at all points */

/* write results to a SAS data set and visualize the function by using a heat map */
M = xy || f`;
create Heatmap from M[c={'x' 'y' 'f'}];  append from M;  close;

/***************************************************************/
/******************************************************************/
/* Can we eliminate the loop and use only matrix computations?
   Yes, but it is not efficient! 
   Revise the "EvalQuad" function so that X can be a matrix,
   where each column is a vector at which to evaluate the 
   quadratic form. */

/* Evaluate the quadratic form at each column of X and return a row vector. */
start EvalQuadMat(X, Q, L, const=0);
   M = 0.5 * X`*Q*X;
   /* M is an n x n matrix. 
      Extract diagonal elements of M: x1`*Q*x1, x2`*Q*x2, etc */
   diag = vecdiag(M);     
   /* diag is a column vector. Return a row vector. */
   return diag` + L`*x + const;
finish;

/* test performance: compare the matrix method 
   (which requires computing unneeded cells of the form xi`*Q*xj) 
   and the loop method. 
*/
call randseed(1);
k = 10000;
p = 200;
Q = toeplitz(p:1) / p;
L = T(1:p)/p;
X = j(p, k);
call randgen(X, "Uniform");
t0 = time();
   f0 = EvalQuadMat(X, Q, L, const);
tMat = time() - t0;
t0 = time();
   f = EvalQuadVec(X, Q, L, const);
tLoop = time() - t0;
print tMat tLoop (max(abs(f-f0)))[L="MaxDiff"]; 
/******************************************************************/
/******************************************************************/
QUIT;

data optimal;
   xx=0.64; yy=0.42;  /* optional: add the optimal (x,y) value */
run;
data All;  set Heatmap optimal; run;

title "Heat Map of Quadratic Function";
proc sgplot data=All;
   heatmapparm x=x y=y colorresponse=f / colormodel= (WHITE CYAN YELLOW RED BLACK);
   scatter x=xx y=yy / markerattrs=(symbol=StarFilled);
   xaxis offsetmin=0 offsetmax=0;
   yaxis offsetmin=0 offsetmax=0;
run;

/*********************************************************/
/* "Evaluate a function on a linear subspace"
   https://blogs.sas.com/content/iml/2019/11/27/function-on-linear-subspace.html
*/
PROC IML;
/* Evaluate the quadratic function at each column of X and return a row vector. */
start Func(XY);
   x = XY[1,]; y = XY[2,];
   return (9*x##2  + x#y + 4*y##2) - 12*x - 4*y + 6;
finish;

x0 = {0, -1};          /* evaluate polynomial at this point */
d  = {3, 1};           /* vector that determines direction */
u  = d / norm(d);      /* unit vector (direction) */

t = do(-2, 2, 0.1);    /* parameter values */
v = x0 + t @ u;        /* evenly spaced points along the linear subspace */
f = Func(v);           /* evaluate the 2-D function along the 1-D subspace */
title "Quadratic Form Evaluated on Linear Subspace";
call series(t, f) grid={x y};

tMin = 0.92;           /* t* = parameter near the minimum */
vMin = x0 + tMin * u;  /* corresponding (x(t*), y(t*)) */
fMin = Func(vMin);     /* f(x(t*), y(t*)) */
print tMin vMin fMin;

/* graph the 2-D and restricted function */
x = do(-2,    2,   0.1);    nx = ncol(x);
y = do(-1.7, -0.3, 0.1);    ny = ncol(y);
xy = expandgrid(x, y);
z = Func(xy`);

M = xy || z`;
create Heatmap from M[c={'x' 'y' 'z'}];
append from M;
close;

/* axis and tick marks */
t = do(-1.5, 2, 0.5);
vx = x0 + t @ u;
v = t` || vx`;
create ParamLine from v[c={'t' 'xt' 'yt'}];
append from v;
close;
QUIT;

data All;
set Heatmap ParamLine;
textX = xt + 0.05;    /* adjust position of tick values */
textY = yt - 0.05;
run;

title "Quadratic Function and Linear Subspace";
proc sgplot data=All noautolegend;
heatmapparm x=x y=y colorresponse=z / 
      colormodel= (WHITE CYAN YELLOW RED BLACK);
series  x=xt y=yt / markers markerattrs=(symbol=X);
text x=textX y=textY text=t / rotate=35 position=center textattrs=(size=12) strip;
xaxis offsetmin=0 offsetmax=0;
yaxis offsetmin=0 offsetmax=0;
gradlegend;
run;

/*********************************************************************/
/*********************************************************************/

/* Appendix: Verify that Q can be larger than 2x2 
   and you can still use the matrix operations */
proc iml;
/* Evaluate the quadratic form at each column of X and return a row vector. */
start EvalQuadVec(X, Q, L, const=0);
   f = j(1, ncol(X), .);
   do i = 1 to ncol(X);
      v = X[,i];
      f[i] = 0.5 * v`*Q*v + L`*v + const;
   end;
   return f;
finish;

Q = { 0.08 -0.05 -0.05 -0.05,
     -0.05  0.16 -0.02 -0.02,
     -0.05 -0.02  0.35  0.06,
     -0.05 -0.02  0.06  0.35}; 
c = j(nrow(Q), 1, 0);
v = do(0, 1, 0.5);
x = expandgrid(v, v, v, v);
f = EvalQuadVec(x`, Q, c);
M = x || f`;
*print M[c={x1 x2 x3 x4 f}];
QUIT;
