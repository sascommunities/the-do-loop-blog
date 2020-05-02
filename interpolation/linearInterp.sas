/* SAS program to accompany the article 
   "Linear interpolation in SAS"
   by Rick Wicklin, published 04MAY2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/05/04/linear-interpolation-sas.html

   This program shows how to perform linear interpolation in SAS.
   The primary (recommended) method is to use PROC IML.
   If you do not have PROC IML, I show how you can use
   PROC EXPAND in SAS/ETS or PROC TRANSREG in SAS/STAT software,
   but each of these procedures have shortcomings that are 
   discussed in the article.

   An older (less efficient) version of linear interpolation is available at
   https://blogs.sas.com/content/iml/2012/03/16/linear-interpolation-in-sas.html
*/

/* Example dats for 1-D interpolation */
data Points;  /* these points define the model */
input x y;
datalines;
0  1
1  3
4  5
5  4
7  6
8  3
10 3
;

data Score; /* these points are to be interpolated */
input t @@;
datalines;
2 -1 4.8 0 0.5 1 9 5.3 7.1 10.5 9
;

/* Linear interpolation based on the values (x1,y1), (x2,y2),....
   The X  values must be nonmissing and in increasing order: x1 < x2 < ... < xn
   The values of the t vector are linearly interpolated.

   The algorithm is vectorized, but the basics are:
   1. Find the first interval [x_i, x_{i+1}] such that x_i <= t <= x_{i+1}
   2. Let xL = x_i and xR = x_{i+1}
      Let yL = y_i and yR = y_{i+1}
   3. Let s = (t - xL) / (xR - xL) be the proportion of the interval to the left of t
   4. p = (1 - s)#yL + s#yR is the linear interpolation at t
*/
proc iml;
start LinInterp(x, y, _t);
   d = dif(x, 1, 1);                     /* check that x[i+1] > x[i] */
   if any(d<=0) then stop "ERROR: x values must be nonmissing and strictly increasing.";
   idx = loc(_t>=min(x) && _t<=max(x));  /* check for valid scoring values */
   if ncol(idx)=0 then stop "ERROR: No values of t are inside the range of x.";

   p = j(nrow(_t)*ncol(_t), 1, .);     /* allocate output (prediction) vector */
   t = _t[idx];                        /* subset t values inside range(x) */
   k = bin(t, x);                      /* find interval [x_i, x_{i+1}] that contains s */
   xL = x[k];   yL = y[k];             /* find (xL, yL) and (xR, yR) */
   xR = x[k+1]; yR = y[k+1];
   f = (t - xL) / (xR - xL);           /* f = fraction of interval [xL, xR] */
   p[idx] = (1 - f)#yL + f#yR;        /* interpolate between yL and yR */
   return( p );
finish;
store module=LinInterp;
QUIT;

/* example of linear interpolation in SAS */
proc iml;
load module=LinInterp;
use Points; read all var {'x' 'y'}; close;
use Score; read all var 't'; close;

pred = LinInterp(x, y, t);
create PRED var {'t' 'pred'}; append; close;
QUIT;

/* Visualize: concatenate data and predicted (interpolated) values */
data All;
set Points Pred;
run;

title "Linear Interpolation";
title2 "No Extrapolation";
proc sgplot data=All noautolegend;
   *dropline x=t y=Pred / dropto=both;  /* optional */
   series x=x y=y;
   scatter x=x y=y / markerattrs=(symbol=CircleFilled size=12)
                    name="data" legendlabel="Data";
   scatter x=t y=Pred / markerattrs=(symbol=asterisk size=12 color=red)
                    name="interp" legendlabel="Interpolated Values";
   fringe t / lineattrs=(color=red thickness=2)
                    name="score" legendlabel="Values to Score";
   xaxis grid values=(0 to 10) valueshint label="X";
   yaxis grid label="Y" offsetmin=0.05;
   keylegend "data" "score" "interp";
run;

/* ------------------------------- */
/* Demonstrate linear interpolation by using PROC TRANSREG or PROC EXPAND:
   SAS Usage Note: https://support.sas.com/kb/24/560.html
*/
data Combine;
set Points Score(rename=(t=x));
run;

proc sort data=Combine ; by x ;  run;
proc print noobs;run;
/* NOTE: PROC EXPAND does not allow duplicate values for X.
   It gives an ERROR for the COMBINE data. There are 
   programming ways to get rid of the dupliate values, 
   but I will manually remove them for this example
*/
data Combine2;   /* delete dup X */
input x y;
datalines;
-1.0 . 
0.0 1 
0.5 . 
1.0 3 
2.0 . 
4.0 5 
4.8 . 
5.0 4 
5.3 . 
7.0 6 
7.1 . 
8.0 3 
9.0 . 
10.0 3 
10.5 . 
;
proc expand data=Combine2 out=LinInterp;
   convert y=Pred / method=join;
   id x;
run;

title "Linear Interpolation by using PROC EXPAND";
proc sgplot data=LinInterp;
series x=x y=y;
scatter x=x y=Pred;
run;


/* You can also perform linear interpolation by using PROC TRANSREG.
   PROC TRANSREG requires knowing the number of data points,
   not counting the scoring data. Count only obs for which 
   X and Y are both nonmissing. */
data _null_;
set Combine END=EOF;
n + ^cmiss(x,y);
if EOF then 
   call symput('NumIntKnots', n-2);
run;
%put &=NumIntKnots;

proc transreg data=Combine;
   model identity(y)=spline(x / degree=1 nknots=&NumIntKnots);
   output out=LinInterp2 predicted;
run;

title "Linear Interpolation by using PROC TRANSREG";
title "Notice Extrapolation";
proc sgplot data=LinInterp2 noautolegend;
   series x=x y=y / markers;
   scatter x=x y=Py / markerattrs=(symbol=X size=12);
   xaxis grid;
   yaxis grid;
run;



/*******************************************/
/* investigate performance */
/*******************************************/
%let nPts = 1000;
%let nInterp = 1e6;
data APoints;
call streaminit(123);
do i = 1 to &nPts;
   x = rand("Uniform", 0, 100);
   y = x + round(rand("Normal"), 0.01);
   output;
end;
keep x y;
run;

proc sort data=APoints; by x; run;

data AScore; /* these points are to be interpolated */
call streaminit(54321);
do i = 1 to &nInterp;
   t = rand("Uniform", -1, 101);
   output;
end;
keep t;
run;

proc iml;
load module=LinInterp;
/* example of linear interpolation in SAS */
use APoints; read all var {'x' 'y'}; close;
use AScore; read all var 't'; close;

t0 = time();
pred = LinInterp(x, y, t);
t = time() - t0;
print "Time for &NPts data points and &nInterp scoring points", t;
QUIT;

