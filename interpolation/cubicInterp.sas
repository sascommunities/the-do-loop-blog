/* SAS program to accompany the article 
   "Cubic interpolation in SAS"
   by Rick Wicklin, published 06MAY2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/05/06/cubic-interpolation-sas.html

   This program shows how to perform cubic interpolation in SAS
   by using PROC IML.

   For linear regression, see the article 
   https://blogs.sas.com/content/iml/2020/05/04/linear-interpolation-sas.html
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


/*****************/
/* Cubic interpolating spline where the second derivative of the 
   interpolating curve is zero at both endpoints. 

   The interpolation is based on the values (x1,y1), (x2,y2),....
   The X  values must be nonmissing and in increasing order: x1 < x2 < ... < xn
   The values of the t vector are interpolated.
*/
proc iml;
start CubicInterp(x, y, t);
   d = dif(x, 1, 1);                     /* check that x[i+1] > x[i] */
   if any(d<=0) then stop "ERROR: x values must be nonmissing and strictly increasing.";
   idx = loc(t>=min(x) && t<=max(x));    /* check for valid scoring values */
   if ncol(idx)=0 then stop "ERROR: No values of t are inside the range of x.";

   /* fit the cubic model to the data */
   call splinec(splPred, coeff, endSlopes, x||y) smooth=0 type="zero";

   p = j(nrow(t)*ncol(t), 1, .);       /* allocate output (prediction) vector */
   call sortndx(ndx, colvec(t));       /* get sort index for t */
   sort_t = t[ndx];                    /* sorted version of t */
   sort_pred = splinev(coeff, sort_t); /* evaluate model at points of t */
   p[ndx] = sort_pred[,2];             /* "unsort" by using the inverse sort index */
   return( p );
finish;
store module=CubicInterp;
QUIT;

/* example of linear interpolation in SAS */
proc iml;
load module=CubicInterp;
use Points; read all var {'x' 'y'}; close;
use Score; read all var 't'; close;

pred = CubicInterp(x, y, t);
create PRED var {'t' 'pred'}; append; close;

/* OPTIONAL: Visualize the interpolant by drawing the entire curve */
s = do(min(x), max(x), 0.05);
curve = CubicInterp(x, y, s);
create CURVE var {'s' 'curve'}; append; close;
QUIT;

/* Visualize: concatenate the data and the predicted (interpolated) 
   values. Optionally, concatenate the curve. */
data All;
set Points Pred Curve;
run;

title "Cubic Spline Interpolation";
proc sgplot data=All noautolegend;
   series x=s y=curve;
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
load module=CubicInterp;
/* example of linear interpolation in SAS */
use APoints; read all var {'x' 'y'}; close;
use AScore; read all var 't'; close;

t0 = time();
pred = CubicInterp(x, y, t);
t = time() - t0;
print "Time for &NPts data points and &nInterp scoring points", t;
QUIT;
