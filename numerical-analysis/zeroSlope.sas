/* SAS program to accompany the article 
   "Find points where a regression curve has zero slope"
   by Rick Wicklin, published 13MAY2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/05/13/regression-curve-zero-slope.html ‎

   This program shows how to find local minima and maxima for a regression
   curve by finding points at which the slope of the curve is zero.
   Examples are shown for 
   - PROC LOESS, which uses the SCORE statement to score a model
   - PROC GAMPL, which uses the missing value trick to score a model
   - A previous article
     https://blogs.sas.com/content/iml/2018/07/05/derivatives-nonparametric-regression.html
     shows how to use the STORE statement and PROC PLM for 
     parametric regression models, such as PROC GLM.   

   This example uses the Sashelp.ENSO data, which are distributed as 
   part of SAS/STAT software.
*/


/* 1. For clarity, create VIEW where x is the independent 
      variable and y is the response variable 
*/
data Have / view=Have;
set Sashelp.Enso(rename=(Month=x Pressure=y));
keep x y;
run;

/* 2. Put min and max into macro variables.
      Create a grid of points at which to evaluate the 
      regression curve and estimate derivatives 
*/
proc sql noprint;
  select min(x), max(x) into :min_x, :max_x 
  from Have;
quit;

data Grid;
dx = (&max_x - &min_x)/201;    /* choose the step size wisely */
do x = &min_x to &max_x by dx;   
   output;
end;
drop dx;
run;


/***************************************************************/
/* PROC LOESS supports a SCORE statement for scoring the model */
/***************************************************************/

/* 3. Score the model on the grid */
ods select none;
proc loess data=Have plots=none;
   model y = x;
   score data=Grid;  /* PROC LOESS does not support an OUT= option */
   /* Most procedures support an OUT= option to save the scored values.
      PROC LOESS displays the scored values in a table, so use ODS to
      save the table to an output data set */
   ods output ScoreResults=ScoreOut;
run;
ods select all;

/* 4. Compute slope by using finite difference formula. */
data Deriv0;
set ScoreOut;
Slope = dif(p_y) / dif(x);      /* (f(x) - f(x-dx)) / dx */
/* save previous values of x, y, and slope */
xPrev = lag(x);  yPrev = lag(p_y); SlopePrev = lag(Slope);
if n(SlopePrev) AND sign(SlopePrev) ^= sign(Slope) then do;
   /* The slope changes sign between this obs and the previous.
      Assuming linearity on the interval, find (t, f(t))
      where slope is exactly zero */
   t0 = xPrev - SlopePrev * (x - xPrev)/(Slope - SlopePrev); 
   /* use linear interpolation to find the corresponding y value:
      f(t) ~ y0 + (y1-y0)/(x1-x0) * (t - x0)       */
   f_t0 = yPrev + (yPrev - p_y)/(x - xPrev) * (t0 - xPrev);
   if sign(SlopePrev) > 0 then _Type_ = "Max";
   else _Type_ = "Min";
   output; 
end;
keep t0 f_t0 Slope _Type_;
label f_t0 = "f(t0)";
run;

proc print data=Deriv0 label;
run;

/* 5. (Optional) Visualize the points where the curve has local extrema */
data Combine;
merge Have                           /* data   : (x, y)     */
      ScoreOut(rename=(x=t p_y=p_t)) /* curve  : (t, p_t)   */
      Deriv0;                        /* extrema: (t0, f_t0) */
run;

title "Loess Smoother";
title2 "Red Markers Indicate Zero Slope for Smoother";
proc sgplot data=Combine noautolegend;
   scatter x=x y=y;
   series x=t y=p_t / lineattrs=GraphData2;
   scatter x=t0 y=f_t0 / markerattrs=(symbol=circlefilled color=red);
   yaxis grid; 
run;

/***************************************************************/
/* PROC GAMPL requires using the MISSING VALUE TRICK:
   https://blogs.sas.com/content/iml/2014/02/17/the-missing-value-trick-for-scoring-a-regression-model.html
*/
/***************************************************************/

/* 3. Score the model on the grid */
data Both;
set Have Grid(in=s);
score = s;              /* indicate that these are scoring points */
run;

ods select none;
proc gampl data=Both plots=none;
   model y = Spline(x);
   output out=ScoreOut pred=p_y;
   ID x score;  /* use the ID statement to get other vars in the output data set */
run;
ods select all;

/* 4. Compute slope by using finite difference formula. */
data Deriv0;
set ScoreOut(where=(score=1));  /* only use grid points */
Slope = dif(p_y) / dif(x);      /* (f(x) - f(x-dx)) / dx */
xPrev = lag(x);  yPrev = lag(p_y);  SlopePrev = lag(Slope);
if n(SlopePrev) AND sign(SlopePrev) ^= sign(Slope) then do;
   /* The slope changes sign between this obs and the previous.
      Assuming linearity on the interval, find (t, f(t))
      where slope is exactly zero */
   t0 = xPrev - SlopePrev * (x - xPrev)/(Slope - SlopePrev); 
   /* use linear interpolation to find the corresponding y value:
      f(t) ~ y0 + (y1-y0)/(x1-x0) * (t - x0)       */
   f_t0 = yPrev + (yPrev - p_y)/(x - xPrev) * (t0 - xPrev);
   if sign(SlopePrev) > 0 then _Type_ = "Max";
   else _Type_ = "Min";
   output; 
end;
keep t0 f_t0 _Type_;
label f_t0 = "f(t0)";
run;

proc print data=Deriv0 label;
run;

/* 5. (Optional) Visualize the points where the curve has local extrema */
data Combine;
merge Have                           /* data   : (x, y)     */
      ScoreOut(where=(score=1) 
               rename=(x=t p_y=p_t)) /* curve  : (t, p_t)   */
      Deriv0;                        /* extrema: (t0, f_t0) */
run;

title "GAMPL Smoother";
title2 "Red Markers Indicate Zero Slope for Smoother";
proc sgplot data=Combine noautolegend;
   scatter x=x y=y;
   series x=t y=p_t / lineattrs=GraphData2;
   scatter x=t0 y=f_t0 / markerattrs=(symbol=circlefilled color=red);
   yaxis grid; 
run;

