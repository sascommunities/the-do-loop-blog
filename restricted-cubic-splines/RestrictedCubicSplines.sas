/* SAS program to accompany the articles
   "An easier way to perform regression with restricted cubic splines in SAS"
   by Rick Wicklin, published 18FEB2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/02/18/regression-restricted-cubic-splines-sas.html
AND 
   "Regression with restricted cubic splines in SAS"
   by Rick Wicklin, published 19APR2017 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2017/04/19/restricted-cubic-splines-sas.html

   This program shows how to 
   - Use the EFFECT statement with the SPLINE option to generate spline effects
   - Specify the spline basis, the number of knots, and the placement of the knots
   - Reproduce the results of the %RCSPLINE macro (Harrell, 2009)

   The data are the X=Weight and Y=mpg_city variables in the Sashelp.Cars data.
*/

/* create (X,Y) data from the Sashelp.Cars data. Sort by X for easy graphing. */
data Have;
   set sashelp.cars;
   rename mpg_city = Y  weight = X  model = ID;
run;
proc sort data=Have;  by X;  run;


/* Method 1: The EFFECT statement is directly supported by many procedures:
   GLIMMIX, GLMSELECT, LOGISTIC, PHREG, PLS, QUANTREG, ROBUSTREG, ....
   You can generate spline effects and use them as explanatory variables.
   The restricted cubic splines are generated by using 
   the NATURALCUBIC BASIS=TPF options. The placement of the knots are
   controlled by using the KNOTMETHOD= option. */

/* In SAS/STAT 15.1 (SAS 9.4M6) you can use the 
   KNOTMETHOD=PERCENTILELIST to specify the location of internal
   knots in terms of percentiles of the data. 
   The following statements fit data by using restricted cubic splines. */
ods select ANOVA ParameterEstimates SplineKnots;
proc glmselect data=Have;
   effect spl = spline(X/ details naturalcubic basis=tpf(noint)
             /* NEW in SAS/STAT 15.1 (SAS 9.4M6): PERCENTILELIST 
                provides an easy way to use Harrell's suggested knot placement */
             knotmethod=percentilelist(5 27.5 50 72.5 95));
   model Y = spl / selection=none;       /* fit model by using spline effects */
   output out=SplineOut predicted=Fit;   /* output predicted values */
quit;

title "Restricted Cubic Spline Regression";
title2 "Five Knots Placed at Percentiles of the DATA";
proc sgplot data=SplineOut noautolegend;
   scatter x=X y=Y;
   series x=X y=Fit / lineattrs=(thickness=3 color=red);
run;


/* The remainder of this program is for 
   "Regression with restricted cubic splines in SAS"
   https://blogs.sas.com/content/iml/2017/04/19/restricted-cubic-splines-sas.html
 */
ods select ANOVA ParameterEstimates SplineKnots;
/* fit data by using restricted cubic splines */
proc glmselect data=Have;
   effect spl = spline(X/ details naturalcubic basis=tpf(noint)                 
                knotmethod=percentiles(5)
             /*
                knotmethod=equal(5)
                knotmethod=rangefractions(0.05 0.275 0.50 0.725 0.95)
                knotmethod=list(2513,3174,3474.5,3903,5000)
                knotmethod=percentilelist(5 27.5 50 72.5 95)
             */
                                );
   model Y = spl / selection=none;       /* fit model by using spline effects */
   output out=SplineOut predicted=Fit;   /* output predicted values for graphing */
quit;

title "Restricted Cubic Spline Regression";
title2 "Five Knots Placed at Percentiles";
proc sgplot data=SplineOut noautolegend;
   scatter x=X y=Y;
   series x=X y=Fit / lineattrs=(thickness=3 color=red);
run;


/* Harrell (2010) recommends knots placed at the percentiles (5 27.5 50 72.5 95).
   In SAS/STAT 15.1 (SAS 9.4M6) you can use the PERCENTILELIST option:
       knotmethod=percentilelist(5 27.5 50 72.5 95) 
   Prior to SAS 9.4M6, you can use PROC UNIVARIATE to generate the percentiles,
   then use the KNOTMETHOD=LIST option to specify the know locations, as shown below
*/
proc univariate data=Have noprint;
   var X;
   output out=percentiles pctlpts=5 27.5 50 72.5 95 pctlpre=p_;
run;
data _null_;
   set percentiles;
   call symput('pctls',catx(',', of _numeric_));
run;
%put &=pctls;
/* you can use this macro on the EFFECT statement:
  knotmethod=list(&pctls)
*/

/*************************************************************/
/* How does the regression curve depend on knot placement?   */
/* Compare knot locations for four knot placement methods:   */
/* EQUAL, PECENTILES RANGEFRACTIONS, and LIST/PERCENTILELIST */
/*************************************************************/

*ods select SplineKnots(persist);
ods exclude all;
proc glmselect data=Have;
   effect spl = spline(X / details naturalcubic basis=tpf(noint)                 
                knotmethod=equal(5));
   model Y = spl / selection=none; 
   output out=Equal predicted=Fit; 
   ods output SplineKnots=KEqual;
quit;
proc glmselect data=Have;
   effect spl = spline(X / details naturalcubic basis=tpf(noint)                 
                knotmethod=percentiles(5));
   model Y = spl / selection=none;  
   output out=Percentiles predicted=Fit;
   ods output SplineKnots=KPercentiles;
quit;
proc glmselect data=Have;
   effect spl = spline(X / details naturalcubic basis=tpf(noint)                 
                knotmethod=rangefractions(0.05, 0.275, 0.50, 0.725, 0.95));
   model Y = spl / selection=none;
   output out=Range predicted=Fit; 
   ods output SplineKnots=KRange;
quit;
proc glmselect data=Have;
   effect spl = spline(X / details naturalcubic basis=tpf(noint)                 
                knotmethod=list(&pctls) /* OR knotmethod=list(2513,3174,3474.5,3903,5000) */
                /* SAS 9.4M6:  knotmethod=percentilelist(5 27.5 50 72.5 95) */
                );
   model Y = spl / selection=none;  
   output out=PctlList predicted=Fit;
   ods output SplineKnots=KPList;
quit;
ods select all;

data KnotLocations;
merge KEqual(rename=X=Equal)
      KPercentiles(rename=X=Percentiles)
      KRange(rename=X=RangeFractions)
      KPList(rename=X=PercentileList);
drop ConstructedEffect Control;
run;
title "Knot Locations for Four Knot Placement Methods";
proc print data=KnotLocations noobs; run;

/* Merge the fitted curves for each knot-placement method.
   Create a graph that overlays the four curves */
data AllFit;
length Method $14;
set Equal(keep=X Y Fit)
    Percentiles(keep=X Y Fit)
    Range(keep=X Y Fit)
    PctlList(keep=X Y Fit) indsname=source;
Method = propcase(scan(source,2,'.'));
label Fit= "Predicted Y";
if Method = "Range" then Method="RangeFract";
if Method = "Pctllist" then Method="PercentileList";
run;

title2 "Fitted Curves from Knots Placed By Various Schemes";
proc sgplot data=AllFit;
   series x=X y=Fit / group=method;
run;

/* Show the know locations on the graph by using a fringe plot */
data AllKnots;
length Method $14;
set KEqual(in=E)
    KPercentiles(in=P)
    KRange(in=R)
    KPList(in=L);
if      E then Method = "Equal";
else if P then Method = "Percentiles";
else if R then Method = "RangeFract";
else if L then Method = "PercentileList";
drop ConstructedEffect Control;
run;

data All;
set AllFit AllKnots(rename=X=KnotLocation);
run;

title "Knot Locations By Various Schemes";
proc sgplot data=All;
   where method NOT in ("Percentiles");/* Percentiles and List almost the same */
   series x=X y=Fit / group=method lineattrs=(thickness=2);
   fringe KnotLocation / group=method lineattrs=(thickness=3);
run;


/*********************************************************************/
/* END blog post */
/*********************************************************************/



/*********************************************************************/
/* Compare the knot generation by PROC GLMSELECT to Harrell's macro  */
/*********************************************************************/

/* METHOD 2: %RCSpline macro */
 /*MACRO RCSPLINE

   For a given variable named X and from 3-10 knot locations,
   generates SAS assignment statements to compute k-2 components
   of cubic spline function restricted to be linear before the
   first knot and after the last knot, where k is the number of
   knots given.  These component variables are named c1, c2, ...
   ck-2, where c is the first 7 letters of X.

   Usage:

   DATA; ....
   %RCSPLINE(x,knot1,knot2,...,norm=)   e.g. %RCSPLINE(x,-1.4,0,2,8)

        norm=0 : no normalization of constructed variables
        norm=1 : divide by cube of difference in last 2 knots
                 makes all variables unitless
        norm=2 : (default) divide by square of difference in outer knots
                 makes all variables in original units of x

   Reference:

   Devlin TF, Weeks BJ (1986): Spline functions for logistic regression
   modeling. Proc Eleventh Annual SAS Users Group International.
   Cary NC: SAS Institute, Inc., pp. 646-51.

   Author  : Frank E. Harrell Jr.
             Clinical Biostatistics, Duke University Medical Center
   Date    : 10 Apr 88
   Mod     : 22 Feb 91 - normalized as in S function rcspline.eval
             06 May 91 - added norm, with default= 22 Feb 91
             10 May 91 - fixed bug re precedence of <>

                                                                      */
%MACRO RCSPLINE(x,knot1,knot2,knot3,knot4,knot5,knot6,knot7,
                  knot8,knot9,knot10, norm=2);
%LOCAL j v7 k tk tk1 t k1 k2;
%LET v7=&x; %IF %LENGTH(&v7)=8 %THEN %LET v7=%SUBSTR(&v7,1,7);
  %*Get no. knots, last knot, next to last knot;
    %DO k=1 %TO 10;
    %IF %QUOTE(&&knot&k)=  %THEN %GOTO nomorek;
    %END;
%LET k=11;
%nomorek: %LET k=%EVAL(&k-1); %LET k1=%EVAL(&k-1); %LET k2=%EVAL(&k-2);
%IF &k<3 %THEN %PUT ERROR: <3 KNOTS GIVEN.  NO SPLINE VARIABLES CREATED.;
%ELSE %DO;
 %LET tk=&&knot&k;
 %LET tk1=&&knot&k1;
 DROP _kd_; _kd_=
 %IF &norm=0 %THEN 1;
 %ELSE %IF &norm=1 %THEN &tk - &tk1;
 %ELSE (&tk - &knot1)**.666666666666; ;
    %DO j=1 %TO &k2;
    %LET t=&&knot&j;
    &v7&j=max((&x-&t)/_kd_,0)**3+((&tk1-&t)*max((&x-&tk)/_kd_,0)**3
        -(&tk-&t)*max((&x-&tk1)/_kd_,0)**3)/(&tk-&tk1)%STR(;);
    %END;
 %END;
%MEND;

/*********************************************************************/
/*********************************************************************/

/* Old way: 
   1. Call %RCSPLINE macro from DATA step to generate splines
   2. Use PROC GLM (or REG) to predict response from the spline bases
*/
data RCS;
   set Have;
   %rcspline(X, &pctls);
run;

title 'Basis functions from %RCSPLINE';
proc sort data=RCS out=Splines; by X; run;
proc sgplot data=Splines;
   series x=X y=X;
   series x=X y=X1;
   series x=X y=X2;
   series x=X y=X3;
run;

/* after the splines are generated, you can use PROC GLM to predict the response */
ods select OverallANOVA ParameterEstimates;
proc glm data=RCS;
   model Y = X:;
   output out=RCSOut predicted=RCSpred;
quit;

/* Modern way: Use PROC GLMSELECT to generate splines.
   PROC GLMSELECT can fit the curves as well, but here we write
   out the spline basis and use PROC REG (or GLM) to fit the curve */
proc glmselect data=Have outdesign(addinputvars fullmodel)=SplineBasis;
   effect spl = spline(X / naturalcubic basis=tpf(noint)
                knotmethod=list(&pctls)  ); /* restricted cubic spline */
   model Y = spl / selection=none;  
quit;
title "Basis functions from GLMSELECT";
proc sgplot data=SplineBasis noautolegend;
   series x=X y=Intercept;
   series x=X y=spl_1;
   series x=X y=spl_2;
   series x=X y=spl_3;
   series x=X y=spl_4;
run;

proc reg data=SplineBasis plots=none;
   model Y = spl_:;
   output out=SASOut p=SASPred;
run;

/* Overlay the two fits: %RSCSPLINE vs OUTDESIGN= from GLMSELECT */
title "Natural Cubic Spline";
title2 "RCSSPLINE Macro vs PROC GLMSELECT";
data all;
merge SASout(keep=SASPred)
      RCSOut(keep=RCSPred X:);
diff = SASPred - RCSPred;
d1 = divide(spl_2, X1);
d2 = divide(spl_3, X2);
d3 = divide(spl_4, X3);
run;
/* show that the two methods are the same to numerical precision */
proc means data=all;
var diff;
run;

title "Overlay Restricted Cubic Splines (5 knots)";
proc sgplot data=all noautolegend;
   series x=X y=SASpred / lineattrs=(thickness=2) ;
   series x=X y=RCSpred / lineattrs=(thickness=2) ;
run;