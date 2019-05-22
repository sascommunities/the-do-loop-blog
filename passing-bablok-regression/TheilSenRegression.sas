/* SAS program to accompany the article
   "The Theil-Sen robust estimator for simple linear regression"
   by Rick Wicklin, published 27MAY2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/05/27/theil-sen-robust-regression.html

   This program shows how to implement a simple (one variable) robust 
   regression model due to Theil (1950) and Sen (1968). 
   Data is from Sen (7 obs), but 2 outliers were added.

   Long history of robust statistics:
   Not just "modern" method of regression. Univariate stats include
   trimmed means (idea as early as 1820), Winsorized means (Winsor died in 1951),
   and Tueky's contaminated distributions.

   Theil-Sen robust estimation of slope (Theil, 1950; Sen 1968) 
   https://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator
   Wikipedia claims "the most popular nonparametric technique for estimating a linear trend"
   Unbiased estimator of slope. Breakdown value 29.3%
   Theil (1950) and Sen (1968)
   https://www.pacificclimate.org/~wernera/zyp/Sen%201968%20JASA.pdf
*/
proc iml;
XY = {1  9,  2 15,  3 19, 4 20, 10 45,  12 55, 18 78, /* 7 data points used by Sen (1968) */
     12.5 30,  4.5 50};                               /* 2 outliers (not considered by Sen) */

/* Theil uses all "N choose 2" combinations of slopes of segments.
   Assume that the first coordinates (X) are distinct */
c = allcomb(nrow(XY), 2);        /* all "N choose 2" combinations of pairs */
Pt1 = XY[c[,1],];                /* first point of line segments */
Pt2 = XY[c[,2],];                /* second point of line segments */
slope = (Pt1[,2] - Pt2[,2]) / (Pt1[,1] - Pt2[,1]); /* Careful! Assumes x1 ^= x2 */
m = median(slope);
b = median( XY[,2] - m*XY[,1] );  /* median(y-mx) */
print (b||m)[c={'Intercept' 'Slope'} L="Method=Theil Pairs=All"];

lp = round(0 || b || m,  1e-3);
title "Theil Estimate of Robust Regression Line";
call scatter(XY[,1], XY[,2]) grid={x y} lineparm=lp;  /* line through (0,b) with slope m */

title "Distribution of Slopes of Line Segments Connecting Pairs";
call histogram(slope) rebin={-50, 5} other="refline 3.96875/axis=x;";

/* compare to least squares linear regression */
X = j(nrow(XY),1,1) || XY[,1];
Y = XY[,2];
b = solve( X`*X, X`*Y );
print b;
lp = round(0 || b`, 1e-3);
title "Least Squares Regression Line";
call scatter(XY[,1], XY[,2]) grid={x y} lineparm=lp;  /* line through (0,b) with slope m */

create MyData from XY[c={x y}]; append from XY; close;
x1 = Pt1[,1]; x2 = Pt2[,1]; y1 = Pt1[,2]; y2 = Pt2[,2]; 
create Pts var {x1 y1 x2 y2}; append; close;
quit;

data All; set MyData Pts; run;

title "Line Segments Connecting Pairs of Points";
proc sgplot data=All noautolegend;
   vector x=x1 y=y1 / xorigin=x2 yorigin=y2 noarrowheads;
   scatter x=x y=y / markerattrs=(symbol=CircleFilled color=darkred size=12);
   xaxis grid; yaxis grid;
run;

/***************************/
proc iml;
/* Sen (1968) handles repeated X coords by using only pairs with distinct X */
XY = {1  9,  2 15,  3 19, 4 20, 10 45,  12 55, 18 78,
     18 30,  4 50};  /* last two obs are outliers not considered by Sen */
c = allcomb(nrow(XY), 2);        /* all "N choose 2" combinations of pairs */
Pt1 = XY[c[,1],];                /* first point of line segments */
Pt2 = XY[c[,2],];                /* second point of line segments */
idx = loc(Pt1[,1]-Pt2[,1]^=0);   /* find pairs with same X value */
Pt1 = Pt1[idx,];                 /* keep only pairs with different X values */
Pt2 = Pt2[idx,];

slope = (Pt1[,2] - Pt2[,2]) / (Pt1[,1] - Pt2[,1]);
m = median(slope);
b = median( XY[,2] - m*XY[,1] );  /* median(y-mx) */
print (b||m)[c={'Intercept' 'Slope'} L="Method=Sen Pairs=All"];


/* of course, "N choose 2" = N*(N-1)/2 grows quadratically with N. 
   For N=15, you cn easily generate all 105 pairs.
   However, for N = 1000, you would have to compute almost 500,000 slopes. 
   An alternative is to use random combinations of points for slopes. 
*/
numPairs = 100;
call randseed(54321,1);
c = rancomb(nrow(XY), 2, numPairs); /* random combinations */
Pt1 = XY[c[,1],]; Pt2 = XY[c[,2],];
idx = loc(Pt1[,1]-Pt2[,1]^=0);   /* exclude any pairs with same X value */
Pt1 = Pt1[idx,];  Pt2 = Pt2[idx,];

slope = (Pt1[,2] - Pt2[,2]) / (Pt1[,1] - Pt2[,1]);
m = median(slope);
b = median( XY[,2] - m*XY[,1] );  /* median(y-mx) */
print (b||m)[c={'Intercept' 'Slope'} L="Method=Sen Pairs=200"];


proc iml;
/* Return (intercept, slope) for Theil-Sen robust estimate of a regression line.
   XY is N x 2 matrix. The other arguments are:
   METHOD: 
      If method="Theil" and a pair of points have the same X coordinate, 
         assign a large positive value instead of +Infinity and a large negative 
         value instead of -Infinity. 
      If method="Sen", omit any pairs of points that have the same first coordinate. 
   NUMPAIRS:
      If numPairs="All", generate all "N choose 2" combinations of the N points.
      If numPairs=K (integer), generate K random pairs of points. 
*/
start TheilSenEst(XY, method="SEN", numPairs="ALL");
   Infinity = 1e99;             /* big value for slope instead of +/- infinity */
   if type(numPairs)='N' then
      c = rancomb(nrow(XY), 2, numPairs);  /* random combinations of pairs */
   else if upcase(numPairs)="ALL" then 
      c = allcomb(nrow(XY), 2);            /* all "N choose 2" combinations of pairs */
   else stop "ERROR: The numPairs option must be 'ALL' or a postive integer";

   Pt1 = XY[c[,1],];                       /* first points for slopes */
   Pt2 = XY[c[,2],];                       /* second points for slopes */
   dy = Pt1[,2] - Pt2[,2];                 /* change in Y */
   dx = Pt1[,1] - Pt2[,1];                 /* change in X */ 
   idx = loc( dx ^= 0 );  
   if upcase(method) = "SEN" then do;      /* exclude pairs with same X value */
      slope = dy[idx] / dx[idx];           /* slopes of line segments */
   end;
   else do;                        /* assign big slopes for pairs with same X value */
      slope = j(nrow(Pt1), 1, .);  /* if slope calculation is 0/0, assign missing */
      /* Case 1: x1 ^= x2. Do the usual slope computation */
      slope[idx] = dy[idx] / dx[idx];
      /* Case 2: x1 = x2. Assign +Infinity if sign(y1-y2) > 0, else assign -Infinity */
      jdx = loc( dx = 0 & sign(dy)>0 );
      if ncol(jdx)>0 then 
         slope[jdx] = Infinity;
      jdx = loc( dx = 0 & sign(dy)<0 );
      if ncol(jdx)>0 then 
         slope[jdx] = -Infinity;
   end;
   *print c Pt1 Pt2 slope dy dx;
   m = median(slope);
   b = median( XY[,2] - m*XY[,1] );  /* median(y-mx) */
   return( b || m );
finish;

/* Test all four calls */
XY = {1  9,  2 15,  3 19, 4 20, 10 45,  12 55, 18 78,
     18 30,  4 50};  /* last two obs are outliers not considered by Sen */

ods layout gridded columns=2 advance=table;
   est = TheilSenEst(XY, "Theil", "All");
   print est[c={'Intercept' 'Slope'} L="Method=Theil; Pairs=All"];

   est = TheilSenEst(XY, "Sen", "All");
   print est[c={'Intercept' 'Slope'} L="Method=Sen; Pairs=All"];

   call randseed(123, 1);
   est = TheilSenEst(XY, "Theil", 200);
   print est[c={'Intercept' 'Slope'} L="Method=Theil; Pairs=200"];

   call randseed(123, 1);
   est = TheilSenEst(XY, "Sen", 200);
   print est[c={'Intercept' 'Slope'} L="Method=Sen; Pairs=200"];
ods layout end;

QUIT;
