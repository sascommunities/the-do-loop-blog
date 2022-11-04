/* SAS program to accompany the article 
   "The area of the convex hull of random points"
   by Rick Wicklin, published 07NOV2022 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2022/11/07/area-random-convex-hull.html

   This program shows Buchta's (1984) formula, which gives
   the expected value of the area of a convex hull of 
   n random points in the unit square.

   The program runs a Monte Carlo simulation to estimate
   other properties of the area statistic for random points,
   such as the median and a prediction interval.
*/

/* Expected area of the convex hull of random sample of size n in the unit square.
   Result by Buchta (1984, in German), as reported by 
   https://mathoverflow.net/questions/93099/area-enclosed-by-the-convex-hull-of-a-set-of-random-points
*/
data EArea;
do n = 3 to 100;
   s = 0;
   do k = 1 to n+1;
      s + (1 - 1/2**k) / k;
   end;
   EArea = 1 - 8/3 /(n+1) * (s - 2**(-n-1) / (n+1));
   output;
end;
keep n EArea;
run;

/* print (n, E(Area(n))) as fractions */
proc print data=EArea noobs; 
format EArea fract32.;
where n <= 6;
run;

/* print selected values of (n, E(Area(n))) */
proc print data=EArea noobs; 
where n in (3,4,6,8,10,12,13) or mod(n,10)=0;
run;

ods graphics/reset;
title "Expected Area of Convex Hull";
proc sgplot data=EArea;
   series x=n y=EArea;
   xaxis grid label="Number of Random Points in Unit Square";
   yaxis grid label="Expected Area" min=0 max=1;
run;


/* define function for the area of a polygon. See
   https://blogs.sas.com/content/iml/2022/11/02/area-perimeter-convex-hull.html
*/

proc iml;
/* PolyArea: Return the area of a simple polygon.
   P is an N x 2 matrix of (x,y) values for vertices. N > 2.
   This function uses Gauss' "shoelace formula," which is also 
   known as the "triangle formula."  
   See https://en.wikipedia.org/wiki/Shoelace_formula  */
start PolyArea(P);
   lagIdx = 2:nrow(P) || 1;   /* vertices of the lag; appends first vertex to the end of P */
   xi   = P[ ,1];       yi = P[ ,2];
   xip1 = xi[lagIdx]; yip1 = yi[lagIdx]; 
   area = 0.5 * sum(xi#yip1 - xip1#yi);
   return( area );
finish;
store module=(PolyArea);
QUIT;



proc iml;
/* load the PolyArea function or copy it from 
   https://blogs.sas.com/content/iml/2022/11/02/area-perimeter-convex-hull.html
*/
load module=(PolyArea);

/* compute the area of the convex hull of n pts in 2-D */
start CHArea(Pts);
   indices = cvexhull(Pts);            /* compute convex hull of the N points */
   hullIdx = indices[loc(indices>0)];  /* the positive indices are on the CH */
   P = Pts[hullIdx, ];                 /* extract rows to form CH polygon */
   return PolyArea(P);                 /* return the area */
finish;

/* generate distribution of Area(n=4) */
call randseed(1234);
B = 1E4;                               /* number of Monte Carlo samples */
n = 4;
X = randfun((B*n)//2, "Uniform");   /* matrix of uniform variates in [0,1]x[0,1] */
Area = j(B, 1, .);                  /* allocate array for areas */
do i = 1 to B;
   beginIdx = (i-1)*n + 1;             /* first obs for this sample */
   ptIdx    = beginIdx:(beginIdx+n-1); /* range of obs for this sample */
   Pts = X[ptIdx, ];                   /* use these N random points */
   Area[i] = CHArea(Pts);              /* find the area */
end;

meanArea = mean(Area);
medianArea = median(Area);
call qntl(CL95, Area, {0.025 0.975});
print n meanArea medianArea (CL95`)[c={LCL UCL}];

title "Distribution of Areas of Convex Hull";
title2 "n=4; B = 10,000";
refStmt = "refline 0.15278 / axis=x lineattrs=(color=red) label='E(Area)';";
call histogram(Area) other=refStmt;


/* repeat for N=50 */
B = 1E4;                               /* number of Monte Carlo samples */
n = 50;
X = randfun((B*n)//2, "Uniform");   /* matrix of uniform variates in [0,1]x[0,1] */
Area = j(B, 1, .);                  /* allocate array for areas */
do i = 1 to B;
   beginIdx = (i-1)*n + 1;             /* first obs for this sample */
   ptIdx    = beginIdx:(beginIdx+n-1); /* range of obs for this sample */
   Pts = X[ptIdx, ];                   /* use these N random points */
   Area[i] = CHArea(Pts);              /* find the area */
end;

meanArea = mean(Area);
medianArea = median(Area);
call qntl(CL95, Area, {0.025 0.975});
print n meanArea medianArea (CL95`)[c={LCL UCL}];

title "Distribution of Areas of Convex Hull";
title2 "n=50; B = 10,000";
refStmt = "refline " + strip(char(meanArea)) + 
          " / axis=x lineattrs=(color=red) label='E(Area)';";
call histogram(Area) other=refStmt;
QUIT;

/* run Monte Carlo simulation for many values of n. 
   Estimate properties such as median and prediction interval.
   Write results to data set and visualize 
*/
proc iml;
load module=(PolyArea);
/* compute the area of the convex hull of n pts in 2-D */
start CHArea(Pts);
   indices = cvexhull(Pts);            /* compute convex hull of the N points */
   hullIdx = indices[loc(indices>0)];  /* the positive indices are on the CH */
   P = Pts[hullIdx, ];                 /* extract rows to form CH polygon */
   return PolyArea(P);                 /* return the area */
finish;

/* repeat Monte Carlo simulation for many values of N. Write results to data set and visualize */
create ConvexHull var {'n' MeanArea LCL UCL};

numPts = 3:100;
NSim = 1E4;                            /* number of Monte Carlo samples */
do j = 1 to ncol(numPts);                 /* number of random pts for each convex hull */
   N = numPts[j];
   X = randfun((NSim*N)//2, "Uniform");   /* matrix of uniform variates in [0,1]x[0,1] */
   Area = j(NSim, 1, .);                  /* allocate array for areas */
   LCL = j(NSim, 1, .);                   /* allocate arrays for 95% lower/upper CL */
   UCL = j(NSim, 1, .);
   do i = 1 to NSim;
      beginIdx = (i-1)*N + 1;             /* first obs for this sample */
      ptIdx    = beginIdx:(beginIdx+N-1); /* range of obs for this sample */
      Pts = X[ptIdx, ];                   /* use these N random points */
      Area[i] = CHArea(Pts);
   end;

   meanArea = mean(Area);
   call qntl(CL95, Area, {0.025 0.975});
   LCL = CL95[1];
   UCL = CL95[2];
   append;
end;
close ConvexHull;
QUIT;

title "Monte Carlo Estimate of Expected Area of Convex Hull";
proc sgplot data=ConvexHull noautolegend;
   band x=n lower=LCL upper=UCL / legendlabel="95% Prediction";
   series x=n y=meanArea / legendlabel="Expected Area";
   keylegend / location=inside position=E opaque across=1;
   xaxis grid label="Number of Random Points in Unit Square";
   yaxis grid label="Expected Area" min=0 max=1;
run;

