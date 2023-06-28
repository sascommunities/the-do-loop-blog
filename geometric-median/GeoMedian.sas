/* SAS program to accompany the article 
   "Compute the geometric median in SAS"
   by Rick Wicklin, published 28JUN2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/06/28/geometric-median-sas.html

   This program shows how to compute the geometric median for N points
   in k-dimensional space. The article discusses two methods:

   1. In SAS Viya LTS 2022.09, PROC OPTMODEL supports computing 
      geometric medians 
      https://go.documentation.sas.com/doc/en/sasstudiocdc/v_037/pgmsascdc/casmopt/casmopt_conicsolver_examples02.htm
   2. In SAS IML, you can implement the Ostresh (1978) variation of the 
      Weiszfeld (1937) algorithm for the geometric median of 
      k-dimensional points:
      Ostresh, L. M. (1978) "On the Convergence of a Class of 
         Iterative Methods for Solving the Weber Location Problem"
         Operations Research, 26(4), pp. 597-609
         https://www.jstor.org/stable/169721
*/

/******** GEOMETRIC MEDIAN OF N-D POINTS ************************/

/******** PROC OPTMODEL *****************************************/
/* In SAS Viya LTS 2022.09, PROC OPTMODEL supports computing 
   geometric medians 
   https://go.documentation.sas.com/doc/en/sasstudiocdc/v_037/pgmsascdc/casmopt/casmopt_conicsolver_examples02.htm
*/
/*
%let n = 7;
%let k = 2;
proc optmodel;
   var x, y;                            * (x, y) is the geometric median;
   var d{1..&n} >= 0;                   * d[i] = dist from geometric median to each sample point;
   minimize f = sum{i in 1..&n}d[i];    * minimize the sum of distances;

   * define the sample points ;
   num samples{1..&n, 1..&k} = [0 0, 2 0, 1 1, 1 3, 0 2, 3 2, 6 3 ];
   * calculate each distance ;
   con distance{i in 1..&n}: soc(d[i], (x-samples[i,1]) (y-samples[i,2]));
   expand;

   solve with conic;
   print x y;                           * print optimal solution;
quit;
*/


/******** PROC IML **********************************************/
/* Implement the Ostresh (1978) variation of the Weiszfeld (1937) 
   algorithm for the geometric median of k-dimensional points.
   Ostresh, L. M. (1978) "On the Convergence of a Class of 
      Iterative Methods for Solving the Weber Location Problem"
      Operations Research, 26(4), pp. 597-609
      https://www.jstor.org/stable/169721
/****************************************************************/

/* The SAS IML program uses the following user-defined functions:
   IndexRowNoMatch: To avoid a technical problem in Weiszfeld's 
      original algorithm, this function checks whether the estimate 
      of the geometric median is one of the data points. 
      It returns the indices to the data points that are different 
      from the estimate.

   SumDist: The objective function. It returns the sum of the 
      distance from an arbitrary point to the data points.</li>

   UpdateGeoMedian: Implement one step in Weiszfeld's algorithm. 
      Given a point, the function returns another point that is 
      "downhill" from the input point in the direction of the 
      gradient vector.

   GeoMedian: The main driver routine. The function uses the relative 
      change of the objective function to determine convergence, but
      you could implement more sophisticated convergence criteria.
      It uses the centroid of the data points as an initial guess. 
      It iterates the initial guess until convergence and returns the 
      result.
*/
proc iml;
/* for a matrix, X, return the rows of X that are not all zero */
start IndexRowNoMatch0(X);
   count = (X=0)[,+];          /* count number of 0s in each row */
   idx = loc(count<ncol(X));   /* at least one elements in row is not 0 */
   return idx;
finish;
/* For a matrix, X, return the rows of X that are not equal to
   a target vector. By default, the vector is the zero
   vector j(1, ncol(X), 0). */
start IndexRowNoMatch(X, target=);
   if IsSkipped(target) then 
      return IndexRowNoMatch0(X);
   else 
      return IndexRowNoMatch0(X-target);
finish;

/* Return the sum of weighted distances from m to rows of P.
   The rows of P are reference points. The point m is an 
   arbitrary point.
*/
start SumDist(m, P, wt=j(nrow(P),1,1));
   dist = distance(P, m);    /* dist from each row to m */
   return sum(wt # dist);    /* sum of weighted distances */
finish;

/* Implement one step in Weiszfeld's iterative algorithm to
   find the geometric median. This function assumes that 
   no row of P equals m. */
start UpdateGeoMedian(m, P, wt);
   dist = distance(P, m);    /* dist from each row to m */
   w = wt / dist;            /* scaled weights (avoid zeros in dist) */
   mNew = w`*P / sum(w);     /* new estimate for geometric median */
   return mNew;
finish;

/* Implement Ostresh's modification of Weiszfeld's iterative 
   algorithm to find the geometric median. 
   Input
      P: The rows of P specify the reference points.
      PrintIter: Set PrintIter=1 to display the iteration history.
      wt: (optional) Specify a column vector of weights if you
          want a weighted median 

   Output: An estimate of the geometric median.
*/
start GeoMedian(P, PrintIter=0, wt=j(nrow(P),1,1));
   MaxIt = 100;     /* maximum number of iterations */
   FCrit = 1E-6;    /* criterion: stop when |fPrev - f| < FCrit */
   k = ncol(P);     /* dimension of points */
   n = nrow(P);     /* number of points points */
   m = (wt`*P)[:,] / sum(wt);  /* initial guess is weighted centroid */

   IterHist = j(MaxIt+1, k+1, .);
   IterHist[1,] = m || sum(distance(P,m));
   converged = 0;
   do i = 1 to MaxIt until(converged);
      mPrev = IterHist[i, 1:k];        /* previous guess for m */
      fPrev = IterHist[i, k+1];        /* previous sum of distances */
      idx = IndexRowNoMatch(P, m);     /* is the guess a data point? */
      if ncol(idx)=n then
         m = UpdateGeoMedian(m, P, wt);           /* no, use all pts */
      else
         m = UpdateGeoMedian(m, P[idx,], wt[idx,]); /* yes, exclude pt */
      f = SumDist(m, P, wt);           /* evaluate objective function */
      fDiff = abs(fPrev - f) / f;      /* relative change */
      converged = (fDiff < FCrit);     /* has algorithm converged? */
      IterHist[i+1,] = m || f;         
   end;
   if PrintIter then do;
      IterHist = T(0:i) || IterHist[1:i+1,];
      colNames = 'Iter' || ('x1':strip('x'+char(k))) || 'Sum Dist';
      print IterHist[c=colNames L="Iteration History"];
   end;
   return m;
finish;
store module=(IndexRowNoMatch0 IndexRowNoMatch SumDist 
              UpdateGeoMedian GeoMedian);
QUIT;


/* load the computational routines and compute the geometric median */
proc iml;
reset wide;
load module=_all_;

/* N=7 points in 2-D */
P = {0 0, 2 0, 1 1, 1 3, 0 2, 3 2, 5 3};
gm = GeoMedian(P, 1);       /* estimate the geometric median */
sumDist = SumDist(gm, P);   /* get the sum of distances */
print gm[L='Geometric Median'],
      sumDist[L='Sum of Distances'];

/* perform a computation for the weighted median */
wt = {1,1,1,2,1,2,2};
gmWt = GeoMedian(P, 0, wt);    /* estimate weighted geometric median */
sumDist = SumDist(gmWt, P, wt);/* weighted sum of distances */
print gmWt[L='Weighted Geo Median'],
      sumDist[L='Weighted Sum of Distances'];



/********* VISUALIZATION SECTION (optional) *********************/
/* visualize the objective function by creating a heat map 
   for a grid of points in the bounding box (BBox) */

/* Define a helper function:
   LinSpace returns a vector, v, that has N points, where
   v[1]=Min, v[N]=Max, and the elements are evenly spaced.
   Example; t = LinSpace(0, 1, 7);
*/
start LinSpace(min, max, N);
   dt = (max - min)/(N-1);
   return do(min, max, dt);
finish;

xy = expandgrid( LinSpace(P[><,1], P[<>,1], 41),
                 LinSpace(P[><,2], P[<>,2], 25) );
f = j(nrow(xy), 1, .);
do i = 1 to nrow(xy);
   f[i] = SumDist(xy[i,], P);
end;
create SS from xy f[c={'x' 'y' 'SumDist'}];
append from xy f;
close;

create TT from P[c={'px' 'py'}];
append from P;
close;

create GG from gm[c={'gx' 'gy'}];
append from gm;
close;

submit;
data Have;
   set SS TT GG;
run;
ods graphics / width=400px height=250px;
title "Sum of Distances from (x,y) to Points";
%let Spectral7 = (CX3288BD CX99D594 CXE6F598 CXFFFFBF CXFEE08B CXFC8D59 CXD53E4F );

proc sgplot data=Have noautolegend aspect=0.5;
   scatter x=x y=y / colorresponse=SumDist name="scat"
           markerattrs=(symbol=SquareFilled size=9) colormodel=&Spectral7;
   scatter x=px y=py / markerattrs=(symbol=Circle size=12 color=black);
   scatter x=gx y=gy / markerattrs=(symbol=X size=12 color=black);
   gradlegend "scat" / title="Sum of Distances";
   xaxis grid label='x';
   yaxis grid label='y';
run;
endsubmit;

/* now do WEIGHTED distance */
wt = {1,1,1,2,1,2,2};
xy = expandgrid( LinSpace(P[><,1], P[<>,1], 41),
                 LinSpace(P[><,2], P[<>,2], 25) );
f = j(nrow(xy), 1, .);
do i = 1 to nrow(xy);
   f[i] = SumDist(xy[i,], P, wt);
end;
create SS from xy f[c={'x' 'y' 'SumWtDist'}];
append from xy f;
close;

create TT from P wt[c={'px' 'py' 'Weight'}];
append from P wt;
close;

create GG from gmWt[c={'gx' 'gy'}];
append from gmWt;
close;

submit;
data Have;
set SS TT GG;
run;

ods graphics / width=400px height=250px;
title "Sum of Weighted Distances from (x,y) to Points";
%let Spectral7 = (CX3288BD CX99D594 CXE6F598 CXFFFFBF CXFEE08B CXFC8D59 CXD53E4F );

proc sgplot data=Have noautolegend aspect=0.5;
   scatter x=x y=y / colorresponse=SumWtDist name="scat"
           markerattrs=(symbol=SquareFilled size=9) colormodel=&Spectral7;
   bubble x=px y=py size=Weight / bradiusmax=10;
   scatter x=gx y=gy / markerattrs=(symbol=X size=12 color=black);
   gradlegend "scat" / title="Sum of Weighted Distances";
   xaxis grid label='x';
   yaxis grid label='y';
run;
endsubmit;

/***** show that the same code works for arbitrary dimensions */
/* 3-D points */
P = {0.8 -0.2 0, 
     0.3  1   0,
     0    0   0,
     0.3  1   1};
gm = GeoMedian(P);          /* estimate the geometric median */
sumDist = SumDist(gm, P);   /* get the sum of distances */
print gm[L='Geometric Median'],
      sumDist[L='Sum of Distances'];
QUIT;

/* Examine the performance.
   Create data set with k variables and N observations. */
%let N = 100000;
%let k = 10;
data Sim;
array x[&k];
call streaminit(123);
do i = 1 to &N;
   do j = 1 to &k;
      x[j] = rand("Normal");
   end;
   output;
end;
drop i j;
run;

/* test the performance on large data */
proc iml;
load module=_all_;

* read the data;
use Sim; read all var _NUM_ into X; close;

t0 = time();
gm = GeoMedian(X);    /* compute the geometric median */
time = time() - t0;
print time[F=5.3];
QUIT;


/* how does the performance depend on N, the number of observations? */
proc iml;
load module=_all_;

call randseed(12345);
k = 10;
size = T( {10000 50000} || do(1e5, 1e6, 1e5) );
Y = j(max(size), k);
call randgen(Y, "Normal");
time = j(nrow(size), 1, .);
do i = 1 to nrow(size);
   X = Y[1:size[i], ];
   t0 = time();
   gm = GeoMedian(X);    /* compute the geometric median */
   time[i] = time() - t0;
end;

print size time;

title "Time to Compute Geometric Median of N Points";
call series(size, time) grid={x y} label={"Num Points", "Time (s)"};

QUIT;
