/* SAS program to accompany the article 
   "Compute 2-D cumulative sums and ogives"
   by Rick Wicklin, published 30JUN2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/06/30/2d-cumulative-sums-ogives.html

   This program defines the CUSUM2D function, which computes the 
   2-D cumulative sum of a matrix.
   It also defines the HIST2D and OGIVE2D functions.
   An ogive is a histogram-based estimate of the CDF.
*/

/*************************************************************/
/* Optional: review the 1-D computation of a histogram 
   of counts and an ogive. 
   For explanation of GSCALE, BIN, and TABULATE to count obs
   in bins, see 
   https://blogs.sas.com/content/iml/2018/12/05/histogram-table-of-counts.html
   For a way to create an ogive from PROC UNIVARIATE, see
   https://blogs.sas.com/content/iml/2016/09/26/create-ogive-sas.html
*/
options ls=120;
proc iml;
/* Hist : compute a 1-D histogram
   Input: X     is an (n x 1) column vector of data 
          _cutX is a (1 x k) row vector of cutpoints 
   Output: The counts of X in the bins of a histogram defined by _cutX
*/ 
start Hist(x, _cutX);
   cutX = rowvec(_cutX);           /* make sure input is row vector */
   b = bin(x, cutX);               /* find bin for each obs */
   call tabulate(bins, freq, b);   /* count how many obs in each bin */
   /* in case not every bin has an obs, assign counts to bins */
   Count = j(1, ncol(cutX)-1, 0);  /* allocate  */
   Count[bins] = Freq;             /* fill nonzero counts */
   return Count;
finish;

use sashelp.Cars; read all var "MPG_City" into x; close;

/* get (start, stop, increment) for a histogram
  (or define your own regular grid) */
call gscale(s, x, 11);
print s;

cutX = do(s[1], s[2], s[3]);    /* use "nice" cutpoints from GSCALE */
Count = Hist(x, cutX);    /* fill nonzero counts */
binLabels = char(cutX);         /* use left endpoint as labels for bins */
print Count[c=binLabels];

/* You can use these values to create an ogive.
   Just append 0 on the left and compute the cumulative proportions */
rtPts = cutX || (s[2]+s[3]);              /* ogive computed at right point */
ECDF = 0 || (cusum(Count) / sum(Count));  /* add 0 as leftmost value */
title "Ogive Estimate of CDF";
call series(rtPts, ECDF) option="markers" grid={x y};
QUIT;

/*************************************************************/
/* START THE 2-D OGIVE PROGRAM HERE */
/*************************************************************/

proc iml;

/* compute the 2-D CUSUM:
   1. First, compute CR = cusum of counts across rows
   2. Then compute C = cusum of CR down columns
*/
start Cusum2D(X);
   C = j(nrow(X), ncol(X), .);  /* allocate result matrix */
   do i = 1 to nrow(X);
      C[i,] = cusum( X[i,] );   /* cusum across each row */
   end;
   do j = 1 to ncol(X);
      C[,j] = cusum( C[,j] );   /* cusum down columns */
   end;
   return C;
finish;
store module=(Cusum2D);

Count = { 7  6 1 3,
          6  3 5 2,
          1  2 4 1};
Cu = Cusum2D(Count);
print Count, Cu;


/* HIST2D : Given data and two vectors of endpoints for bins, count 
   the number of observations in each 2-D bin. Return the count.
   https://blogs.sas.com/content/iml/2013/07/17/2d-binning.html
*/
start Hist2D(u, cutX, cutY);
   bX = bin(u[,1], cutX);   /* bins in X direction: 1,2,...,kx */
   bY = bin(u[,2], cutY);   /* bins in Y direction: 1,2,...,ky */
   bin = bX + (ncol(cutX)-1)*(bY-1);     /* bins 1,2,...,kx*ky */
   call tabulate(binNumber, Freq, bin);    /* count in each bin */
   Count = j(ncol(cutY)-1, ncol(cutX)-1, 0);  /* allocate  */
   Count[binNumber] = Freq;              /* fill nonzero counts */
   return(Count);
finish;
store module=(Hist2D);

/* OGIVE2D : Given a matrix of counts, such as from HIST2D,
   that represent counts in the bins defined by the vectors
   cutX and cutY, return a list that contains
   L$'Ogive' : the matrix of cumulative proportions
   L$'Xpts'  : the X cut points with the point x_{N_x}+dx appended
   L$'YPts'  : the Y cut points with the point x_{N_y}+dy appended
*/
start Ogive2D(Count, cutX, cutY);
   Cu = Cusum2D(Count);                         /* cumulative counts */
   /* To get estimate of CDF, pad Cu with zeros on left and top */
   ogive = j(nrow(Cu)+1, ncol(Cu)+1, 0);
   rIdx = 2:nrow(ogive);
   cIdx = 2:ncol(ogive);
   ogive[rIdx, cIdx] = Cu / sum(Count);         /* cumulative proportions */
   /* adjust the cut points by appending another endpoint */
   Nx=ncol(cutX); dx= cutX[nX]-cutX[nX-1]; XX=cutX[nX]+dx;
   Ny=ncol(cutY); dy= cutY[nY]-cutY[nY-1]; YY=cutY[nY]+dy;
   rtX = cutX || XX;                            /* right-side endpoint */
   rtY = cutY || YY;
   L = [#"Ogive"=ogive,                         /* return a list of results */
        #"Xpts"=rtX,   
        #"Ypts"=rtY ];
   return L;
finish;
store module=(Ogive2D);
QUIT;

/* demonstrate how to use the HIST2D and OGIVE2D functions */
proc iml;
/* load computational modules */
load module=(Hist2D Cusum2D Ogive2D);

/* read bivariate data */
use sashelp.cars;
   read all var {MPG_City Weight} into Z;
close;

* define cut points in X direction (or use the GSCALE function);
cutX = do(10, 60, 5);        
cutY = do(1500, 7500, 500);
Count = Hist2D(Z, cutX, cutY);   /* matrix of counts in each 2-D bin */
L = Ogive2D(Count, cutX, cutY);  /* return value is a list */
ogive = L$'Ogive';               /* matrix of cumulative proportions */
print ogive[F=Best5. r=(L$'YPts') c=(L$'XPts')];

/* ALTERNATIVE:
call gscale(xScale, Z[,1], 11);
call gscale(yScale, Z[,2], 11);
print xScale yScale;

cutX = do(xScale[1], xScale[2], xScale[3]);  * cut points in X direction;
cutY = do(yScale[1], yScale[2], yScale[3]);  * cut points in Y direction;
*/

/* A quick-and-dirty viz:
   In a matrix, the rows increase down the matrix. But in a 
   plot, the Y axis points up, so reverse the rows to visualize */
CDFEst = ogive[nrow(ogive):1, ];  * reverse rows to plot;
call HeatmapCont(CDFEst) title="Ogive for 2-D Histogram" displayoutlines=0 
     colorramp="ThreeColor" range={0 1};

/* A better heat map: write ogive and labels to data set */
Z = ogive;
col = col(Z);  CNAMES = cutX[col];
row = row(Z);  RNAMES = cutY[row];
create _Heatmap var {RNAMES CNAMES Z}; append; close _Heatmap;
QUIT;


title "Ogive for 2-D Histogram";
proc sgplot data=_Heatmap;
   label CNAMES='MPG_City' RNAMES='Weight' Z='Cumulative Probability';
   heatmapparm x=CNAMES y=RNAMES colorresponse=Z; 
run;
title;
