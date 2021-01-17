/* SAS program to accompany the article:
   "The stationary block bootstrap in SAS"
   by Rick Wicklin, published 17JAN2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/01/20/stationary-bootstrap-sas.html

   This program shows how to perform the stationary block bootstrap 
   technique for time series.
*/

/*********** THE DATA *************/
/* A time series requires block resampling of residuals, which are
   a stationary series. Case resampling for ordinary least squares regression
   https://blogs.sas.com/content/iml/2018/10/29/bootstrap-regression-residual-resampling.html 
   can't be applied to time series because the error terms in a time series
   analysis are not assumed to be independent.

   Sashelp.Air has 144 months of data. Twelve months is a good block length,
   but that would result in 12 blocks of length 12.
   For clarity, let's drop the first year of data, which results
   in 11 blocks of length 12.
*/

/******************************/
/* STATIONARY BLOCK BOOTSTRAP */
/******************************/
/* Use AUTOREG to decompose series as Y = Predicted + Residuals
   Similar to Getting Started example in PROC AUTOREG */
data Air;
   set Sashelp.Air;
   if Date >= '01JAN1950'd;
   Time = _N_;
run;
proc autoreg data=Air;
   AR12: model Air = Time / nlag=12;
   output out=OutReg pm=Pred rm=Resid;  /* mean prediction and residuals */
   ods select FinalModel.ParameterEstimates ARParameterEstimates;
run;


/* The stationary block bootstrap:
   Use a random starting point and a random block length.
   The random block length is chosen from a Geom(p) 
   distribution, which has the interpretation:
   - probability (1-p) of using the next residual from the current block
   - probability p of jumping to a new starting position/block
*/
%let p = (1/12);       /* expected block length is 1/p */ 
proc iml;

/* generate indices for a block of random length L~Geom(p)
   chosen from the sequence 1:n */
start RandBlockIdx(n, p);
   beginIdx = randfun(1, "Integer", n); /* random start: 1:n */
   L = randfun(1, "Geom", p);           /* random length */
   endIdx = ((beginIdx+L-1) >< n);      /* truncate into 1:n */
   return( beginIdx:endIdx );           /* return sequence of indices in i:n */
finish;

/* generate indices for one stationary bootstrap resample */
start StationaryIdx(n, p);
   bootIdx = j(n,1,.);
   s = 1;                       /* starting index for first block */
   do while( s <= n );
      idx = RandBlockIdx(n, p); /* get block of random length */
      if ncol(idx) > n-s then   /* truncate if block too long */
         idx = idx[, 1:n-s+1];
      L = ncol(idx);            /* length of this block */
      jdx = s:s+L-1;      
      bootIdx[jdx] = idx;       /* assign this block */
      s = s + L;                /* starting index for next block */
   end;
   return bootIdx;
finish;

use OutReg;
   read all var {'Time' 'Pred' 'Resid'};
close;

call randseed(1234);
n = nrow(Pred);                 /* length of series */
prob = &p;                      /* probability of jumping to a new block */

/*---- Visualization Example: generate one bootstrap resample ----*/
idx = StationaryIdx(n, prob);
YBoot = Pred + Resid[idx];

/* visualize the places where a sequence is not ascending:
https://blogs.sas.com/content/iml/2013/08/28/finite-diff-estimate-maxi.html */
blocksIdx = loc(dif(idx) <= 0); /* where is the sequence not ascending? */
b = Time[blocksIdx];
f = T(blocksIdx) // n;
print f (dif(f)) (mean(dif(f)));

title "One Bootstrap Resample";
title2 "Stationary Block Bootstrap";
refs = "refline " + char(b) + " / axis=x;";
call series(Time, YBoot) other=refs;
/*---- End of visualization example ----*/

/* A complete stationary block bootstrap repeats this process B times
   and usually writes the resamples to a SAS data set. */
B = 1000;
SampleID = j(n,1,.);
create BootOut var {'SampleID' 'Time' 'YBoot'};  /* create outside of loop */
do i = 1 to B;
   SampleId[,] = i;   /* fill array. See https://blogs.sas.com/content/iml/2013/02/18/empty-subscript.html */
   idx = StationaryIdx(n, prob);
   YBoot = Pred + Resid[idx];
   append;
end;
close BootOut;
QUIT;

/* Analyze the bootstrap samples by using a BY statement. See
   https://blogs.sas.com/content/iml/2012/07/18/simulation-in-sas-the-slow-way-or-the-by-way.html
*/
proc autoreg data=BootOut plots=none outest=BootEst noprint;
   by SampleID;
   AR12: model YBoot = Time / nlag=12;
run;

title "Bootstrap Estimates of Mean and 90% CIs";
proc means data=BootEst mean stddev P5 P95;
   var Intercept Time _A:;
run;

title "Distribution of Parameter Estimates";
proc sgplot data=BootEst;
   scatter x=Intercept y=Time;
   ellipse x=Intercept y=Time;
   xaxis grid; yaxis grid;
   refline 77.5402 / axis=x;
   refline 2.7956  / axis=y;
run;

