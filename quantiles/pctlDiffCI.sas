/* SAS program to accompany the article 
   "How to estimate the difference between percentiles"
   by Rick Wicklin, published 03JUN2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/06/01/difference-between-percentiles.html

   This program shows how to use quantile regression in SAS to estimate
   the difference between a percentile in two groups. Given a sample of 
   size Nx from group 'A' and a sample of size Ny from group 'B', the 
   problem is to test whether the p_th percentile is different in the 
   two groups and to construct a confidence interval for the difference.

   This article is motivated by the 2019 paper:
   "Percentile-based grain size distribution analysis tools (GSDtools) – 
      estimating con?dence limits and hypothesis tests 
      for comparing two samples"
   by Brett C. Eaton, R. Dan Moore, and Lucy G. MacKenzie
   Earth Surf. Dynam., 7, 789–806, 2019 
   https://doi.org/10.5194/esurf-7-789-2019 

   In this program, data are simulated from lognormal distributions.
*/

/* Prelude and background:
   LN dist has median exp(mu). Use mu1=0.64 and mu2=0.53.
   Compute the 50th and 84th percentile for the LN distrib. 

   Why the 0.84 quantile? 
   Recall from the 68–95–99.7 rule that the area under the 
   standard normal curve between +/-1 standard deviation 
   from the mean is 0.68. Thus the quantile for 0.5+0.34=0.84 
   corresponds to mean + 1 stdDev.
*/
proc iml;
med1 = exp(0.64);
med2 = exp(0.53);
print med1 med2;

med1 = quantile("Lognormal", 0.5, 0.64, 1);
med2 = quantile("Lognormal", 0.5, 0.53, 1);
print med1 med2 (med1-med2)[L='Diff'];

P84_1 = quantile("Lognormal", 0.84, 0.64, 1);
P84_2 = quantile("Lognormal", 0.84, 0.53, 1);
print P84_1 P84_2 (P84_1-P84_2)[L='Diff'];
QUIT;


/* 1. Simulate data from LN(0.64, 1) and LN(0.53, 1) */
data Grains;
call streaminit(12345);
Site = 'A';
do i = 1 to 502;  /* simulate a 2 extra because 2 are out of range */
   Diameter = rand("Lognormal", 0.64, 1);
   Diameter = round(Diameter, 0.1);
   if 0 < Diameter < 35 then output;
end;
Site = 'B';
do i = 1 to 300;
   Diameter = rand("Lognormal", 0.53, 1);
   Diameter = round(Diameter, 0.1);
   if 0 < Diameter < 35 then output;
end;
drop i;
run;

/* 2. Compute sample median and 84th pctl. 
      Also create comparative histogram */
proc univariate data=Grains;
   class Site;
   histogram Diameter / endpoints=(0 to 28);
   ods select histogram;
   output out=Pctls pctlpre=P pctlpts=50,84; 
run; 
 
/* Difference between medians = 0.45
   Difference between P84 = 0.4 */
proc print data=Pctls noobs; run;

/* 3. Use PROC QUANTREG to estimate the difference between medians.
      Default CI is CI=RANK */
ods select ParameterEstimates Estimates;
proc quantreg data=Grains;
   class Site;
   model Diameter = Site / quantile=0.5 0.84;
   estimate 'Diff in Pctl' Site 1 -1 / CL;
run;

/* Technical Remark:
The median of Site 'A' has a nonunique solution.  
Any number between 0.4 and 0.5 is an optimal estimate for 
Site A for median regression in term of minimizing sum 
of check losses, so that the median of Site A diameter is in 
1.6+(0.4,0.5)=(2,2.1).   
The default simplex algorithm outputs 0.4, while if you use the 
ALGORITHM=SMOOTH option, the estimate will be 0.425.  
*/
/*
proc iml;
use Grains;
read all var {'Site' 'Diameter'};
close;
x = Diameter[ loc(Site='A') ];
med = do(1.5, 2.4, 0.001);
loss = j(1, ncol(med), .);
do i = 1 to ncol(med);
   loss[i] = sum(abs(x-med[i]));
end;
title "Median Is Not Unique";
title2 "Any value in [2, 2.1] Minimizes L1 Difference";
call series(med, loss);
QUIT;
*/

/* 4. QUANTREG supports multiple choices for estimating the CI.
   See if changing the method changes the inference.
   CI methods are:
   CI=Rank                  [default for small data]
   CI=Sparsity(BF)
   CI=Sparsity(HS)/IID
   CI=Resampling(NREP=5000) [default for large data]
*/
ods select none;
ods output Estimates(PERSIST) = CI;
proc quantreg data=Grains CI=RANK;
   class Site;
   model Diameter = Site / quantile=0.5 0.84;
   estimate 'Rank' Site 1 -1 / CL;
run;
proc quantreg data=Grains CI=Sparsity(BF);
   class Site;
   model Diameter = Site / quantile=0.5 0.84;
   estimate 'Sparsity(BF)' Site 1 -1 / CL;
run;
proc quantreg data=Grains CI=Sparsity(HS)/IID;
   class Site;
   model Diameter = Site / quantile=0.5 0.84;
   estimate 'Sparsity(HS) / IID' Site 1 -1 / CL;
run;
proc quantreg data=Grains CI=Resampling(NREP=5000);
   class Site;
   model Diameter = Site / quantile=0.5 0.84;
   estimate 'Resampling(NREP=5000)' Site 1 -1 / CL seed=12345;
run;
ods select all;  /* turn off PERSIST; restore normal output */

proc sort data=CI; by Quantile Label; run;

/* create graph and table for the results */
ods graphics / width=480px height=360px;
title "Difference Between Quantiles";
proc sgpanel data=CI noautolegend;
   panelby Quantile / columns=1;
   scatter x=Estimate y=Label/ xerrorlower=Lower xerrorupper=Upper;
   refline 0 / axis=x;
   label Estimate = "Estimate of Difference"
         Label = "CI= Method";
run;

proc print data=CI noobs;
   var Quantile Label Estimate Lower Upper;
run;

/************************************/
/* 5. Eaton, Moore, and MacKenzie (2019) used bootstrap methods to 
   estimate the difference between percentiles. Show how to use
   the bootstrap in SAS/IML. For each quantile (0.5 and 0.84),
   print the CI by using the percentile
   method and display the bootstrap distribution of the difference.
*/
proc iml;
call randseed(12345);
/* bootstrap estimate */
use Grains;
read all var {'Site' 'Diameter'};
close;
x = Diameter[ loc(Site='A') ];
y = Diameter[ loc(Site='B') ];
Nx = nrow(x); Ny = nrow(y);

/* naive point estimates of P50 and P84 */
p = {0.5, 0.84};        
call qntl(qxEst, x, p);
call qntl(qyEst, y, p);
diffEst = qxEst - qyEst;
print p[L="quantile"] qxEst qyEst diffEst;

/* do the bootstrap (next 5 lines) */
B = 5000;                    /* number of bootstrap samples */
BootX = sample(x, Nx // B);  /* each row is a bootstrap sample */
BootY = sample(y, Ny // B);
call qntl(qx, BootX`, p);    /* quantiles of each sample */
call qntl(qy, BootY`, p);
D = qx - qy;                 /* each row is a difference of quantiles */

/* bootstrap estimates of P50 and P84 */
diffBoot = T( mean(D`) );
print p[L="quantile"] diffBoot;

/* compute 95% bootstrap CIs by using the percentile method */
call qntl(CI50, T(D[1,]), {0.025, 0.975});
call qntl(CI84, T(D[2,]), {0.025, 0.975});
print CI50[r={"Lower" "Upper"}] , CI84[r={"Lower" "Upper"}];

/* display the results in a table */
Estimate = p || diffBoot || (CI50`//CI84`);
print Estimate[c={'Quantile' 'DiffBoot' 'Lower' 'Upper'} L='Bootstrap Estimate'];

/* 6. Graph the bootstrap distributions and overlay the bootstrap estimate of
   the difference and the CI */
title "Difference of Medians in Bootstrap Samples";
r = rowcat( char(diffBoot[1] || CI50`) );
refLineStmt = "refline " + r +"/axis=x label=('Mean' 'Lower' 'Upper');";
call histogram(D[1,]) other=refLineStmt;

title "Difference of P84 in Bootstrap Samples";
r = rowcat( char(diffBoot[2] || CI84`) );
refLineStmt = "refline " + r +"/axis=x label=('Mean' 'Lower' 'Upper');";
call histogram(D[2,]) other=refLineStmt;

QUIT;

