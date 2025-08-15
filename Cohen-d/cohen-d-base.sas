/* SAS program to accompany the articles by Rick Wicklin on The DO Loop:
   "Cohen's D statistic in SAS" (11AUG2025)
   https://blogs.sas.com/content/iml/2025/08/11/cohens-d-statistic-in-sas.html
   and
   "Confidence intervals for Cohen's d statistic in SAS" (18AUG2025)
   https://blogs.sas.com/content/iml/2025/08/18/cohens-d-ci-in-sas.html

   This program shows how to compute estimates of the effect size for a 
   standardized mean difference (SMD) in a two-sample independent design.
   The estimates are 
       Cohen's d
       Hedges' g
   For each statistic, we compute a standard error and a 
   confidence interval for the SMD.
 
   The best overview of effect size statistics is
     Goulet-Pelletier & Cousineau (2018)​ 
     https://www.tqmp.org/RegularArticles/vol14-4/p242/p242.pdf 
   When you read the paper, you will learn:
     1. There are many different estimates for the SMD.
     2. There are many different formulas for the standard error.
     3.  There are several ways to estimate a confidence interval.
   In the paper, G-P & C state:
     1. "Hedges’ g should always be preferred over Cohen’s d" (p. 245) and
        "Hedges g is superior to the other estimators" (p. 255).
     2. The best estimate of SE is "True SE."
     3. CIs obtained from noncentral t distribution are exact for Hedges' g
        if the true noncentrality parameter is known.
        "The noncentral CI is without a doubt the most reliable method especially 
         when n is smaller than 20." G-P & C (p. 251).
        "The noncentral method is superior to the central method" (p. 255).
   This program shows how to compute both Cohen's d and Hedges' g, as well as the 
   recommended formulas for SE and CI.It uses only Base SAS procedures:
   PROC TTEST and the DATA step.
*/


/* data from PROC TTEST documentation: 
   Golf scores for students in a physical education class */
data scores;
   input Gender $ Score @@;
   datalines;
f 75  f 76  f 80  f 77  f 80  f 77  f 73
m 82  m 80  m 85  m 81  m 78  m 83  m 82 m 76 m 81
;

/* run t test for two independent samples; 
   save the statistics n1,n2,m1,m2,PooledSD */
%let DSName = scores;
%let Group = Gender;
%let Y = Score;

proc ttest data=&DSName plots=none;
   class &Group;
   var &Y;
   ods select Statistics;
   ods output Statistics=TTestOut;
run;

/* read statisitcs n1,n2,m1,m2,PooledSD; use to form Cohen's d and Hedges' g */
data EffectSize;
retain n1 m1 n2 m2 PooledSD;
keep n1 m1 n2 m2 PooledSD d g J;
label d="Cohen's d" g="Hedges' g" J="Hedges' correction factor";
set TTestOut end=EOF;
if _N_=1 then do;
   n1 = N; m1 = Mean;
end;
if _N_=2 then do;
   n2 = N; m2 = Mean;
end;
if _N_=3 then 
   PooledSD = StdDev;
if EOF then do;
   /* Cohen's d, which is a biased estimator */
   d = (m1 - m2) / PooledSD;   
   /* J is Hedge's correction factor. See Cousineau & Goulet-Pelletier (2020, p. 244) */
   df = n1 + n2 - 2;
   J =  exp( lgamma(df/2) - lgamma((df-1)/2) - log(sqrt(df/2)) ); 
   g = J * d;                  /* Hedge's g, which is an unbiased estimator */
   output;
end;
run;

title "Effect Sizes";
title2 "Two Independent Samples";
proc print label noobs data=EffectSize;
run;

/* The formulas for standard errors are in Table 3 (p. 246) of
   Goulet-Pelletiera and Cousineau (2018).
   https://www.tqmp.org/RegularArticles/vol14-4/p242/p242.pdf 

   The formulas depend on whether you use d or g
   to estimate the population standardized effect size. The symbol 'delta'
   is the population effect size, which we estimate by using d or g.
   
   Standard error of statistic for standardized effect size:
   Var     Label               Reference
   SE      "True SE"           Hedges (1981), p. 111, Eqn. 6b 
   SE_mle  "MLE Approx SE" 
   SE_h    "Hedges' Approx SE" 

   The best estimate of SE is "True SE". */
data SEEffectSize;
length Statistic $10;
keep Statistic Estimate SE SE_mle SE_h delta_lower delta_upper;
label Estimate="Estimate" SE="True SE" SE_mle="MLE Approx SE" SE_h="Hedges' Approx SE"
      delta_lower="Lower 95% CL" delta_upper="Upper 95% CL";
set EffectSize;
do i = 1 to 2;
   if i=1 then do;
      Statistic="Hedges' g"; Estimate = g;    /* g is a better estimate of effect size */
   end;
   else do;
      Statistic="Cohen's d"; Estimate = d;    /* but you can use d, if you want */
   end;

   /* the harmonic mean of n1 and n2 is  2*n1*n2/(n1+n2) */
   delta = Estimate;          /* use either d or g as the estimate */
   w = 2/harmean(n1,n2);      /* 2/H(n1,n2) = (n1+n2)/(n1*n2) */
   df = n1 + n2 - 2;          /* degrees of freedom */
   /* The best estimate of standard error is the "true formula" */
   var = df/(df-2) * (w + delta**2) - (delta/J)**2;
   SE = sqrt( var );           

   /* The following formulas are approximations. They are included because some
      references use them */
      /* Hedges' Approx SE: Hedges (1981, p. 117) approx for N > 50 */
      var_h = w + delta**2/(2*df);   
      SE_h   = sqrt( var_h );
      /* MLE Approx SE: Hedges & Olkin (1985), Eqn.11, p.82.
         This is used in the Handbook of Research Synthesis (p. 238)
         https://stats.stackexchange.com/questions/8487/how-do-you-calculate-confidence-intervals-for-cohens-d */
      factor = (df+2)/df;       /* MLE correction factor */
      SE_mle = sqrt( var_h * factor );
  /* END approx formula section */

  /* CIs obtained from noncentral t distribution; this is exact for Hedges' g. */
  /* The assumptions for the CI are 
      1. the data are sampled independently
      2. the populations are normal
      3. groups have the same variance (eg, homogeneity of variance)
   */
   lambda = delta / sqrt(w);  /* noncentrality parameter for t distribution (2 indep groups) */
   alpha = 0.05;
   t_lower = quantile("t",   alpha/2, df, lambda);  /* noncentral t,   alpha/2 quantile */
   t_upper = quantile("t", 1-alpha/2, df, lambda);  /* noncentral t, 1-alpha/2 quantile */
   /* transform CI for t into CI for delta. Note: This actually reduces to t*sqrt(w) */
   delta_lower = t_lower / (lambda/delta);          /* LCL for effect size, delta */
   delta_upper = t_upper / (lambda/delta);          /* UCL for effect size, delta */
   output;
end;
run;

title "Comparison of Standard Error Estimates";
title2 "Two Independent Samples";
proc print label noobs data=SEEffectSize;
   var Statistic Estimate SE SE_mle SE_h;
run;

title "Effect Size and Confidence Interval";
title2 "Two Independent Samples";
proc print label noobs data=SEEffectSize;
   var Statistic Estimate SE  delta_lower delta_upper;
run;

