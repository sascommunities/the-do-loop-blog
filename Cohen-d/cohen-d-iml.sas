
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
   recommended formulas for SE and CI. It implements the formulas
   in PROC IML.
*/

/* data from PROC TTEST documentation: 
   Golf scores for students in a physical education class */
data scores;
   input Gender $ Score @@;
   datalines;
f 75  f 76  f 80  f 77  f 80  f 77  f 73
m 82  m 80  m 85  m 81  m 78  m 83  m 82 m 76 m 81
;

/* computation of Cohen's d and Hedges' g in SAS IML.
   Helper Functions:
   run Indep2SampleStats(...) : get descriptive statistics mean1, mean2, n1, n2, df, pooled std dev
   delta = EffectSize(...)    : estimate effect size by using Cohen's d or Hedges' g
   SE =  EffectSize_SE(delta,...) : estimate standard error of estimate of effect size
   CI =  EffectSize_SE(delta,...) : estimate confidence interval for effect size
*/
proc iml;
/* get descriptive statistics mean1, mean2, n1, n2, df, pooled std dev */
start Indep2SampleStats(m1, m2, n1, n2, df, s_p, /* output stats */
                        Group, Y);               /* input stats  */
   u = unique(Group);
   x1 = Y[ loc(Group=u[1]) ];   /* reference group */
   n1 = countn(x1); m1 = mean(x1); s1 = std(x1);
   x2 = Y[ loc(Group=u[2]) ];   /* comparison group */
   n2 = countn(x2); m2 = mean(x2); s2 = std(x2);
   df = n1 + n2 - 2;
   s_p = sqrt( ((n1-1)*s1##2 + (n2-1)*s2##2) / df );  /* pooled SD */
finish;

/* estimate the effect size. Return Cohen's d if STAT="D", otherwise return Hedges' g */
start EffectSize(m1, m2, n1, n2, s_p, stat="G"); /* input statistics */
   d = (m1 - m2) / s_p;        /* Cohen's d, which is a biased estimator */
   if upcase(stat)="D" then 
      return( d );
   /* J is Hedge's correction factor. See Cousineau & Goulet-Pelletier (2020, p. 244) */
   df = n1 + n2 - 2;           /* degrees of freedom */
   J =  exp( lgamma(df/2) - lgamma((df-1)/2) - log(sqrt(df/2)) ); 
   g = J * d;                  /* Hedge's g, which is an unbiased estimator */
   return( g );
finish;

/* estimate standard error of estimate of effect size */
start EffectSize_SE(delta, n1, n2);
   /* the harmonic mean of n1 and n2 is  2*n1*n2/(n1+n2) */
   w = 2/harmean(n1//n2);     /* 2/H(n1,n2) = (n1+n2)/(n1*n2) */
   df = n1 + n2 - 2;          /* degrees of freedom */
   J =  exp( lgamma(df/2) - lgamma((df-1)/2) - log(sqrt(df/2)) ); 
   /* The best estimate of standard error is the "true formula" */
   var = df/(df-2) * (w + delta**2) - (delta/J)**2;
   return( sqrt(var) );
finish;

/* estimate confidence interval for effect size */
start EffectSize_CI(delta, n1, n2, alpha=0.05);
   /* the harmonic mean of n1 and n2 is  2*n1*n2/(n1+n2) */
   w = 2/harmean(n1//n2);     /* 2/H(n1,n2) = (n1+n2)/(n1*n2) */
   df = n1 + n2 - 2;          /* degrees of freedom */
   /* CIs obtained from noncentral t distribution; this is exact for Hedges' g */
   lambda = delta / sqrt(w);  /* noncentrality parameter for t distribution (2 indep groups) */
   t_lower = quantile("t",   alpha/2, df, lambda);  /*   alpha/2 quantile */
   t_upper = quantile("t", 1-alpha/2, df, lambda);  /* 1-alpha/2 quantile */
   /* transform CI for t into CI for delta */
   delta_lower = t_lower / (lambda/delta);
   delta_upper = t_upper / (lambda/delta);
   return( delta_lower || delta_upper );
finish;

/* test the functions on the golf score data */
DSName = "Scores";
ClassName = "Gender";
YName = "Score";
use (DSName);
   read all var ClassName into Group;
   read all var YName into Y;
close;

run Indep2SampleStats(m1, m2, n1, n2, df, s_p, /* output stats */
                      Group, Y);               /* input stats  */
*print n1 n2 m1 m2 df s_p;

/* Hedges' g statistic */
g    = EffectSize(m1, m2, n1, n2, s_p);  /* default is Hedges' g */
g_SE = EffectSize_SE(g, n1, n2);
g_CI = EffectSize_CI(g, n1, n2);
g_results = n1 || n2 || m1 || m2 || g || g_SE || g_CI;
lbl = {'n1' 'n2' 'mean1' 'mean2' 'Estimate' 'SE' 'LCL' 'UCL'};
print g_results[c=lbl L="Hedges' g for Two Independent Samples"];

/* Cohen's d statistic */
d    = EffectSize(m1, m2, n1, n2, s_p, "D");
d_SE = EffectSize_SE(d, n1, n2);
d_CI = EffectSize_CI(d, n1, n2);
d_results = n1 || n2 || m1 || m2 || d || d_SE || d_CI;
print d_results[c=lbl L="Cohen's d for Two Independent Samples"];

/* Notice that the estimate is not centered in the CI */
Left_width = g - g_CI[1];
Right_width = g_CI[2] - g;
print Left_width Right_width;


/* for fun, also look at the CI for delta +/- SE*t_df.
   This gives a symmetric CI.
*/
   /* CIs obtained from central t distribution */
   alpha = 0.05;
   t_lower = quantile("t",   alpha/2, df);  /*   alpha/2 quantile */
   t_upper = quantile("t", 1-alpha/2, df);  /* 1-alpha/2 quantile */
   print t_lower t_upper;
   CI_d_SE = (d + t_lower*d_SE) || (d + t_upper*d_SE);
   print CI_d_SE;
   CI_g_SE = (g + t_lower*g_SE) || (g + t_upper*g_SE);
   print CI_g_SE;
   
