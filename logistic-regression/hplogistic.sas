/* SAS program to accompany the article 
   "Jump-start PROC LOGISTIC by using parameter estimates from PROC HPLOGISTIC"
   by Rick Wicklin, published 26JUN2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/06/26/logistic-estimates-from-hplogistic.html

   This program shows how to
   1. Simulate data from a logistic regression model.
   2. Use PROC HPLOGISTIC to fit a model.
   3. Examine the time spent in various portions of the regression analysis.
   4. Use the parameter estimates from PROC HPLOGISTIC in PROC LOGISTIC,
      thus eliminating the need for the optimization step of the 
      logistic regression.

   For information about the simulation program, see
   https://blogs.sas.com/content/iml/2019/01/28/simulate-data-regression-categorical-continuous.html
   https://blogs.sas.com/content/iml/2014/06/25/simulate-logistic-data.html
*/

/* Program based on Simulating Data with SAS, Chapter 11 (Wicklin, 2013, p. 208-209) */
%let N = 500000;                 /* 1a. number of observations in simulated data */
%let numCont = 30;               /* number of continuous explanatory variables */
%let numClass = 4;               /* number of categorical explanatory variables */
%let numLevels = 3;              /* (optional) number of levels for each categorical variable */

/*
   1. Simulate data from a logistic regression model.
*/
data SimLogi; 
  call streaminit(12345); 
  /* 1b. Use macros to define arrays of variables */
  array x[&numCont]  x1-x&numCont;   /* continuous variables named x1, x2, x3, ... */
  array c[&numClass] c1-c&numClass;  /* CLASS variables named c1, c2, ... */
  /* the following statement initializes an array that contains the number of levels
     for each CLASS variable. You can hard-code different values such as (2, 3, 3, 2, 5) */
  array numLevels[&numClass] _temporary_ (&numClass * &numLevels);  
 
  do k = 1 to &N;                 /* for each observation ... */
    /* 2. Simulate value for each explanatory variable */ 
    do i = 1 to &numCont;         /* simulate independent continuous variables */
       x[i] = round(rand("Normal"), 0.001);
    end; 
    do i = 1 to &numClass;        /* simulate values 1, 2, ..., &numLevels with equal prob */
       c[i] = rand("Integer", numLevels[i]);        /* the "Integer" distribution requires SAS 9.4M5 */
    end; 
 
    /* 3. Simulate response as a function of certain explanatory variables */
    eta = -6 - 3*x[1] +0.7*x[2] - 2*x[&numCont] +    /* coefficients for continuous effects */
       -3*(c[1]=1) - 3.01*(c[1]=2) + 5*c[&numClass];  /* coefficients for categorical effects */
    mu = logistic(eta);                             /* transform by inverse logit */
    y = rand("Bernoulli", mu);                      /* random binary response with probability mu */
    output;  
  end;  
  drop i k;
run;    

/*
   2. Use PROC HPLOGISTIC to fit a model.
   3. Examine the time spent in various portions of the regression analysis.
*/
ods select PerformanceInfo IterHistory Timing;
proc hplogistic data=simLogi OUTEST;
   class c1-c&numClass;
   model y(event='1') = x1-x&numCont c1-c&numClass;
   ods output ParameterEstimates=PE;
   performance details;
run;

/*
   4. Use the parameter estimates from PROC HPLOGISTIC in PROC LOGISTIC,
      thus eliminating the need for the optimization step of the 
      logistic regression.
*/
/* create INEST= data set from ParameterEstimates */
proc transpose data=PE out=inest(type=EST) label=_TYPE_;
   label Estimate=PARMS;
   var Estimate;
   id ParmName;
run;

/*
   Re-examine and compare the time spent in various portions of the regression analysis.
*/
ods select PerformanceInfo IterHistory Timing ParameterEstimates;
proc hplogistic data=simLogi INEST=inest MAXITER=0;
   class c1-c&numClass;
   model y(event='1') = x1-x&numCont c1-c&numClass;
   performance details;
run;

/* Now use PROC LOGISTIC. Initialize with parameter estimates from HPLOGISTIC.
   Request statistics that are not supported by PROC HPLOGISTIC.
*/
proc logistic data=simLogi INEST=inest 
   plots(only MAXPOINTS=NONE)=oddsratio(range=clip);
   class c1-c&numClass;
   model y(event='1') = x1-x&numCont c1-c&numClass / maxiter=0 covb;
   oddsratio c1;
run;
