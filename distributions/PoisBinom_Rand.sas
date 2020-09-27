/* SAS program to accompany the article 
   "The Poisson-binomial distribution"
   by Rick Wicklin, published 28Sep2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/09/28/poisson-binomial-distribution.html

   This program shows how to simulate random samples from the
   Poisson-binomial distribution. The P-B distribution is a generalization
   of the binomial distribution, so first the program
   simulates data from the binomial distribution in two different ways.
   The simulation code is presented by using the SAS DATA step and by using SAS/IML software.
*/

title; footnote;

/* Generate a random sample from the binomial distribution */
/* The Easy Way: call rand("Binom", p, N) */
%let m = 1000;              /* number of obs in random sample */
data BinomSample(keep=i k);
call streaminit(123);
p = 0.8; N = 10;            /* p = prob of success; N = num trials */
label k="Number of Successes";
do i = 1 to &m;
   k = rand("Binom", p, N); /* k = number of successes in N trials */
   output;
end;
run;

/*
proc means data=BinomSample N Mean Var;
   var k;
run;
*/

ods graphics/ reset;
title "Binomial Sample";
footnote J=L "Sample Size = &m";
proc sgplot data=BinomSample;
   vbar k;
   yaxis grid;
run;

/* The Alternative Way: Make N calls to rand("Bernoulli", p) */
data BinomSample2(keep=i k);
call streaminit(123);
p = 0.8; N = 10;
label k="Number of Successes";
do i = 1 to &m;
   k = 0;                      /* initialize k=0 */
   do j = 1 to N;              /* accumulate N Bernoulli variables */
      k = k + rand("Bernoulli", p);
   end;
   output;
end;
run;

title "Binomial Sample 2";
proc sgplot data=BinomSample2;
   vbar k;
   yaxis grid;
run;

/* Now generate a random sample from the Poisson-binomial distribution 
   that has 10 parameters: p1, p2, p3, ..., p10 */
data PoisBinomSample(keep=i k);
/* p[j] = probability of success for the j_th trial, i=1,2,...,10 */
array p[10] (0.2 0.2 0.3 0.3 0.4 0.6 0.7 0.8 0.8 0.9);  /* N = dim(p) */

/* if desired, use these arrays of parameters to see how the shapes changes:
array p[10] (0.2 0.2 0.3 0.3 0.4 0.1 0.1 0.1 0.1 0.2);
array p[10] (0.9 0.8 0.9 0.9 0.9 0.9 0.7 0.8 0.8 0.9);
*/

call streaminit(123);
label k="Number of Successes";
do i = 1 to &m;
   k = 0;                      /* initialize k=0 */
   do j = 1 to dim(p);         /* accumulate N Bernoulli variables */
      k = k + rand("Bernoulli", p[j]);
   end;
   output;
end;
run;

title "Poisson-Binomial Sample";
proc sgplot data=PoisBinomSample;
   vbar k;
   yaxis grid;
run;

/*
proc means data=PoisBinomSample N Mean Var;
   var k;
run;
*/

options ps=32000 ls=128;
proc iml;
/* Simulate from the Poisson-Binomial distribution. 
   Input: m = number of observations in sample
          p = column vector of N probabilities. The probability 
              of success for the i_th Bernoulli trial is p_i.
   Output: m realizations of X ~ PoisBinom(p) 
*/
start RandPoisBinom(m, _p);
   p = colvec(_p);
   b = j(nrow(p), m);     /* each column is a binary indicator var */
   call randgen(b, "Bernoulli", p); 
   return( b[+, ] );      /* return numSuccesses = sum each column */
finish;

/* The Poisson-binomial has N parameters: p1, p2, ..., pN */
p = {0.2 0.2 0.3 0.3 0.4 0.6 0.7 0.8 0.8 0.9};  /* 10 trials */
call randseed(123);
S = RandPoisBinom(&m, p);
call bar(S) grid="y" label="Number of Successes";

/* mean and variance of distribution (expected values) */
mu = sum(p);
var = sum( (1-p)#p );
/* sample estimates of mean and variance */
XBar = mean(S`);
s2 = var(S`);
Moments = (mu//XBar) || (var//s2);
print Moments[f=5.2 r={"Distribution" "Sample Estimate"} c={"Mean" "Variance"}];

QUIT;
title; footnote;
