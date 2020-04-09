/* SAS program to accompany the article 
   "The geometric distribution in SAS"
   by Rick Wicklin, published 06Apr2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/04/06/geometric-distribution-sas.html

   This program shows how to compute and use the four essential functions
   for working with the geometric distribution in SAS.
   The four functions compute the PDF, CDF, quantile, and generate random variates.
*/

/* 
For background on the geometric distribution, see
https://en.wikipedia.org/wiki/Geometric_distribution

Note that there are two standard definitions of the geometric distribution.
Wicklin (2013) says:
 For the geometric distribution, PDF("Geometric",t,0.5) computes the 
 probability of t FAILURES, t=0,1,2,...  Use PDF("Geometric",t-1,0.5) 
 to compute the number of TOSSES until heads appears, t=1,2,3,.... 

Wicklin (2013) also says:
The two definitions are interchangeable  
because if $X$ is a geometric random variable that 
satisfies the first definition, then $X-1$ is a random 
variable that satisfies the second definition. 
SAS uses both definitions. The first definition is 
used by the RAND and RANDGEN functions, whereas (regrettably) 
the PDF function uses the second definition.  
Notice that a random variable that obeys the first definition 
takes on positive values; for the second definition, the variable 
is nonnegative.
*/

/* PCD and CDF functions.
1/2  = 0.5: the probability of "heads" when you toss a coin
1/6  = 0.166: the probability of rolling a 6 with a six-sided die
1/13 = 0.077: the probability of drawing an ace from a shuffled deck of 52 cards
*/
data Geometric;
do p = 0.5, 0.166, 0.077;
   /* We want the probability that the event occurs on the n_th trial.
      This is the same as the probability of n-1 failure before the event. */
   do n = 1 to 16;                     
      pdf = pdf("Geometric", n-1, p);  /* probability that event occurs on the n_th trial */
      cdf = cdf("Geometric", n-1, p);  /* probability that event happens on or before n_th Trial*/
      cdf2 = 1 - (1-p)**n;             /* explicit formula */
      output;
   end;
end;
run;

ods graphics / reset;

title "Probability Mass Function for Geom(p)";
proc sgplot data=Geometric;
   scatter x=n y=pdf / group=p markerattrs=(symbol=CircleFilled);
   series  x=n y=pdf / group=p lineattrs=(color=Gray);
   label n="Number of Trials" pdf="Probability of Event";
   xaxis grid integer values=(0 to 16 by 2) valueshint;
run;

title "Cumulative Distribution for Geom(p)";
title2 "Probability That Event Occurs on or before n_th Trial";
proc sgplot data=Geometric;
   step x=n y=cdf / group=p curvelabel markers markerattrs=(symbol=CircleFilled);
   label n="Number of Trials" cdf="Cumulative Probability";
   xaxis grid integer values=(0 to 16 by 2) valueshint;
   yaxis grid;
run;

/* The CDF is a step function, so not inverse. The quantile function 
   is also a step function. 
*/
data Quantile;
prob = 0.166;
do p = 0.01 to 0.99 by 0.005;
   n = quantile("Geometric", p, prob); * number of failures until event;
   n = n+1;                            * number of trials until event;
   output;
end;
run;

title "Quantiles for Geom(0.166)";
proc sgplot data=Quantile;
   scatter x=p y=n / markerattrs=(symbol=CircleFilled);
   label n="Number of Trials" p="Cumulative Probability";
   xaxis grid;
   yaxis grid;
run;

/* Generate random draw from Geometric distribution */
%let numSamples = 20;           *sample of 20 people;
data RandGeom;
call streaminit(12345);         * random seed for covid-19;
p = 0.166;
do ID = 1 to &numSamples;
   x = rand("Geometric", p);    *number of trials until event occurs;
   /* if you want the number of failures before the event: x=x-1 */
   Zero = 0;
   output;
end;
run;
/* Check correctness against population values 
proc means data=RandGeom min max mean stddev median;
var x;
run;
*/
proc sort data=RandGeom;
by x;
run;

data Stats;
p = 0.166;
PopMean = 1/p;
PopMedian = ceil( -1/log2(1-p) );
call symputx("PopMean", PopMean);
call symputx("PopMedian", PopMedian);
run;
proc print;run;

/* use CLIPCAP amd MAX= option to prevent long bars */
title "&numSamples Random Rolls of a Six-side Die";
title2 "Random Sample of Geom(0.166)";
proc sgplot data=RandGeom;
   highlow y=ID low=Zero high=x / type=bar highlabel=x CLIPCAP; 
   refline &PopMean / axis=x label="Mean" labelpos=min; 
   refline &PopMedian / axis=x label="Population Median" labelpos=min; 
   xaxis grid label="Number of Trials until Event" minor          MAX=10;  
   yaxis type=discrete discreteorder=data reverse valueattrs=(size=7);
run;
