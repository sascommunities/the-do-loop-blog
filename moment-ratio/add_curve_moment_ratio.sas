/* SAS program to accompany the article 
   "How to add a curve to a moment-ratio diagram"
   by Rick Wicklin, published 23JUN2025 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2025/06/23/add-curve-moment-ratio.html 

   This program shows how to add a skewness-kurtosis curve for the Weibull
   distribution to a simple moment-ratio diagram. It builds on the M-R diagram 
   in Appendix E of Simulating Data with SAS (Wicklin, 2013).
   https://support.sas.com/content/dam/SAS/support/en/books/simulating-data-with-sas/65378_Appendix_E_Constructing_a_Moment_Ratio_Diagram_in_SAS.pdf
*/

/* BEFORE RUNNING THIS PROGRAM:
   Download the SAS code to create a moment-ratio diagram from GitHub:
   https://github.com/sascommunities/the-do-loop-blog/blob/master/moment-ratio/define_moment_ratio.sas
   Download the file to some location. When you %include that file, it 
     1. Creates an annotation data set named "Anno"
     2. Defines a macro %PlotMRDiagram(DS, annoDS, Transparency=0)
        that you can use to overlay sample estimates of skew and kurt on a moment-ratio diagram.
*/
%include "define_moment_ratio.sas";   /* specify a full path, such as "C:/<path>/define_moment_ratio.sas" */

/* Display a simple moment-ratio diagram and overlay a (skew,kurt) value for a data sample */
data SK;
Skew = 1.0704; Kurt = 1.3983;
run;

title "Moment-Ratio Diagram for Sample";
%PlotMRDiagram(SK, Anno);


/* If you look at the file, you'll see that the M-R diagram is 
   drawn by using an annotation data set. The data set is the concatenation of 
   several pieces:

   data Anno;
      set Beta Boundary LogNormal GammaCurve MRPoints;
   run;

   If we want to add more curves, we can mimic the existing code and create a new annotation set
   for whatever distribution we want to support.
*/

/***************************************/

/* The M-R diagram is very simple. It does not contains a curve for the Weibul(k) distributions.
   How do you add a Weibull curve to the M-R diagram? 

   I show how to do it in Appendix E: Constructing a Moment-Ratio Diagram in SAS
   (Wicklin, 2013)
*/


/* First, write SAS versions of the formulas for skewness and kurtosis from Wikipedia:
   https://en.wikipedia.org/wiki/Weibull_distribution
*/

%macro SkewExKurt(k);
   mu = Gamma(1 + 1/k);                         /* mean (lambda=1) */
   var = (Gamma(1 + 2/k) - Gamma(1 + 1/k)**2);  /* variance */
   sigma = sqrt(var);                           /* std dev */
   skew = ( Gamma(1 + 3/k) - 3*mu*var - mu**3 ) / sigma**3;
   x1 = skew;                                   /* x1=skew; y1=excess kurt */
   y1 = ( Gamma(1 + 4/k) - 4*mu*sigma**3*skew - 6*mu**2*var - mu**4 ) / var**2;
   y1 = y1 - 3;                                 /* excess kurtosis */
%mend;

/* Compute the skew-kurt curve for the Weibull distribution for various values of k.
   Experiments show that k < 0.9 results in the (skew,kurt) curve being too extreme.
   The range k in [0.9, 1000] seems to trace nearly the full curve.
*/
data WeibullCurve1;
keep k x1 y1;
label x1="Skewness" y1="Excess Kurtosis";
k0 = 0.9;
k_infin = 1000;
do k = k0 to 10 by 0.2, 50, 100, k_infin; 
   %SkewExKurt(k);
   output;
end;
run;

title "Skewness and Excess Kurtosis Curve for Weibull Distribution";
title2 "k in [0.9, 1000]";
proc sgplot data=WeibullCurve1;
scatter x=x1 y=y1 / datalabel=k;
   xaxis grid;
   yaxis grid reverse;
run;


/* After you can successfully draw the curve, you'll want to incorporate these points 
   into the previous M-R diagram. To do this, you must create an annotation data set
   and use the POLYLINE and POLYCONT statements to draw the curve.  See the examples in 
   define_moment_ratio.sas
*/
data WeibullCurve;
length function $12 Label $24 Curve $12 LineColor $20 ;
retain DrawSpace "DataValue"
       LineColor "Brown"
       Curve     "Weibull";
drop k0 k_infin k mu var sigma skew;
k0 = 0.9; 
k_infin = 1000;
k = k0;
%SkewExKurt(k);
function = "POLYLINE"; Anchor="     ";
output;
function = "POLYCONT";
do k = 1 to 1.9 by 0.1, 
       2 to 5 by 0.25, 
       6 to 12,
       15, 20, 25, 50, 100, k_infin;
   %SkewExKurt(k);
   output;
end;
label = "Weibull"; function = "TEXT"; Anchor="Left ";
output;
run;


/* Add the new curve to the existing set of curves and regions: */
data Anno;
  set Beta Boundary LogNormal GammaCurve WeibullCurve MRPoints;
run;
/* Draw the new moment-ratio diagram with the Weibull curve */
title "Moment-Ratio Diagram with Weibull Curve";
%PlotMRDiagram(SK, Anno);


/* How variable are the skewness and kurtosis for a random sample from the Weibull distribution?
   Simulate 100 samples of size N=250 from a Weibull(k=1.5) distribution.
   You can set lambda=1 because the skew and kurt do not depend on lambda.
*/
%let N = 200;
%let NumSamples = 100;
data Weibull(keep=x SampleID);
call streaminit(12345);
do SampleID = 1 to &NumSamples;
   do i = 1 to &N;
      x = rand("Weibull", 1.5, 1);
      output;
   end;
end;
run;

/* compute the skewness and kurtosis of each sample */
proc means data=Weibull noprint;
  by SampleID;
  var x;
  output out=MomentsWeibull skew=Skew kurt=Kurt;
run;

/* Display a moment-ratio diagram with a scatter plot in the background 
   of the sample skewness and kurtosis for simulated samples. */
title "Moment-Ratio Diagram and Sample Moments for Weibull(k=1.5)";
%PlotMRDiagram(MomentsWeibull, Anno);
