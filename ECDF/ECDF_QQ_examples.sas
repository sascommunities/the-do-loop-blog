/* SAS program to accompany articles 
   by Rick Wicklin, published on The DO Loop blog:
   https://blogs.sas.com/content/iml/2026/05/26/create-ecdf.html
   https://blogs.sas.com/content/iml/2026/06/08/confidence-bands-ecdf.html
   https://blogs.sas.com/content/iml/2026/06/15/confidence-band-qqplot.html

   Before running this program, you should have defined the 
   three functions for creating ECDFs and QQ plots. Run the program in
   the file ECDF_QQ.sas, which is included in the same directory as this file.

   The data are from the SAS documentation for PROC UNIVARIATE.
   The data set contains 50 measurements of the breaking strength of 
   50 random fiber-optic cords.
*/

data Cord;
   label Strength="Breaking Strength (psi)";
   input Strength @@;
datalines;
6.94 6.97 7.11 6.95 7.12 6.70 7.13 7.34 6.90 6.83
7.06 6.89 7.28 6.93 7.05 7.00 7.04 7.21 7.08 7.01
7.05 7.11 7.03 6.98 7.04 7.08 6.87 6.81 7.11 6.74
6.95 7.05 6.98 6.94 7.06 7.12 7.19 7.12 7.01 6.84
6.91 6.89 7.23 6.98 6.93 6.83 6.99 7.00 6.97 7.01
;
 
title 'Cumulative Distribution Function of Breaking Strength';
proc univariate data=Cord noprint;
   cdfplot Strength / odstitle = title;
run;

/*****************************************************/

/* Test the ECDF function */
proc iml;
/* In SAS Viya, the ECDF function is built into SAS IML and does not need to be loaded. */
load module=(ECDF);

use Cord;
read all var "Strength" into x;
close;
 
/* estimate the probability that the breaking strength is less than 
   6.8, 6.9, 7.0, and 7.1 psi. Note that some of these values 
   are not observed data values.
*/
t = {6.8, 6.9, 7.0, 7.1};
ecdf_t = ECDF(x, t);
print t ecdf_t;

/* typically, an ECDF plot evaluates the ECDF at all data points */
call sort(x);
ecdf = ECDF(x);
create ECDF var {"x" "ECDF"}; append; close;
QUIT;
 
/* Graph a step function. See
   https://blogs.sas.com/content/iml/2016/09/06/graph-step-function-sas.html */
title "Empirical CDF";
proc sgplot data=ECDF noautolegend;
   label x="Breaking Strength (psi)" ECDF="Cumulative Proportion";
   step x=x y=ecdf;
   fringe x;
   xaxis grid label="x" offsetmin=0.05 offsetmax=0.05;
   yaxis grid min=0 offsetmin=0.03;
run;

/*****************************************************/
/* Test the ECDF_KSCL function */
proc iml;
load module=(ECDF ECDF_KSCL);
 
/* read the data for the ECDF */
use Cord;  read all var "Strength" into x;  close;
/* evaluate the ECDF at all data points */
call sort(x);
ecdf = ECDF(x);
CL = ECDF_KSCL(ecdf);        /* 95% confidence bands for the ECDF */
*CL = ECDF_KSCL(ecdf, 0.90); /* get bands for any confidence levels in [0.75, 0.99] */
CL_Lower = CL[,1];
CL_Upper = CL[,2];
create ECDF_bands var {"x" "ECDF" "CL_Lower" "CL_Upper" };
  append;
close;
QUIT;
 
/* Graph the step functions by using PROC SGPLOT */
title "Empirical CDF and 95% Confidence Bands";
proc sgplot data=ECDF_bands noautolegend;
   label x="Breaking Strength (psi)" ECDF="Cumulative Proportion";
   step x=X y=ECDF / lineattrs=(thickness=2);
   step x=X y=CL_Lower / lineattrs=(color=gray);
   step x=X y=CL_Upper / lineattrs=(color=gray);
   xaxis grid label="x";
   yaxis grid min=0 offsetmin=0.03;
run;


/*****************************************************/
/* Test the QQ_KSCL function */
proc iml;
/* In SAS Viya, the ECDF function is built into SAS IML and does not need to be loaded. */
load module=(ECDF ECDF_KSCL QQ_KSCL);  
 
/* Read and sort the data */
use Cord;  read all var "Strength" into x;  close;
M = QQ_KSCL(x);

/* Create a dataset for plotting */
create QQ_Bands from M[c={"x" "q" "Q_Lower" "Q_Upper"}];
append from M;
close;
QUIT;

title "Normal Q-Q Plot with 95% Confidence Bands";
proc sgplot data=QQ_Bands noautolegend;
   /* Draw the confidence bands */
   series x=Q_Lower y=x / lineattrs=(color=gray);
   series x=Q_Upper y=x / lineattrs=(color=gray);
   
   /* Overlay the scatter plot of the quantiles */
   scatter x=q y=x;
   xaxis label="Theoretical Normal Quantiles" grid;
   yaxis label="Sample Quantiles (Strength)" grid;
run;
