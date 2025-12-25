/* SAS program to accompany the articles 
   "Nine ways to visualize a continuous univariate distribution in SAS"
   "Nine ways to visualize a continuous univariate distribution in SAS - Part 2"
   by Rick Wicklin, published 07JAN2026 and 14JAN2026 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2026/01/07/9-dataviz-distrib-sas-part1.html
   https://blogs.sas.com/content/iml/2026/01/14/9-dataviz-distrib-sas-part2.html

   This program shows how to visualize univariate data distributions in SAS.

   Inspired by "Nine ways to visualize a distribution" by Nicola Rennie
   https://www.linkedin.com/posts/nicola-rennie_tidytuesday-rstats-dataviz-share-7404110322237669376-iRsT/
*/

/* you can change the DSName and varName macro values to run this code
   on your own data */
%let DSName = sashelp.heart;
%let varName = Cholesterol;

proc means ndec=1 data=&DSName N NMiss Mean Std Min Q1 Median Q3 Max;
   var &varName;
run;

/* Part 1: Truly univariate visualization methods that take one variable */

ods graphics / width=300px height=200px;
title "Histogram";
/* Alternative: PROC UNIVARIATE + HISTOGRAM statement */
proc sgplot data=&DSName;
   histogram &varName;  * vertical scale is percent;
run;

title "Kernel Density Plot";
proc sgplot data=&DSName;
   density &varName / type=kernel scale=percent; * default vertical scale is density ;
run;


proc freq data=&DSName order=freq;
where &varName> 160;
tables &varName / maxlevels=20;
run;

title "Box Plot";
/* Alternative: PROC BOXPLOT */
proc sgplot data=&DSName;
   hbox &varName /  lineattrs=(thickness=2) nofill;
run;

title "Fringe Plot";
proc sgplot data=&DSName noautolegend;
   histogram &varName;
   fringe &varName / transparency=0.8;
   yaxis offsetmin=0.08;
run;


/* ECDF plot */
proc univariate data=&DSName noprint;
   var &varName;
   cdfplot &varName / vref=(10 25 50 75 90) odstitle="ECDF Plot";
run;

ods exclude all;
/* use PROC UNIVARIATE to create the CDFPLOT, but use SGPLOT to visualize 
   by using DROPLINE statement instead of REFLINE statement */
proc univariate data=&DSName noprint;
   var &varName;
   cdfplot &varName / vref=(10 25 50 75 90) statref=P 10 Q1 Q2 Q3 P 90;
   ods output cdfplot=ecdf;
run;
ods exclude none;
title "ECDF Plot";
proc sgplot data=ECDF;
   series x=ECDFX y=ECDFY;
   dropline x=HRefValue y=VRefValue / dropto=both;
   yaxis values=(0 to 100 by 10);
   xaxis label="&varName";
run;

/* Part 2: To get dot plots and heatmaps, you need to add a second variable
   that has a constant value */

/* The OneCategory data set includes a categorical variable with one level.
   This is useful for using the SCATTER and HEATMAP statements, 
   which require two variables.
*/
data OneCategory;
set &DSName(keep=&varName);
_ID = 1;
run;


/* a gradient is one strip of a lasagna plot */
title "Gradient Plot";
proc sgplot data=OneCategory;
   heatmap  x=&varName y=_ID / discretey colorstat=pct colormodel=TwoColorRamp;
   yaxis display=(nolabel noticks novalues);
   xaxis label="&varName";
run;

title "Strip Plot";
proc sgplot data=OneCategory noautolegend;
   scatter x=&varName y=_ID  / jitter transparency=0.9 markerattrs=(symbol=CircleFilled);
   yaxis type=discrete display=(nolabel noticks novalues);
run;

title "Strip Plot and Box Plot Overlay";
proc sgplot data=OneCategory noautolegend;
   scatter x=&varName y=_ID  / jitter transparency=0.9 markerattrs=(symbol=CircleFilled);
   yaxis type=discrete display=(nolabel noticks novalues);
   hbox &varName / category=_ID lineattrs=(thickness=2) medianattrs=(thickness=3) nofill nooutliers;
run;

title "Swarm Plot";
proc sgplot data=OneCategory noautolegend;
   scatter x=&varName y=_ID  / jitter=uniform transparency=0.9;
   yaxis type=discrete display=(nolabel noticks novalues);
run;

ods graphics / reset;

/* Others: 
   Violin plot: essentially a kernel density plot, often used when there are multiple groups
   Raincloud plot: Combine a kernel density plot and a strip plot
       (Very similar to the histogram + Fringe plot)
*/
