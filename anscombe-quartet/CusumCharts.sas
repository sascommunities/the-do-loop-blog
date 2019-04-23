/* SAS program to accompany the article 
   "A CUSUM test for autregressive models"
   by Rick Wicklin, published 24APR2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/04/24/cusum-test-autregressive-models.html

   This program shows how to use PROC AUTOREG in SAS to generate CUSUM
   charts for recursive residuals. You can use the CUSUM charts as a
   diagnostic goodness-of-fit tool for a misspeccified model.

   Data from "Create your own version of Anscombe's quartet: 
   Dissimilar data that have similar statistics" by Rick Wicklin
   https://blogs.sas.com/content/iml/2019/04/17/create-version-of-anscombes-quartet.html
*/
data Anscombe2;
input X Y1 Y2 ;
datalines;
4.0 5.26 3.42619 
4.2 6.17 3.74636 
4.4 6.02 4.05710 
4.6 4.75 4.35843 
4.8 6.94 4.65034 
5.0 4.27 4.93282 
5.2 5.46 5.20589 
5.4 6.74 5.46955 
5.6 5.87 5.72378 
5.8 7.13 5.96859 
6.0 5.85 6.20399 
6.2 6.39 6.42996 
6.4 5.05 6.64652 
6.6 5.82 6.85366 
6.8 5.95 7.05138 
7.0 6.77 7.23968 
7.2 6.38 7.41857 
7.4 6.92 7.58803 
7.6 6.38 7.74808 
7.8 7.15 7.89870 
8.0 7.91 8.03991 
8.2 8.58 8.17170 
8.4 8.86 8.29407 
8.6 5.62 8.40702 
8.8 6.64 8.51056 
9.0 8.13 8.60467 
9.2 9.03 8.68937 
9.4 8.06 8.76465 
9.6 7.78 8.83050 
9.8 8.37 8.88694 
10.0 8.74 8.93397 
10.2 7.54 8.97157 
10.4 9.57 8.99975 
10.6 8.38 9.01852 
10.8 7.09 9.02786 
11.0 6.57 9.02779 
11.2 9.20 9.01830 
11.4 9.52 8.99939 
11.6 9.85 8.97106 
11.8 7.90 8.93332 
12.0 10.08 8.88615 
12.2 8.60 8.82957 
12.4 8.04 8.76356 
12.6 8.84 8.68814 
12.8 10.34 8.60330 
13.0 7.53 8.50904 
13.2 9.17 8.40536 
13.4 9.75 8.29227 
13.6 8.76 8.16975 
13.8 10.69 8.03782 
14.0 10.39 7.89647 
;

/**********************************************************/
/* PROC AUTOREG models a time series with autocorrelation */
/**********************************************************/

ods graphics on;
proc autoreg data=Anscombe2;
  Linear: model y1 = x;     /* Y1 is linear. Model is oorrectly specified. */
  output out=CusumLinear cusum=cusum cusumub=upper cusumlb=lower 
                         r=Residual recres=RecursiveResid;
run;
proc autoreg data=Anscombe2;
  Quadratic: model y2 = x;     /* Y1 is quadratic. Model is misspecified. */
  output out=CusumQuad cusum=cusum cusumub=upper cusumlb=lower 
                       r=Residual recres=RecursiveResid;
run;

title "Cumulative Recursive Residuals"
title2 "Linear Data (left) and Quadratic Data (right)";

ods graphics / width=400px height=300px;
ods layout gridded columns=2 advance=table;
 proc sgplot data=CusumLinear noautolegend;
    band x=x lower=lower upper=upper;
    series x=x y=cusum / break markers;
    refline 0  /axis=y noclip;
    xaxis grid; yaxis grid;
 run;
 proc sgplot data=CusumQuad noautolegend;
    band x=x lower=lower upper=upper;
    series x=x y=cusum / break markers;
    refline 0  /axis=y noclip;
    xaxis grid; yaxis grid;
 run;
ods layout end;

