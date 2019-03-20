/* SAS program to accompany the article
   "The geometry of multivariate versus univariate outliers"
   by Rick Wicklin, published 25MAR2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/03/25/geometry-multivariate-univariate-outliers.html

   This program computes and visualizes bivariate normal data 
   and two kinds of outliers: those that are extreme in each 
   coordinates and those that are not. The data are created in 
   standardized coordinates (where the geometry is easy) 
   and then transformed to data coordinates by using the Cholesky matrix
   https://blogs.sas.com/content/iml/2012/02/08/use-the-cholesky-transformation-to-correlate-and-uncorrelate-variables.html
*/

title;
ods graphics / width=400px height=400px attrpriority=none;
proc iml;
call randseed(123);
N = 200;

/* generate standardized uncorrelated data */
/* covariance matrix is product of diagonal (scaling) and correlation matrices */
x = randfun(N//2, "Normal");
R = {1   0.9,
     0.9 1  };
sd = {4 9};
Sigma = corr2cov(R, sd);
print Sigma;

/* The Cholesky root of Sigma is the matrix that 
   linearly scales and correlates the data */
U = root(Sigma);

/********************************/
/* add a few unusual points ("outliers") */
pi = constant('pi');
extra = (2.5*cos(0*pi/12) || 2.5*sin(0*pi/12)) //  /* MD approx 2.5 */
          (5*cos(7*pi/12) ||   5*sin(7*pi/12));    /* MD approx 5   */
extra = round(extra, 0.01);
v = x // extra;                /* concatenate MVN data and outliers */
w = v*U;                       /* transform from stdized to data coords */
call scatter(w[,1], w[,2]) grid={x y};
 
/****************************************************/
/* So that we can color the outliers, write orig data and 
   outliers to separate data sets */
/* write orig data */
w1 = x*U;           /* transform from stdized to data coords */
center = mean(w1);  std = std(w1);  cov = cov(w1);
z = (w1 - center) / std;
MD = MAHALANOBIS(w1, center, cov); 

Y = w1 || z || MD;
Label = j(N, 1, ' '); 
varNames = {'x1' 'x2' 'z1' 'z2' 'MD'};
create MVNorm from Y[r=Label colname=varNames];
append from Y[r=Label];
close;

/* write outliers */
w2 = extra*U;      /* transform from stdized to data coords */
z = (w2 - center) / std;
MD = MAHALANOBIS(w2, center, cov); 

Y = w2 || z || MD;
Label = {A, B}; 
varNames = {'xx' 'yy' 'z1' 'z2' 'MD'}; /* different var names */
create Label from Y[r=Label colname=varNames];
append from Y[r=Label];
close;

QUIT;

/* concatenate data sets in block form */
data All;
set MVNorm Label;
run;

title "Bivariate Normal Data";
title2 "'Unusual' Observations Added";
proc sgplot data=All noautolegend;
   scatter x=x1 y=x2;
   scatter x=xx y=yy / markerattrs=(symbol=CircleFilled) datalabel=Label;
   xaxis grid values=(-12 to 12 by 4) valueshint max=12; 
   yaxis grid values=(-30 to 30 by 10) valueshint max=30;
run;

proc print data=Label label noobs;
   var Label xx yy z1 z2 MD;
   label xx="x1" yy="x2";
run;

/* Print probabilities of obs with large MD for various distances.
   In 2-D, P(MD < t) = 1 - exp(-t**2 / 2) */
data Prob;
input t @@;
pLT = 1 - exp(-t**2/2);   /* prob MD < t */
pGT = 1 - pLT;            /* prob MD > t */
datalines;
2 2.5 3 5
run;

proc print data=Prob; 
format pGT E10.;
run;
