/* SAS program to accompany the article
   "How to simulate multivariate outliers"
   by Rick Wicklin, published 27MAR2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/03/27/simulate-multivariate-outliers.html

   This program shows how to create outliers in multivariate
   normal data that are a specified Mahalanobis distance from the 
   mean of the distribution,
*/

title;
ods graphics / width=400px height=400px attrpriority=none;
proc iml;
/* generate standardized uncorrelated data */
call randseed(123);
N = 200;
x = randfun(N//2, "Normal"); /* 2-dimensional data. X ~ MVN(0, I(2)) */

/* covariance matrix is product of diagonal (scaling) and correlation matrices
   See https://blogs.sas.com/content/iml/2010/12/10/converting-between-correlation-and-covariance-matrices.html
*/
R = {1   0.9,             /* correlation matrix */
     0.9 1  };
D = {4 9};                /* standard deviations */
Sigma = corr2cov(R, D);   /* Covariance matrix Sigma = D*R*D */
print Sigma;

/* U is the Cholesky transformation that scales and correlates the data */
U = root(Sigma);

/* add a few unusual points (outliers) */
pi = constant('pi');
t = T(0:11) * pi / 6;            /* evenly spaced points in [0, 2p) */
outliers = 5#cos(t) || 5#sin(t); /* evenly spaced on circle r=5 */
v = x // outliers;               /* concatenate MVN data and outliers */
w = v*U;                         /* transform from stdized to data coords */

labl = j(N,1," ") // T('0':'11');
Type = j(N,1,"Normal") // j(nrow(t),1,"Outlier");
title "Markers on Circle r=5";
title2 "Evenly spaced angles t[i] = i*pi/6";
call scatter(w[,1], w[,2]) grid={x y} datalabel=labl group=Type;

center = {0 0};
MD = MAHALANOBIS(w, center, Sigma);
obsNum = 1:nrow(w);
title "Mahalanobis Distance of MVN Data and Outliers";
call scatter(obsNum, MD) grid={x y} group=Type;


/********************************/

/* you can also add points at different distances */
t = {0, 3, 10} * pi / 6;       /* angular positions */
MD = {5, 6, 10};               /* MD for outlier */
outliers = MD#cos(t) ||  MD#sin(t);
v = x // outliers;             /* concatenate MVN data and outliers */
w = v*U;                       /* transform from stdized to data coords */

labl = j(N,1," ") // {'p1','p2','p3'};
Type = j(N,1,"Normal") // j(nrow(t),1,"Outlier");
call scatter(w[,1], w[,2]) grid={x y} datalabel=labl group=Type;
 

/********************************/

/* General technique to generate random outliers at distance r */
d = 2;                              /* use any value of d */
MD = {5, 6, 10};                    /* MD for outlier */
Y = randfun(nrow(MD)//d, "Normal"); /* d-dimensional Y ~ MVN(0, I(d)) */
outliers = MD # Y / sqrt(Y[,##]);   /* surface of spheres of radii MD */

v = x // outliers;             /* concatenate MVN data and outliers */
w = v*U;                       /* transform from stdized to data coords */

labl = j(N,1," ") // {'p1','p2','p3'};
Type = j(N,1,"Normal") // j(nrow(t),1,"Outlier");
title "MVN Data with Random Outliers";
title2 "Outliers at MD = {5, 6, 10}";
call scatter(w[,1], w[,2]) grid={x y} datalabel=labl group=Type;
