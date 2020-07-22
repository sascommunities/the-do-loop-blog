/* SAS program to accompany the article 
   "Fit a multivariate Gaussian mixture model by using the expectation-maximization (EM) algorithm"
   by Rick Wicklin, published 23JUL2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/07/22/fit-multivariate-gaussian-mixture-em-algorithm.html ‎

   This program shows how to fit a Gaussian mixture model by using 
   model based clustering with "hard thresholding," as discussed in 
   the documentation for PROC MBC (SAS Visual Statistics) and in 
   Kessler, D. (2019) "Introducing the MBC Procedure for Model-Based Clustering"
   https://www.sas.com/content/dam/SAS/support/en/sas-global-forum-proceedings/2019/3016-2019.pdf
*/

/* initial exploration of the data */
title "Iris Data and Real Clusters";
title2 "And 95% Prediction Ellipse, assuming MVN";
proc sgplot data=Sashelp.Iris;
   ellipse x=PetalWidth y=SepalWidth / group=Species outline fill transparency=0.6;
   scatter x=PetalWidth y=SepalWidth / group=Species  transparency=0.2
                       markerattrs=(symbol=CircleFilled size=10) jitter;
   xaxis grid; yaxis grid;
run;

/* standardize and use k-means clustering (k=3) for initial guess */
proc stdize data=Sashelp.Iris out=StdIris method=std;
   var Sepal: Petal:;
run;

proc fastclus data=StdIris out=Clust maxclusters=3 maxiter=100 random=123;
   var Sepal: Petal:;
run;

data Iris;
merge Sashelp.Iris Clust(keep=Cluster);
/* for consistency with the Species order, remap the cluster numbers  */
if Cluster=1 then Cluster=2;
else if Cluster=2 then Cluster=1;
run;

title "k-Means Clustering of Iris Data";
proc sgplot data=Iris;
   ellipse x=PetalWidth y=SepalWidth / group=Cluster;
   scatter x=PetalWidth y=SepalWidth / group=Cluster transparency=0.2
                       markerattrs=(symbol=CircleFilled size=10) jitter;
   xaxis grid; yaxis grid;
run;

/****************************************************************/
/* HELPER FUNCTIONS */
/* Store two modules that will be called during the EM algorithm:
   The LogPdfMVN module is explained in 
   https://blogs.sas.com/content/iml/2020/07/15/multivariate-normal-log-likelihood.html

   The MLEstMVN module is explained in 
   https://blogs.sas.com/content/iml/2020/07/20/list-multivariate-statistics.html
*/
proc iml;
/* This function returns the log-PDF for a MVN(mu, Sigma) density at each row of X.
   The output is a vector with the same number of rows as X. */
start LogPdfMVN(X, mu, Sigma);
   d = ncol(X);
   log2pi = log( 2*constant('pi') );
   logdet = logabsdet(Sigma)[1];             /* sign of det(Sigma) is '+' */
   MDsq = mahalanobis(X, mu, Sigma)##2;      /* (x-mu)`*inv(Sigma)*(x-mu) */
   Y = -0.5 *( MDsq + d*log2pi + logdet );   /* log-PDF for each obs. Sum for LL */
   return( Y );
finish;

/* GIVEN Y and a classification variable Group that indicate group membership,
   return a list that contains:
   1. Number of observations in each group
   2. ML estimates of within-group mean
   3. ML estimate of within-group covariance 
   The k_th row of each matrix contains the statistics for the k_th group.
*/
start MLEstMVN(Y, Group);
   d = ncol(Y);
   u = unique(Group);
   G = ncol(u);
   ns    = J(G,   1, 0);           /* k_th row is number of obs in k_th group */
   means = J(G,   d, 0);           /* k_th row is mean of k_th group */
   covs  = J(G, d*d, 0) ;          /* k_th row is cov of k_th group */

   do k = 1 to G;
      X = Y[loc(Group=u[k]), ];    /* obs in k_th group */
      n = nrow(X);                 /* number of obs in this group */
      C = (n-1)/n * cov(X);        /* ML estimate of COV does not use the Bessel correction */
      ns[k]      = n;             
      means[k, ] = mean(X);        /* ML estimate of mean */
      covs[k, ]  = rowvec( C );    /* store COV in k_th row */
   end;
   outList = [#'n'=ns, #'mean'=means, #'cov'=Covs, #'group'=u]; /* pack into list */
   return (outList);
finish;
store module=(LogPdfMVN MLEstMVN);
QUIT;

/****************************************************************/
/* Demonstrate the EM algorithm:
   1. Assign observations to groups by k-means clustering.
   2. Compute group estimates and complete data LL (CD LL).
   3. Reassign the observations to groups based on the CD LL.
   4. Go to 2 and continue until the CD LL barely changes.

   This implementation uses hard clustering. 
*/
proc iml;
load module=(LogPdfMVN MLEstMVN);  /* you need to STORE these modules */

/* 1. Read data. Initialize 'Cluster' assignments from PROC FASTCLUS */
use Iris;
varNames = {'SepalLength' 'SepalWidth' 'PetalLength' 'PetalWidth'};
read all var varNames into X;
read all var {'Cluster'} into Group;  /* from PROC FASTCLUS */
close;

nobs = nrow(X); d = ncol(X); G = ncol(unique(Group));
prevCDLL = -constant('BIG');    /* set to large negative number */
converged = 0;                  /* iterate until converged=1 */
eps = 1e-5;                     /* convergence criterion */
iterHist = j(100, 3, .);        /* monitor EM iteration */
LL = J(nobs, G, 0);             /* store the LL for each group */

/* EM algorithm: Solve the M and E subproblems until convergence */
do nIter = 1 to 100 while(^converged);
   /* 2. M Step: Given groups, find MLE for n, mean, and cov within each group */
   L = MLEstMVN(X, Group);      /* function returns a list */
   ns = L$'n';       tau = ns / sum(ns);
   means = L$'mean'; covs = L$'cov';  

   /* 3. E Step: Given MLE estimates, compute the LL for membership 
         in each group. Use LL to update group membership. */
   do k = 1 to G;
      LL[ ,k] = log(tau[k]) + LogPDFMVN(X, means[k,], shape(covs[k,],d));
   end;
   Group = LL[ ,<:>];               /* predicted group has maximum LL */

   /* 4. The complete data LL is the sum of log(tau[k]*PDF[,k]).
         For "hard" clustering, Z matrix is 0/1 indicator matrix.
         DESIGN function: https://blogs.sas.com/content/iml/2016/03/02/dummy-variables-sasiml.html
   */
   Z = design(Group);               /* get dummy variables for Group */
   CDLL = sum(Z # LL);              /* sum of LL weighted by group membership */
   /* compute the relative change in CD LL. Has algorithm converged? */
   relDiff = abs( (prevCDLL-CDLL) / prevCDLL );
   converged = ( relDiff < eps );               

   /* monitor convergence; if no convergence, iterate */
   prevCDLL = CDLL;
   iterHist[nIter,] = nIter || CDLL || relDiff;
end;

/* remove unused rows and print EM iteration history */
iterHist = iterHist[ loc(iterHist[,2]^=.), ];  
print iterHist[c={'Iter' 'CD LL' 'relDiff'}];

/* print final parameter estimates for Gaussian mixture */
GroupNames = strip(char(1:G));
rows = repeat(T(1:d), 1, d);  cols = repeat(1:d, d, 1);
SNames = compress('S[' + char(rows) + ',' + char(cols) + ']');
print tau[r=GroupNames F=6.2],
      means[r=GroupNames c=varNames F=6.2],
      covs[r=GroupNames c=(rowvec(SNames)) F=6.2];

create FinalClusters var "Group";
append; close;
QUIT;

/* visualize the data and the assigned clusters */
data PlotEM;
merge Iris FinalClusters;
diff = (Cluster ^= Group);  /* did obs change from initial group? */
run;
proc sort data=PlotEM; by Group; run;  /* for graphing */

title "Clusters Determined by EM";
proc sgplot data=PlotEM;
   ellipse x=PetalWidth y=SepalWidth / group=Group;
   scatter x=PetalWidth y=SepalWidth / group=Group transparency=0.2
                       markerattrs=(symbol=CircleFilled size=10) jitter;
   xaxis grid; yaxis grid;
run;

proc sgscatter data=PlotEM;
   matrix SepalLength SepalWidth PetalLength PetalWidth / group=Group;
run;

/* how many obs changed groups? Which ones? */
proc means data=PlotEM sum; var diff; run;
title "Observations that Changed Clusters";
ods graphics / attrpriority=none;
proc sgplot data=PlotEM;
   styleattrs datasymbols=(circle X);
   ellipse x=PetalWidth y=SepalWidth / group=Group;
   scatter x=PetalWidth y=SepalWidth / group=diff;
   xaxis grid; yaxis grid;
run;
ods graphics / attrpriority=color;

/**********************************************************/
