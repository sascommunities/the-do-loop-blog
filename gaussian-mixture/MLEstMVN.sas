/* SAS program to accompany the article 
   "Compute within-group multivariate statistics and store them in a list"
   by Rick Wicklin, published 20JUL2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/07/20/list-multivariate-statistics.html

   This program shows how to use a list to return several 
   multivariate statistics for several groups. 
   A SAS/IML module computes several statistics for G groups.
   There are two ways to return the statistics:
   1. If feasible, you can return a list that contains matrices.
      In each matrix, the k_th row contains the data for the k_th group.
   2. You can return a list that contains G sublists.
      Each sublist contains the statistics for the corresponding group.
*/
ods graphics/attrpriority=color;
title "Iris Data, Colored by Species";
title2 "95% Prediction Ellipse, assuming MVN";
proc sgplot data=Sashelp.Iris;
   ellipse x=PetalWidth y=SepalWidth / group=Species outline fill transparency=0.6;
   scatter x=PetalWidth y=SepalWidth / group=Species transparency=0.2
                       markerattrs=(symbol=CircleFilled size=10) jitter;
   xaxis grid; yaxis grid;
run;

proc iml;
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

/* Example: read the iris data  */
varNames = {'PetalWidth' 'SepalWidth'};
use Sashelp.Iris;
read all var varNames into X;
read all var {'Species'};
close;

L = MLEstMVN(X, Species);            /* call function; return named list of results */
ns = L$'n';                          /* extract obs numbers          */
means = L$'mean';                    /* extract mean vectors         */
Covs = L$'cov';                      /* extract covariance matrix    */
lbl = L$'group';
print ns[r=lbl c='n'], 
      means[r=lbl c={'m1' 'm2'}], 
      covs[r=lbl c={C11 C12 C21 C22}];

Setosa_Cov = shape(covs[1,], 2);

/* Alternative: Return a list of G sublists. Each sublist has only 
   the statistics relevant to the corresponding group.    */
start MLEst2(Y, Group);
   u = unique(Group);
   G = ncol(u);
   outList = ListCreate(G);

   do k = 1 to G;
      X = Y[loc(Group=u[k]), ];    /* obs in k_th group */
      n = nrow(X);                 /* number of obs in k_th group */
      SL = [#'n'=n, 
            #'mean'=mean(X),        /* ML estimate of mean */
            #'cov'=(n-1)/n*cov(X)]; /* ML estimate of COV */
      call ListSetItem(outList, k, SL);   /* set sublist item */
   end;
   call ListSetName(outList, 1:G, u);
   return (outList);
finish;

L2 = MLEst2(X, Species);
package load ListUtil;  /* load a package to print lists */
call Struct(L2);        /* or: call ListPrint(L2);       */

SetosaList = L2$'Setosa';    /* get the list of statistics for the 'Setosa' group */
m = SetosaList$'mean';       /* get the mean vector */
C = SetosaList$'cov';        /* get the MLE covariance matrix */
print m[c=varNames], C[c=varNames r=varNames];

QUIT;

