title;
ods graphics / reset;

/* SAS program to accompany the article 
   "Simulate correlated variables by using the Iman-Conover transformation"
   by Rick Wicklin, published 14JUN2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/06/14/simulate-iman-conover-transformation.html ‎

   This program shows how to define and call the Iman-Conover transformation, 
   which induces a specified correlation structure onto an arbitrary set of
   multivariate data.

   Iman and Conover, 1982, "A distribution-free approach to inducing rank 
   correlation among input variables", 
   URL: https://www.tandfonline.com/doi/abs/10.1080/03610918208812265
*/

%let N = 100; /* sample size */
/* independently simulate variables from known distributions */
data SimIndep;
call streaminit(12345);
do i = 1 to &N;
   Normal    = rand("Normal");
   Lognormal = rand("Lognormal");
   Expon     = rand("Expo");
   Uniform   = rand("Uniform");
   output;
end;
drop i;
run;

proc corr data=SimIndep Spearman noprob plots=matrix(hist);
   var Normal Lognormal Expon Uniform;
run;

/* define a store a SAS/IML function that computes the Iman-Conover transformation */
proc iml;
/* Input: X is a data matrix with k columns
          C is a (k x k) rank correlation matrix
   Output: A matrix W. The i_th column of W is a permutation of the i_th
          columns of X. The rank correlation of W is close to the 
          specified correlation matrix, C.
*/          
start ImanConoverTransform(X, C);
   N = nrow(X);
   S = J(N, ncol(X));
   /* T1: Create normal scores of each column */
   do i = 1 to ncol(X);
      ranks = ranktie(X[,i], "mean");          /* tied ranks */
      S[,i] = quantile("Normal", ranks/(N+1)); /* van der Waerden scores */
   end;
   /* T2: apply two linear transformations to the scores */
   CS = corr(S);        /* correlation of scores */
   Q = root(CS);        /* Cholesky root of correlation of scores */
   P = root(C);         /* Cholesky root of target correlation */
   T = solve(Q,P);      /* same as  T = inv(Q) * P; */
   Y = S*T;             /* transform scores: Y has rank corr close to target C */

   /* T3: Permute or reorder data in the columns of X to have the same ranks as Y */
   W = X; 
   do i = 1 to ncol(Y);
      rank = rank(Y[,i]);          /* use ranks as subscripts, so no tied ranks */
      tmp = W[,i]; call sort(tmp); /* sort column by ranks */
      W[,i] = tmp[rank];           /* reorder the column of X by the ranks of M */
   end;
   return( W );
finish;
store module=(ImanConoverTransform);  /* store definition for later use */
quit;

/* Use Iman-Conover method (1982) to generate MV data with known marginals
   and approximate rank correlation. */
proc iml;
load module=(ImanConoverTransform);  /* load the function definition */

/* Step 1: Read in the data */
varNames = {'Normal' 'Lognormal' 'Expon' 'Uniform'};
use SimIndep; read all var varNames into X; close;

/* Step 2: specify target rank correlation */
/*    X1    X2    X3    X4                 */
C = { 1.00  0.75 -0.70  0,           /* X1 */
      0.75  1.00 -0.95  0,           /* X2 */
     -0.70 -0.95  1.00 -0.3,         /* X3 */
      0     0    -0.3   1.0};        /* X4 */

/* Step 3: Strategically reorder the columns of X. 
           The new ordering has the specified rank corr (approx)
           and preserves the original marginal distributions. */
W = ImanConoverTransform(X, C);
RankCorr = corr(W, "Spearman");
print RankCorr[format=6.3 r=varNames c=varNames];

/* write new data to a SAS data set */
newvarNames = {'newNormal' 'newLognormal' 'newExpon' 'newUniform'};
create SimCorr from W[c=newvarNames];  append from W;  close;
quit;

proc corr data=SimCorr Spearman noprob plots=matrix(hist);
   var newNormal newLognormal newExpon newUniform;
run;



/**************************************/
/* does it work for categorical vars? Kind of, but not really.
   Multinomial distribution has limitations on the feasible correlations.
*/
data CatIndep;
call streaminit(12345);
do i = 1 to &N;
   C1 = rand("Bern", 0.6);
   C2 = rand("Table", 0.3, 0.1, 0.2, 0.1, 0.3);
   C3 = rand("Table", 0.5, 0.3, 0.2);
   C4 = rand("Table", 0.1, 0.4, 0.5); 
   output;
end;
drop i;
run;

proc iml;
load module=(ImanConoverTransform);
/* Step 1: Read in the data */
varNames = 'C1':'C4';
use CatIndep; read all var varNames into X; close;

/* Step 2: specify target rank correlation to induce */
C = { 1.00  0.75 -0.70  0,
      0.75  1.00 -0.95  0,
     -0.70 -0.95  1.00 -0.3,
      0     0    -0.3   1.0};

/* Step 3: Strategically reorder the columns of X. 
           The new ordering has the specified rank corr (approx)
           and preserves the original marginal distributions. */
Y = ImanConoverTransform(X, C);
RankCorr = corr(Y, "Spearman");
print RankCorr[format=6.3];
quit;
/**************************************/
