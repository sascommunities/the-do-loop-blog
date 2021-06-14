/* SAS program to accompany the article 
   "The geometry of the Iman-Conover transformation"
   by Rick Wicklin, published 16JUN2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/06/16/geometry-iman-conover-transformation.html

   This program shows the geometry ot the Iman-Conover transformation, 
   which induces a specified correlation structure onto an arbitrary set of
   multivariate data. The Iman-Conover transformation is defined at
   https://blogs.sas.com/content/iml/2021/06/14/simulate-iman-conover-transformation.html‎

   REFERENCE: 
   Iman and Conover, 1982, "A distribution-free approach to inducing rank 
   correlation among input variables", 
   URL: https://www.tandfonline.com/doi/abs/10.1080/03610918208812265
*/

/*******************************/
/* How does Iman-Conover work? */
/*******************************/
ods graphics / reset;
title;

/* Let's look at the geometry in terms of 2-D Cholesky transforms. */
/* Start with continuous data that is not normal. */


/* simulate data for example from known distributions */
/*
%let N = 1000;
data SimData;
call streaminit(12345);
do i = 1 to &N;
   X = rand("Lognormal", 0.5);
   Y     = 0.3*LogNormal + sqrt(rand("Expo", 2));
   if X<20 & Y< 20 then 
      output;
end;
drop i;
run;

proc corr data=SimData spearman 
         plots(maxpoints=none)=matrix(hist) ;
var LogNormal Expon;
run;
*/

/* ... or use real data */
data Orig;
set sashelp.cars(rename=(EngineSize=X MPG_Highway=Y));
keep X Y;
label X= Y=;
run;

proc corr data=Orig spearman nosimple noprob 
         plots(maxpoints=none)=matrix(hist) ;
   var X Y;
run;

proc iml;
varNames = {"X" "Y"};
use Orig; read all var varNames into X; close;
corrX = corr(X, "Spearman")[2];

/* T1: Create normal scores of each column 
   Columns of S have exactly the same ranks as columns of X.
   ==> Spearman corr is the same    */
N = nrow(X);
S = J(N, ncol(X));
do i = 1 to ncol(X);
   ranks = ranktie(X[,i], "mean");          /* tied ranks */
   S[,i] = quantile("Normal", ranks/(N+1)); /* van der Waerden scores */
end;

corrScores = corr(S, "Spearman")[2];
print corrX, corrScores;

/* T2: use the inverse Cholesky transform of the scores to uncorrelate */
CS = corr(S);        /* correlation of scores */
Q = root(CS);        /* Cholesky root of correlation of scores */
Z = S*inv(Q);
corrZ = corr(Z, "Spearman")[2];
print corrZ;

/* T3: use the Cholesky transform of the target matrix to induce correlation */
/* define the target correlation */
C = {1    0.6,
     0.6 1  };
P = root(C);         /* Cholesky root of target correlation */
Y = Z*P;
corrY = corr(Y, "Spearman")[2];
print corrY;

/* these operations have changed the marginal distributions
   of the variables. */

/* T4: Permute or reorder data in the columns of X to have the same ranks as Y */
W = X; 
do i = 1 to ncol(Y);
   rank = rank(Y[,i]);          /* use ranks as subscripts, so no tied ranks */
   tmp = W[,i]; call sort(tmp); /* sort column by ranks */
   W[,i] = tmp[rank];           /* reorder the column of X by the ranks of Y */
end;

corrW = corr(W, "Spearman")[2];
print corrW;

/* Iman and Conover assume no duplicate values. Look at the number of unique values
   in each column of the matrices: */
nuX = ncol(unique(X[,1])) || ncol(unique(X[,2]));
nuS = ncol(unique(S[,1])) || ncol(unique(S[,2]));
nuZ = ncol(unique(Z[,1])) || ncol(unique(Z[,2]));
nuY = ncol(unique(Y[,1])) || ncol(unique(Y[,2]));
nuW = ncol(unique(W[,1])) || ncol(unique(W[,2]));
numUnique = nuX // nuS // nuZ // nuY // nuW;
print numUnique[r={X S Z Y W} c={"1st Col" "2nd Col"}];

/************************************************/
/* VISUALIZE THE GEOMETRY OF THE TRANSFORMATION */
/************************************************/

/* write all variables to data set for visualization */
M = X || S || Z || Y || W;
name = {'X1' 'X2' 'Score1' 'Score2' 'Z1' 'Z2' 'Y1' 'Y2' 'W1' 'W2'};
create ImanConover from M[c=name];
append from M;
close;
QUIT;

/* create all graphs for the geometry diagram */
ods graphics / width=200px height=200px;
title "Original Data";
title2 "Rank Correlation = -0.77";
proc sgplot data=ImanConover;
   scatter x=X1 y=X2;
   xaxis grid; yaxis grid;
run;

ods graphics / width=160px height=200px;
title "Normal Scores";
title2 "Rank Correlation = -0.77";
proc sgplot data=ImanConover aspect=1;
   scatter x=Score1 y=Score2;
   refline 0 / axis=x; refline 0 / axis=y;
   xaxis grid; yaxis grid;
run;

title "Uncorrelated Scores";
title2 "Rank Correlation ~ 0";
proc sgplot data=ImanConover aspect=1;
   scatter x=Z1 y=Z2;
   refline 0 / axis=x; refline 0 / axis=y;
   xaxis grid; yaxis grid;
run;

title "Induce New Correlation";
title2 "Rank Correlation ~ 0.6";
proc sgplot data=ImanConover aspect=1;
   scatter x=Y1 y=Y2;
   refline 0 / axis=x; refline 0 / axis=y;
   xaxis grid; yaxis grid;
run;

ods graphics / width=200px height=200px;
title "Restore Marginal Distributions";
title2 "Rank Correlation ~ 0.6";
proc sgplot data=ImanConover;
   scatter x=W1 y=W2;
   xaxis grid; yaxis grid;
run;

title;
*ods graphics / reset;
