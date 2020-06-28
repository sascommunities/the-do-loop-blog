/* SAS program to accompany the article 
   "Pooled, within-group, and between-group covariance matrices"
   by Rick Wicklin, published 01JUL2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/07/01/pooled-covariance-between-group.html

   This program shows how to compute the 
   within-group, pooled, and between-group 
   CSSCP and convariance matrices for a multivariate analysis.

   If you assume is that the covariance is 
   constant across different groups in the data,
   the pooled covariance is an estimate of the common covariance.
   It is a weighted average of the sample within-group covariances,
   where larger groups are weighted more heavily than smaller groups.

   Data is Fisher's Iris data.
*/

/* Visualize all pairs of dimensions */
ods graphics/ width=600px height=600px;
title "68% Prediction Ellipses for Iris Data";
proc sgscatter data=sashelp.iris;
  matrix SepalLength SepalWidth PetalLength PetalWidth / group=species;

run;

/* Visualize two pairs of variables. Add prediction ellipse for groups */
ods graphics/ width=300px height=300px;
proc sgplot data=Sashelp.Iris;
   scatter x=SepalLength y=SepalWidth / group=Species transparency=0.5;
   ellipse x=SepalLength y=SepalWidth / group=Species alpha=0.32 lineattrs=(thickness=2);
run;
proc sgplot data=Sashelp.Iris;
   scatter x=PetalLength y=PetalWidth / group=Species transparency=0.5;
   ellipse x=PetalLength y=PetalWidth / group=Species alpha=0.32 lineattrs=(thickness=2);
run;

/* Compute within-group, pooled, and between-group covariance (and CSSCP and Correlation) */
proc discrim data=sashelp.iris method=normal pool=yes outstat=Cov noprint;
   class Species;
   var SepalLength SepalWidth PetalLength PetalWidth;
run;

proc print data=Cov noobs; 
   where _TYPE_ = "PCOV";
   format _numeric_ 6.2;
   var _TYPE_ _NAME_ Sepal: Petal:;
run;

proc print data=Cov noobs; 
   where _TYPE_ = "COV" and Species^=" ";
   format _numeric_ 6.2;
run;

proc print data=Cov noobs; 
   where _TYPE_ = "BCOV";
   format _numeric_ 6.2;
run;

/* Reproduce the computation in SAS/IML */

/* Compute a pooled covariance matrix when observations
   belong to k groups with sizes n1, n2, ..., nk
   where n1+n2+...+nk = N
*/
proc iml;
varNames = {'SepalLength' 'SepalWidth' 'PetalLength' 'PetalWidth'};
use Sashelp.iris;
   read all var varNames into Z;
   read all var "Species" into Group;
close;

/* assume complete cases, otherwise remove rows with missing values */
N   = nrow(Z);

/* You can compute the within-group covariance by 
   computing the covariance for the observations in each group/
*/
u = unique(Group);
k = ncol(u);                /* number of groups */
p   = ncol(varNames);       /* number of variables */
M = j(p, p, 0);             /* sum of within-group CSCCP matrices */
do i = 1 to k;
   idx = loc(Group = u[i]);
   X = Z[idx,];             /* extract obs for i_th group */
   n_i = nrow(X);           /* n_i = size of i_th group */
   S   = cov(X);            /* within-group cov */
   /* accumulate the weighted sum of within-group covariances */
   M = M + (n_i-1) * S;     /* (n_i-1)*S is centered X`*X */
end;

/* The pooled covariance is an average of the 
   within-class SSCP matrices.
*/
Sp = M / (N-k);
print Sp[L="Pooled Cov" c=varNames r=VarNames format=6.2];

/* The between-class CSSCP is the difference between total 
   CSSCP and the sum of the within-group CSSCPs.
   SAS doc for PROC DISCRIM defines the between-class 
   covariance matrix as 
   the between-class SSCP matrix divided by N*(k-1)/k, 
   where N is the number of observations and k 
   is the number of classes.
*/
/* the total covariance matrix ignores the groups */
C = (N-1)*cov(Z);    
BCSSCP = C - M;   /* between = Full - Sum(Within) */
BCov = BCSSCP * k/( N*(k-1) );
print BCov[L="Between Cov" c=varNames r=VarNames format=6.2];

QUIT;
/******************************************/

/* Use SAS/IML to create 68% prediction ellipses for the pooled covariance 
   See https://blogs.sas.com/content/iml/2014/07/23/prediction-ellipses-from-covariance.html
*/
proc iml;
start PredEllipseFromCov(m, S, n, ConfLevel=0.95, nPts=128);
/* Compute prediction ellipses centered at m from the 2x2 covariance matrix S.
   The return matrix is a matrix with three columns.
   Col 1: The X coordinate of the ellipse for the confidence level.
   Col 2: The Y coordinate of the ellipse for the confidence level.
   Col 3: The confidence level. Use this column to extract a single ellipse.
 
   Input parameters:
   m           a 1x2 vector that specifies the center of the ellipse. 
               This can be the sample mean or median.
   S           a 2x2 symmetric positive definite matrix. This can be the 
               sample covariance matrix or a robust estimate of the covariance.
   n           The number of nonmissing observations in the data. 
   ConfLevel   a 1 x k vector of (1-alpha) confidence levels that determine the
               ellipses. Each entry is 0 < ConfLevel[i] < 1.  For example,
               0.95 produces the 95% confidence ellipse for alpha=0.05.
   nPts       the number of points used to draw the ellipse. The default is 0.95.
*/
 
   /* parameterize standard ellipse in coords of eigenvectors */
   call eigen(lambda, evec, S);   /* eigenvectors are columns of evec */
   t = 2*constant("Pi") * (0:nPts-1) / nPts;
   xy  = sqrt(lambda[1])*cos(t) // sqrt(lambda[2])*sin(t);
   stdEllipse = T( evec * xy );   /* transpose for convenience */
 
   /* scale the ellipse; see documentation for PROC CORR */
   c = 2 * (n-1)/n * (n+1)/(n-2);          /* adjust for finite sample size */
   p = rowvec(ConfLevel);                  /* row vector of confidence levels */
   F = sqrt(c * quantile("F", p, 2, n-2)); /* p = 1-alpha */
 
   ellipse = j(ncol(p)*nPts, 3);  /* 3 cols: x y p */
   startIdx = 1;                  /* starting index for next ellipse */
   do i = 1 to ncol(p);           /* scale and translate */
      idx = startIdx:startIdx+nPts-1;
      ellipse[idx, 1:2] = F[i] # stdEllipse + m; 
      ellipse[idx, 3] = p[i];     /* use confidence level as ID */
      startIdx = startIdx + nPts;
   end;           
   return( ellipse );
finish;

use Cov where (_TYPE_ = "PCOV");
read all var _NUM_ into COV[c=varNames];
close;
*print COV[r=varNames c=varNames];

use Cov where (_TYPE_ = "MEAN" & Species=" ");
read all var _NUM_ into grandmean;
close;

mean = grandmean[,1:2];
C = Cov[1:2, 1:2];  /* sepal length and width */
/* The N is supposed to be the number of observations in the group.
   The pooled covariance is for all obs ....
   but maybe we need to subtract off k-1 for DF?
   I am NOT SURE WHAT value of N to use: 150? 148? 147?
*/
N = 148;
E = PredEllipseFromCov(mean, C, N, 0.68); /* 68% ellipse */
Name = j(nrow(E), 1, "Pooled");
create Ellipse from E[r=Name c={"xx" "yy" "ConfLevel"}];
append from E[r=Name];
close Ellipse; 

data Ell;
set Sashelp.Iris Ellipse;
run;
ods graphics / width=400px height=400px;
title2 "Shaded Region Shows Pooled Covariance";
proc sgplot data=Ell;
   polygon x=xx y=yy id=Name / fill;
   scatter x=SepalLength y=SepalWidth / group=Species transparency=0.5;
   ellipse x=SepalLength y=SepalWidth / group=Species alpha=0.32 lineattrs=(thickness=2);
run;

