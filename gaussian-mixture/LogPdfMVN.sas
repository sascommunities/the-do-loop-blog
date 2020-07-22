/* SAS program to accompany the article 
   "How to evaluate the multivariate normal log-likelihood"
   by Rick Wicklin, published 15JUL2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/07/15/multivariate-normal-log-likelihood.html

   This program shows how to efficiently evaluate the log-PDF of the 
   multivariate normal distribution. You can use the same function to 
   compute the log lkelihood: just evaluate the function on a data
   matrix and sum the log-PDF values.
*/
proc sgscatter data=sashelp.iris; 
matrix SepalLength--PetalLength /group=Species diagonal=(histogram kernel);
run;

/* log likelihood for each observation of MVN data in SAS/IML */
proc corr data=Sashelp.Iris COV plots=Scatter;
where species="Setosa";
var SepalLength SepalWidth;
run;

title "Iris Data and 95% Prediction Ellipse";
title2 "Assuming Multivariate Normality";
proc sgplot data=Sashelp.Iris noautolegend;
where species="Setosa";
scatter x=SepalLength y=SepalWidth / jitter;
ellipse x=SepalLength y=SepalWidth;
run;

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

/* read the iris data for the Setosa species */
use Sashelp.Iris where(species='Setosa');
read all var {'SepalLength' 'SepalWidth'} into X;
close;

n = nrow(X);           /* assume no missing values */
m = mean(X);           /* maximum likelihood estimate of mu */
S = (n-1)/n * cov(X);  /* maximum likelihood estimate of Sigma */
/* evaluate the log likelihood for each observation */
LL = LogPdfMVN(X, m, S);

   d = ncol(X);
   log2pi = log( 2*constant('pi') );
   logdet = logabsdet(S)[1];             /* sign of det(Sigma) is '+' */
   maxLL = -0.5 *(d*log2pi + logdet );  /* loglikelihood for each obs */
   print maxLL;


/* What if we use "wrong" parameters? For example, if the covariance 
   matrix indicates negative correlation? */
m2 = {45 30};
S2 = {12 -10,  -10 14};
LL_Wrong = LogPdfMVN(X, m2, S2);


/* If you want the total log likelihood, compute sum(LL) over all obs */
TotalLL = sum(LL);
TotalLL_Wrong = sum(LL_Wrong);
print TotalLL TotalLL_Wrong;


/* the largest the LL can be is when x=mu. For that value of x, the LL is
   -0.5 *(d*log2pi + logdet)
*/
   d = ncol(X);
   log2pi = log( 2*constant('pi') );
   logdet = logabsdet(S)[1];
   LargestLL = -0.5 *( d*log2pi + logdet ); 
   print LargestLL;


A = X || LL || LL_Wrong;
create MVN from A[c={'x1' 'x2' 'LL' 'LL_Wrong'}];
append from A;
close;


proc means data=MVN Min Max;
var LL LL_Wrong;
run;
proc sql;
select floor(min(min(LL), min(LL_Wrong))), 
       ceil(max(max(LL), max(LL_Wrong))) into :min, :max 
from MVN; 
quit; 

/* add fake  min and max values so that each colorramp will use the same scale. */
data Fake;
LL = &max; LL_Wrong = &max; output;
LL = &min; LL_Wrong = &min; output;
run;


data Plot;
set MVN Fake;
run;
/*
proc iml;
p = palette('SPECTRAL', 5);
print p;
*/
title "Log likelihood of Data";
title2 "Correct Model Parameters";
proc sgplot data=Plot;
scatter x=x1 y=x2 / colorresponse=LL
                    colormodel=(CXD7191C CXFDAE61 CXFFFFBF CXABDDA4 CX2B83BA )
                    markerattrs=(symbol=CircleFilled size=10) jitter;
xaxis grid; yaxis grid;
label x1 = 'SepalLength' x2 = 'SepalWidth';
run;

title "Log likelihood of Data";
title2 "Incorrect Model Parameters";
proc sgplot data=Plot;
scatter x=x1 y=x2 / colorresponse=LL_Wrong 
                    colormodel=(CXD7191C CXFDAE61 CXFFFFBF CXABDDA4 CX2B83BA )
                    markerattrs=(symbol=CircleFilled size=10) jitter;
xaxis grid; yaxis grid;
label x1 = 'SepalLength' x2 = 'SepalWidth';
run;

