/* An introduction to copulas */

/* SAS program to accompany the article 
   "An introduction to simulating correlated data by using copulas"
   by Rick Wicklin, published 05JUL2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/07/05/introduction-copulas.html ‎

   This program explores the geometry of a copula. You can use copulas
   to simulate correlated multivariate data with specified marginal 
   distributions. The copula captures the correlation structure,
   independent from the marginal distributions.

   This article and example are based on Chapter 9 of 
   Simulating Data with SAS  (Wicklin, 2013)
   https://support.sas.com/en/books/authors/rick-wicklin.html
*/

proc iml;
N = 1e4;
call randseed(12345);
/* 1. Z ~ MVN(0, Sigma) */
Sigma = {1.0  0.6,
         0.6  1.0};
Z = RandNormal(N, {0,0}, Sigma);          /* Z ~ MVN(0, Sigma) */

/* 2. transform marginal variables to U(0,1) */
U = cdf("Normal", Z);                     /* U_i are correlated U(0,1) variates */

/* 3. construct the marginals however you wish */
gamma = quantile("Gamma", U[,1], 4);         /* gamma ~ Gamma(alpha=4)   */
LN = quantile("LogNormal", U[,2], 0.5, 0.8); /* LN ~ LogNormal(0.5, 0.8) */
X = gamma || LN;

/* the rank correlations for Z and U and X are exactly the same */
rankCorrZ = corr(Z, "Spearman")[2]; 
rankCorrU = corr(U, "Spearman")[2]; 
rankCorrX = corr(X, "Spearman")[2]; 
print rankCorrZ rankCorrU;
print rankCorrZ rankCorrX;

/* Not true for Pearson corr. 
   If Z~MVN(0,Sigma), corr(X) might not be close to Sigma,
   where X=(X1,X2,...,Xm) and X_i = F_i^{-1}(Phi(Z_i)) */
rhoZ = corr(Z, "Pearson")[2];
rhoX = corr(X, "Pearson")[2];
print rhoZ rhoX;

/* write simulated data to a SAS data set for visualization */
Q = Z||U||X;
labels = {Z1 Z2 U1 U2 X1 X2};
create CorrData from Q[c=labels];
append from Q;
close CorrData;


/* even though corr(X) ^= Sigma, you can often choose a target
   correlation, such as 0.6, and then choose Sigma so that corr(X)
   has the target correlation. */
/* re-run the example with new correlation 0.642 */
Sigma = {1.0    0.642,
         0.642  1.0};
newZ = RandNormal(N, {0,0}, Sigma);
newU = cdf("Normal", newZ);           /* columns of U are U(0,1) variates */
gamma = quantile("Gamma", newU[,1], 4);      /* gamma ~ Gamma(alpha=4) */
expo = quantile("Expo", newU[,2]);           /* expo ~ Exp(1)          */
newX = gamma || expo;

rhoZ = corr(newZ, "Pearson")[2];
rhoX = corr(newX, "Pearson")[2];
print rhoZ rhoX;
QUIT;

/* visualize the final data */
proc corr data=CorrData plots(maxpoints=none)=matrix(histogram);
   var x1 x2;     
run;

/* Create a scatter plot with marginal histograms, See
   https://blogs.sas.com/content/iml/2011/05/20/how-to-create-a-scatter-plot-with-marginal-histograms-in-sas.html
*/
proc template;
  define statgraph scatterhist;
  dynamic XVAR YVAR TITLE TRANSVAL;
  begingraph / designwidth=600px designheight=400px;
    entrytitle TITLE;
    layout lattice / rows=2 columns=2 rowweights=(.2 .8) columnweights=(.8 .2)
                     rowdatarange=union columndatarange=union rowgutter=0 columngutter=0;
    /* histogram at X2 axis position */
    layout overlay / walldisplay=(fill)
                     xaxisopts=(display=none) yaxisopts=(display=none offsetmin=0);
        histogram XVAR / binaxis=false;
    endlayout;
    /* upper right cell */  
    layout overlay; 
      entry " ";
    endlayout;
    /* scatter plot */
    layout overlay;         
      scatterplot y=YVAR x=XVAR / datatransparency=TRANSVAL markerattrs=(symbol=circlefilled);
    endlayout;
    /* histogram at Y2 axis position */
    layout overlay / walldisplay=(fill)
                     xaxisopts=(display=none offsetmin=0) yaxisopts=(display=none);
      histogram YVAR / orient=horizontal binaxis=false;
    endlayout;
    endlayout;    
  endgraph;
  end;
run;

/* visualize each step in the three-step process */
title;
ods graphics / width=500 height=500;
proc sgrender data=CorrData template=scatterhist;  
  dynamic XVAR="Z1" YVAR="Z2"  TRANSVAL=0.95
          TITLE="Multivariate Normal Data";
run;

proc sgrender data=CorrData template=scatterhist;  
  dynamic XVAR="U1" YVAR="U2"  TRANSVAL=0.95
          TITLE="The Copula Data";
run;

proc sgrender data=CorrData template=scatterhist;  
  dynamic XVAR="X1" YVAR="X2"  TRANSVAL=0.95
          TITLE="Gamma-Lognormal Data";
run;

/* visualize the copula, which is the 2-D cumulative distribution 
   for the uniform marginals. See
   https://blogs.sas.com/content/iml/2021/06/28/bivariate-cdf-sas.html
*/
proc kde data=CorrData;
   bivar U1 U2 / plots=NONE
         CDF out=KDEOut(rename=(distribution=kerCDF Value1=U1 Value2=U2));
run;

ods graphics / width=500px height=450px;
title "Visualization of the Gaussian Copula";
proc sgplot data=KDEOut aspect=1;
   label kerCDF="Bivariate CDF" U1="U1" U2="U2";
   heatmapparm x=U1 y=U2 colorresponse=kerCDF; 
run;



/*************************************/
/* Now do a four-dimensional example */
/*************************************/


proc iml;
N = 1e3;
call randseed(12345);
/* 1. Z ~ MVN(0, Sigma) */
Sigma = { 1.00  0.64 -0.70  0,
          0.64  1.00 -0.95  0,
         -0.70 -0.95  1.00 -0.25,
          0     0    -0.25   1.0};
Z = RandNormal(N, j(1,ncol(Sigma),0), Sigma);    /* Z ~ MVN(0, Sigma) */

/* 2. transform marginal variables to U(0,1) */
U = cdf("Normal", Z);           /* U_i are correlated U(0,1) variates */

/* visualize the copula */
create Copula from U[c=('U1':'U4')];  
append from U;  close;

/* 3. construct the marginals however you wish */
gamma = quantile("Gamma",  U[,1], 4);        /* ~ Gamma(alpha=4)      */
LN    = quantile("LogN",   U[,2], 0.5, 0.8); /* ~ LogNormal(0.5, 0.8) */
expo  = quantile("Expo",   U[,3], 1.5);      /* ~ Expo(1.5)           */
iGauss= quantile("IGauss", U[,4], 1, 0.8);   /* ~ InvGauss(1, 0.8)    */
X = gamma || LN || expo || iGauss;

/* the rank correlations for Z and X are exactly the same */
/*
rankCorrZ = corr(Z, "Spearman"); 
rankCorrX = corr(X, "Spearman"); 
print rankCorrZ rankCorrX;
*/
/* write to SAS data set */
create MVData from X[c={"Gamma" "LogN" "Expo" "IGauss"}];  
append from X;  close;
quit;


title "Visualization of a Multivariate Simulated Data";
proc sgscatter data=MVData;
   matrix Gamma LogN Expo IGauss / transparency=0.85 diagonal=(histogram);
run;

title "Visualization of a Copula";
proc sgscatter data=Copula;
   matrix U1 U2 U3 U4 / transparency=0.75 diagonal=(histogram);
run;

