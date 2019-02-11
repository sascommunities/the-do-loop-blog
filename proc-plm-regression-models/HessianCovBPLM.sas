/* SAS program to accompany the articles:

   "4 reasons to use PROC PLM for linear regression models in SAS"
   by Rick Wicklin, published 11FEB2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/02/11/proc-plm-regression-models-sas.html

   This program shows how to us PROC PLM to score a regression model, visualize a model,
   compute estimates and hypothesis tests, and display statistics.

ALSO

   "3 ways to obtain the Hessian at the MLE solution for a regression model"
   by Rick Wicklin, published 13FEB2019 on The DO Loop blog:
    https://blogs.sas.com/content/iml/2019/02/13/hessian-maximum-likelihood-covariance.html

   This program shows three ways to get the Hessian matrix for an MLE solution.

   Data from the PROC LOGISTIC documentation.
*/

Data Neuralgia;
   input Treatment $ Sex $ Age Duration Pain $ @@;
   PainY = (Pain^='Yes');
   datalines;
P F 68  1 No  B M 74 16 No  P F 67 30 No  P M 66 26 Yes B F 67 28 No  B F 77 16 No
A F 71 12 No  B F 72 50 No  B F 76  9 Yes A M 71 17 Yes A F 63 27 No  A F 69 18 Yes
B F 66 12 No  A M 62 42 No  P F 64  1 Yes A F 64 17 No  P M 74  4 No  A F 72 25 No
P M 70  1 Yes B M 66 19 No  B M 59 29 No  A F 64 30 No  A M 70 28 No  A M 69  1 No
B F 78  1 No  P M 83  1 Yes B F 69 42 No  B M 75 30 Yes P M 77 29 Yes P F 79 20 Yes
A M 70 12 No  A F 69 12 No  B F 65 14 No  B M 70  1 No  B M 67 23 No  A M 76 25 Yes
P M 78 12 Yes B M 77  1 Yes B F 69 24 No  P M 66  4 Yes P F 65 29 No  P M 60 26 Yes
A M 78 15 Yes B M 75 21 Yes A F 67 11 No  P F 72 27 No  P F 70 13 Yes A M 75  6 Yes
B F 65  7 No  P F 68 27 Yes P M 68 11 Yes P M 67 17 Yes B M 70 22 No  A M 65 15 No
P F 67  1 Yes A M 67 10 No  P F 72 11 Yes A F 74  1 No  B M 80 21 Yes A F 69  3 No
;

/* store model into item store */
title 'Logistic Model on Neuralgia';
proc logistic data=Neuralgia;
   class Sex Treatment;
   model Pain(Event='Yes')= Sex Age Duration Treatment;
   store PainModel / label='Neuralgia Study';  /* or use mylib.PaimModel for permanent storage */
   ods select ParameterEstimates;
run;

/* 1.Use PLM to score new obs */
data NewPatients;
   input Treatment $ Sex $ Age Duration;
   datalines;
B F 63  5 
B F 79 16 
A M 74 12 
;

proc plm restore=PainModel;
   score data=NewPatients out=NewScore predicted LCLM UCLM / ilink;
run;

proc print data=NewScore;
run;

/* 2. Use PROC PLM to create an effect plot */
proc plm restore=PainModel;
   effectplot slicefit(x=Age sliceby=Treatment plotby=Sex);
run;

/* 3. Use PROC PLM to create contrasts and estimates */
proc plm restore=PainModel;
   /* 'Exponentiated' column is odds ratio between treatment or gender */
   estimate 'Pairwise A vs B' Treatment 1 -1 / exp CL; 
run;

/* 4. Use PROC PLM to show statistics or the original program */
proc plm restore=PainModel;
   show Parameters COVB Program;
run;


/***************************************************/
/*    
  3 ways to obtain the Hessian at the MLE solution for a regression model
*/
/***************************************************/

/* 1. use PROC PLM, you can get the Hessian matrix evaluated at the MLE */
proc plm restore=PainModel;
   show Hessian CovB;
   ods output Cov=CovB;
run;

/* Why isn't the Hessian part of the PROC LOGISTIC output? It is, in a way,
   because you can request the COVB matrix.
   The connection between Hessian and COVB is that 
   the inverse of the Hessian is an estimator of the asymptotic covariance matrix. 
   Therefore, you can estimate the standard errors of the parameter estimates 
   by using the square root of the diagonal elements of the COVB matrix.
*/

/* 2. Show that you can get the Hessian as the inverse of the COVB matrix
      and vice versa */

proc iml;
use CovB nobs p;                         /* read number of obs (p) */
   cols = "Col1":("Col"+strip(char(p))); /* variable names are Col1 - Colp */
   read all var cols into Cov;           /* read COVB matrix */
   read all var "Parameter";             /* read names of parameters */
close;

/* Hessian and covariance matrices are inverses */
Hessian = inv(Cov);
print Hessian[r=Parameter c=Parameter F=BestD8.4];

v = eigval(Hessian); /* show Hessian is positive definite */
print v;

/* incidentally, stderr are the sqrt of diagonal elements */
stderr = sqrt(vecdiag(Cov)); 
*print stderr;
quit;

/* 3. Use PROC NLMIXED to set up and solve a custom MLE problem.
      Use sme data and model. Because NLMIXED has no CLASS statement,
      get design matrix from PROC LOGISTIC and ass a numerica binary response. */
/* output design matrix and EFFECT parameterization */
proc logistic data=Neuralgia outdesign=Design outdesignonly;
   class Pain Sex Treatment;
   model Pain(Event='Yes')= Sex Age Duration Treatment; /* use NOFIT option for design only */
run;
/* PROC NLMIXED required a numeric response */
data Design;
   set Design;
   PainY = (Pain='Yes');  
run;

ods exclude IterHistory;
proc nlmixed data=Design COV HESS;
   parms b0 -18 bSexF bAge bDuration bTreatmentA bTreatmentB 0;
   eta    = b0 + bSexF*SexF + bAge*Age + bDuration*Duration +
                 bTreatmentA*TreatmentA + bTreatmentB*TreatmentB;
   p = logistic(eta);       /* or 1-p to predict the other category */
   model PainY ~ binary(p);
run;


/* 4. You can do the same thing with other procedures
      and other parameterizations */
proc genmod data=Neuralgia;
   class Sex Treatment / param=EFFECT;
   model PainY= Sex Age Duration Treatment / COVB dist=bin;
   store GMPainModel / label='Neuralgia Study';  /* or use mylib.PaimModel for permanent storage */
   ods select COV;
run;

proc plm restore=GMPainModel;
   show Hessian CovB;
   ods output Cov=CovB;
   ods exclude StoreInfo ClassLevels;
run;

