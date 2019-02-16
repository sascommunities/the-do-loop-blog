/* SAS program to accompany the articles
   "Calibration plots in SAS"
   by Rick Wicklin, published 14MAY2018 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2018/05/14/calibration-plots-in-sas.html

AND 
   "Decile calibration plots in SAS"
   by Rick Wicklin, published 16MAY2018 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2018/05/16/decile-calibration-plots-sas.html

AND
   "An easier way to create a calibration plot in SAS"
   by Rick Wicklin, published 20FEB2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/02/20/easier-calibration-plot-sas.html

   This program shows how to 
   - Manually create a calibration plot (with a loess fit) for a 
     logistic model with a binary response.
   - Manually create a decile calibration plot for a logistic model 
     with a binary response.
   - Automatically generate a calibration plot by using the PLOTS=CALIBRATION 
     option on the PROC LOGISTIC statement.

   The data are simulated according to a paper by 
   Hosmer, Hosmer, Le Cessie, and Lemeshow (1997)

REFERENCES:
   Graphical methods of assessing calibration, based on 
   Hosmer, Hosmer, Le Cessie, and Lemeshow (1997)
   "A comparison of goodness-of-fit tests for the logistic regression model," 
   https://www2.stat.duke.edu/~zo2/dropbox/goflogistic.pdf
   and referenced in 
   Austin and Steyerberg (2013) 
   "Graphical assessment of internal and external calibration of logistic regression 
    models by using loess smoothers"
   https://doi.org/10.1002/sim.5941 

   Get quadratic coefficients from 
   https://www2.stat.duke.edu/~zo2/dropbox/goflogistic.pdf
   as implemented in My SAS Files/blog/CalibrationParams.sas
*/

/* Step 1: Simulate data for a quadratic logistic regression model */
/* simulate sample of size N=500 that follows a quadratic model from Hosmer et al. (1997) */
%let N = 500;                              /* sample size */ 
data LogiSim(keep= Y x);
call streaminit(1234);              /* set the random number seed */
do i = 1 to &N;                      /* for each observation in the sample  */
   /* Hosmer et al. chose x ~ U(-3,3). Use rand("Uniform",-3,3) for modern versions of SAS */
   x = -3 + 6*rand("Uniform");
   eta = -2.337 + 0.8569*x + 0.3011*x**2;  /* quadratic model used by Hosmer et al. */
   mu = logistic(eta);
   Y = rand("Bernoulli", mu);
   output;
end;
run;

/******************************************************/
/* Third blog post:
   "An easier way to create a calibration plot in SAS"
   https://blogs.sas.com/content/iml/2019/02/20/easier-calibration-plot-sas.html

   This section requires SAS/STAT 15.1 (SAS 9.4M6).
*/
/******************************************************/

/*----------------------------------------------------*/
/* NEW in SAS/STAT 15.1 (SAS 9.4M6):                  */
/* Create calibration plot by using the               */
/* PLOTS=CALIBRATION option in PROC LOGISTIC          */
/*----------------------------------------------------*/


/* NEW in SAS/STAT 15.1 (SAS 9.4M6): PLOTS=CALIBRATION option in PROC LOGISTIC */
ods graphics / width=600px height=600px;
ods select CalibrationPlot;
title "Calibration Plot for a Quadratic Model";
title2 "Created by PROC LOGISTIC";
proc logistic data=LogiSim plots=calibration(CLM ShowObs);
   model y(Event='1') = x x*x / GOF;      /* New in 15.1: More goodness-of-fit statistics */
   output out=LogiOut predicted=PredProb; /* optional: output predicted probability */
run;

/*----------------------------------------------------*/
/* Traditional method: Use output from PROC LOGIST    */
/* and PROC SGPLOT to create calibration plot with    */
/* a loess smoother.                                  */
/*----------------------------------------------------*/

proc sort data=LogiOut;  by PredProb;  run;  /* Then sort */

/* Then create the calibration plot */
ods graphics / reset;
title "Calibration Plot for a Correct Model";
proc sgplot data=LogiOut noautolegend aspect=1;
   loess x=PredProb y=y / interpolation=cubic clm nomarkers; /* optional: SMOOTH=0.75 */
   lineparm x=0 y=0 slope=1 / lineattrs=(color=grey pattern=dash);
   yaxis grid; xaxis grid;
run;

/* calibration plots for polytomous response */

/* Derr (2013) "Ordinal Response Modeling with the LOGISTIC Procedure"
   https://support.sas.com/resources/papers/proceedings13/446-2013.pdf
   The following data, from McCullagh and Nelder (1989, p. 179), contain a 
   measure of the severity of pneumoconiosis (black lung disease) in coal miners
   and the number of years of exposure.
*/
data Coal; 
input Severity $ @@; 
do i=1 to 8; 
   input Exposure freq @@; 
   log10Exposure=log10(Exposure); 
   output; 
end; 
datalines; 
Normal   5.8 98 15 51 21.5 34 27.5 35 33.5 32 39.5 23 46 12 51.5 4 
Moderate 5.8  0 15  2 21.5  6 27.5  5 33.5 10 39.5  7 46  6 51.5 2 
Severe   5.8  0 15  1 21.5  3 27.5  8 33.5  9 39.5  8 46 10 51.5 5 
;

title 'Severity of Black Lung vs Log10(Years Exposure)';
proc logistic data=Coal rorder=data plots=Calibration(CLM);
   freq freq; 
   model Severity(descending)=log10Exposure; 
   effectplot / noobs individual;
run; 


/******************************************************/
/* First blog post:
   "Calibration plots in SAS"
   https://blogs.sas.com/content/iml/2018/05/14/calibration-plots-in-sas.html
*/
/******************************************************/

/*----------------------------------------------------*/
/* A calibration curve for a misspecified model       */
/*----------------------------------------------------*/

/* Use PROC LOGISTIC and output the predicted probabilities.
   Intentionally MISSPECIFY the model as linear. */
proc logistic data=LogiSim plots=effect;
   model Y(event='1') = x;
   output out=LogiOut predicted=PredProb; /* save predicted probabilities in data set */
   ods exclude InfluencePlots;
run;

proc sort data=LogiOut;  by PredProb;  run;

/* let the data choose the smoothing parameter */
title "Calibration Plot for Misspecified Model";
title2 "True Model Is Quadratic; Fit Is Linear";
proc sgplot data=LogiOut noautolegend aspect=1;
   loess x=PredProb y=y / interpolation=cubic clm;   /* smoothing value by AICC (0.657) */
   lineparm x=0 y=0 slope=1 / lineattrs=(color=grey pattern=dash);
run;

/* 0.75 = smoothing value used by Austin and Steyerberg */
proc sgplot data=LogiOut noautolegend;
   loess x=PredProb y=y / smooth=0.75 interpolation=cubic clm; /* or use NOMARKERS */
   lineparm x=0 y=0 slope=1 / lineattrs=(color=grey pattern=dash);
run;

/*----------------------------------------------------*/
/* A calibration plot for a correctly specified model */
/*----------------------------------------------------*/

proc logistic data=LogiSim noprint;
   model Y(event='1') = x x*x;  /* fit a quadratic model */
   output out=LogiOut2 predicted=PredProb;
run;

proc sort data=LogiOut2;  by PredProb;  run;

/* let the data choose the smoothing parameter */
title "Calibration Plot for a Correct Model";
proc sgplot data=LogiOut2 noautolegend aspect=1;
   loess x=PredProb y=y / interpolation=cubic clm nomarkers;
   lineparm x=0 y=0 slope=1 / lineattrs=(color=grey pattern=dash);
   yaxis grid; xaxis grid;
run;

proc loess data=LogiOut2 plots=fit;
model y = PredProb / interp=cubic clm;
run;


/******************************************************/
/* Second blog post:
   "Decile calibration plots in SAS"
   https://blogs.sas.com/content/iml/2018/05/16/decile-calibration-plots-sas.html
*/
/******************************************************/

/*----------------------------------------------------*/
/* A decile calibration curve for a misspecified model*/
/*----------------------------------------------------*/

/* Use PROC LOGISTIC and output the predicted probabilities.
   Intentionally MISSPECIFY the model as linear. */
proc logistic data=LogiSim noprint;
   model Y(event='1') = x;
   output out=LogiOut predicted=PredProb; /* save predicted probabilities in data set */
run;

/* Now use decile method from Hosmer et al. First identify deciles of the predicted prob. */
proc rank data=LogiOut out=LogiDecile groups=10;
   var PredProb;
   ranks Decile;
run;

/* Then compute the mean predicted prob and the empirical proportions (and CI) for each decile */
proc means data=LogiDecile noprint;
   class Decile;
   types Decile;
   var y PredProb;
   output out=LogiDecileOut mean=yMean PredProbMean
          lclm=yLower uclm=yUpper;
run;

proc  print data=LogiDecileOut; run;

title "Calibration Plot for Misspecified Model";
title2 "True Model Is Quadratic; Fit Is Linear";
proc sgplot data=LogiDecileOut noautolegend aspect=1;
   lineparm x=0 y=0 slope=1 / lineattrs=(color=grey pattern=dash);
   *loess x=PredProbMean y=yMean;  /* if you want a smoother based on deciles */
   series x=PredProbMean y=yMean; /* if you to connect the deciles */
   scatter x=PredProbMean y=yMean / yerrorlower=yLower yerrorupper=yUpper;
   yaxis label="Observed Probability of Outcome";
   xaxis label="Predicted Probability of Outcome";
run;

/*----------------------------------------------------*/
/* A decile calibration curve for a correctly specified model */
/*----------------------------------------------------*/

proc logistic data=LogiSim noprint;
   model Y(event='1') = x x*x;  /* fit a quadratic model */
   output out=LogiOut2 predicted=PredProb;
run;

/* Now use decile method from Hosmer et al. */
proc rank data=LogiOut2 out=LogiDecile2 groups=10;
   var PredProb;
   ranks Decile;
run;

proc means data=LogiDecile2 noprint;
   class Decile;
   types Decile;
   var y PredProb;
   output out=LogiDecileOut2 mean=yMean PredProbMean
          lclm=yLower uclm=yUpper;
run;

/* let the data choose the smoothing parameter */
title "Calibration Plot for a Correct Model";
proc sgplot data=LogiDecileOut2 noautolegend aspect=1;
   lineparm x=0 y=0 slope=1 / lineattrs=(color=grey pattern=dash);
   series x=PredProbMean y=yMean; /* if you want to connect the deciles */
   scatter x=PredProbMean y=yMean / yerrorlower=yLower yerrorupper=yUpper;
   yaxis label="Observed Probability of Outcome";
   xaxis label="Predicted Probability of Outcome";
run;
