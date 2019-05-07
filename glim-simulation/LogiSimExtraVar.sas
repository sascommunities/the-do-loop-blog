/* SAS program to accompany the article
   "How to simulate data from a generalized linear model"
   by Rick Wicklin, published 14MAY2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/05/06/simulate-glim.html
*/

/* Example from _Simulating Data with SAS_, p. 226--229 
   See also 
   https://blogs.sas.com/content/iml/2014/06/25/simulate-logistic-data.html

   Logistic model with parameters (-2.7, -0.03, 0.07) */

/*********************************/

/* simulate 100 samples from two models:
   SimType=Correct: the correct way to simulate a logistic model
   SimType=Wrong: a wrong way to simulate a logistic model. This 
       simulation adds a random normal effect to the linear predictor.
*/  
%let numSim = 100;
%let Std = 0.8;
data Heart;
call streaminit(12345);
set Sashelp.Heart(keep=Cholesterol Systolic);
x1 = Cholesterol;
x2 = Systolic;
do SampleID = 1 to &numSim;
   N = rand("Normal", 0, &Std);
   U = rand("Uniform");

   /* The correct simulation: use linear predictor w/o error */
   SimType = "Correct";              
   eta = -2.7 - 0.03*x1 + 0.07*x2;   /* linear predictor */
   mu = logistic(eta);               /* transform by inverse logit */
   y = (U < mu);                     /* random binary response with probability mu */
   output;

   /* The wrong simulation: add error to the linear predictor */
   SimType = "Wrong";                
   eta = eta + N;                    /* Wrong: linear predictor plus error */
   mu = logistic(eta);               /* transform by inverse logit */
   y = (U < mu);                     /* random binary response with probability mu */
   output;
end;
run;

proc sort data=Heart;
by SimType SampleID;
run;

/* Macros to disable ODS. See
   https://blogs.sas.com/content/iml/2013/05/24/turn-off-ods-for-simulations.html 
*/
%macro ODSOff(); /* Call prior to BY-group processing */
ods graphics off; ods exclude all; ods noresults;
options nonotes;
%mend;
 
%macro ODSOn(); /* Call after BY-group processing */
ods graphics on; ods exclude none; ods results;
options notes;
%mend;

/* fit both sets of data to the logistic model Y = X1 X2 */
%ODSOff
proc logistic data=Heart;
   by SimType SampleID;
   model y(Event='1') = x1 x2;
   ods output ParameterEstimates = PE;
run;
%ODSOn

/* convert parameter estimates from long to wide */
data PEWide;
keep SimType SampleID Intercept x1 x2;
retain Intercept x1 x2;
set PE;
if mod(_N_, 3)=1      then Intercept = Estimate;
else if mod(_N_, 3)=2 then x1 = Estimate;
else                       x2 = Estimate;
if mod(_N_,3)=0;
run;

proc sgscatter data=PEWide;
matrix Intercept x1 x2 / group=SimType transparency=0.5;
run;


/* create clibration plots for correct model and misspecified model */
/*
title "Calibration Plot";
ods graphics / LOESSMAXOBS=6000;
proc logistic data=Heart plots(maxpoints=none only)=calibration;
   where SampleID = 1;
   by SimType;
   model y(Event='1') = x1 x2 / GOF;      * New in 15.1;
   ods select CalibrationPlot;
run;
*/
