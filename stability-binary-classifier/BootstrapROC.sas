/* SAS program to accompany the article
   "Discrimination, accuracy, and stability in binary classifiers"
   by Rick Wicklin, published 08MAY2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/05/08/stability-binary-classifier.html

   This analyses illustrates many of the ideas in
   D. Ling, "Measuring Model Stability", Proceedings of the SAS Global Forum 2019 Conference,
   and presented at the conference.
   https://www.sas.com/content/dam/SAS/support/en/sas-global-forum-proceedings/2019/3568-2019.pdf
*/

/*  Data and example from "Comparing Receiver Operating Characteristic Curves"
    in the SAS/STAT documentation for PROC LOGISTIC.

    Used in "Create and compare ROC curves for any predictive model"
    https://blogs.sas.com/content/iml/2018/11/14/compare-roc-curves-sas.html
*/
data roc;
   input alb tp totscore popind @@;
   totscore = 10 - totscore;
   datalines;
3.0 5.8 10 0   3.2 6.3  5 1   3.9 6.8  3 1   2.8 4.8  6 0
3.2 5.8  3 1   0.9 4.0  5 0   2.5 5.7  8 0   1.6 5.6  5 1
3.8 5.7  5 1   3.7 6.7  6 1   3.2 5.4  4 1   3.8 6.6  6 1
4.1 6.6  5 1   3.6 5.7  5 1   4.3 7.0  4 1   3.6 6.7  4 0
2.3 4.4  6 1   4.2 7.6  4 0   4.0 6.6  6 0   3.5 5.8  6 1
3.8 6.8  7 1   3.0 4.7  8 0   4.5 7.4  5 1   3.7 7.4  5 1
3.1 6.6  6 1   4.1 8.2  6 1   4.3 7.0  5 1   4.3 6.5  4 1
3.2 5.1  5 1   2.6 4.7  6 1   3.3 6.8  6 0   1.7 4.0  7 0
3.7 6.1  5 1   3.3 6.3  7 1   4.2 7.7  6 1   3.5 6.2  5 1
2.9 5.7  9 0   2.1 4.8  7 1   2.8 6.2  8 0   4.0 7.0  7 1
3.3 5.7  6 1   3.7 6.9  5 1   3.6 6.6  5 1
;

proc logistic data=roc plots=roc(id=prob);
   LogisticModel: model popind(event='0') = alb tp totscore;
   output out=LogiOut predicted=LogiPred;       /* output predicted value, to be used later */
run;

data ExpertPred;
   input ExpertPred @@;
   datalines;
0.95 0.2  0.05 0.3  0.1  0.6  0.8  0.5 
0.1  0.25 0.1  0.2  0.05 0.1  0.05 0.1 
0.4  0.1  0.2  0.25 0.4  0.7  0.1  0.1 
0.3  0.2  0.1  0.05 0.1  0.4  0.4  0.7
0.2  0.4  0.1  0.1  0.9  0.7  0.8  0.25
0.3  0.1  0.1 
;
data Survival;
   merge LogiOut ExpertPred;
run;
 
/* overlay two or more ROC curves by using variables of predicted values */
proc logistic data=Survival alpha=0.1;
   model popind(event='0') = LogiPred ExpertPred / nofit;
   roc 'Expert'   pred=ExpertPred;
   roc 'Logistic' pred=LogiPred;
   ods select ROCOverlay ROCAssociation;
run;

data Prob;
  x = 0.73;
  xbar = 0.8162;
  s = 0.0801;
  /* x is how many std errors away from the point estimate? */
  z = abs(xbar - x) / s;
  /* probability of being that far or farther from the mean */
  p = 2*cdf("Normal", -z);
run;
proc print noobs; run;


/************************************************/
/* Bootstrap estimates */
/************************************************/


%let NumSamples = 5000;       /* number of bootstrap resamples */
/* 2. Generate many bootstrap samples */
proc surveyselect data=Survival NOPRINT seed=12345
     out=BootOut
     method=urs              /* resample with replacement */
     samprate=1              /* each bootstrap sample has N observations */
     /* OUTHITS                 option to suppress the frequency var */
     reps=&NumSamples;       /* generate NumSamples bootstrap resamples */
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

%ODSOff
proc logistic data=BootOut alpha=0.1;
   by Replicate;
   freq NumberHits;
   model popind(event='0') = LogiPred ExpertPred / nofit;
   roc 'Expert'   pred=ExpertPred;
   roc 'Logistic' pred=LogiPred;
   ods output ROCAssociation=ROCEst;
run;
%ODSOn

/* The data set ROCEst contains the bootstrap distribution 
   of the Area statistic. */
proc means data=ROCEst Mean Std P5 P95 P10 P90 P25 P75 ndec=4;
   class ROCModel;
   var Area;
run; 

/* output custom percentiles
https://blogs.sas.com/content/iml/2013/10/23/percentiles-in-a-tabular-format.html
*/
proc univariate data=ROCEst noprint;
   class ROCModel;
   var Area;
   histogram Area;
   output out=WidePctls pctlpre=P_ pctlpts=2.5 97.5 mean=Mean Std=Std; 
run; 

proc print data=WidePctls noobs label;
   format Mean Std P_2_5 P_97_5 6.4;
   label Mean="BootMean" Std="BootStdErr" P_2_5="95% Lower CL" P_97_5="95% Upper CL";
run;


ods graphics / height=481px width=640px;
title "Bootstrap Distributions of Area Under the Curve for Two Models";
proc sgpanel data=ROCEst;
  panelby ROCModel / rows=2 layout=rowlattice novarname;
  histogram Area;
  rowaxis grid; colaxis values=(0.4 to 1 by 0.1) grid;
run;

