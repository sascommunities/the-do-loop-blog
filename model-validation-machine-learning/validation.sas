/* 1. Simulate data from a cubic polynomial regression model.
      nTrain = 50; nValidate = 200
*/
data Have;
length Type $10.;
call streaminit(54321);
do i = 1 to 250;
   if i <= 50 then Type = "Train";
   else            Type = "Validate";
   x = rand("uniform", -3, 3);
   /* 2 - 1.105 x - 0.2 x^2 + 0.5 x^3 */
   y = 2 + 0.5*x*(x+1.3)*(x-1.7) + rand("Normal");
   output;
end;
run;

/* Visualize the training and validation data. */
title "Training and Validation Data";
title2 "True Model Is Cubic Polynomial";
proc sgplot data=Have;
   scatter x=x y=y / group=Type grouporder=data;
   xaxis grid;
   yaxis grid;
run;


/* 2. You can use the EFFECT statement to define a POLYNOMIAL effect 
   of degreed=d. See
   https://blogs.sas.com/content/iml/2017/09/07/polynomial-effects-regression-sas.html
*/
%let Degree = 3;   
proc glmselect data=Have;
   effect poly = polynomial(x / degree=&Degree);              /* model is polynomial of specified degree */
   partition rolevar=Type(train="Train" validate="Validate"); /* specify training/validation observations */
   model y = poly / selection=NONE;                           /* fit model on training data */
   ods select FitStatistics ParameterEstimates;
   *output out=glmout P=pred R=resid;
run;

/* Of course, in reality, we don't know the true model! 
   Let's fit many polynomial models of different degrees and 
   choose the best one according
   (A) A classical information criterion such as AICC or SBC
   (B) The average square error of the predicted model 
       when evaluated on the validation data.
*/
%MACRO DoPolyFit(MaxDegree);
proc datasets noprint nowarn;   /* delete and data sets with prefix 'FitStats' */
   delete FitStats:;
quit;
options nonotes;
%DO Degree = 1 %TO &MaxDegree;
   /* use POLYNOMIAL effect to fit polynomial of degree=D */
   title "Fit Polynomial of Degree = &Degree.";
   proc glmselect data=Have;
      partition rolevar=Type(train="Train" validate="Validate");
      effect poly = polynomial(x / degree=&Degree);
      model y = poly / selection=NONE;
      ods output FitStatistics = Stats;
      ods select FitStatistics ParameterEstimates;
   run;
   /* add Degree variable to data set */
   data FitStats&Degree.;
      Degree = &Degree.;   set Stats(drop=cValue1);
   run;
%END;
options notes;
/* concatenate all data sets into one data set. Rename some variables. */
data PolyFit;
   set FitStats:;
   rename Label1 = Statistic nValue1 = Value;
run;
%MEND;

%DoPolyFit(7);

/* 3. Visualize the goodness of fit as a function of the degree of 
      the polynomial models. Plot ASE on training and validation 
      data sets versus the degree of the polynomial. Note that 
      the minimum ASE on the validation data occurs when d=3. */
title "Fit Polynomial Models to Data";
title2 "nTrain = 50; nValidation = 200";
proc sgplot data=PolyFit;
   where Statistic in ('ASE (Train)' 'ASE (Validate)');
   series x=Degree y=Value / markers group=Statistic;
   yaxis grid type=log;
   xaxis grid;
run;


/* Compare with the classical method: AIC and SBC also choose d=3 */
title "Fit Statistics for Polynomial Models";
title2 "Sample Size = 50";
proc sgpanel data=PolyFit;
   where Statistic in ('AIC' 'AICC' 'SBC');
   panelBy Statistic / columns=1 uniscale=column onepanel sort=data;
   series x=Degree y=Value / markers;
   colaxis grid;
   rowaxis grid;
run;


/* PROC GLMSELECT can actually automate this process by using 
   variable selection techniques. Use validation data to choose
   effects to enter and leave the model. Effects chosen from 
   Intercept, x, x**2, ..., x**7 */
proc glmselect data=Have seed=1 plots=(ASEPlot Coefficients);
   effect poly = polynomial(x / degree=7);
   model y = poly / selection= stepwise(choose=validate select=validate);
   partition rolevar=Type(train="Train" validate="Validate");
run;
