/* SAS program to accompany the articles
   "Simulate data from a logistic regression model: How the intercept 
    parameter affects the probability of the event"
   by Rick Wicklin, published 16JAN2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/01/16/simulate-logistic-intercept.html

   and
 
   "Visualize how parameters in a binary logistic regression model 
    affect the probability of the event"
   by Rick Wicklin, published 18JAN2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/01/18/visualize-logistic-probability.html

   This program shows how to simulate data from the binary regression model
   where the probability of the event (Y=1) is given by 
   Pr(Y=1) = mu = logistic(eta;) = 1 / (1 + exp(-eta))
   where the linear predictor is 
   eta = beta0 + beta1*X
   and the explanatory variable, X, has a known distribution.
*/

ods graphics / reset;

/* graph the logistic function */
data Logi;
do x = -5 to 5 by 0.1;
   y = logistic(x);
   output;
end;
run;
title "Logistic Function";
proc sgplot data=Logi;
   series x=x y=y;
   xaxis grid;
   yaxis grid label="logistic(x)";
run;

/* A generalized linear model consists of 
 - a linear predictor 
 - a monotonic mapping between the mean of the data and the linear predictor
 - a response distribution 
 The distribution of the explanatory variables affects the results.
*/
/* Simulate X ~ N(0,1) */
data Explanatory;
call streaminit(12345);
do i = 1 to 1000;
   x = rand("Normal",  0, 1);   /* ~ N(0,1) */
   output;
end;
drop i;
run;

/*******************************************************/
/* Run simulations that show that Pr(Y=1) depends on the Intercept */
/*******************************************************/

/* Simulate data from one-variable logistic model.
   X ~ N(0,1) and Y ~ Bern(eta)
   where eta = logistic(Intercept + Slope*X).
   Display frequency table for Y and a graph of Y and Prob(Y=1) vs eta */
%macro SimLogistic(Intercept, Slope);
   title "Logistic Model for Simulated Data";
   title2 "Intercept=&Intercept; Slope=&Slope";
   data Sim1;
      call streaminit(54321);
      set Explanatory;
      eta = &Intercept + &Slope*x;   /* eta = linear predictor */
      mu = logistic(eta);            /* mu = Prob(Y=1) */
      Y = rand("Bernoulli", mu);     /* simulate binary response */
   run;

   proc freq data=Sim1;
      tables Y / nocum;
   run;
   proc sort data=Sim1; by eta; run;

   proc sgplot data=Sim1;
      label Y="Observed" mu="P(Y=1)" eta="Linear Predictor";
      scatter x=eta y=Y / transparency=0.8 markerattrs=(symbol=CircleFilled);
      series x=eta y=mu;
      xaxis grid values = (-16 to 16) valueshint;
      yaxis grid label="Probability";
      keylegend / location=inside position=E opaque;
   run;
%mend;

%SimLogistic(5, 2);
%SimLogistic(1, 2);
%SimLogistic(-1, 2);

/*******************************************************/
/* VISUALIZE HOW P(Y=1) DEPENDS ON SLOPE AND INTERCEPT */
/*******************************************************/
/* Simulate X ~ N(0,1) */
data Explanatory;
call streaminit(12345);
do i = 1 to 1000;
   x = rand("Normal",  0, 1);   /* ~ N(0,1) */
   output;
end;
drop i;
run;

/* as Intercept changes, how does P(Y=1) change when Slope=2? */
%let Slope = 2;
data SimLogistic;
call streaminit(54321);
set Explanatory;
do Intercept = -3 to 3 by 0.2;
   do nSim = 1 to 100;               /* Optional: param est better for large samples */
      eta = Intercept + &Slope*x;    /* eta = linear predictor */
      mu = logistic(eta);            /* mu = Prob(Y=1) */
      Y = rand("Bernoulli", mu);     /* simulate binary response */
      output;
   end;
end;
run;

proc sort data=SimLogistic; by Intercept; run;
ods select none;
proc freq data=SimLogistic;
   by Intercept;
   tables Y / nocum;
   ods output OneWayFreqs=FreqOut;
run;
ods select all;

title "Percent of Y=1 in Logistic Model";
title2 "Slope=&Slope";
footnote J=L "X ~ N(0,1)";

proc sgplot data=FreqOut(where=(Y=1));
   series x=Intercept y=percent;
   xaxis grid;
   yaxis grid label="Percent of Y=1";
run;


/* Run two-parameter simulation study where 
   slope and intercept are varied systematically on a grid */
data SimLogistic2;
call streaminit(54321);
set Explanatory;
do Slope = 0 to 4 by 0.2;
do Intercept = -3 to 3 by 0.2;
   do nSim = 1 to 50;                  /* optional: the parameter estimate are better for larger samples */
   eta = Intercept + Slope*x;          /* eta = linear predictor */
   mu = logistic(eta);                 /* transform by inverse logit */
   Y = rand("Bernoulli", mu);          /* simulate binary response */
   output;
   end;
end;
end;
run;

/* Monte Carlo estimate of Pr(Y=1) for each (Int,Slope) pair */
proc sort data=SimLogistic2; by Intercept Slope; run;
ods select none;
proc freq data=SimLogistic2;
   by Intercept Slope;
   tables Y / nocum;
   ods output OneWayFreqs=FreqOut2;
run;
ods select all;


/* Create template for a contour plot
   https://blogs.sas.com/content/iml/2012/07/02/create-a-contour-plot-in-sas.html
*/
proc template;
define statgraph ContourPlotParm;
dynamic _X _Y _Z _TITLE _FOOTNOTE;
begingraph;
   entrytitle _TITLE;
   entryfootnote halign=left _FOOTNOTE;
   layout overlay;
      contourplotparm x=_X y=_Y z=_Z /
        contourtype=fill nhint=12  colormodel=twocolorramp name="Contour";
      continuouslegend "Contour" / title=_Z;
   endlayout;
endgraph;
end;
run;

/* render the Monte Carlo estimates as a contour plot */
title; footnote;
proc sgrender data=FreqOut2 template=ContourPlotParm;
where Y=1;
dynamic _TITLE="Percent of Y=1 in a One-Variable Logistic Model"
        _FOOTNOTE="X ~ N(0, 1)"
        _X="Intercept" _Y="Slope" _Z="Percent";
run;

/* print (Int, Slope) pairs for which Pr(Y=1) is close to target value */
%let Target = 70;
proc print data=FreqOut2;
   where Y=1 and %sysevalf(&Target-1) <= Percent <= %sysevalf(&Target+1);
   var Intercept Slope Y Percent;
run;

/* The previous simulation used X ~ N(0,1), which is a symmetric distribution.
   Simulate exponential variates, standardize, and rerun the study */
/* Simulate X ~ Exp(1.5) */
data Expo;
call streaminit(12345);
do i = 1 to 1000;
   x = rand("Expon", 1.5);   /* ~ Exp(1.5) */
   output;
end;
drop i;
run;
/* standardize data some mean=0 and SD=1 */
proc stdize data=Expo method=Std out=Explanatory;
   var x;
run;

proc sgplot data=Explanatory;
histogram x;
run;

/* repeat the simulation study for the new explanatory variable */
data SimLogistic3;
call streaminit(54321);
set Explanatory;
do Slope = 0 to 4 by 0.2;
do Intercept = -3 to 3 by 0.2;
   do nSim = 1 to 50;                  /* optional: the parameter estimate are better for larger samples */
   eta = Intercept + Slope*x;          /* eta = linear predictor */
   mu = logistic(eta);                 /* transform by inverse logit */
   Y = rand("Bernoulli", mu);          /* simulate binary response */
   output;
   end;
end;
end;
run;

/* Monte Carlo estimate of Pr(Y=1) for each (Int,Slope) pair */
proc sort data=SimLogistic3; by Intercept Slope; run;
ods select none;
proc freq data=SimLogistic3;
   by Intercept Slope;
   tables Y / nocum;
   ods output OneWayFreqs=FreqOut3;
run;
ods select all;

/* render the Monte Carlo estimates as a contour plot */
proc sgrender data=Freqout3 template=ContourPlotParm;
where Y=1;
dynamic _TITLE="Percent of Y=1 in a One-Variable Logistic Model"
        _FOOTNOTE="X ~ Exponential"
        _X="Intercept" _Y="Slope" _Z="Percent";
run;

/* print (Int, Slope) pairs for which Pr(Y=1) is close to target value */
%let Target = 70;
proc print data=FreqOut3;
   where Y=1 and %sysevalf(&Target-1) <= Percent <= %sysevalf(&Target+1);
   var Intercept Slope Y Percent;
run;

