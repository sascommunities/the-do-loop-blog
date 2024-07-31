/* SAS program to accompany the article 
   "Fit, simulate, fit: How models can collapse after generations of recursive fitting"
   by Rick Wicklin, published 31JUL2024 on The DO Loop blog:
    https://blogs.sas.com/content/iml/2024/07/31/synthetic-data-model-collapse.html

   This program shows that models trained recursively from synthetic data 
   can experience parameter drift and model collapse, which is a decrease in the 
   variance of the model.

   The programs are intended to illustrate the ideas related to 
   model collapse in synthetic text generation:
   Gibney, 2024, "AI models fed AI-generated data quickly spew nonsense"
       https://www.nature.com/articles/d41586-024-02420-7
   which reports are the article by 
   Shumailov, I., Shumaylov, Z., Zhao, Y. et al. (2024)
       "AI models collapse when trained on recursively generated data,"
       Nature 631, 755–759. 
       https://doi.org/10.1038/s41586-024-07566-y
*/

/* clear macros and set define some useful macros from 
   https://blogs.sas.com/content/iml/2013/05/24/turn-off-ods-for-simulations.html
*/
%symdel N / nowarn;
%symdel SeedVal / nowarn;
%macro ODSOff(); /* Call prior to macro loop */
ods graphics off;
ods exclude all;
ods noresults;
options nonotes;
%mend;

%macro ODSOn(); /* Call after BY-group processing */
ods graphics on;
ods exclude none;
ods results;
options notes;
%mend;

/***************************/
/* Start with normally distributed data.
   At each generation, fit a model to the data, then simulate 
   more data from the model. Repeat many times.
   The following simulations take about 30-45 seconds to complete.
*/
%let N = 200;
%let SeedVal = 123;

%let mu0 = 125;
%let sigma0 = 30;
data Param0;
mu = &mu0;
sigma = &sigma0;
run;

%macro Synthetic(i, seed);
%global N;
data Model;
set Param%eval(&i-1);
call streaminit(&seed);
do i = 1 to &N;
   x = rand("Normal", mu, sigma);
   x = round(x, 1e-5);
   output;
end;
keep x;
run;

proc means data=Model noprint;
  var x;
  output out=Param&i.(keep=mu sigma) mean=mu std=sigma;
run;
data Param&i.;
   set Param&i.;
   ID = &i;
   label ID = "Generation";
run;
%mend;


%macro RunSynthetic(Niter, seed);
%if &seed=0 %then %do;
   %do i = 1 %to &NIter;
      %Synthetic(&i, 0);
   %end;
%end;
%else %do;
%do i = 1 %to &NIter;
   %Synthetic(&i, %eval(&seed+&i));
%end;
%end;
%mend;

%let Iters = 400;
%let SeedVal = 12345;
%ODSOff;
%RunSynthetic(&Iters, &SeedVal);

data AllParams;
set Param:;
Lower = mu - sigma;
Upper = mu + sigma;
run;
proc sort data=AllParams; by ID; run;

proc delete lib=work data=Param1-Param&Iters.;
run;
%ODSOn;

title "Evolution of Parameter Estimates";
title2 "Each model is trained on the preceding set of synthetic data";
title3 "Initial Seed = &SeedVal";
proc sgplot data=AllParams;
   scatter x=ID y=mu / yerrorlower=Lower yerrorupper=Upper;
   refline &mu0;
   yaxis grid label="Mean";
   xaxis label="Generation";
   WHERE mod(ID,2)=0; /* graph is very crowded; show only even generations */
run;

/******************************************************/
/* repeat with a different set of random number seeds */
%let SeedVal = 2468;
%ODSOff;
%RunSynthetic(&Iters, &SeedVal);

data AllParams;
set Param:;
Lower = mu - sigma;
Upper = mu + sigma;
run;
proc sort data=AllParams; by ID; run;

proc delete lib=work data=Param1-Param&Iters.;
run;
%ODSOn;

title "Evolution of Parameter Estimates";
title2 "Each model is trained on the preceding set of synthetic data";
title3 "Initial Seed = &SeedVal";
proc sgplot data=AllParams;
   scatter x=ID y=mu / yerrorlower=Lower yerrorupper=Upper;
   refline &mu0;
   yaxis grid label="Mean";
   xaxis label="Generation";
   WHERE mod(ID,2)=0; /* graph is very crowded; show only even generations */
run;

/*********************************/

/* Simulate lognormal data. See
   https://blogs.sas.com/content/iml/2017/05/10/simulate-lognormal-data-sas.html
*/
/* These data were simulated from a lognormal model with
   theta = 10;       * threshold parameter;
   zeta = 2;         * scale or log-location parameter;
   sigma = 0.5;      * shape or log-scale parameter;
*/
%let N = 200;          /* sample size */
data Have;
input x @@;
datalines;
18.97 20.85 25.34 17.41 18.25 17.01 13.92 16.03 16.71 16.11 
22.38 26.13 15.36 25.3 17.54 20.45 24.67 12.87 20.39 21.96 
22.45 13.9 20.86 15.98 14.79 18.59 17.82 21.55 18.69 22.48 
18.28 18.31 14.97 17.41 14.51 16.95 21.52 14.73 16.9 20.81 
18.63 18.2 15.31 20.06 12.96 33.35 16.34 15.78 15.08 18.4 
14.78 13.26 15.49 16.84 20.97 19.87 17.47 16.32 11.95 16.15 
25.22 20.15 17.59 14.72 18.35 21.05 15.12 19 21.43 12.23 
22.5 19.1 22.65 16.25 14.48 18.49 19.34 13.55 13.91 15.43 
18.12 13.99 21.25 15.42 16.39 17.97 15.17 15.97 18.2 15.29 
18.42 26 17.88 13.31 14.01 14.02 11.38 13.21 26.44 17.41 
17.59 21.3 22.78 18.21 20.91 15.08 15.82 18.26 14.32 16.51 
15.23 17.88 15.51 18.46 18.11 16.24 21.12 26.43 22.95 18.67 
19.63 16.32 16.27 18.5 15.35 13.22 14.1 15.48 20.39 16.31 
14.28 18.35 10.36 24.15 15.71 12.9 19.48 17.99 21.49 13.79 
15.69 19.69 16.82 19.16 15.88 13.12 14.84 18.02 15.86 18.98 
12.85 12.67 14.34 23.46 11.46 21.61 14.52 15.9 18.66 18.21 
21.68 26.86 18.85 24.15 24.34 19.97 17.43 19.67 16.32 16 
22.26 21.45 18.37 15.79 21.88 16.75 23.38 19.73 19.4 27.98 
24.82 17.99 16.52 18.1 15.81 17.43 16.98 16.87 20.42 26.21 
29.21 25.79 18.88 14.92 31.02 17.12 19.36 24.38 13.87 11.72 
;

title "Original Data and Lognormal Model";
proc univariate data=Have;
histogram x / lognormal(threshold=est) odstitle=title;
ods select Moments ParameterEstimates Histogram;
ods output ParameterEstimates = PELong;
run;

/* transpose from long to wide */
data PE0;
retain theta zeta sigma;
set PELong end=EOF;
if Symbol="Theta" then theta=Estimate;
if Symbol="Zeta"  then  zeta=Estimate;
if Symbol="Sigma" then sigma=Estimate;
if EOF then output;
keep theta zeta sigma;
run;

/* show synthetic data simulated from the model */
data Model;
call streaminit(54321);
set PE0;
do i = 1 to &N;
   Y = rand("Normal", zeta, sigma);
   X = theta + exp(Y);
   output;
end;
keep X;
run;

title "First Generation Lognormal Model";
proc univariate data=Model;
histogram x / lognormal(threshold=est) odstitle=title;
ods select Moments ParameterEstimates Histogram;
ods output ParameterEstimates = PELong;
run;

/**************************/

%macro SyntheticLN(i, seed);
%global N;
%put  &=N, &=i, &=seed;
ods select NONE;
data Model;
set PE%eval(&i-1);
call streaminit(&seed);
do i = 1 to &N;
   Y = rand("Normal", zeta, sigma);
   X = theta + exp(Y);
   output;
end;
keep X;
run;

proc univariate data=Model;
histogram x / lognormal(threshold=est);
ods output ParameterEstimates = PELong;
run;

/* transpose from long to wide */
data PE&i.;
retain theta zeta sigma;
set PELong end=EOF;
if Symbol="Theta" then theta=Estimate;
if Symbol="Zeta"  then  zeta=Estimate;
if Symbol="Sigma" then sigma=Estimate;
if EOF then output;
keep theta zeta sigma;
run;

data PE&i.;
   set PE&i.;
   ID = &i;
   label ID = "Generation";
run;
ods select ALL;
%mend;

/*
%SyntheticLN(1 ,124);
proc print data=PE1 noobs; run;
%RunSyntheticLN(1, 123);
*/

%macro RunSyntheticLN(Niter, seed);
%if &seed=0 %then %do;
   %do i = 1 %to &NIter;
      %SyntheticLN(&i, 0);
   %end;
%end;
%else %do;
%do i = 1 %to &NIter;
   %SyntheticLN(&i, %eval(&seed+&i));
%end;
%end;
%mend;

data Model;
set PE0;
call streaminit(124);
do i = 1 to &N;
   Y = rand("Normal", zeta, sigma);
   X = theta + exp(Y);
   output;
end;
keep X;
run;

%let Iters = 100;
%let SeedVal = 123;
%ODSOff;
%RunSyntheticLN(&Iters, &SeedVal);

data AllParams;
set PE:;
run;
proc sort data=AllParams; by ID; run;

proc delete lib=work data=PE1-PE&Iters.;
run;
%ODSOn;

title "Evolution of Parameter Estimates in a Lognormal Model";
title2 "Each model is trained on the preceding set of synthetic data";
title3 "Initial Seed = &SeedVal";
proc sgplot data=AllParams;
series x=id y=theta;
series x=id y=zeta;
series x=id y=sigma;
yaxis grid label="Mean";
xaxis label="Generation";
run;

/*
title "Simulated Data and Model after &Iters Generations";
proc univariate data=Model;
histogram x / lognormal(threshold=est) odstitle=title;
ods select Moments ParameterEstimates Histogram;
run;
*/
