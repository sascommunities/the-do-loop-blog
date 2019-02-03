/* SAS program to accompany the article 
   "Model selection with PROC GLMSELECT"
   by Rick Wicklin, published 04FEB2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/02/04/model-selection-glmselect.html

   This program shows how to use PROC GLMSELECT to build models 
   from a set of 8 monomial effects. The sequence of models are built on 
   training data by adding or removing effects that minimize the SBC criterion.
   The final model is chosen to the one that minimizes the ASE on the validation
   data. 
*/

/* Simulate data from a cubic polynomial regression model */
data Have;
call streaminit(54321);
do i = 1 to 250;
   x = rand("uniform", -3, 3);
   /* 2 - 1.105 x - 0.2 x^2 + 0.5 x^3 */
   y = 2 + 0.5*x*(x+1.3)*(x-1.7) + rand("Normal");
   output;
end;
run;

/* PROC GLMSELECT can actually automate this process by using 
   variable selection techniques. Use validation data to choose
   effects to enter and leave the model. Effects chosen from 
   Intercept, x, x**2, ..., x**7 */
ods trace off;
ods exclude ANOVA FitStatistics;
title "Select Model from 8 Effects";
proc glmselect data=Have seed=1 
               plots(startStep=1)=(ASEPlot Coefficients(unpack));
   effect poly = polynomial(x / degree=7);    /* generate monomial effects: x, x^2, ..., x^7 */
   partition fraction(validate=0.4);          /* use 40% of data for validation */
   model y = poly / selection=stepwise(select=SBC choose=validate)
                details=steps(Candidates ParameterEstimates); /* OPTIONAL */
run;


/* Run again, but this time save parameter estimates so that we can 
   visualize the sequence of chosen models */
title "Select Model from 8 Effects";
proc glmselect data=Have seed=1;
   effect poly = polynomial(x / degree=7);    /* generate monomial effects: x, x^2, ..., x^7 */
   partition fraction(validate=0.4);          /* use 40% of data for validation */
   model y = poly / selection=stepwise(select=SBC choose=validate stop=NONE)
                details=steps(ParameterEstimates); 
   ods output ParameterEstimates=PE;
   output out=GLSOut;
run;

proc print data=PE;
where Step=4; 
var Step Effect Estimate;
run;

/* visualize the model building process. Read the ParameterEstimates 
   for each step and overlay the series of models on the data. */

proc iml;
/* generate n evenly spaced points (a linearly spaced vector) in the interval [a,b] 
   https://blogs.sas.com/content/iml/2018/09/17/linearly-spaced-vectors-in-sas.html */
start linspace(a, b, numPts=100);
   n = floor(numPts);               /* if n is not an integer, truncate */
   if n < 1 then return( {} );      /* return empty matrix */
   else if n=1 then return( b );    /* return upper endpoint */
   return( do(a, b, (b-a)/(n-1)) ); /* return n equally spaced points */
finish;

/* Given a character matrix, this module returns a string with 
   the N elements separated by spaces or a user-specified delimiter. 
   https://blogs.sas.com/content/iml/2015/07/27/vector-to-string.html */
start MakeString(vv, delim=" ");
   v = rowvec(vv);
   if type(v)='N' then v = putn(v,"BEST6."); /* convert numbers to char */
   N = ncol(v);
   s = "";
   do i = 1 to N-1;
      s = s + strip(v[i]) + delim;
   end;
   return ( s + strip(v[N]) );
finish;

/* fit model where estimates are {"Intercept", "x", "x^2",...,"x^d"} 
   Input: t = x values on which to evaluate model
          Effect (col vec) = character vector identify polynomial effects 
          b (col vec) = parameter estimates for each effect
*/
start FitModel(t, Effect, b);
   f = 0#t;
   do j = 1 to nrow(b);
       ef = Effect[j];
      if ef = "Intercept" then f = f + b[j];
      else do;
         caretLoc = findc(ef, "^");  /* get degree */
         if caretLoc=0 then d=1;
         else d = num( substr(ef, caretLoc+1) );
         f = f + b[j]*t##d;
      end;
   end;
   return f;
finish;

use PE;
read all var {"Effect" "Estimate"};
read all var "Step" into stp;
close;
use GLSOut where (_ROLE_="TRAIN");
read all var {x y};
close;

MaxStep = max(stp);
numPts = 21;
L = max(nrow(x), numPts);
Step = j(L, 1, "       ");
Model = j(L, 1, BlankStr(256));
create Models var {"Step" "x" "y" "t" "f" "Model"};
t = linspace(min(x), max(x), numPts);
do k = 0 to maxStep;
   idx = loc(stp=k);
   f = FitModel(t, Effect[idx], Estimate[idx]);
   Step[,] = strip(char(k));
   Model[,] = MakeString(Effect[idx]);
   append;
end;
/* final model has Step=. */
   idx = loc(stp=.);
   f = FitModel(t, Effect[idx], Estimate[idx]);
   Step[,] = "Final";
   Model[,] = MakeString(Effect[idx]);
   append;
close Models;
quit;

/* Use BY group processessing to generate a plot for each model */
ods graphics / width=400px height=300px;
title "Build Model from Polynomial Effects";
proc sgplot data=Models uniform=scale noautolegend;
by Step notsorted Model;
scatter x=x y=y;
series x=t y=f;
xaxis min=-3 max=3 grid; yaxis min=-10 max=10 grid;
run;
title2;

/* use animation to overlay the models in one animated gif */
ods graphics / imagefmt=GIF width=4in height=3in;     /* each image is 4in x 3in GIF */
options papersize=('4 in', '3 in')                    /* set size for images */
        nodate nonumber                               /* do not show date, time, or frame number */
        animduration=1.0 animloop=yes noanimoverlay   /* animation details */
        printerpath=gif animation=start;              /* start recording images to GIF */
ods printer file='C:\AnimGif\GLMSelect\Anim.gif';  /* images saved into animated GIF */
 
ods html select none;                                 /* suppress screen output */
title "Build Model from Polynomial Effects";
proc sgplot data=Models uniform=scale noautolegend;
by Step notsorted Model;
scatter x=x y=y;
series x=t y=f;
xaxis grid; yaxis grid;
xaxis min=-3 max=3 grid; yaxis min=-10 max=10 grid;
run;

/* trick: generate the final frame again so it lasts twice as long 
   when displaying the animation */
proc sgplot data=Models uniform=scale noautolegend;
by Step notsorted Model;
where Step="Final";
scatter x=x y=y;
series x=t y=f;
xaxis min=-3 max=3 grid; yaxis min=-10 max=10 grid;
run;

ods html select all;                                  /* restore screen output */
 
options printerpath=gif animation=stop;               /* stop recording images */
ods printer close;                                    /* close the animated GIF file */

