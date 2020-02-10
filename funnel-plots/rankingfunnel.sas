/* SAS program to accompany the article 
   "Funnel plots: An alternative to ranking"
   by Rick Wicklin, published 15APR2011 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2011/04/15/funnel-plots-an-alternative-to-ranking.html
   
   SAS program to create a funnel plot of group means. 
   The "funnel" limits are computed according to 
   the usual CI for means:  
      mean(x) +/- z_alpha s / sqrt(n)
   where 
   mean(x) is the overall mean
   z_alpha is the 1-alpha quantile for the normal distribution
   s       is the overall sample standard deviation
   n       is the number of observations

   This computation (but not the data) appears in the following paper:
   D. Spiegelhalter (2005), _Statistics in Medicine_, 
   http://medicine.cf.ac.uk/media/filer_public/2010/09/24/spiegelhalter_stats_in_med_funnel_plots.pdf
*/
title "Temperatures of Car Roofs";
title2 "71 Degrees (F) in the Shade";

/*
Data from Clark Andersen, posted on
http://www.stat.columbia.edu/~cook/movabletype/archives/2011/03/a_possible_reso.html
*/
data CarTemps;
/** MISSOVER option ==> read a missing
    value when beyond the record **/
infile datalines missover;
input color $8. Temperature @;

/** read until we get a missing value **/
do while (^missing(Temperature));
  output;
  input Temperature @;
end;
datalines;
black    145.5 154.0 143.0 126.5 128.0 131.5 128.5 141.5
blue     134.5 123.0 113.5 126.0 139.0 125.5 143.5
burgundy 131.0 141.0 141.5 122.0
gray     124.0 141.5 128.0 127.5 132.0 130.0
green    149.5 130.5 114.0 129.5
red      140.5 128.5 124.5 119.5 132.0 126.0
silver   118.5 112.5 96.5 105.5 111.0 103.5
tan      127.0 103.5 121.0 113.0
white    103.5 97.5 93.0 99.0 91.5 104.5 99.5
;
run;

title; title2;
proc iml;
use CarTemps;
read all var {Color Temperature};
close CarTemps;

/**********************************************************************/
/** IMPORTANT: The VAR function was introduced in SAS/IML 9.22.      **/
/** For SAS/IML versions prior to 9.22, uncomment the following code **/
/**********************************************************************/
/**
start Var(x);
   mean = x[:,];
   countn = j(1, ncol(x));
   do i = 1 to ncol(x);
      countn[i] = sum(x[,i]^=.);
   end;
   var = (x-mean)[##,] / (countn-1);
return ( var );
finish;
**/

/** for each car color, compute the mean and
    standard error of the mean **/
u = unique(Color); /** unique colors **/
p = ncol(u); /** how many colors? **/
mean = j(p,1); sem = j(p,1); n = j(p,1);
do i = 1 to p;
   idx = loc(Color=u[i]);
   n[i] = ncol(idx);
   T = Temperature[idx];
   mean[i] = T[:];
   sem[i] = sqrt(var(T)/n[i]); 
end;

/** print means and stderrs in descending order **/
call sortndx(jdx, mean, 1, 1);
/** reorder all the values **/
color = u[jdx];   n = n[jdx];
mean = mean[jdx]; sem = sem[jdx];
print color n mean[format=5.1] sem[format=4.1];
/** save to data set **/
create Temp var {color n mean sem};
append;
close Temp;



/** compute overall mean and variance **/
y = Temperature[:];
s = sqrt(var(Temperature));
print y s[label="StdDev"];

/** confidence limits for a range of sample sizes **/
n = T( do(3, 8.5, 0.1) );
p = {0.001 0.025 0.975 0.999}; /** lower/upper limits **/
z = quantile("normal", p);
/** compute 56 x 4 matrix that contains confidence limits
    for n = 3 to 8.5 by 0.1 **/
limits = y + s/sqrt(n)*z;

/** write to data set **/
d = n || limits;
varNames = {"N" "L998" "L95" "U95" "U998"};
create Bounds from d[c=varNames];
append from d;
close Bounds;
quit;

data All;
set Temp Bounds;
run;

proc sgplot data=All;
scatter x=N y=Mean /datalabel=Color;
refline 123.4 / axis=y;
band x=N lower=L95 upper=U95 / nofill 
   name="CI95" LegendLabel="95% limits"
   lineattrs=(color=gray);
band x=N lower=L998 upper=U998 / nofill 
   name="CI998" LegendLabel="99.8% limits"
   lineattrs=(pattern=dash color=gray);
xaxis label="Number of Cars";
yaxis label="Average Temperature";
keylegend "CI998" "CI95" /
   location=inside position=TopRight
   across=1;
run;

/* compare with other approaches:
   Analysis of Means plot (ANOM)
   and mean with standard errors */
/*
proc anom data=CarTemps;
   xchart Temperature*Color;
   label Temperature = 'Mean Temperature (F)';
   label Color  = 'Car Color';
run;

proc sgplot data=CarTemps noautolegend;
  dot Color / response=temperature stat=mean limitstat=clm;
  refline 123.4 / axis=x lineattrs=(pattern=dash);
  xaxis label="Mean Temperature with 95% Confidence Intervals";
  yaxis discreteorder=data label="Color"; 
run;

proc glm data=CarTemps order=freq;
   class Color;
   model Temperature = Color;
   lsmeans Color / diff=anom adjust=nelson;
run;

*/

