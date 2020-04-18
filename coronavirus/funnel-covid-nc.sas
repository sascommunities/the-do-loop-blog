/* SAS program to accompany the article 
   "Use a funnel plot to visualize rates: The case fatality rate for COVID-19 
   in North Carolina counties"
   by Rick Wicklin, published 20APR2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/04/20/funnel-plot-covid-nc.html

   This program shows how to create a funnel plot for the case fatality rates
   for COVID-19 in North Carolina counties.
   Data from the New York Times, accessed 16Apr2020 from 
   https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv

   To learn more about the case fatality rate, see 
   https://newslit.org/updates/case-fatality-rate-vs-mortality-rate/
*/

/* Read 15Apr2020 data for county-level confirmed cases and deaths */
*filename USData 'c:/<PATH-TO-DATA>/us-COVID-15APR2020.csv';
filename USData 'c:/Users/frwick/OneDrive - SAS/Documents/My SAS Files/Blog/us-COVID-15APR2020.csv';

data CovidUSCounty;
infile USData delimiter = ',' MISSOVER DSD lrecl=32767;
informat Date anydtdte10.;
length County $16;
length State $14;
format Date DATE9. FIPS Z5.;
input Date County $ State $ FIPS Cases Deaths;
if cases > 0 then
   Proportion = Deaths / Cases;  /* observed proportion (case fatality rate) */
else delete;
run;

/* subset to NC and label counties with largest number of cases */
data CovidNC;
set CovidUSCounty;
if State="North Carolina";
if cases>=100 then 
   Label = substr(County, 1,11);
else Label = "        ";
run;

proc means data=CovidNC P10 P25 P50 P75 P90 P95 P99 Max ndec=0;
var Cases;
run;

/*******************************/
%let DSName = CovidNC;
%let Events = Deaths;
%let Trials = Cases;

/* compute overal rate in NC */
proc sql noprint;                              
 select  sum(&Events) / sum(&Trials) format=best6. into :AvgProp 
 from &DSName;
quit;
%put &=&AvgProp;

proc print data=CovidNC ;
where Cases >= 100;
run;

/* Isolate the computation of the the control limits for the funnel plot 
   so it can be reused in multiple programs.
*/
proc iml;
/* Funnel plot for proportion:
   https://blogs.sas.com/content/iml/2011/11/23/funnel-plots-for-proportions.html

   Compute binomial limits as a function of the nunber of trials.
   Adjust the binomial quantiles according to method in Appendix A.1.1
   of Spiegelhalter (2005) p. 1197.

   These limits do not depend on the data, although we assume that theta
   (often the average proportion from the data) is the probability parameter.

Parameters:
   theta (scalar): empirical probability of event = sum(Events) / sum(Trials)
   nTrials (column vector): sequence of trials at which to evaluate the binomial limits.
             For example, nTrials = T(1:100)
   prob (row vector): probabilities at which to compute quantiles. 
             For example, prob = {0.025 0.975} for 95% (~ 2 StdDev) limits
                          prob = {0.001 0.025 0.975 0.999} for 2 and 3 StdDev limits
*/
start BinomialFunnel(theta, nTrials, prob);
   n = round(colvec(nTrials));
   n = n <> 1;         /* must be integer > 0 */
   p = rowvec(prob);
   /* compute columns of matrix, one for each CL */
   limits = j(nrow(n), ncol(p), .);
   do i = 1 to ncol(p);
      r = quantile("binom", p[i], theta, n);
      numer = cdf("binom", r  , theta, n) - p[i];
      denom = cdf("binom", r  , theta, n) - 
              cdf("binom", r-1, theta, n);
      alpha = numer / denom;
      limits[,i] = (r-alpha)/n;
   end;
   /* we can't have negative limits */
   idx = loc(limits < 0);
   if ncol(idx)>0 then limits[idx]=0;

   /* return these control limits */
   return ( n || limits );
finish;
store module=BinomialFunnel;
QUIT;


/* Compute the funnels for the COVID-19 data */
proc iml;
use CovidNC;
read all var "&Events" into Events;
read all var "&Trials" into Trials;
close;

/* compute overall proportion */
theta = sum(Events) / sum(Trials); 

/* compute confidence limits for range of sample size */
load module=BinomialFunnel;
nTrials = T(do(4,200,10)) // T( do(200, round(max(Trials)+20, 50), 50) );
prob = {0.025 0.975};
results = BinomialFunnel(theta, nTrials, prob);
labels = {"N" "L2sd" "U2sd"};
create Limits from results[colname=labels];
   append from results;
close;
QUIT;

/* concatenate the data and the limits for the funnel plot */
data FunnelNC;
set CovidNC(keep=County Cases Deaths Proportion Label) 
    Limits; 
run;

/* create funnel plot with data tips. See
   https://blogs.sas.com/content/iml/2011/12/09/creating-tooltips-for-scatter-plots-with-proc-sgplot.html
*/
ods graphics / imagemap=ON TIPMAX=5000; 
title "NC Counties";
title "Case Fatility Rate: Deaths per Confirmed Cases";
title2 "Counties with 5 or More Cases";
footnote J=L "NYT data from 15Apr2020";

proc sgplot data=FunnelNC noautolegend;
where Cases=. OR Cases >= 5;
format Cases comma7.;
band x=N lower=L2sd upper=U2sd /
     fill fillattrs=GraphConfidence outline lineattrs=GraphPrediction 
     transparency=0.5 legendlabel="95% limits" name="band95"
     curvelabelloc=outside curvelabelupper="Upper Limit" curvelabellower="Lower Limit";
scatter x=Cases  y=Proportion / legendlabel="County Proportion" 
        datalabel=Label tip=(County Cases Deaths);
refline &AvgProp / axis=y label="NC Proportion = &AvgProp" 
                   SplitChar="/" name="ref" 
                   lineattrs=GraphData2 label;
xaxis grid ;
yaxis grid Label="Case Fatality Rate"; 
run;

ods graphics / imagemap=OFF;
