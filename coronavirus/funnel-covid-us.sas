/* SAS program to accompany the article 
   "Visualize the case fatality rate for COVID-19 in US counties"
   by Rick Wicklin, published 20APR2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/04/20/funnel-plot-covid-us.html

   This program shows how to create a funnel plot for the case fatality rates
   for COVID-19 in US counties.
   Data from the New York Times, accessed 16Apr2020 from 
   https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv

   To learn more about the case fatality rate, see 
   https://newslit.org/updates/case-fatality-rate-vs-mortality-rate/
*/

/* Read 15Apr2020 data for county-level confirmed cases and deaths. The
   data are available at
   https://raw.githubusercontent.com/sascommunities/the-do-loop-blog/master/coronavirus/us-COVID-15APR2020.csv
   Fetch the file from the web site 
*/
filename USData temp;
proc http
 url="https://raw.githubusercontent.com/sascommunities/the-do-loop-blog/master/coronavirus/us-COVID-15APR2020.csv"
 method="GET"
 out=USData;
run;

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

/* label counties with the most cases */
data CovidUS;
set CovidUSCounty;
if cases>10000 then 
   Label = substr(County, 1,9);
else Label = "         ";
run;

/* look at percentiles of the cases. Highly skewed! */
proc means data=CovidUS P10 P25 P50 P75 P90 P95 P99 Max ndec=0;
   var Cases;
run;

/* print the counties with the most cases */
proc sort data=CovidUsCounty(where=(cases>10000)) out=Top;
   by descending Cases;
run;

proc print data=Top noobs;
   format Cases Deaths COMMA7. Proportion 5.3;
   var County State Cases Deaths Proportion;
run;

/*******************************************/
/* Create a funnel plot. Based on 
   https://blogs.sas.com/content/iml/2011/11/23/funnel-plots-for-proportions.html
   https://blogs.sas.com/content/iml/2018/11/26/funnel-plot-immunization-rates.html
*/
%let DSName = CovidUS;
%let Events = Deaths;
%let Trials = Cases;

/* overall proportion */
proc sql noprint;                              
 select  sum(&Events) / sum(&Trials) format=best6. into :AvgProp 
 from &DSName;
quit;
%put &=&AvgProp;

/* Compute the funnels for the COVID-19 data */
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

proc iml;
use &DSName;
read all var "&Events" into Events;
read all var "&Trials" into Trials;
close;

/* compute overall proportion */
theta = sum(Events) / sum(Trials); 

/* compute confidence limits for range of sample size */
load module=BinomialFunnel;
minN = min(Trials); maxN = max(Trials);
dx = (maxN - minN) / 500;
nTrials = T(do(20,200,10)) // T( do(200, maxN, dx) );
prob = {0.025 0.975};
results = BinomialFunnel(theta, nTrials, prob);
labels = {"N" "L2sd" "U2sd"};
create Limits from results[colname=labels];
   append from results;
close;
QUIT;

/* append control limits */
data FunnelUS; 
set &DSName(keep=State County Cases Deaths Proportion Label) 
    Limits; 
run;

/* create the funnel plot */
title "Case Fatility Rate: Deaths per Confirmed Cases";
title2 "US Counties with More Than 20 Cases";
footnote J=L "NYT data from 15Apr2020";

ods graphics / imagemap=ON TIPMAX=5000; /* enable data tips */
proc sgplot data=FunnelUS noautolegend;
where Cases=. OR Cases > 20;
format Cases comma7.;
band x=N lower=L2sd upper=U2sd / /* nofill lineattrs=(color=gray thickness=2)  */
     fill fillattrs=GraphConfidence outline lineattrs=GraphPrediction 
     transparency=0.5 legendlabel="95% limits" name="band95"
     curvelabelloc=outside
     curvelabelupper="Upper Limit" curvelabellower="Lower Limit";
scatter x=Cases  y=Proportion / legendlabel="County Proportion" 
        datalabel=Label datalabelpos=Top tip=(County State Cases Deaths Proportion);
refline &AvgProp / axis=y label="US Proportion = &AvgProp" 
                   SplitChar="/" name="ref" 
                   lineattrs=GraphData2 label;
xaxis grid offsetmin=0.01 offsetmax=0.06
      /* for a LOG scale axis, add the following:
         type=log logbase=10 label="Cases (log scale)" */
      ; 
yaxis grid Label="Case Fatality Rate"; 
run;

/* Zoom in on a range of Case values */
data FunnelZoom;
set FunnelUS;
if (20 < Cases < 5000) or (10< N < 5000);
if County = "Macomb" and State = "Michigan" then Label = "Macomb, MI";
if County = "King" and State = "Washington" then Label = "King, WA";
if County = "Hartford" and State = "Connecticut" then Label = "Hartford, CT";
if County = "Harris" and State = "Texas" then Label = "Harris, TX";
run;

title "Case Fatility Rate: Deaths per Confirmed Cases";
title2 "US Counties with More Than 20 and Less Than 5000 Cases";
footnote J=L "NYT data from 15Apr2020";

proc sgplot data=FunnelZoom noautolegend;
format Cases comma7.;
band x=N lower=L2sd upper=U2sd / /* nofill lineattrs=(color=gray thickness=2)  */
     fill fillattrs=GraphConfidence outline lineattrs=GraphPrediction 
     transparency=0.5 legendlabel="95% limits" name="band95"
     curvelabelloc=outside
     curvelabelupper="Upper Limit" curvelabellower="Lower Limit";
scatter x=Cases  y=Proportion / legendlabel="County Proportion" 
        datalabel=Label datalabelpos=Right tip=(County State Cases Deaths Proportion);
refline &AvgProp / axis=y label="US Proportion = &AvgProp" 
                   SplitChar="/" name="ref" 
                   lineattrs=GraphData2 label;
xaxis grid offsetmax=0.08;
yaxis grid Label="Case Fatality Rate"; 
run;

ods graphics / imagemap=OFF;
