/* SAS program to accompany the article 
   "What does 'flatten the curve' mean? To which curve does it apply?"
   by Rick Wicklin, published 06MAY2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/05/06/what-flatten-the-curve-means.html

   This program plots the new cases of confirmed coronavirus by day
   for the US until 03May2020.
   The data are from:
   www.cdc.gov/coronavirus/2019-ncov/cases-updates/cases-in-us.html
*/
/***************************************************/
/* read code by using 
   https://gitlab.sas.com/covid-19/data/jhu-data/-/blob/master/git-nyt-covid.sas
   which was written by Chris Hemedinger.

   Run this part just once. It creates a temp folder
   and clones the NYT GitHub repo into it          */
/***************************************************/
options dlcreatedir;
%let repoPath = %sysfunc(getoption(WORK))/covid-19-data;
libname repo "&repoPath.";
libname repo clear;
 
/* Fetch latest data from GitHub */
/* Note that GIT functions are in SAS 9.4m6 and SAS Viya 3.5 */
data _null_;
 rc = gitfn_clone(
   "https://github.com/nytimes/covid-19-data",
   "&repoPath."
 );
 put rc=;
run;
/***************************************************/

data us_states_nyt;
 length date 8 state $ 30 fips $ 6 cases 8 deaths 8;
 format date date9.;
 infile "&repoPath./us-states.csv" 
   dlm=',' missover dsd firstobs=2;
 input date : yymmdd10.
       state
       fips
       cases
       deaths;
run;

/* calculate delta of new cases and new deaths each day */
/* At state level */
proc sort data=us_states_nyt
 out=prepfordif;
by state date;
run;

data us_states_newperday
  (rename=(cases=totalcases deaths=totaldeaths));
 set prepfordif;
 by state date;
 retain newcases 0 newdeaths 0;
 newcases = dif(cases);
 newdeaths = dif(deaths);
 if first.state then do;
   newcases=cases; newdeaths=deaths;
 end;
run;

/*********************************************/
/* Data download and summary completed. Now
   compute the overall new cases by day and 
   cumulative curves.                        */
/*********************************************/
proc means data=us_states_newperday noprint;
class Date;
var newcases;
output out=USNewCases(where=(_type_=1)) sum=Cases;
run;
/*
proc sgplot data=USNewcases;
series x=Date y=Cases;
run;
*/
data USCases;
set USNewCases;
/* 7-day moving average */
Roll7 = mean(Cases, lag (Cases), lag2(Cases), lag3(Cases),
                   lag4(Cases), lag5(Cases), lag6(Cases));
Cumul + Cases;
run;

ods graphics /reset;
title "New Confirmed Cases versus Time (US)";
title2 "7-Day Rolling Average";
footnote J=L H=8pt "Source: github.com/nytimes/covid-19-data";
proc sgplot data=USCases noautolegend;
where Date > '16Feb2020'd;
vbar Date / response=Cases;
vline Date / response=Roll7 lineattrs=GraphData2(thickness=3)
             name="roll" legendlabel="7-Day Average";
xaxis type=Time display=(nolabel);
yaxis grid values=(0 to 40000 by 10000) label="Number of New Cases";
keylegend "roll" / location=inside position=NE;
run;

title "Cumulative Cases versus Time (US)";
proc sgplot data=USCases noautolegend;
where Date > '16Feb2020'd and cumul>0;
series x=Date y=Cumul;
xaxis type=Time display=(nolabel);
yaxis grid label="Cumulative Cases" /* type=log logbase=10 */;
run;


proc means data=USCases;
where Date>'12Apr2020'd;
var Cases;
run;
