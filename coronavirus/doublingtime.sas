/* SAS program to accompany the article 
   "Estimates of doubling time for exponential growth"
   by Rick Wicklin, published 01APR2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/04/01/estimate-doubling-time-exponential-growth.html

   This program shows how to estimate the slope of a log(Y) vs Time curve
   at the last time point. The method applies linear regression to the 
   most recent data to estimate the slope of the curve. The doubling time is
   DT = log(2) / m, where m is the estimate of teh slope.
*/

/* Data downloaded on 28MAR2020 from 
   https://data.humdata.org/dataset/novel-coronavirus-2019-ncov-cases
   and the  file 
   time_series_covid19_confirmed_global_narrow.csv
*/

/* Confirmed ronavirus cases: 3MAR2020 to 27MAR2020 */
data Italy;
length Region $15;
Region="Italy";
input Cumul @@; 
Day = _N_ - 1;
datalines;
2502 3089 3858 4636 5883 7375 9172 10149 12462 12462 17660 21157 
24747 27980 31506 35713 41035 47021 53578 59138 63927 69176 74386 80589 86498 
;
data US;
length Region $15;
Region="US";   
input Cumul @@; 
Day = _N_ - 1;
datalines;
118 149 217 262 402 518 583 959 1281 1663 2179 2727 3499 4632 
6421 7783 13677 19100 25489 33276 43847 53740 65778 83836 101657 
;
data Canada;
length Region $15;
Region="Canada";
input Cumul @@; 
Day = _N_ - 1;
datalines;
1 1 2 2 3 4 4 4 8 9 17 17 24 50 74 94 121 139 181 219 628 1013 1342 1632 2024 
;
data SouthKorea;
length Region $15;
Region="S. Korea";
input Cumul @@; 
Day = _N_ - 1;
datalines;
5186 5621 6088 6593 7041 7314 7478 7513 7755 7869 7979 8086 
8162 8236 8320 8413 8565 8652 8799 8961 8961 9037 9137 9241 9332 
;

/* concatenate the country data */
data Virus;
set Italy US Canada SouthKorea;
log10Cumul = log10(Cumul);
label Cumul = "Cumulative Cases"
      Day = "Days Since 03Mar2020";
run;

proc contents data=Virus short varnum; run;

/* Create a cumulative curve plot and use a logarithmic axis for the cumulatve counts.
   That is, plot log(Cumul) vs Time.
*/
title "Cumulative Counts (log scale)";
footnote H=0.8 J=L "Data from https://data.humdata.org/dataset/novel-coronavirus-2019-ncov-cases";
footnote2 H=0.8 J=L "time_series_covid19_confirmed_global_narrow.csv";
proc sgplot data=Virus;
  where Cumul > 0;
  series x=Day y=Cumul / group=Region curvelabel;
  xaxis grid;
  yaxis type=LOG logbase=10 grid
        values=(100 500 1000 5000 10000 50000 100000) valueshint;
run;
footnote;

/* Use most recent days (Day 20, 21, 22, 23, and 24) to estimate the
   slope of each curve on the last day */
%let MaxDay = 24;
proc reg data=Virus outest=Est noprint;
  where Day >= %eval(&MaxDay-4);  * previous 5 days;
  by Region notsorted;
  model log10Cumul = Day;
quit;

/* FOr exponential growth, the doubling time is DT = log(2) / m, 
   where m is the estimate of the slope of the log(Cumul)-vs-Time curve */
data DoublingTime;
set Est(rename=(Day=Slope));
label Slope = "Slope of Line" DoublingTime = "Doubling Time";
DoublingTime = log10(2) / Slope;
keep Region Intercept Slope DoublingTime;
run;

proc print data=DoublingTime label noobs;
  format Slope DoublingTime 6.3;
  var Region Slope DoublingTime;
run;

/* You can visualize the doubling time by adding an arrow to the end of each curve */
/* Base of arrow is the last data point */
data LastDataPoint;
set Virus;
by Region notsorted;
if last.Region;
run;

/* Arrow is a vector from the base at (t0, Y0) to (t0+DT, 2*Y0) */
data Vectors;
merge LastDataPoint DoublingTime;
keep Region t0 y0 t1 y1 Day;
t0 = &MaxDay; 
y0 = Cumul;
t1 = t0 + DoublingTime;
y1 = 2*y0;
Day = t0;
run;

/* add in a block variable to show prediction region */
data All;
set Virus Vectors(in=V);
if V then Block = "Predicted";
else      Block = "Observed ";
run;

/* show the graph with arrows and doubling times */
title "Confirmed Coronavirus Cases (Log Scale)";
title2 "Estimates of Doubling Time from Linear Trends";
proc sgplot data=All noautolegend;
  where Region ^= "S. Korea";
  block x=Day block=Block / blocklabel=Block FILLTYPE=ALTERNATE VALUEVALIGN=BOTTOM
              fillattrs=(color=White) altfillattrs=(color=CXE0E0E0);
  series x=Day y=Cumul / group=Region curvelabel curvelabelpos=START;
  scatter x=t0 y=y0 / group=region markerattrs=(symbol=CircleFilled);
  vector x=t1 y=y1 / xorigin=t0 yorigin=y0 group=Region lineattrs=(pattern=dash) name="vec";
  xaxis grid label = "Day Since 03Mar2020";
  yaxis type=LOG logbase=10 grid label="Cumulative Cases"
        values=(100 500 1000 5000 10000 50000 100000 20000) valueshint;
run;
