/* SAS program to accompany the article 
   "Create a deviation plot to visualize values relative to a baseline"
   by Rick Wicklin, published 02MAR2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/03/02/deviation-plot-baseline.html

   This program shows how to create a deviation-from-baseline plot 
   (Or "deviation plot," for short). The deviation plot is similar to a waterfall chart
   https://blogs.sas.com/content/iml/2015/04/20/waterfall-plot.html

   The temperature data for Raleigh is available from 
   Data from https://waterdata.usgs.gov/nwis/uv?site_no=02079490
   The blood glucose data is artificial.

   This article was inspired by 
   "3 Steps to Building an Air Temperature Circle Graph"
   by Mike Drutar 
   https://communities.sas.com/t5/SAS-Communities-Library/3-Steps-to-Building-an-Air-Temperature-Circle-Graph/ta-p/620899
*/

/* 0. Read the data */
data Series;
informat Date date.;
format Date Date.;
input Date y @@;
label y = "Blood Glucose (mg/dL)";
datalines;
01SEP19 100 02SEP19  96 03SEP19  86 04SEP19  93 05SEP19 105 06SEP19 106 07SEP19 123 
08SEP19 121 09SEP19 115 10SEP19 108 11SEP19  94 12SEP19  96 13SEP19  95 14SEP19 120
15SEP19 112 16SEP19 104 17SEP19  97 18SEP19 101 19SEP19 108 20SEP19 108 21SEP19 117 
22SEP19 103 23SEP19 109 24SEP19  97 25SEP19  93 26SEP19 100 27SEP19  98 28SEP19 122 
29SEP19 116 30SEP19  99 01OCT19 102 02OCT19  99 03OCT19  95 04OCT19  99 05OCT19 116 
06OCT19 109 07OCT19 106 08OCT19  94 09OCT19 104 10OCT19 112 11OCT19 119 12OCT19 111 
13OCT19 104 14OCT19 101 15OCT19  99 16OCT19  92 17OCT19 101 18OCT19 115 19OCT19 109 
20OCT19  98 21OCT19  91 22OCT19  92 23OCT19 100 24OCT19 109 25OCT19 102 26OCT19 117 
27OCT19 106 28OCT19  98 29OCT19  98 30OCT19  95 31OCT19  97 01NOV19 129 02NOV19 120 
03NOV19 117 04NOV19   . 05NOV19 101 06NOV19 105 07NOV19 105 08NOV19 106 09NOV19 118 
10NOV19 109 11NOV19 102 12NOV19  98 13NOV19  97 14NOV19 .   15NOV19  92 16NOV19 114 
17NOV19 107 18NOV19  98 19NOV19  91 20NOV19  97 21NOV19 109 22NOV19  98 23NOV19 95 
24NOV19  95 25NOV19  94 26NOV19   . 27NOV19  98 28NOV19 115 29NOV19 123 30NOV19 114 
01DEC19 104 02DEC19  96 03DEC19  97 04DEC19 100 05DEC19  94 06DEC19  93 07DEC19 105 
08DEC19   . 09DEC19  88 10DEC19  84 11DEC19 101 12DEC19 122 13DEC19 114 14DEC19 108 
15DEC19 103 16DEC19  88 17DEC19  74 18DEC19  92 19DEC19 110 20DEC19 118 21DEC19 106 
22DEC19 100 23DEC19 106 24DEC19 107 25DEC19 116 26DEC19 113 27DEC19 113 28DEC19 117 
29DEC19 101 30DEC19  96 31DEC19 101  
;

/* 1. Compute the deviation and encode the data as 'High' or 'Low' by using the reference value */
%let RefValue = 100;

data Center;
set Series;
if (y > &RefValue) then Group="High";
else Group="Low";
Low  = min(y, &RefValue);    /* lower end of highlow bar */
High = max(y, &RefValue);    /* upper end of highlow bar */
run;

/* 2. Create a discrete attribute map that maps values to colors */
data DAttrs;
length c FillColor LineColor $16.;
ID = "HighLow";
Value = "High"; c="DarkRed";  FillColor=c; LineColor=c; output;
Value = "Low";  c="DarkBlue"; FillColor=c; LineColor=c; output;
run;

/* 3. Use a HIGHLOW plot to graph the deviations from the reference value */
title "Deviations from Reference Value (&RefValue)";
title2 "Morning Fasting Blood Glucose";
ods graphics / width=600px height=400px;
proc sgplot data=Center DATTRMAP=DAttrs noautolegend;
   highlow x=Date low=Low high=High / group=Group ATTRID=HighLow;
   refline &RefValue / axis=y;
   yaxis grid label="Blood Glucose Level";
run;


/***************************************************************/
/* Do the same for a time series for a full year               */
/* of daily temperature averages.                               */
/***************************************************************/


/* 0. Read in raw data. Use PROC MEANS to calculate daily average temperature */
title;
%inc "RDUTemp.sas";
proc means data=RDUTemp noprint;
   by Date;
   var Temperature;
   output out=RDUDailyAvg mean=Temp;
run;

/* 1. Compute the deviation and encode the data as 'High' or 'Low' by using the reference value */
%let RefValue = 65;

data Center;
length Group $10;
set RDUDailyAvg;
if (Temp > &RefValue) then Group="Warm";
else Group="Cold";
Low  = min(Temp, &RefValue);    /* lower end of highlow bar */
High = max(Temp, &RefValue);    /* upper end of highlow bar */
run;

/* 2. Create a discrete attribute map that maps values to colors */
data DAttrs;
length c FillColor LineColor $16.;
ID = "WarmCold";              /* the ID value */
Value = "Warm"; c="DarkRed";  FillColor=c; LineColor=c; output;
Value = "Cold"; c="DarkBlue"; FillColor=c; LineColor=c; output;
run; 

/* 3. Use a HIGHLOW plot to graph the deviations from the reference value */
title "Deviations from Reference Value (&RefValue)";
title2 "Daily Averages in Raleigh, NC";
ods graphics / width=900px height=400px;
proc sgplot data=Center DATTRMAP=DAttrs;
   highlow x=Date low=Low high=High / group=Group ATTRID=WarmCold;
   refline &RefValue / axis=y;
   keylegend / location=inside title="" across=1 sortorder=reverseauto;
   yaxis grid label="Daily Average Temperature (F)";
run;

