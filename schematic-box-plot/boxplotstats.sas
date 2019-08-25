/*************************************************************/
/* SAS program to accompany the article 
   "Annotate features of a schematic box plot in SGPLOT"
   by Rick Wicklin, published 28AUG2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/08/28/schematic-box-plot.html

   This program shows how to overlay text and curves lines that explain
   the main statistical features of a schematic box plot. The features are.
   The lower/upper fences, the lower/upper whiskers, the quartiles (Q1, 
   median, and Q3), and the mean. Also mark outliers that exceed the 
   lower/upper fences.  

   The steps of teh construction are:

	1. Plan the graph.
	2. Define the data.
	3. Compute values for the whiskers, outliers, and mean value.
	4. Compute values for quartiles and fences.
	5. Create an annotation data set for lines, arrows, and text.
*/

/* Step 1: Put text for outliers, whiskers, and mean on the right. Put
           text for quartiles and fences on the left.
*/
/* Step 2. Create example data */
data Have(keep= Group x);
call streaminit(1);
Group=1;         /* required for PROC BOXPLOT */
do i = 1 to 40;
   x = rand("Normal", 0, 1.5);  /* normal data */
   output;
end;
x = 6; output;   /* upper outlier      */
x = 10; output;  /* far upper outliers */
x = -4; output;  /* lower outlier      */
run;

/*  Alternative data sets:

data Have;
call streaminit(1);
Group=1;
do i = 1 to 40;
   x = rand("logNormal", 0, 1.5);
   output;
end;
run;

data Have(keep=x Group);
   set Sashelp.Cars(rename=(MPG_City=x));
   label x = "MPG_City";
   Group = 1;
run;
*/

/* Step 3: Compute values for box plot features:
   Mean, lower/upper whiskers, outliers
*/
proc boxplot data=Have;
   plot x*Group / boxstyle=schematic outbox=outbox;
run;

proc print data=outbox;
  where _Type_  in ('MIN', 'MAX' 'MEAN' 
             'HIGH' 'LOW' 'FARHIGH' 'FARLOW' 'HIWHISKR' 'LOWHISKR');
  var _Type_ _Value_;
run;

/* Step 3a. Solve potential problem: 
   If there are no lower outliers, the output data set does not
   contain an entry for LOWHISKR. Must look at data and assign LOWHISKR=MIN.
   Similarly, if there are no upper outliers, the output data set does not
   contain an entry for HIWHISKR. Must look at data and assign HIWHISKR=MAX.
   See
   https://blogs.sas.com/content/iml/2019/08/26/conditionally-append-to-sas-data-set.html
*/
data outBox2;
drop HighFound LowFound Min Max;   /* temporary variables */
retain HighFound LowFound 0        /* binary indicator variables */
       Min Max .;                  /* VALUE of 'Min' and 'Max' obs */
length _Type_ $12;
set outBox(keep=_Type_ _Value_ 
    where=(_Type_  in ('MIN' 'MEAN' 'MAX' 
          'HIGH' 'LOW' 'FARHIGH' 'FARLOW' 'HIWHISKR' 'LOWHISKR')))
    END=EOF;                       /* EOF is temporary indicator variable */
label _Type_= _Value_=;

/* remember the Min and Max values. Set value to missing so MIN/MAX
   do not appear on the graph. */
if      _Type_ = 'MIN' then do;  Min=_VALUE_;  _VALUE_=.; end;
else if _Type_ = 'MAX' then do;  Max=_VALUE_;  _VALUE_=.; end;
else if _Type_ = 'MEAN'     then _Type_ = 'Mean';
else if _Type_ = 'FARHIGH'  then _Type_ = 'Far Outlier';
else if _Type_ = 'FARLOW'   then _Type_ = 'Far Outlier';
else if _Type_ = 'HIGH'     then _Type_ = 'Outlier';
else if _Type_ = 'LOW'      then _Type_ = 'Outlier';
else if _Type_ = 'LOWHISKR' then do;
   LowFound = 1; _Type_ = 'Low Whisker';  /* Low found; no need to create */
end;
else if _Type_ = 'HIWHISKR' then do;
   HighFound = 1; _Type_ = 'High Whisker'; /* High found; no need to create */
end;
output;                      /* need OUTPUT because of EOF processing */
/* end-of-file processing: conditionally append new observations */
if EOF then do;
   if ^LowFound then do;     /* Low value not found. Add it. */
      _Type_ = 'Low Whisker'; _Value_ = Min; output;
   end;
   if ^HighFound then do;    /* High value not found. Add it. */
      _Type_ = 'High Whisker'; _Value_ = Max; output;
   end;
end;
run;

proc print data=outbox2;run;

data Schematic1;
   merge outBox2 Have;
run;

/* axis table with box plot:
https://blogs.sas.com/content/graphicallyspeaking/2015/12/23/box-plot-with-stat-table-and-markers/
*/
/* How are we doing? Display box plot with axis table on right */
ods graphics/ width=300px height=400px;
title "Schematic Box Plot: First Draft";
proc sgplot data=Schematic1;
   vbox x;
   yaxistable _Type_  / y=_Value_ nolabel valueattrs=(size=10) location=inside;
   yaxis display=none;            /* OPTIONAL: Suppress Y ticks and values */
   xaxis offsetmin=0 offsetmax=0;
run;


/* Step 4: Use PROC MEANS and a DATA step to compute the quantile stats */
proc means data=Have Q1 Median Q3 noprint;
   var x;
   output out=Q Q1=Q1 Median=Median Q3=Q3;
run;

option validvarname=ANY;    /* permit 'Upper Fence'n to be a var name */
data IQR;
set Q;
IQR = Q3 - Q1;
'Upper Fence'n = Q3 + 1.5*IQR;
'Lower Fence'n = Q1 - 1.5*IQR;
drop _TYPE_ _FREQ_ IQR;
run;

proc transpose data=IQR out=IQR2(rename=(COL1=Value2)) name=Stat;
run;

proc print data=IQR2; run;

data Schematic2;
   merge Schematic1 IQR2(keep=Stat Value2);
run;

/* How are we doing? Display box plot with two axis tables */
title "Schematic Box Plot: Second Draft";
proc sgplot data=Schematic2;
   vbox x;
   yaxistable _Type_  / y=_Value_ nolabel valueattrs=(size=10) location=inside;
   yaxistable Stat  / y=Value2 nolabel valueattrs=(size=10) location=inside
                      position=left valuejustify=right;
   yaxis display=none;            /* OPTIONAL: Suppress Y ticks and values */
   xaxis offsetmin=0 offsetmax=0;
run;

/* Step 5: Create annotation data set */
data anno; 
retain Function 'Line' 
       x1Space x2Space 'WallPercent  '
       y1Space y2Space 'DataValue    '
       LinePattern 2               /* short dashed line */
       x1 35 
       x2 65 
       y1 y2 0;
set IQR2(where=(upcase(Stat) contains 'FENCE'));
y1 = Value2; y2 = Value2;
run;

title "Schematic Box Plot: Final Version";
proc sgplot data=Schematic2 sganno=anno;
   vbox x;
   yaxistable _Type_  / y=_Value_ nolabel valueattrs=(size=10) location=inside;
   yaxistable Stat  / y=Value2 nolabel valueattrs=(size=10) location=inside
                      position=left valuejustify=right;
   yaxis display=none;            /* OPTIONAL: Suppress Y ticks and values */
   xaxis offsetmin=0.2 offsetmax=0;
run;
