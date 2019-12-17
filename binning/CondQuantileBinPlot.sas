/* SAS program to accompany the article 
   "Create a conditional quantile bin plot in SAS"
   by Rick Wicklin, published 18DEC2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/12/18/conditional-quantile-bin-plot-sas.html

   This program creates a variation on the quantile bin plot, as shown in 
   https://blogs.sas.com/content/iml/2014/09/24/quantile-bin-plot.html

   In the conditional quantile bin plot, each cell has approximately the 
   same number of observations. However, the cells do not form a regular grid.
*/

/* Use PROC RANK to get unconditional quantiles for Y. 
   Use PROC RANK with a BY statement to get conditional quantiles for X 
*/
%let XVar = MomWtGain; /* name of X variable */
%let YVar = Weight;    /* name of Y variable */
%let NumXDiv = 3;      /* number of subdivisions (quantiles) for the X variable */
%let NumYDiv = 4;      /* number of subdivisions (quantiles) for the Y variable */
 
/* create the example data */
data Have;
set sashelp.bweight(obs=5000 keep=&XVar &YVar);
if ^cmiss(&XVar, &YVar);   /* only use complete cases */
run;


/* optional: get quantiles of Y for manually adding reflines to plot */
proc stdize data=Have PctlMtd=ORD_STAT outstat=StdLongPctls
           pctlpts=25 50 75;
var &YVar;
run;
 proc print data=StdLongPctls noobs;
where _type_ =: 'P';
run;

/*************************************/ 
/* group the Y values into quantiles */
proc rank data=Have out=RankY groups=&NumYDiv ties=high;
   var &YVar;
   ranks rank_&YVar;
run;
proc sort data=RankY;  by rank_&YVar;  run;
 
/* Conditioned on each Y quantile, group the X values into quantiles */
proc rank data=RankY out=RankBoth groups=&NumXDiv ties=high;
   by rank_&YVar;
   var &XVar;
   ranks rank_&XVar;
run;
proc sort data=RankBoth;  by rank_&YVar rank_&XVar;  run;
 
/* Trick: Use the CATT fuction to create a group for each unique (RankY, RankX) combination */
data Points;
set RankBoth;
Group = catt(rank_&YVar, rank_&XVar);
run;
 
proc sgplot data=Points;
   scatter x=&XVar y=&YVar / group=Group markerattrs=(symbol=CircleFilled) transparency=0.75;
   /* These REFLINESs are hard-coded for these data. They are the unconditional Y quantiles. */
   refline 3061 3402 3716 / axis=Y;
run;


/*************************************/

/* Create polygons from ranks. Each polygon is a rectangle of the form
   [min(X), max(X)] x [min(Y), max(Y)] 
*/
/* generate min/max for each group. Results are in wide form. */
proc means data=Points noprint;
   by Group;
   var &XVar &YVar;
   output out=Limits(drop=_Type_) min= max= / autoname;  
run;
 
proc print data=Limits noobs;
run;

proc means data=Limits Min Max;
   var _FREQ_;
run;

/* generate polygons from min/max of groups. Write data in long form. */
data Rectangles;
set Limits;
xx = &XVar._min; yy = &YVar._min; output;
xx = &XVar._max; yy = &YVar._min; output;
xx = &XVar._max; yy = &YVar._max; output;
xx = &XVar._min; yy = &YVar._max; output;
xx = &XVar._min; yy = &YVar._min; output;   /* repeat the first obs to close the polygon */
drop &XVar._min &XVar._max &YVar._min &YVar._max;
run;
 
/* combine the data and the polygons. Plot the results */
data All;
set Points Rectangles;
label _FREQ_ = "Frequency";
run;

/************************************/

title "Conditional Quantile Bin Plot";
proc sgplot data=All noautolegend;
   scatter x=&XVar y=&YVar / markerattrs=(symbol=CircleFilled) transparency=0.9;
   polygon x=xx y=yy ID=Group;  /* plot in the foreground, no fill */
run;

proc sgplot data=All noautolegend;
   polygon x=xx y=yy ID=Group / fill outline group=Group; /* background, filled */
   scatter x=&XVar y=&YVar / markerattrs=(symbol=CircleFilled) transparency=0.9;
run;

proc sgplot data=All;
   polygon x=xx y=yy ID=Group / fill outline colorresponse=_FREQ_; /* background, filled */
   scatter x=&XVar y=&YVar / markerattrs=(symbol=CircleFilled) transparency=0.9;
   gradlegend;
run;
title;
