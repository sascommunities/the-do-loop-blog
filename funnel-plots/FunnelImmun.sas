/* SAS program to accompany the article 
   "A funnel plot for immunization rates"
   by Rick Wicklin, published 26NOV2018 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2018/11/26/funnel-plot-immunization-rates.html

   This program shows how to create a funnel plot for proportions.
   Data from Robert Allison's blog post
   https://blogs.sas.com/content/sastraining/2018/11/20/immunization-rates-in-north-carolina-schools/
*/


/**************************************************/
/* 1. Reproduce Robert Allison's data cleaning    */
/**************************************************/
%let name=nc_kindergarten_immunization;
filename odsout '.';

/*
From chicken pox article:
https://www.washingtonpost.com/nation/2018/11/19/anti-vaccination-stronghold-nc-hit-with-states-worst-chickenpox-outbreak-decades/?utm_term=.1307778bc9fd

'data' link from the article:
http://mediad.publicbroadcasting.net/p/wcqs/files/201808/2017-2018_kindergarten_immunization_reporting_data_by_school__003_.pdf
*/

data school_data1 county_data (drop = school_type) state_data (drop = school_type);
infile "C:\Users\frwick\Documents\My SAS Files\Blog\&name..txt" lrecl=200 pad;
input whole_line $ 1-200;
format percent_not_immunized percent7.1;
percent_not_immunized=.;
percent_not_immunized=input(scan(whole_line,-1,' '),percent7.1);
format enrollment comma8.0;
enrollment=.;
enrollment=input(scan(whole_line,-4,' '),comma8.0);
length school_type $10;
school_type=scan(whole_line,-5,' ');
length county_name $50;
county_name=scan(whole_line,-6,' ');
if county_name='Hanover' then county_name='New Hanover';
statecode='NC';
if index(whole_line,'Statewide Total')^=0 then output state_data;
else if index(whole_line,'Total')^=0 then output county_data;
else output school_data1;
run;

data school_data; set school_data1;
length school_name $100;
school_name=substr(whole_line,1,index(whole_line,trim(left(county_name))||' '||trim(left(school_type)))-2);
school_name=propcase(school_name);
/* 
It claimed 930 in the pdf.
This is unlikely, since a Google search says there are only 609 total students in their K-5.
(this throws off their % value, and county total, and probably county stats as well)
*/
if school_name not in ('China Grove Elementary School') then output; 
run;

/* merge in numeric county number */
proc sql noprint;

create table county_data as
select unique county_data.*, us_counties_attr.county
from county_data left join mapsgfk.us_counties_attr
on county_data.statecode=us_counties_attr.statecode
 and county_data.county_name=us_counties_attr.idname;

quit; run;

proc sort data=county_data out=county_data;
by descending percent_not_immunized county_name;
run;
data county_data; set county_data;
order=_n_;
run;

/* merge in the county data & order with the individual school dataset */
proc sql noprint;
create table school_data as
select school_data.*, county_data.order, county_data.percent_not_immunized as county_not_immunized
from school_data left join county_data
on school_data.county_name=county_data.county_name;
quit; run;

/***********************************************/
/* 2. Rick's analysis starts here              */
/***********************************************/

/* create funnel plot with 99.8% upper CL for proportions
   See  https://blogs.sas.com/content/iml/2011/11/23/funnel-plots-for-proportions.html
*/
data Schools;
set School_Data;
Trials = enrollment;
Events = round(Percent_not_immunized * enrollment); /* num at risk */
drop whole_line statecode order;
run;

%let DSName = Schools;

proc iml;
use &DSName;
read all var {"Events" "Trials"};
close;

/* 1. Compute observed proportions */
Proportion = Events/ Trials;
*call histogram(Proportion);

/* 2. compute overall proportion */
theta = sum(Events) / sum(Trials); 
Expected = theta * Trials; 
*print (sum(Events))[L="NumEvents"] (sum(Trials))[L="NumTrials"] theta;

/* write theta to macro to use as reference line */
call symputx("AvgProp", theta);

/* write results to data set */
results = Proportion || Expected;
labels = {"Proportion" "Expected"};
create Stats from results[colname=labels];
append from results;
close;

/* 3. compute confidence limits for range of sample size */
/*    Notice that these limits do not depend on the data,
      although we use theta = average proportion as a "target" */
n = T( (1:4) || do(5, 300, 5) );
p = {0.001 0.025 0.975 0.999}; /* lower/upper limits */
/* compute matrix with four columns, one for each CL */
limits = j(nrow(n), ncol(p));
do i = 1 to ncol(p);
   r = quantile("binom", p[i], theta, n);
   /* adjust quantiles according to method in Appendix A.1.1
      Spiegelhalter (2005) p. 1197 */
   numer = cdf("binom", r  , theta, n) - p[i];
   denom = cdf("binom", r  , theta, n) - 
           cdf("binom", r-1, theta, n);
   alpha = numer / denom;
   limits[,i] = (r-alpha)/n;
end;

/* write these control limits to data set */
results = n || limits;
labels = {"N" "L3sd" "L2sd" "U2sd" "U3sd"};
create Limits from results[colname=labels];
append from results;
close;
quit;

%put &=AvgProp;

/* merge data with proportions */
data &DSName.1; 
merge &DSName Stats; 
if School_type^='Federal';  /* simplify analysis dy discarding one school */
run;
proc sort data=&DSName.1; by descending School_type; run;

/* append control limits */
data Funnel; set &DSName.1 Limits; run;

/* https://blogs.sas.com/content/iml/2011/12/09/creating-tooltips-for-scatter-plots-with-proc-sgplot.html */
ods graphics on / imagemap tipmax=2500 attrpriority=none;

/* 4. Plot proportions versus sample size. Overlay control limits */
title "Proportion of NC Kindergarten Students without Required Immunizations";
footnote1 J=L "Based on Robert Allison's 2018 analysis:";
footnote2 J=L "blogs.sas.com/content/sastraining/2018/11/20/immunization-rates-in-north-carolina-schools/";
proc sgplot data=Funnel noautolegend;
styleattrs datasymbols=(CircleFilled TriangleFilled TriangleRightFilled);
scatter x=Trials y=Proportion / group=school_type name="Scat"
        tip=(School_Name School_Type County_Name Events Trials Proportion);
series x=N y=U3sd / lineattrs=(color=lightpink thickness=3)
     curvelabel="99.8% limits" curvelabelloc=outside;
refline &AvgProp / axis=y label="NC Avg";
yaxis grid values=(0 to 1 by 0.1); xaxis grid;
keylegend "Scat" / exclude=(" ");
run;

/**************************************************/
/* use WHERE statement to make plots for counties */
/**************************************************/

title2 "Mecklenburg County";
proc sgplot data=Funnel noautolegend;
where County_Name in ("Mecklenburg" " ");
styleattrs datasymbols=(CircleFilled TriangleFilled TriangleRightFilled);
scatter x=Trials y=Proportion / group=school_type name="Scat"
        tip=(School_Name School_Type County_Name Events Trials Proportion);
series x=N y=U3sd / lineattrs=(color=lightpink thickness=3)
     curvelabel="99.8% limits" curvelabelloc=outside;
refline &AvgProp / axis=y label="NC Avg";
yaxis grid values=(0 to 1 by 0.1); xaxis grid;
keylegend "Scat" / exclude=(" ");
run;

title2 "Wake County";
proc sgplot data=Funnel noautolegend;
where County_Name in ("Wake" " ");
styleattrs datasymbols=(CircleFilled TriangleFilled TriangleRightFilled);
scatter x=Trials y=Proportion / group=school_type name="Scat"
        tip=(School_Name School_Type County_Name Events Trials Proportion);
series x=N y=U3sd / lineattrs=(color=lightpink thickness=3)
     curvelabel="99.8% limits" curvelabelloc=outside;
refline &AvgProp / axis=y label="NC Avg";
yaxis grid values=(0 to 1 by 0.1); xaxis grid;
keylegend "Scat" / exclude=(" ");
run;

/************************************************/
/* create data for Public/Private/Charter panel */
/************************************************/
data Limits3; 
length School_type $10;
School_type = 'Public';
set Limits Limits(in=ds2) Limits(in=ds3);
if ds2 then School_type = 'Private';
if ds3 then School_type = 'Charter';
run;

data Funnel3; 
set &DSName.1(where=(School_type = 'Public')) Limits3(where=(School_type = 'Public'))
    &DSName.1(where=(School_type = 'Private')) Limits3(where=(School_type = 'Private'))
    &DSName.1(where=(School_type = 'Charter')) Limits3(where=(School_type = 'Charter'));
run;

proc sgpanel data=Funnel3 noautolegend;
panelby school_type / sort=descending;
styleattrs datasymbols=(CircleFilled TriangleFilled TriangleRightFilled);
scatter x=Trials y=Proportion / group=school_type name="Scat"
        tip=(School_Name School_Type County_Name Events Trials Proportion);
series x=N y=U3sd / lineattrs=(color=lightpink thickness=3);
refline &AvgProp / axis=y;
rowaxis grid values=(0 to 1 by 0.1); colaxis grid;
run;
