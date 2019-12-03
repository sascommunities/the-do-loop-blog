/* SAS program to accompany the article 
   "Longitudinal data: The response-profile model"
   by Rick Wicklin, published 03DEC2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/12/03/longitudinal-data-response-profile-model.html
   -and-
   "Longitudinal data: The mixed model"
   by Rick Wicklin, published 05DEC2019:
   https://blogs.sas.com/content/iml/2019/12/05/longitudinal-data-mixed-model.html

   This program shows how to analyze longitudinal data by following the 
   the article 
   "A Primer in Longitudinal Data Analysis"
   G. Fitzmaurice and C. Ravichandran (2008), Circulation, 118(19), p. 2005-2010.

   You can download the TLC data from the web site for the book 
   Applied Longitudinal Analysis (2011, 2nd Ed) 
   G. Fitzmaurice, N. Laird, and J. Ware.
   https://content.sph.harvard.edu/fitzmaur/ala2e/
*/

/* Read data into SAS, and create a dataset (TLC) in long format */

data tlc;
   input id Treatment $ lead0 lead1 lead4 lead6;
   y=lead0; Time=0; output;
   y=lead1; Time=1; output;
   y=lead4; Time=4; output;
   y=lead6; Time=6; output;
   drop lead0 lead1 lead4 lead6;
   label y = "Blood Lead Level (mcg/dL)";
datalines;
                1      P     30.8    26.9    25.8    23.8
                2      A     26.5    14.8    19.5    21.0
                3      A     25.8    23.0    19.1    23.2
                4      P     24.7    24.5    22.0    22.5
                5      A     20.4     2.8     3.2     9.4
                6      A     20.4     5.4     4.5    11.9
                7      P     28.6    20.8    19.2    18.4
                8      P     33.7    31.6    28.5    25.1
                9      P     19.7    14.9    15.3    14.7
               10      P     31.1    31.2    29.2    30.1
               11      P     19.8    17.5    20.5    27.5
               12      A     24.8    23.1    24.6    30.9
               13      P     21.4    26.3    19.5    19.0
               14      A     27.9     6.3    18.5    16.3
               15      P     21.1    20.3    18.4    20.8
               16      P     20.6    23.9    19.0    17.0
               17      P     24.0    16.7    21.7    20.3
               18      P     37.6    33.7    34.4    31.4
               19      A     35.3    25.5    26.3    30.3
               20      A     28.6    15.8    22.9    25.9
               21      P     31.9    27.9    27.3    34.2
               22      A     29.6    15.8    23.7    23.4
               23      A     21.5     6.5     7.1    16.0
               24      P     26.2    26.8    25.3    24.8
               25      A     21.8    12.0    16.8    19.2
               26      A     23.0     4.2     4.0    16.2
               27      A     22.2    11.5     9.5    14.5
               28      P     20.5    21.1    17.4    21.1
               29      A     25.0     3.9    12.8    12.7
               30      P     33.3    26.2    34.0    28.2
               31      A     26.0    21.4    21.0    22.4
               32      A     19.7    13.2    14.6    11.6
               33      P     27.9    21.6    23.6    27.7
               34      P     24.7    21.2    22.9    21.9
               35      P     28.8    26.4    23.8    22.0
               36      A     29.6    17.5    21.0    24.2
               37      P     32.0    30.2    30.2    27.5
               38      P     21.8    19.3    16.4    17.6
               39      A     24.4    16.4    11.6    16.6
               40      A     33.7    14.9    14.5    63.9
               41      P     24.9    20.9    22.2    19.8
               42      P     19.8    18.9    18.9    15.5
               43      A     26.7     6.4     5.1    15.1
               44      A     26.8    20.4    19.3    23.8
               45      A     20.2    10.6     9.0    16.0
               46      P     35.4    30.4    26.5    28.1
               47      P     25.3    23.9    22.2    27.2
               48      A     20.2    17.5    17.4    18.6
               49      A     24.5    10.0    15.6    15.2
               50      P     20.3    21.0    16.7    13.5
               51      P     20.4    17.2    15.9    17.7
               52      P     24.1    20.1    17.9    18.7
               53      A     27.1    14.9    18.1    21.3
               54      A     34.7    39.0    28.8    34.7
               55      P     28.5    32.6    27.5    22.8
               56      P     26.6    22.4    21.8    21.0
               57      A     24.5     5.1     8.2    23.6
               58      P     20.5    17.5    19.6    18.4
               59      P     25.2    25.1    23.4    22.2
               60      P     34.7    39.5    38.6    43.3
               61      P     30.3    29.4    33.1    28.4
               62      P     26.6    25.3    25.1    27.9
               63      P     20.7    19.3    21.9    21.8
               64      A     27.7     4.0     4.2    11.7
               65      A     24.3    24.3    18.4    27.8
               66      A     36.6    23.3    40.4    39.3
               67      P     28.9    28.9    32.8    31.8
               68      A     34.0    10.7    12.6    21.2
               69      A     32.6    19.0    16.3    18.6
               70      A     29.2     9.2     8.3    18.4
               71      A     26.4    15.3    24.6    32.4
               72      A     21.8    10.6    14.4    18.7
               73      P     27.2    28.5    35.0    30.5
               74      P     22.4    22.0    19.1    18.7
               75      P     32.5    25.1    27.8    27.3
               76      P     24.9    23.6    21.2    21.1
               77      P     24.6    25.0    21.7    23.9
               78      P     23.1    20.9    21.7    19.9
               79      A     21.1     5.6     7.3    12.3
               80      P     25.8    21.9    23.6    24.8
               81      P     30.0    27.6    24.0    23.7
               82      A     22.1    21.0     8.6    24.6
               83      P     20.0    22.7    21.2    20.5
               84      P     38.1    40.8    38.0    32.7
               85      A     28.9    12.5    16.7    22.2
               86      P     25.1    28.1    27.5    24.8
               87      A     19.8    11.6    13.0    23.1
               88      P     22.1    21.1    21.5    20.6
               89      A     23.5     7.9    12.4    18.9
               90      A     29.1    16.8    15.1    18.8
               91      A     30.3     3.5     3.0    11.5
               92      P     25.4    24.3    22.7    20.1
               93      A     30.6    28.2    27.0    25.5
               94      A     22.4     7.1    17.2    18.7
               95      A     31.2    10.8    19.8    22.2
               96      A     31.4     3.9     7.0    17.8
               97      A     41.1    15.1    10.9    27.1
               98      A     29.4    22.1    25.3     4.1
               99      A     21.9     7.6    10.8    13.0
              100      A     20.7     8.1    25.7    12.3
;

ods graphics /reset;

/* Create discrete data map for attributes:  
   https://blogs.sas.com/content/graphicallyspeaking/2016/09/13/legend-order-and-group-attributes/ */
data Order;
   input Value $ n;
   retain ID 'Treat' Show 'AttrMap';
   FillStyle        = cats('GraphData', n);
   LineStyle        = cats('GraphData', n);
   MarkerStyle      = cats('GraphData', n);
   datalines;
A 2
P 1
;

/* Use PROC GLM to perform analysis of response profiles on the TLC data */

/* Plot means vs discrete time periods:  
   https://blogs.sas.com/content/iml/2019/10/07/graph-means-vs-time.html
*/
/* For information about the LEGENDITEM statement, see
   https://blogs.sas.com/content/iml/2018/02/12/merged-legends.html 
*/

title "Response Profiles of Children Exposed to Lead";
/* Spaghetti plot of the 100 response profiles */
proc sgplot data=tlc dattrmap=Order;
   series x=Time y=Y / group=ID groupLC=Treatment break lineattrs=(pattern=solid)
                       attrid=Treat;
   /* LEGENDITEM is a SAS 9.4M5 feature. Delete the following statements in older versions of SAS */
   legenditem type=line name="P" / label="Placebo (P)" lineattrs=GraphData1; 
   legenditem type=line name="A" / label="Succimer (A)" lineattrs=GraphData2; 
   keylegend "A" "P";
   xaxis values=(0 1 4 6) grid;
run;

proc sgplot data=tlc dattrmap=Order;
   vline Time / response=y group=Treatment stat=mean limitstat=stderr markers attrid=Treat;
run;
/*
proc sgplot data=tlc;
vbox y / category=Time group=Treatment connect=mean groupdisplay=cluster 
            connect=mean clusterwidth=0.35;
run;
*/

/* Do the two Treatments differ in their mean change from baseline?
   With PROC GLM, you can answer yes or no, but you do not know HOW they differ.
*/
/* NOTE: The ORDER of the variable are important for the graph
   if you use  plots=(intplot(clm)) */
ods trace off;
proc glm data=tlc;
   class Time(ref='0') Treatment(ref='P');
   model y = Treatment Time Treatment*Time / solution;
   *ods select ModelANOVA 'Type III Model Anova' ParameterEstimates IntPlot;
   output out=GLMOut Predicted=Pred;
quit;

/* Look at predictions for each individual.
   They are identical! */
proc sort data=GLMOut(where=(ID in (1,2,4,5,6,7))) out=GLMSubset;
   by Treatment ID;
run;

title2 "Fixed Effect Model";
proc sgpanel data=GLMSubset dattrmap=Order;
   panelby Treatment ID;
   scatter x=Time y=y / group=Treatment attrid=Treat;
   series x=Time y=Pred / group=Treatment attrid=Treat;
run;

/**************************************************/
/* You can model the intercept of individuals as a random variable. 
   In more sophisticated analyses, you can use random effects to 
   model clusters such as classrooms or hospitals.
*/
%macro SortAndPlot(DSName);
proc sort data=&DSName;
   by descending Treatment ID Time;
run;

proc sgplot data=&DSName dattrmap=Order;
   series x=Time y=Pred / group=ID groupLC=Treatment break lineattrs=(pattern=solid)
                       attrid=Treat;
   legenditem type=line name="P" / label="Placebo (P)" lineattrs=GraphData1; 
   legenditem type=line name="A" / label="Succimer (A)" lineattrs=GraphData2; 
   keylegend "A" "P";
   xaxis values=(0 1 4 6) grid;
   yaxis label="Predicted Blood Lead Level (mcg/dL)";
run;
%mend;

/* Make "discrete time" (t) to use in REPEATED statement.
   Make spline effect with knot at t=1. */
data TLCView / view=TLCView;
set tlc;
t = Time;       /* discrete copy of time */
T1 = ifn(Time<1, 0, Time - 1);  /* knot at Time=1 for PWL analysis */
run;

/* Repeat the response-profile analysis, but 
   use the RANDOM statement to add random intercept for each subject */
proc mixed data=TLCView;
   class id Time(ref='0') Treatment(ref='P');
   model y = Treatment Time Treatment*Time / s chisq outpred=MixedOut;
   repeated Time / type=un subject=id r;   /* measurements are repeated for subjects */
   random intercept / subject=id;          /* each subject gets its own intercept */
run;
 
title "Predicted Individual Growth Curves";
title2 "Random Intercept Model";
%SortAndPlot(MixedOut);

/* Model time as continuous and use a quadratic model in Time. 
   For more about quadratic growth models, see
   https://support.sas.com/resources/papers/proceedings/proceedings/sugi27/p253-27.pdf */
proc mixed data=TLCView;
   class id t(ref='0') Treatment(ref='P');
   model y = Treatment Time Time*Time Treatment*Time / s outpred=MixedOutQuad;
   repeated t / type=un subject=id r;      /* measurements are repeated for subjects */
   random intercept / subject=id;          /* each subject gets its own intercept */
run;
 
title2 "Random Intercept; Quadratic in Time";
%SortAndPlot(MixedOutQuad);

/* Piecewise linear (PWL) model with knot at Time=1.
   For more about PWL models, see Hwang (2015) 
   "Hands-on Tutorial for Piecewise Linear Mixed-effects Models Using SAS PROC MIXED"
   https://www.lexjansen.com/pharmasug-cn/2015/ST/PharmaSUG-China-2015-ST08.pdf    */
proc mixed data=TLCView;
   class id t(ref='0') Treatment(ref='P');
   model y = Treatment Time T1 Treatment*Time Treatment*T1 / s outpred=MixedOutPWL;
   repeated t / type=un subject=id r;      /* measurements are repeated for subjects */
   random intercept / subject=id;          /* each subject gets its own intercept */
run;
 
title2 "Random Intercept; Piecewise Linear Model";
%SortAndPlot(MixedOutPWL);
