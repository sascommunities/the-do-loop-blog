/* SAS program to accompany the article 
   "Visualize collinearity diagnostics"
   by Rick Wicklin, published 17FEB2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/02/17/visualize-collinearity-diagnostics.html

   This program shows how to create a heat map that visualizes
   the collinearity diagnostics that are produced by PROC REG in SAS.

   Main steps of the program:
   1. Use 
      ods output CollinDiag = Collin;
      to save the Collinearity Diagnostics table to a data set.
   2. Use PROC FORMAT to define a format. The format converts
      the table values into discrete values:
      A. The condition indices (which are in the range [1, infinity])
         ==> condLow, condMed, condHigh, condExtreme
      B. The proportion of variance values (which are in the range [0, 1))
         ==> propLow, propMed, propHigh
   3. Convert the data set from wide form to long form.
   4. Create a discrete attribute map that maps categories to colors.
   5. Create a discrete heat map that visualizes the collinearity
      diagnostics. Overlay (rounded) values for the condition indices
      and the important (relatively large) values of the proportion of variance. 

   Data from the PROC REG documentation.
*/
ods graphics/reset;

data fitness;
   input Age Weight Oxygen RunTime RestPulse RunPulse MaxPulse @@;
   datalines;
44 89.47 44.609 11.37 62 178 182   40 75.07 45.313 10.07 62 185 185
44 85.84 54.297  8.65 45 156 168   42 68.15 59.571  8.17 40 166 172
38 89.02 49.874  9.22 55 178 180   47 77.45 44.811 11.63 58 176 176
40 75.98 45.681 11.95 70 176 180   43 81.19 49.091 10.85 64 162 170
44 81.42 39.442 13.08 63 174 176   38 81.87 60.055  8.63 48 170 186
44 73.03 50.541 10.13 45 168 168   45 87.66 37.388 14.03 56 186 192
45 66.45 44.754 11.12 51 176 176   47 79.15 47.273 10.60 47 162 164
54 83.12 51.855 10.33 50 166 170   49 81.42 49.156  8.95 44 180 185
51 69.63 40.836 10.95 57 168 172   51 77.91 46.672 10.00 48 162 168
48 91.63 46.774 10.25 48 162 164   49 73.37 50.388 10.08 67 168 168
57 73.37 39.407 12.63 58 174 176   54 79.38 46.080 11.17 62 156 165
52 76.32 45.441  9.63 48 164 166   50 70.87 54.625  8.92 48 146 155
51 67.25 45.118 11.08 48 172 172   54 91.63 39.203 12.88 44 168 172
51 73.71 45.790 10.47 59 186 188   57 59.08 50.545  9.93 49 148 155
49 76.32 48.673  9.40 56 186 188   48 61.24 47.920 11.50 52 170 176
52 82.78 47.467 10.50 53 170 172
;
 
/* Run regression collinearity diagnostics:
   https://blogs.sas.com/content/iml/2020/01/23/collinearity-regression-collin-option.html
*/
proc reg data=fitness plots=none;
   model Oxygen = RunTime Age Weight RunPulse MaxPulse RestPulse / collin;
   ods select ParameterEstimates CollinDiag;
   ods output CollinDiag = Collin;
quit;

/* convert variable ranges into strings for the discrete heat map.
   For proportions, the ranges are determined by the cut points 
   [0.0, 0.4) white
   [0.4, 0.5) pink
   [0.5, 1)   red
   For condition indices, the ranges are determined by the cut points 
   [1, 20)    green
   [20, 30)   yellow
   [30, 100)  orange
   [100, infinity) red
*/
proc format;
value CollinFmt 
  0.0 -<  0.4  = "PropLow"
  0.4 -<  0.5  = "PropMed"
  0.5 -<  1.0  = "PropHigh"
  1.0 -<  20   = "CondLow"
   20 -<  30   = "CondMed"
   30 -<  40   = "CondHigh"
  100 -  1E9   = "CondExtreme";
run;

/* Convert from wide to long for graphing:
   https://blogs.sas.com/content/iml/2019/06/03/graph-wide-data-long-data-sas-sgplot.html
   The variables in the long data set are:
   Number  : The original row number
   VarName : The name of the original (wide) variables
   Val     : The numerical value of VarName for the row
   TextVal : A string value to be overlaid on the heat map.
           : For the condition index, use a 5.0 format
           : For the proportions of variance, use 
           : a blank if Val < 0.4, otherwise use the 
           : percentage round(100*Val)
*/
data CollinLong;
length VarName $32 TextVal $5;
set collin(drop=Model Dependent Eigenvalue);
array v[*] _numeric_;
do j = 2 to dim(v);
   VarName = vname(v[j]);
   Val = v[j];
   if j=2 then TextVal = put(Val, 5.0);
   else do;
      if Val<0.4 then TextVal = ' ';
      else TextVal = put(100*Val, 2.0); *convert proportion to pct;
   end;
   output;
end;
keep Number VarName Val TextVal;
run;

/* create a discrete attribute map:
   https://blogs.sas.com/content/iml/2019/07/15/create-discrete-heat-map-sgplot.html
*/
data Order;                            /* create discrete attribute map */
length Value $15 FillColor $15;
input raw FillColor;
Value = put(raw, CollinFmt.);          /* use format to assign values */
retain ID 'SortOrder'                  /* name of map */
     Show 'AttrMap';                   /* always show all groups in legend */
datalines;
0   White 
0.4 VeryLightRed
0.5 ModerateRed
1   VeryLightGreen 
20  VeryLightYellow
30  LightOrange
100 CXF03B20
;

/* Create a discrete heat map and overlay text for the condition numbers 
   and the important (high) cells of the proportion of variance. */
title "Collinearity Diagnostics";
proc sgplot data=CollinLong dattrmap=Order noautolegend;
  format Val CollinFmt.;
  heatmapparm x=VarName y=Number colorgroup=Val / outline 
              attrid=SortOrder nomissingcolor discretey;
  text x=VarName y=Number text=TextVal / strip textattrs=(size=14);
  refline "ConditionIndex" / axis=x lineattrs=(color=DarkGray)
          discreteoffset=0.5 discretethickness=0.08;
  xaxis display=(nolabel);
run;

