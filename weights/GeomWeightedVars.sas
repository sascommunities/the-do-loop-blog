/* SAS program to accompany the article 
   "Rankings and the geometry of weighted averages"
   by Rick Wicklin, published 16AUG2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/08/16/geometry-weighted-averages.html

   This program demonstrates the geometry of a weighted average of
   variables. The weight vectors determines a line in d-dimensional
   space. The projection of observations onto the line determines the
   ranking for those weights.
*/


proc iml;
Subject = {A, B, C};
M = {10 2,
      7 8,
      9 5};
*print M[r=Subject c={'Test 1' 'Test 2'} L="Test Scores"];

/* make a nice display by using a table */
Tbl = TableCreate("Student", Subject, {"Test1" "Test2"}, M);
call TableSetVarLabel(Tbl, {"Test1" "Test2"}, {"Test 1" "Test 2"});
call TablePrint(Tbl) colheader="Labels" ID="Student" label="Test Scores";

/* save to a data set */
create Pts from M[r=Subject c={'Test1' 'Test2'}];
append from M[r=Subject];
close;

w1 = {1 0};     /* Test1 is all that is important */
score1 = M*w1`;

w2 = {1 1};     /* equal importance = total score */
score2 = M*w2`;

w3 = {2 1};     /* Test1 is twice as important */
score3 = M*w3`;

result = M || score1 || score2 || score3;
*print result[r=Subject c={'Test1' 'Test2' 'Score1' 'Score2' 'Score3'}];

/* make a nice display by adding these scores to the table */
call TableAddVar(Tbl, {'Score1' 'Score2' 'Score3'}, score1 || score2 || score3);
call TableSetVarLabel(Tbl, {'Score1' 'Score2' 'Score3'}, {"Wt Sum 1" "Wt Sum 2" "Wt Sum 3"});
call TablePrint(Tbl) colheader="Labels" ID="Student" label="Weighted Sums of Test Scores";

/* write weight vectors to a data set */
create V1 from w1[c={'vx' 'vy'}]; append from w1; close;
create V2 from w2[c={'vx' 'vy'}]; append from w2; close;
create V3 from w3[c={'vx' 'vy'}]; append from w3; close;
QUIT;

/******************/

ods graphics /width=500px height=500px;
%let LightGray = CXEEEEEE;

/* GRAPH for w1: Weight vector, observations, and projections */
data P1;
label Test1="Var 1" Test2="Var 2";
set Pts V1;
vvx = 10*vx; vvy = 10*vy;
run;
title "Scores for Weight=(1, 0)";
proc sgplot data=P1 aspect=1 noautolegend;
scatter x=Test1 y=Test2 / datalabel=Subject markerattrs=(symbol=CircleFilled);
vector x=vvx y=vvy / noarrowheads lineattrs=(color=Blue pattern=Dash);
vector x=vx y=vy;
refline 7 8 9 10 / axis=x label lineattrs=(color=Gray);
xaxis min=0 grid gridattrs=(color=&LightGray) values=(0 to 10 by 2);
yaxis min=0 grid gridattrs=(color=&LightGray) values=(0 to 10 by 2);
run;


/* GRAPH for w2: Weight vector, observations, and projections */
data Lines2;
array values[4] _temporary_ (12 13 14 15);
array w[2]      _temporary_ (1 1);
do i = 1 to dim(values);
   x0=0; y0 = values[i] / w[2];  /* y intercept */
   lx = values[i] / w[1]; ly=0; 
   value = values[i];
   output;
end;
drop i;
run;
data P2;
label Test1="Var 1" Test2="Var 2";
set Pts V2 Lines2;
vvx = 10*vx; vvy = 10*vy;
run;

title "Scores for Weight=(1, 1)";
proc sgplot data=P2 aspect=1 noautolegend;
scatter x=Test1 y=Test2 / datalabel=Subject markerattrs=(symbol=CircleFilled);
vector x=lx y=ly / xorigin=x0 yorigin=y0 noarrowheads 
            datalabel=value datalabelpos=bottom lineattrs=(color=Gray);
vector x=vvx y=vvy / noarrowheads lineattrs=(color=Blue pattern=Dash);
vector x=vx y=vy;
xaxis min=0 max=15 grid gridattrs=(color=&LightGray) offsetmin=0.0 values=(0 to 14 by 2) valueshint;
yaxis min=0 max=15 grid gridattrs=(color=&LightGray) offsetmin=0.05 values=(0 to 14 by 2) valueshint;
run;


/* GRAPH for w3: Weight vector, observations, and projections */
data Lines3;
array values[4] _temporary_ (20 21 22 23);
array w[2]      _temporary_ (2 1);
do i = 1 to dim(values);
   x0=0; y0 = values[i] / w[2];  /* y intercept */
   lx = values[i] / w[1]; ly=0; 
   value = values[i];
   output;
end;
drop i;
run;
data P3;
label Test1="Var 1" Test2="Var 2";
set Pts V3 Lines3;
vvx = 6*vx; vvy = 6*vy;
run;

title "Scores for Weight=(2, 1)";
proc sgplot data=P3 aspect=1 noautolegend;
scatter x=Test1 y=Test2 / datalabel=Subject markerattrs=(symbol=CircleFilled);
vector x=lx y=ly / xorigin=x0 yorigin=y0 noarrowheads 
            datalabel=value datalabelpos=bottom lineattrs=(color=Gray);
vector x=vvx y=vvy / noarrowheads lineattrs=(color=Blue pattern=Dash);
vector x=vx y=vy;
xaxis min=0 max=12 grid gridattrs=(color=&LightGray) offsetmin=0.0 values=(0 to 12 by 2) ;
yaxis min=0 max=12 grid gridattrs=(color=&LightGray) offsetmin=0.05 values=(0 to 12 by 2);
run;
