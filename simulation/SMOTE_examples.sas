/* SAS program to accompany the article 
   "Implement a SMOTE simulation algorithm in SAS"
   by Rick Wicklin, published 05MAY2025 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2025/05/05/smote-simulation-sas.html

   This program shows how to implement a simple SMOTE algorithm in SAS 
   from "first principals." The algorithm assumes all variables are continuous.

   First, define the helper functions in the file SMOTE_define.sas.
   An easy way to do that is 
   %INCLUDE "<path/>SMOTE_define.sas";
   proc iml;
   load module=(NearestNbr RandRowsAndNbrs SMOTESimContin);
   ...use the functions...
*/

/* Demonstrate SMOTE simulation on a small data set */

data SmallData;
input ID $ x y;
datalines;
X1  0   1
X2  0.6 0.8
X3 -0.1 0.5
X4  0.3 0.4
X5 -0.4 0.1
X6 -1   0 
X7  1   0
;

ods graphics / PUSH width=480px height=280px;
title "Sample Data for SMOTE Algorithm";
proc sgplot data=SmallData aspect=0.5;
   scatter x=x y=y / datalabel=ID datalabelattrs=(size=14) markerattrs=(symbol=CircleFilled size=12);
   xaxis grid;
   yaxis grid;
run;

/* -------------------------------------------------------- */
/* OPTIONAL: Visualize the geometry of the SMOTE simulation */
data lineseg;
length LID $2;
input LID $ lx ly ux uy;
datalines;
P -0.1 0.5 . .
Q  0.3 0.4 . .
U  .   .   0.2 0.425
;
data Example;
set SmallData lineseg;
run;

title "Synthetic Data for SMOTE Algorithm";
proc sgplot data=Example aspect=0.5 noautolegend;
   scatter x=x y=y / datalabel=ID datalabelattrs=(size=14) markerattrs=(symbol=CircleFilled size=12);
   series x=lx y=ly / datalabel=LID datalabelattrs=GraphData2(size=14) markers markerattrs=(size=0) lineattrs=GraphData2;
   scatter x=ux y=uy / markerattrs=(symbol=X size=10) markerattrs=GraphData2;
   xaxis grid;
   yaxis grid;
run;

/* -------------------------------------------------------- */


/* -------------------------------------------------------- */
/* OPTIONAL: visualize the complete graph on these nodes */
proc iml;
use SmallData;
read all var {x y} into Z;
close;

comb = allcomb(nrow(Z), 2);
row = Z[1,];
poly_ID = .;
create CompGraph from poly_ID row [c={'PolyID' 'Px' 'Py'}];
do i = 1 to nrow(comb);
   poly_ID = i;
   row = Z[ comb[i,1], ];
   append from poly_ID row;
   row = Z[ comb[i,2], ];
   append from poly_ID row;
end;
close;
quit;

data All;
set SmallData CompGraph;
run;

title "Complete Graph on Sample Data";
proc sgplot data=All aspect=0.5 noautolegend;
   scatter x=x y=y / datalabel=ID markerattrs=(symbol=CircleFilled size=12);
   polygon ID=PolyID x=Px y=Py;
   xaxis grid;
   yaxis grid;
run;

/* -------------------------------------------------------- */
/* How to call the helper routines for SMOTE simulation */
/* -------------------------------------------------------- */

/* -------------------------------------------------------- */
/* Test 1: call the NearestNbr routine with k=3 */
/* -------------------------------------------------------- */
proc iml;
load module=(NearestNbr RandRowsAndNbrs SMOTESimContin);
use SmallData;
   read all var {'x' 'y'} into X;
close;

k = 3;
run NearestNbr(NN, dist, X, k);
print NN[r=('X1':'X7')];

/* visualize nearest neighbors by connecting NN with lines */
pt = X[1,];
poly_ID = .;
create CompGraph from poly_ID pt [c={'PolyID' 'Px' 'Py'}];
poly_ID = 1;
do i = 1 to nrow(nn);
   pt1 = X[i, ];
   do j = 1 to ncol(nn);
      pt = pt1;
      append from poly_ID pt;
      pt = X[ nn[i,j], ];
      append from poly_ID pt;
      poly_ID = poly_ID + 1;
   end;
end;
close;
QUIT;

data Graph;
set SmallData CompGraph;
run;

title "Graph Connecting 3 Nearest Neighbors";
proc sgplot data=Graph aspect=0.5 noautolegend;
   scatter x=x y=y / datalabel=ID datalabelattrs=(size=14) markerattrs=(symbol=CircleFilled size=12);
   polygon ID=PolyID x=Px y=Py / transparency=0.8;
   xaxis grid;
   yaxis grid;
run;

/* -------------------------------------------------------- */
/* Test 2: call the RandRowsAndNbrs routine with T=20 */
/* -------------------------------------------------------- */

proc iml;
load module=(NearestNbr RandRowsAndNbrs SMOTESimContin);
use SmallData;
   read all var {'x' 'y'} into X;
close;

k = 3;
run NearestNbr(nn, dist, X, k);

call randseed(321);
T = 20;
R = RandRowsAndNbrs(T, NN);
print R[c={"P" "Q"}];

*Note: The following statement complete the SMOTE simulation;
/*  
p1 = X[R[,1], ];
p2 = X[R[,2], ];
u = randfun(T, "Uniform");                    * random variates in U[(0,1) ;
p_new = p1 + u#(p2 - p1);
*/
QUIT;


/* -------------------------------------------------------- */
/* Example: call the SMOTESimContin function with k=3 and T=20 */
/* -------------------------------------------------------- */

proc iml;
load module=(NearestNbr RandRowsAndNbrs SMOTESimContin);
use SmallData;
   read all var {'x' 'y'} into X;
close;

call randseed(321);
k = 3;
T = 20;
X_synth = SMOTESimContin(X, T, k);
print X_synth[c={'sx' 'sy'}];

create synthetic from X_synth[c={'sx' 'sy'}];
append from X_synth;
close;
QUIT;

data Graph;
set SmallData CompGraph synthetic;
run;

title "20 Synthetic Data Points by Using 3 Nearest Neighbors";
proc sgplot data=Graph aspect=0.5 noautolegend;
   scatter x=x y=y / datalabel=ID datalabelattrs=(size=14) markerattrs=(symbol=Circle size=12);
   polygon ID=PolyID x=Px y=Py / transparency=0.85;
   scatter x=sx y=sy / markerattrs=(symbol=X size=10);
   xaxis grid;
   yaxis grid;
run;
ods graphics / POP;

/* -------------------------------------------------------- */
/* simulate from sashelp.iris where Species='Versicolor' */
/* -------------------------------------------------------- */


/* extract one species of iris data */
data Versicolor;
set sashelp.iris(where=(species='Versicolor'));
attrib _all_ label=' ';
run;

proc iml;
load module=(NearestNbr RandRowsAndNbrs SMOTESimContin);
use Versicolor;
read all var _NUM_ into X[c=varNames];
close;

call randseed(321);
k = 5;
T = 100;
X_synth = SMOTESimContin(X, T, k);

create synth_iris from X_synth[c=varNames];
append from X_synth;
close;
QUIT;

/* compare descriptive statistics for the original and synthetic data */
title "50 Iris Data Points (Species='Versicolor')";
proc corr data=Versicolor noprob cov plots=matrix(histogram);
var _numeric_;
run;

title "100 Synthetic Data Points by Using 5 Nearest Neighbors";
proc corr data=synth_iris noprob cov plots=matrix(histogram);
var _numeric_;
run;

/* -------------------------------------------------------- */
/* simulate contin vars from sashelp.BWeight where MomEdLevel=0 and Visit=1 */
/* -------------------------------------------------------- */

ods trace off;
data BWeight;
set sashelp.bweight;
keep CigsPerDay MomAge MomWtGain Weight;
where MomEdLevel=0 and Visit=1;
attrib _all_ label=' ';
run;

proc iml;
load module=(NearestNbr RandRowsAndNbrs SMOTESimContin);
use BWeight;
read all var _NUM_ into X[c=varNames];
close;

call randseed(321);
k = 5;
T = 6000;
X_synth = SMOTESimContin(X, T, k);

create synth_bweight from X_synth[c=varNames];
append from X_synth;
close;
QUIT;

title "BWeight Data";
proc corr data=BWeight noprob cov plots(maxpoints=none)=matrix(histogram);
var _numeric_;
run;

title "6000 Synthetic Data Points by Using 5 Nearest Neighbors";
proc corr data=synth_bweight noprob cov plots(maxpoints=none)=matrix(histogram);
var _numeric_;
run;

/* -------------------------------------------------------- */
/* simulate vars from CRIME data in PROC PRINCOMP documentation */
/* -------------------------------------------------------- */

data Crime;
   input State $1-15 Murder Rape Robbery Assault
         Burglary Larceny Auto_Theft;
   datalines;
Alabama        14.2 25.2  96.8 278.3 1135.5 1881.9 280.7
Alaska         10.8 51.6  96.8 284.0 1331.7 3369.8 753.3
Arizona         9.5 34.2 138.2 312.3 2346.1 4467.4 439.5
Arkansas        8.8 27.6  83.2 203.4  972.6 1862.1 183.4
California     11.5 49.4 287.0 358.0 2139.4 3499.8 663.5
Colorado        6.3 42.0 170.7 292.9 1935.2 3903.2 477.1
Connecticut     4.2 16.8 129.5 131.8 1346.0 2620.7 593.2
Delaware        6.0 24.9 157.0 194.2 1682.6 3678.4 467.0
Florida        10.2 39.6 187.9 449.1 1859.9 3840.5 351.4
Georgia        11.7 31.1 140.5 256.5 1351.1 2170.2 297.9
Hawaii          7.2 25.5 128.0  64.1 1911.5 3920.4 489.4
Idaho           5.5 19.4  39.6 172.5 1050.8 2599.6 237.6
Illinois        9.9 21.8 211.3 209.0 1085.0 2828.5 528.6
Indiana         7.4 26.5 123.2 153.5 1086.2 2498.7 377.4
Iowa            2.3 10.6  41.2  89.8  812.5 2685.1 219.9
Kansas          6.6 22.0 100.7 180.5 1270.4 2739.3 244.3
Kentucky       10.1 19.1  81.1 123.3  872.2 1662.1 245.4
Louisiana      15.5 30.9 142.9 335.5 1165.5 2469.9 337.7
Maine           2.4 13.5  38.7 170.0 1253.1 2350.7 246.9
Maryland        8.0 34.8 292.1 358.9 1400.0 3177.7 428.5
Massachusetts   3.1 20.8 169.1 231.6 1532.2 2311.3 1140.1
Michigan        9.3 38.9 261.9 274.6 1522.7 3159.0 545.5
Minnesota       2.7 19.5  85.9  85.8 1134.7 2559.3 343.1
Mississippi    14.3 19.6  65.7 189.1  915.6 1239.9 144.4
Missouri        9.6 28.3 189.0 233.5 1318.3 2424.2 378.4
Montana         5.4 16.7  39.2 156.8  804.9 2773.2 309.2
Nebraska        3.9 18.1  64.7 112.7  760.0 2316.1 249.1
Nevada         15.8 49.1 323.1 355.0 2453.1 4212.6 559.2
New Hampshire   3.2 10.7  23.2  76.0 1041.7 2343.9 293.4
New Jersey      5.6 21.0 180.4 185.1 1435.8 2774.5 511.5
New Mexico      8.8 39.1 109.6 343.4 1418.7 3008.6 259.5
New York       10.7 29.4 472.6 319.1 1728.0 2782.0 745.8
North Carolina 10.6 17.0  61.3 318.3 1154.1 2037.8 192.1
North Dakota    0.9  9.0  13.3  43.8  446.1 1843.0 144.7
Ohio            7.8 27.3 190.5 181.1 1216.0 2696.8 400.4
Oklahoma        8.6 29.2  73.8 205.0 1288.2 2228.1 326.8
Oregon          4.9 39.9 124.1 286.9 1636.4 3506.1 388.9
Pennsylvania    5.6 19.0 130.3 128.0  877.5 1624.1 333.2
Rhode Island    3.6 10.5  86.5 201.0 1489.5 2844.1 791.4
South Carolina 11.9 33.0 105.9 485.3 1613.6 2342.4 245.1
South Dakota    2.0 13.5  17.9 155.7  570.5 1704.4 147.5
Tennessee      10.1 29.7 145.8 203.9 1259.7 1776.5 314.0
Texas          13.3 33.8 152.4 208.2 1603.1 2988.7 397.6
Utah            3.5 20.3  68.8 147.3 1171.6 3004.6 334.5
Vermont         1.4 15.9  30.8 101.2 1348.2 2201.0 265.2
Virginia        9.0 23.3  92.1 165.7  986.2 2521.2 226.7
Washington      4.3 39.6 106.2 224.8 1605.6 3386.9 360.3
West Virginia   6.0 13.2  42.2  90.9  597.4 1341.7 163.3
Wisconsin       2.8 12.9  52.2  63.7  846.9 2614.2 220.7
Wyoming         5.4 21.9  39.7 173.9  811.6 2772.2 282.0
;

proc iml;
load module=(NearestNbr RandRowsAndNbrs SMOTESimContin);
use Crime;
read all var _NUM_ into X[c=varNames];
close;

call randseed(321);
k = 5;
T = 100;
X_synth = SMOTESimContin(X, T, k);

create synth_crime from X_synth[c=varNames];
append from X_synth;
close;
QUIT;

title "Crime Data";
proc corr data=Crime noprob cov plots(maxpoints=none)=matrix(histogram);
var _numeric_;
run;

title "100 Synthetic Data Points by Using 5 Nearest Neighbors";
proc corr data=synth_crime noprob cov plots(maxpoints=none)=matrix(histogram);
var _numeric_;
run;

/* -------------------------------------------------------- */
/* simulate Longitude and Latitude for S04 satial data, 
   from the PROC LOESS documentation */
/* -------------------------------------------------------- */
data SO4;
   input Latitude Longitude SO4 @@;
   /* Longitudes decrease from west to east in the western hemisphere, 
      so use negative values when graphing longitudes in the western hemisphere. */
   Longitude = -Longitude; 
   format Latitude f4.0;
   format Longitude f4.0;
   format SO4 f4.1;
   label S04 = "Sulfate (g / m^2)";
   datalines;
32.45833  87.24222 1.403 34.28778  85.96889 2.103
33.07139 109.86472 0.299 36.07167 112.15500 0.304
31.95056 112.80000 0.263 33.60500  92.09722 1.950
34.17944  93.09861 2.168 36.08389  92.58694 1.578
36.10056  94.17333 1.708 39.00472 123.08472 0.096
36.56694 118.77722 0.259 41.76583 122.47833 0.065
34.80611 119.01139 0.053 38.53528 121.77500 0.135
37.44139 105.86528 0.247 38.11778 103.31611 0.326
39.99389 105.48000 0.687 39.40306 107.34111 0.225
39.42722 107.37972 0.339 37.19806 108.49028 0.559
40.50750 107.70194 0.250 40.36417 105.58194 0.307
40.53472 106.78000 0.564 37.75139 107.68528 0.557
39.10111 105.09194 0.371 40.80639 104.75472 0.286
29.97472  82.19806 1.279 28.54278  80.64444 1.564
25.39000  80.68000 0.912 30.54806  84.60083 1.243
27.38000  82.28389 0.991 32.14111  81.97139 1.225
33.17778  84.40611 1.580 31.47306  83.53306 0.851
43.46139 113.55472 0.180 43.20528 116.74917 0.103
44.29778 116.06361 0.161 40.05333  88.37194 2.940
41.84139  88.85111 2.090 41.70111  87.99528 3.171
37.71000  89.26889 2.523 38.71000  88.74917 2.317
37.43556  88.67194 3.077 40.93333  90.72306 2.006
40.84000  85.46389 2.725 38.74083  87.48556 3.158
41.63250  87.08778 3.443 40.47528  86.99222 2.775
42.90972  91.47000 1.154 40.96306  93.39250 1.423
37.65111  94.80361 1.863 39.10222  96.60917 0.898
38.67167 100.91639 0.536 37.70472  85.04889 2.693
37.07778  82.99361 2.195 38.11833  83.54694 2.762
36.79056  88.06722 2.377 29.92972  91.71528 1.276
30.81139  90.18083 1.393 44.37389  68.26056 2.268
46.86889  68.01472 1.551 44.10750  70.72889 1.631
45.48917  69.66528 1.369 39.40889  76.99528 2.535
38.91306  76.15250 2.477 41.97583  70.02472 1.619
42.39250  72.34444 2.156 42.38389  71.21472 2.417
45.56083  84.67833 1.701 46.37417  84.74139 1.539
47.10472  88.55139 1.048 42.41028  85.39278 3.107
44.22389  85.81806 2.258 47.53111  93.46861 0.550
47.94639  91.49611 0.563 46.24944  94.49722 0.591
44.23722  95.30056 0.604 32.30667  90.31833 1.614
32.33472  89.16583 1.135 34.00250  89.80000 1.503
38.75361  92.19889 1.814 36.91083  90.31861 2.435
45.56861 107.43750 0.217 48.51028 113.99583 0.387
48.49917 109.79750 0.100 46.48500 112.06472 0.209
41.15306  96.49278 0.743 41.05917 100.74639 0.391
36.13583 115.42556 0.139 41.28528 115.85222 0.075
38.79917 119.25667 0.053 39.00500 114.21583 0.273
43.94306  71.70333 2.391 40.31500  74.85472 2.593
33.22028 108.23472 0.377 35.78167 106.26750 0.315
32.90944 105.47056 0.355 36.04083 106.97139 0.376
36.77889 103.98139 0.326 42.73389  76.65972 3.249
42.29944  79.39639 3.344 43.97306  74.22306 2.322
44.39333  73.85944 2.111 41.35083  74.04861 3.306
43.52611  75.94722 3.948 42.10639  77.53583 2.231
41.99361  74.50361 3.022 36.13250  77.17139 1.857
35.06056  83.43056 2.393 35.69694  80.62250 2.082
35.02583  78.27833 1.729 34.97083  79.52833 1.959
35.72833  78.68028 1.780 47.60139 103.26417 0.354
48.78250  97.75417 0.306 47.12556  99.23694 0.273
39.53139  84.72417 3.828 40.35528  83.06611 3.401
39.79278  81.53111 3.961 40.78222  81.92000 3.349
36.80528  98.20056 0.603 34.98000  97.52139 0.994
36.59083 101.61750 0.444 44.38694 123.62306 0.629
44.63472 123.19000 0.329 45.44917 122.15333 0.716
43.12167 121.05778 0.050 44.21333 122.25333 0.423
43.89944 117.42694 0.071 45.22444 118.51139 0.109
40.78833  77.94583 3.275 41.59778  78.76750 4.336
40.65750  77.93972 3.352 41.32750  74.82028 3.081
33.53944  80.43500 1.456 44.35500  98.29083 0.372
43.94917 101.85833 0.224 35.96139  84.28722 3.579
35.18250  87.19639 2.148 35.66444  83.59028 2.474
35.46778  89.15861 1.811 33.95778 102.77611 0.376
28.46667  97.70694 0.886 29.66139  96.25944 0.934
30.26139 100.55500 0.938 32.37861  94.71167 2.229
31.56056  94.86083 1.472 33.27333  99.21528 0.890
33.39167  97.63972 1.585 37.61861 112.17278 0.237
41.65833 111.89694 0.271 38.99833 110.16528 0.143
41.35750 111.04861 0.172 42.87611  73.16333 2.412
44.52833  72.86889 2.549 38.04056  78.54306 2.478
37.33139  80.55750 1.650 38.52250  78.43583 2.360
47.86000 123.93194 1.144 48.54056 121.44528 0.837
46.83528 122.28667 0.635 46.76056 117.18472 0.255
37.98000  80.95000 2.396 39.08972  79.66222 3.291
45.79639  88.39944 1.054 45.05333  88.37278 1.457
44.66444  89.65222 1.044 43.70194  90.56861 1.309
46.05278  89.65306 1.132 42.57917  88.50056 1.809
45.82278  91.87444 0.984 41.34028 106.19083 0.335
42.73389 108.85000 0.236 42.49472 108.82917 0.313
42.92889 109.78667 0.182 43.22278 109.99111 0.161
43.87333 104.19222 0.306 44.91722 110.42028 0.210
45.07611  72.67556 2.646
;

title "Spatial Observations for SO4 Data";
proc sgplot data=so4;
   scatter x=Longitude y=Latitude;
   xaxis grid values=(-120 to -70 by 10) min=-125 max=-65;
   yaxis grid values=(25 to 50 by 5) min=25 max=50;
run;

proc iml;
load module=(NearestNbr RandRowsAndNbrs SMOTESimContin);

varNames = {'Longitude' 'Latitude'};
use so4;
read all var varNames into X;
close;

/* don't use k too large or you'll get data in Gulf of Mexico, Atlantic ocean,
   and Mexico! */
call randseed(321,1);
k = 30;
T = 3*179;
p_new = SMOTESimContin(X, T, k);

create synthetic from p_new[c=varNames];
append from p_new;
close;

submit T k varNames;
title "&T Synthetic Data Points by Using &k Nearest Neighbors";
proc sgplot data=synthetic;
scatter x=Longitude y=Latitude;
xaxis grid values=(-120 to -70 by 10) min=-125 max=-65;
yaxis grid values=(25 to 50 by 5) min=25 max=50;
run;
endsubmit;
QUIT;
