/* SAS program to accompany the article 
   "Avoid alphabetical order for a categorical axes in a statistical graph"
   by Rick Wicklin, published 13NOV2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/11/13/avoid-alphabetical-order.html

   This program shows why you should avoid using alphabetical ordering for bar charts, or 
   "Why 'Alabama first' is not always the best way to display data"
   Avoiding 'Alabama first' was often advocated by Howard Wainer.

   For a different example see how sorting affects lasagna plots:
   https://blogs.sas.com/content/iml/2016/06/08/lasagna-plot-in-sas.html

   Data from
   https://www.census.gov/data/tables/time-series/demo/popest/2020s-state-total.html
   https://en.wikipedia.org/wiki/List_of_U.S._states_and_territories_by_area

   Lessons: 
   1. For many categories, use tall graphs an horizontal bar charts.
   2. Set FITPOLICY=NONE, then increase graph height and/or decrease label font.
   3. Almost always, use the CATEGORYORDER= option to display the categories in a meaningful order.
   4. If the graph type is not a traditional bar chart, use the COLORBANDS= option to 
      guide the viewer's eyes to the data.
   5. For other category orders, sort data and use CATEGORYORDER=DATA on XAXIS or YAXIS stmts.
*/

data StatePop;
infile datalines delimiter = ","	missover dsd;
length StateName $20;
input StateID StateName Pop2020 PopChange2020 LandArea;
PopM = Pop2020 / 1E6;
label PopM = "Population (Millions)"
      PopChange2020 = "Population Change (2020)"
      LandArea = "Land Area (km^2)";
format PopChange2020 comma7.;
datalines;
1, Alabama,         5024356, 7006   , 135767
2, Alaska,           733378, -455   ,1723337
4, Arizona,         7151507, 28436  , 295234
5, Arkansas,        3011555, 2640   , 137732
6, California,     39538245, -36592 , 423967
8, Colorado,        5773733, 11132  , 269601
9, Connecticut,     3605942, -8580  ,  14357
10, Delaware,        989957, 2157   ,   6446
12, Florida,       21538226, 51376  , 170312
13, Georgia,       10711937, 17891  , 153910
15, Hawaii,         1455273, -4230  ,  28313
16, Idaho,          1839092, 10110  , 216443
17, Illinois,      12812545, -25965 , 149995
18, Indiana,        6785668, 3131   ,  94326
19, Iowa,           3190372, 199    , 145746
20, Kansas,         2937847, 72     , 213100
21, Kentucky,       4505893, 1552   , 104656
22, Louisiana,      4657749, -6085  , 135659
23, Maine,          1362341, 1216   ,  91633
24, Maryland,       6177213, -4008  ,  32131
25, Massachusetts,  7029949, -34220 ,  27336
26, Michigan,      10077325, -7748  , 250487
27, Minnesota,      5706504, 3348   , 225163
28, Mississippi,    2961288, -3147  , 125438
29, Missouri,       6154920, -922   , 180540
30, Montana,        1084197, 2878   , 380831
31, Nebraska,       1961489, 1153   , 200330
32, Nevada,         3104624, 11024  , 286380
33, New Hampshire,  1377518, 1069   ,  24214
34, New Jersey,     9289031, -17342 ,  22591
35, New Mexico,     2117527, 863    , 314917
36, New York,      20201230, -92934 , 141297
37, North Carolina,10439414, 10031  , 139391
38, North Dakota,    779091, 427    , 183108
39, Ohio,          11799374, -1857  , 116098
40, Oklahoma,       3959346, 5566   , 181037
41, Oregon,         4237291, 7504   , 254799
42, Pennsylvania,  13002689, -8249  , 119280
44, Rhode Island,   1097371, -1026  ,   4001
45, South Carolina, 5118429, 13419  ,  82933
46, South Dakota,    886677, 1122   , 199729
47, Tennessee,      6910786, 14833  , 109153
48, Texas,         29145428, 87046  , 695662
49, Utah,           3271614, 12171  , 219882
50, Vermont,         643085, -192   ,  24906
51, Virginia,       8631384, 5087   , 110787
53, Washington,     7705247, 18784  , 184661
54, West Virginia,  1793755, -2335  ,  62756
55, Wisconsin,      5893725, 2546   , 169635
56, Wyoming,         576837, 768    , 253335
;
ods graphics/reset;

/* allow 12-15 px per category, plus top and bottom matter */
ods graphics / push width=480px height=720px;
title "State Populations (2020)";
title2 "Ordered by State Name (Alabama First)";
proc sgplot data=StatePop;
   hbar StateName / response=PopM;
   yaxis display=(nolabel) valueattrs=(size=8) /* small font for values */
                           fitpolicy=none;     /* do not thin labels */
   xaxis grid;
run;

title2 "Ordered by Population";
proc sgplot data=StatePop;
   hbar StateName / response=PopM categoryorder=respdesc; /* or respasc */
   yaxis display=(nolabel) valueattrs=(size=8) /* small font for values */
                           fitpolicy=none;     /* do not thin labels */
   xaxis grid;
run;

title "Change in State Populations (2020)";
proc sgplot data=StatePop;
   hbar StateName / response=PopChange2020 categoryorder=respasc; /* or respdesc */
   yaxis display=(nolabel) colorbands=even     /* faint alternating bands */
                           valueattrs=(size=8) /* small font for values */
                           fitpolicy=none;     /* do not thin labels */
   xaxis grid;
run;

ods graphics / pop;

/* Sort by any variable and display by using DISCRETEORDER=DATA */
proc sort data=StatePop out=State2;
   by PopM;
run;

ods graphics / push width=480px height=640px;
title "State Populations (2020)";
title2 "States Ordered by Population Size";
proc sgplot data=State2;
   hbar StateName / response=PopChange2020;
   yaxis display=(nolabel) discreteorder=data  /* display in data order */
                           colorbands=even     /* faint alternating bands */
                           valueattrs=(size=8) /* small font for values */
                           fitpolicy=none;     /* do not thin labels */
   xaxis grid;
run;
ods graphics / pop;

/* Note: It is often better to use a scatter plot if you want to directly compare two variables */
title "Population vs Land Area";
proc sgplot data=StatePop;
   scatter x=PopChange2020 y=PopM / datalabel=StateName;
   xaxis grid;
   yaxis grid;
run;
