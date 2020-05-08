/* SAS program to accompany the articles
   "What is a geometric mean?"
   by Rick Wicklin, published 30SEP2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/09/30/what-is-a-geometric-mean.html
   and 
   "Compute the geometric mean, geometric standard deviation, and geometric CV in SAS"
   by Rick Wicklin, published 02OCT2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/10/02/geometric-mean-deviation-cv-sas.html ?

   This program shows how to compute the 
   geometric mean (GM),
   geometric standard deviation (GSD), and
   geometric coefficient of variation (GCV) in SAS.

   Data for the yearly return of Gold is from 
   https://www.macrotrends.net/1333/historical-gold-prices-100-year-chart
*/

/*
Yearly return for gold, 2010 - 2018
https://www.macrotrends.net/1333/historical-gold-prices-100-year-chart
*/
/**********************************************/
/* First blog post: What is the geometric mean? */
/**********************************************/

data Gold;
input Year r;
x = 1 + r;
datalines;
2010  0.2774
2011  0.1165
2012  0.0568
2013 -0.2779
2014 -0.0019
2015 -0.1159
2016  0.0863
2017  0.1257
2018 -0.0115
;

data Investment;
retain Value 1000;
set Gold;
Value = Value*x;
run;

proc print data=Investment noobs split='/'; 
label Year = "/Year" r="/r" x = "/1 + r"
      Value = "Value of/Investment";
   var Year r x Value; 
run;

/* several ways to compute the geometric mean. Here
   we use the GEOMMEAN function in Base SAS */
proc transpose data=Investment out=Invest;
   var x;
run;

data Avgx;
set Invest(keep=COL:);
array x[*] COL:;
GM = geomean(of x[*]);
run;

proc print data=Avgx noobs; var GM; run;

/* another alternative is the GEOMEAN function in SAS/IML */
proc iml;
use Gold; read all var "x"; close;
GM = geomean(x);
print GM;
QUIT;

/* the geometric mean is known in finance as the 
   Compound Annual Growth x (CAGR). It tells you 
   the x of return of an investment, assuming
   a constant growth x
*/
data I0;  /* initial investment */
year = 2009; Value = 1000; CD = 1000; output;
run;
data I1;  /* fixed-rate annual compounding at 1.672% */
retain CD 1000;
set Investment;
CD = CD*1.01672;
run;
data I2;  /* combine the data sets and plot */
set I0 I1;
run;
ods graphics / width=640 height=320;
%let Gold = CXDAA520;
title "End-of-Year Value for $1000 Investment";
proc sgplot data=I2;
series x=Year y=CD / curvelabel="CD" curvelabelpos=min markers markers markerattrs=(symbol=CircleFilled);;
series x=Year y=Value / lineattrs=(color=&Gold) curvelabel="Gold" curvelabelpos=min 
       markers markerattrs=(color=&Gold symbol=SquareFilled);
xaxis grid display=(nolabel); 
yaxis grid label="Value" values=(1000 to 1500 by 100) valueshint;
run;



/**********************************************/
/* Second blog post: Compute the geometric mean, 
   geometric standard deviation, and geometric CV in SAS */
/**********************************************/

%let mu = 3;
%let sigma = 0.8;
%let N = 100;

data Have;
call streaminit(12345);
do i = 1 to &N;
   x = round( rand("LogNormal", &mu, &sigma), 0.1);
   logX = log(x);
   output;
end;
run;

proc means data=Have;
var x logx;
run;

title "Geometric Mean of Skewed Positive Data";
proc sgplot data=Have;
histogram x / binwidth=10 binstart=5 showbins;
refline 20.2 / axis=x label="Geometric/Mean" splitchar="/" labelloc=inside lineattrs=GraphData2(thickness=3);
xaxis values=(0 to 140 by 10);
yaxis offsetmax=0.1;
run;
/*
proc univariate data=Have;
var x;
histogram x / endpoints=(0 to 140 by 10) odstitle=title
              kernel lognormal(scale=&mu, shape=&sigma);
run;
*/
/* arithmetic means of log(X). You can exponentiate these values */
proc ttest data=Have dist=normal; 
var logX;
ods select ConfLimits;
run;
/* or directly compute the geometric statistics for X */
proc ttest data=Have dist=lognormal; 
   var x;
   ods select ConfLimits;
run;

/* SURVEYMEANS can do it, but definition of std err is different 
   for survey procs */
/*
proc surveymeans data=Have allgeo;
   var x;
   ods select GeometricMeans;
run;
*/

/* Note: The SAS documentation for the TTEST procedure says
   "For lognormal data, the CV is the natural measure of variability 
   (rather than the standard deviation) because the CV is 
    invariant to multiplication of a lognormal variable by a constant."
*/
/* compute geometric statistics in SAS/IML */
proc iml;
use Have; read all var "x"; close;  /* read in positive data */
GM = geomean(x);               /* built-in GEOMEAN function */
print GM;

/* To estimate the geometric mean and geometric StdDev, compute
   arithmetic estimates of log(X), then EXP transform the results. */
n = nrow(x);
z = log(x);                  /* log-transformed data */
m = mean(z);                 /* arithmetic mean of log(X) */
s = std(z);                  /* arithmetic std dev of log(X) */
GM2 = exp(m);                /* same answer as GEOMEAN function */
GSD = exp(s);                /* geometric std dev */
GCV = sqrt(exp(s**2) - 1);   /* geometric CV */
print GM2 GSD GCV;

/* verify the 64-95 rule */
Upper1 = GM * GSD;
Lower1 = GM / GSD;
within1 = mean(x<=Upper1 & x>=Lower1);
print Lower1 Upper1 within1[F=percent7.2];
Upper2 = GM * (GSD)**2;
Lower2 = GM / (GSD)**2;
within2 = mean(x<=Upper2 & x>=Lower2);
print Lower2 Upper2 within2[F=percent7.2];


/* Note that some researchers use the following formula */
y = log(x/GM)##2;
GSD2 = exp( sqrt(mean(y)) );
*print (GSD // GSD2)[F=16.12];

/* The difference between the two estimates is whether the 
   denominator for StdDev uses Bessel's correction: n vs (n-1) */
GSD3 = exp( sqrt(mean(y)*n/(n-1)) );
*print (GSD // GSD3)[F=16.12];

/* compute CI for mean */
alpha = 0.05;
tCrit = quantile("T", 1 - alpha/2, n-1 );  /* sample mean is t distributed */
SEM = s / sqrt(n);
CImean = m + {-1  1} * tCrit*SEM;          /* mean +/- t_crit * SEM */
CIGM= exp(CImean);                         /* CI for GM = EXP the endpoints for CI mean */
*print SEM (exp(SEM));

/* compute CI for StdDev */
p = (1 - alpha/2) || (alpha/2);
chi2Crit = quantile("chisquare", p, n-1 ); /* variance is chi-sq distributed */
numer = (n-1)*s**2;
CIVar = numer/chi2Crit;                    /* CI for variance */
CISD = sqrt( CIVar );                      /* CI for Std Dev */
CIGSD = exp(CISD);                         /* CI for GSD = EXP endpoints for CI Std Dev */
CIGCV = sqrt(exp(CIVar) - 1);              /* CI for GCV */

/* put arithmetric and geometric stats together */
A = (m || CIMean) // 
    (s || CISD);
print A[r={"Mean" "SD"} c={"Estimate" "Lower" "Upper"} L="Confidence Intervals"];
G = (GM  || CIGM) // 
    (GSD || CIGSD) // 
    (GCV || CIGCV);
print G[r={"GM" "GSD" "GCV"} c={"Estimate" "Lower" "Upper"} L="Confidence Intervals"];
QUIT;

/***********************************************************/
/***********************************************************/
/* Define SAS/IML function that computes geometric 
   mean, stdDev, CV, and CIs in one easy call:
   call GeoStats(x, <alpha=00.05>)), where all(x > 0)
   Run the following program. Then load the module to use it.
/***********************************************************/
/***********************************************************/
proc iml;
start GeoStats(x, alpha=0.05);
if any(x<=0) then stop "ERROR: The data must contains only positive values";

/* To estimate the geometric mean and geometric StdDev, compute
   arithmetic estimates of log(X), then EXP transform the results. */
n = nrow(x);                 /* no missing values (b/c missing < 0) */
z = log(x);                  /* log-transformed data */
m = mean(z);                 /* arithmetic mean of log(X) */
s = std(z);                  /* arithmetic std dev of log(X) */
GM = exp(m);                 /* same answer as GEOMEAN function */
GSD = exp(s);                /* geometric std dev */
GCV = sqrt(exp(s**2) - 1);   /* geometric CV */

/* Note that some researchers use the following formula */
y = log(x/GM)##2;
GSD2 = exp( sqrt(mean(y)) );
/* The difference between the two estimates is whether the 
   denominator for StdDev uses Bessel's correction: n vs (n-1) */
GSD3 = exp( sqrt(mean(y)*n/(n-1)) );

/* compute CI for mean */
tCrit = quantile("T", 1-alpha/2, n-1 );  /* sample mean is t distributed */
SEM = s / sqrt(n);
CImean = m + {-1  1} * tCrit*SEM;        /* mean +/- t_crit * SEM */
CIGM= exp(CImean);                       /* CI for GM = EXP the endpoints for CI mean */

/* compute CI for StdDev */
p = (1 - alpha/2) || (alpha/2);
chi2Crit = quantile("chisquare", p, n-1 ); /* variance is chi-sq distributed */
numer = (n-1)*s**2;
CIVar = numer/chi2Crit;                  /* CI for variance */
CISD = sqrt( CIVar );                    /* CI for Std Dev */
CIGSD = exp(CISD);                       /* CI for GSD = EXP endpoints for CI Std Dev */
CIGCV = sqrt(exp(CIVar) - 1);            /* CI for GCV */

/* put geometric stats together */
G = (GM  || CIGM) // 
    (GSD || CIGSD) // 
    (GCV || CIGCV);
RETURN( G );
finish;
store module=GeoStats;
QUIT;


/* Test the function */
data Have;
input x @@;
datalines;
24.8 47.5 38.6 12.9 68.9  7.5 17.9 46.2 21.2 53.5
17.8 25.3  8.0 13.7 14.0 25.0 16.8 23.9 14.4 24.5 
41.5 65.9 76.0  5.2 10.9 33.2 63.2 26.8 19.7 29.2
36.4 12.9 59.9 21.5  7.0  4.3 32.5 38.6 46.3  9.0
47.8 13.5  7.9 13.9 42.8  4.1 14.2 20.9  8.7 37.8
27.4 10.8 18.9 61.7 37.3 58.1 26.3 18.6 50.0  7.8
26.5  9.0 22.6 34.9 12.5  8.8 19.6 24.4 12.4  4.3
21.0  7.2 37.1 58.5 10.0 18.8 21.0  6.6 35.7 28.7
23.0 52.1 23.8 19.3 10.9 30.9 34.7 11.5  7.3  4.3
137.4 2.8 14.7 30.1 11.7 45.0 11.0 20.4 25.5 15.0 
;

proc iml;
load module=GeoStats;
use Have; read all var "x"; close;  /* read in positive data */

G = geoStats(x);
print G[c={"Estimate" "Lower" "Upper"} L="Confidence Intervals"
        r={"GeoMean" "GeoSD" "GeoCV"} ];
QUIT;
