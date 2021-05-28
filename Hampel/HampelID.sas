/* SAS program to accompany the article 
   "The Hampel identifier: Robust outlier detection in a time series"
   by Rick Wicklin, published 01JUN2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/06/01/hampel-filter-robust-outliers.html

   This program shows how to implement the Hampel identifier and Hampel filter in SAS.
*/


/* Daily temperature of a cow at 6:30 a.m. measured by counting 
   chirps from a telemetric thermometer implanted in the cow. 
   y = (chirps per 5-minute interval) - 800. 
   Peak temperatures are associated with fertility.
   Cows peak about every 15-20 days.
   Source: Velleman and Hoaglin (1981). */
data Cows;
input y @@;
Day + 1;
datalines;
60 70 54 56 70 66 53 95 70 69 56 70 70 60 60 60 50 50 48 59
50 60 70 54 46 57 57 51 51 59 42 46 40 40 54 47 67 50 60 54
55 50 55 54 47 48 54 42 43 62 49 41 45 40 49 46 54 54 60 58
52 47 53 39 55 45 47 41 48 42 45 48 52 49 53
;

/* visualize raw data and runing median */
proc expand data=Cows out=MedOut method=NONE;
   id Day;
   convert y = RunMed / transformout=(cmovmed 7);
run;

ods graphics / height=400px width=600px;

title "Overlay Running Median";
title2 "Window Length = 7";
proc sgplot data=MedOut noautolegend;
   scatter x=Day y=y;
   series x=Day y=RunMed;
   xaxis grid values=(0 to 80 by 10);
   yaxis grid label="Chirps (Body Temperature)";
run;


/* simple example data: Sine curve and four outliers */
data Test;
pi = constant('pi');
do t = 1 to 30;
   y = sin(2*pi*t/30);
   if t in (3 12 13 24) then 
      y = 5;
   output;
end;
run;

/* visualize the classical moving means and standard deviations */
proc expand data=Test out=Classical method=NONE;
   id t;
   convert y = RunMean / transformout=(cmovave 7);
   convert y = RunStd / transformout=(cmovstd 7);
run;

data Classical;
set Classical;
lower = RunMean - 3*RunStd;
upper = RunMean + 3*RunStd;
run;

title "Moving Average and Standard Deviation";
proc sgplot data=Classical;
   band x=t lower=lower upper=upper / name="b" legendlabel="Classical Interval";
   scatter x=t y=y / name="data";
   series x=t y=RunMean / name="c" legendlabel="Running Mean(7)";
   keylegend "b" "c";
run;

/* ----------------------------------------- */
/* extract moving windows, which makes rolling statistics easy */
proc iml;
/* extract the symmetric sequence about x[i], which is x[i-k-1:i+k+1] if
   1 <= i-k-1  and  i+k+1 <= nrow(x).
   Otherwise, pad the sequence. You can pad with a user-defined value 
   (such as missing or zero) or you can repeat the first or last value.
   This function loops over the number of columns and extracts each window.
   INPUT: y is the time series
          k is the radius of the window. The window length is 2*k+1 centered at y[i]
          padVal can be missing or a character value. If character,
                 use y[1] to pad on the left and y[n] to pad on the right.

   See https://blogs.sas.com/content/iml/2021/05/26/running-median-smoother.html
*/
start extractWindows(y, k, padVal="REPEAT");
   n = nrow(y);
   M = j(2*k+1, n, .);
   /* create new series, z, by padding k values to the ends of y */
   if type(padVal)='N' then
      z = repeat(padVal, k) // y // repeat(padVal, k);
   else
      z = repeat(y[1], k) // y // repeat(y[n], k); /* Tukey's padding */
   do i = 1 to n;
      range = i : i+2*k;   /* centered at k+i */
      M[,i] = z[range];
   end;
   return M;
finish;

/* Hampel identifier: Use moving median and MAD to identify outliers.
   See https://blogs.sas.com/content/iml/2021/06/01/hampel-filter-robust-outliers.html
*/
start HampelID(y, k, multiplier=3);
   M = extractWindows(y, k, "REPEAT");
   med = T( median(M) );
   Nmad = T( mad(M, "NMAD") ); 
   idx = loc( abs(y-med) > multiplier*Nmad );
   return idx;
finish;

store module=(extractWindows HampelID);
QUIT;
/* ----------------------------------------- */

/* Test the extractWindows function. See 
   https://blogs.sas.com/content/iml/2021/05/26/running-median-smoother.html */
proc iml;
load module = extractWindows;

/* example data, just to see how the function works */
t = 1:10;
y = T(20:29);
M = extractWindows(y, 3, .);
rowlabls = {'x_{t-3}' 'x_{t-2}' 'x_{t-1}' 'x_t' 'x_{t+1}' 'x_{t+2}' 'x_{t+3}'};
/* numeric colnames supported in SAS/IML 15.1 */
print M[c=t r=rowlabls L="Window of Length=7"];

M = extractWindows(y, 3);     /* Tukey padding: Use y[1] to pad the left and y[n] to pad the right */
print M[c=t r=rowlabls L="Window of Length=7"];

/* now use the function to compare the running median and running trimmed mean */
use Cows; read all var {Day y}; close;

/* median window that contains 7 points */
k = 3;
M = extractWindows(y, k);     /* each column is a window of length 2*k+1 = 7 */
median = median(M);           /* compute the moving median */
mean1 = mean(M, "trim", 1);   /* compute the moving trimmed mean */

create Cows2 var{"Day", "Y", "Median", "Mean1"}; append; close;
QUIT;

title "Overlay Running Median and Trimmed Mean";
title2 "Window Length = 7";
proc sgplot data=Cows2;
   scatter x=Day y=y;
   series x=Day y=Median / legendlabel="Running Median(7)"
                    lineattrs=(thickness=2 color=DarkGreen);
   series x=Day y=Mean1 / legendlabel="Running Trimmed Mean(7)"
                    lineattrs=(thickness=2 color=DarkRed);
   xaxis grid;
   yaxis grid label="Chirps (Body Temperature)";
run;


/*************************************************************/
/* APPENDIX: Discussion of the MAD and 1.4826*MAD statistics */
/* For normally distributed data, 1.4826*MAD is a consistent estimator 
   for the standard deviation. */
proc univariate data=Test robustscale;
  var y;
  ods select robustscale;
run;

proc iml;
call randseed(4321);
x = randfun(10000, "Normal");
sd = std(x);        /* should be close to 1 */
Nmad = mad(x, "NMAD");  /* = 1.4826*mad, which should be close to 1 */
print sd Nmad;
QUIT;
/*************************************************************/


/* Implement the Hampel identifier by using the rolling median and MAD */
proc iml;
/* define or STORE/LOAD the extractWindows function. See
   https://blogs.sas.com/content/iml/2021/05/26/running-median-smoother.html */
load module = extractWindows;

use Test; read all var {'t' 'y'}; close;

M = extractWindows(y, 3);   /* k=3 ==> window length=7 */
med = T( median(M) );       /* column vector of medians */
Nmad = T( mad(M, "NMAD") ); /* column vector of 1.4826*MADs */
/* Hampel identifier: Which obs are outside med +/- 3*NMAD? */
lower = med - 3*Nmad;
upper = med + 3*Nmad;
outlier = ( abs(y-med) > 3*Nmad );

/* Hampel filter: Replace outliers with the medians */
idx = loc(outlier);
newY = y;
if ncol(idx)>0 then 
   newY[idx] = med[idx];

/* visualize the intervals in the Hamepl filter */
create Robust var {'t' 'y' 'med' 'lower' 'upper' 'outlier'};
append;
close;
QUIT;

title "Robust Moving Estimates";
title2 "Hampel Filter Using Median and MAD";
proc sgplot data=Robust;
   band x=t lower=lower upper=upper / name="b" legendlabel="Robust Interval";
   scatter x=t y=y  / name="group" group=outlier;
   series x=t y=med / name="c" legendlabel="Running Median(7)";
   keylegend "b" "c";
run;

/* apply the Hampel identifier to the Cows data */
/* for convenience, use the HampelID function */
proc iml;
/* define or STORE/LOAD the extractWindows function. */
load module = (extractWindows HampelID);

use Cows; read all var {Day y}; close;

idx = HampelID(y, 3);
print idx;
outlier = j(nrow(y), 1, 0); /* most obs are not outliers */
outlier[idx] = 1;           /* these were flagged as potential outliers */

/* visualize the filter's upper and lower limits */
/* robust analysis */
M = extractWindows(y, 3);   /* k=3 ==> window length=7 */
med = T( median(M) );
Nmad = T( mad(M, "NMAD") ); 
lower = med - 3*Nmad;
upper = med + 3*Nmad;
outlier = ( abs(y-med) > 3*Nmad );
create RobustCows var {'Day' 'y' 'med' 'lower' 'upper' 'outlier'};
append;
close;
QUIT;

title "Robust Moving Estimates";
proc sgplot data=RobustCows;
   band x=Day lower=lower upper=upper / name="b" legendlabel="Robust Interval";
   scatter x=Day y=y  / name="group" group=outlier;
   series x=Day y=med / name="c" legendlabel="Running Median(7)";
   keylegend "b" "c";
   xaxis grid values=(0 to 80 by 10) valueshint;
   yaxis grid label="Chirps (Body Temperature)";
run;


