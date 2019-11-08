/* SAS program to accompany the article 
   "Round to even"
   by Rick Wicklin, published 11NOV2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/11/11/round-to-even.html

   This program shows how to use the round-half-to-even method in SAS.
   Data from the documentation for PROC UNIVARIATE
   "Creating a One-Way Comparative Histogram"
   https://bit.ly/36Byroe

   "The effective channel length (in microns) is measured for 1225 field effect 
    transistors."  This data uses only the subset where Lot = 'Lot 1'
*/

/* Create a step plot of the round-half-to-even function */
data Series;
do x = 0 to 6 by 0.05;
   y = rounde(x); 
   group = y;
   output;
end;

data Scatter1;   /* The closed endpoints of the intervale. Use a filled marker */
do x1 = 0.5 to 5.5 by 1;
   y1 = rounde(x1); 
   output;
end;

data Scatter2;   /* The open endpoints of the intervale. Use an unfilled marker */
do x2 = 0.5 to 5.5 by 1;
   y2 = rounde(x2+1) - 1; 
   output;
end;
data RoundToEven;
set Series Scatter1 Scatter2;
run;

Title "The Round-To-Even Function";
proc sgplot data=RoundToEven noautolegend;
   series x=x y=y / group=group lineattrs=GraphFit;
   scatter x=x1 y=y1 / markerattrs=(symbol=circlefilled) FILLEDOUTLINEDMARKERS
 markerfillattrs=(color=cx445694);
   scatter x=x2 y=y2 / markerattrs=(symbol=circlefilled) FILLEDOUTLINEDMARKERS
 markerfillattrs=(color=white);
   yaxis grid label='RoundE(x)';
   xaxis grid label='x';
run;

/* Create a table that compares the ROUND and RONDE (round-to-even) functions */
data Seq;
keep x Round RoundEven;
label Round = "Round(x)" RoundEven="RoundE(x)";
do x = 0 to 3.5 by 0.25;
   Round     = round(x);    /* traditional: half-integers rounded away from 0 */
   RoundEven = rounde(x);   /* round half-integers to the nearest even integer */
   output;
end;
run;

proc print data=Seq noobs label;
run;

/* Data from the documentation for PROC UNIVARIATE
   "Creating a One-Way Comparative Histogram"
   https://bit.ly/36Byroe
*/
/* Lengths of circuit boards for 'Lot 1'. Data from PROC UNIVARIATE example */
data Channel1; 
label Length = 'Length (microns)'; 
input Length @@; 
datalines; 
0.91 1.01 0.95 1.13 1.12 0.86 0.96 1.17 1.36 1.10 
0.98 1.27 1.13 0.92 1.15 1.26 1.14 0.88 1.03 1.00 
0.98 0.94 1.09 0.92 1.10 0.95 1.05 1.05 1.11 1.15 
1.11 0.98 0.78 1.09 0.94 1.05 0.89 1.16 0.88 1.19 
1.01 1.08 1.19 0.94 0.92 1.27 0.90 0.88 1.38 1.02 
1.21 0.99 1.19 1.05 0.95 0.92 1.06 1.03 0.76 1.04 
0.84 0.93 0.96 1.23 1.06 1.20 0.94 0.74 0.80 1.00 
1.00 0.97 0.94 1.05 1.09 1.07 0.79 1.06 0.96 1.00 
1.07 0.98 1.26 0.99 1.12 0.99 1.05 1.12 0.85 0.96 
1.01 1.21 0.95 0.86 0.86 1.19 1.24 0.75 1.03 0.84 
0.93 0.65 0.84 0.92 0.86 1.07 1.05 0.95 1.05 1.10 
0.93 1.03 0.94 0.97 0.95 0.99 0.97 0.81 1.10 0.95 
1.08 0.75 0.97 0.94 0.95 0.85 1.02 1.07 1.11 1.02 
1.16 1.05 0.88 1.29 0.99 0.83 1.06 1.05 1.03 1.15 
1.14 0.91 1.17 1.17 0.70 0.84 1.33 0.95 1.34 1.14 
1.08 1.28 1.14 0.84 1.09 1.12 0.99 1.08 0.95 1.12 
0.70 1.07 1.06 1.24 1.17 0.96 1.00 0.85 1.25 0.84 
1.00 1.01 1.08 1.08 0.75 1.06 1.09 0.86 1.08 1.11 
1.19 1.25 0.94 0.97 0.81 0.91 0.99 0.94 0.94 1.12 
1.10 0.82 1.17 0.89 0.75 0.84 1.10 0.93 1.11 1.07 
1.20 1.22 0.90 1.01 1.17 0.97 1.02 0.77 1.07 0.87 
0.96 1.21 1.04 0.89 0.95 0.86 0.84 1.03 1.26 0.71 
0.75 1.34 1.14 1.09 0.85 0.68 0.79 1.06 1.12 1.10 
0.98 1.11 0.98 0.86 1.15 0.89 1.21 0.79 1.31 0.84 
1.05 0.89 1.01 0.81 1.10 1.04 1.07 1.02 0.91 1.00 
0.98 0.95 0.82 1.11 1.07 1.07 0.89 1.15 1.04 1.19 
1.07 1.03 1.04 0.81 0.93 1.06 0.92 1.29 1.04 0.99 
0.86 0.84 1.22 1.09 1.14 0.76 0.86 0.96 1.07 1.21 
0.93 1.07 1.13 0.95 0.92 1.01 1.01 1.16 1.24 0.82 
1.03 1.31 0.79 1.14 1.15 0.97 0.82 1.10 1.06 1.00 
0.58 0.79 1.08 1.06 1.06 1.14 1.09 1.01 0.87 0.97 
1.08 1.03 1.08 0.92 1.26 0.74 0.98 0.94 0.84 1.07 
1.04 0.93 0.90 0.83 0.91 1.15 1.00 0.86 0.87 0.75 
1.19 1.13 1.14 0.78 1.03 0.95 0.77 0.89 0.90 0.96 
0.95 1.08 1.03 0.98 1.00 0.86 1.16 1.05 0.97 0.91 
1.06 1.12 0.97 1.03 0.88 1.02 0.98 0.95 0.86 1.14 
1.15 1.26 0.84 0.82 0.91 1.02 0.96 0.95 1.17 1.03 
0.89 1.07 1.06 1.09 0.95 0.99 1.38 1.10 0.88 1.05 
1.18 0.93 1.14 0.86 1.24 0.97 1.04 1.12 1.10 1.14 
1.22 0.89 0.93 1.30 0.93 1.00 1.06 0.92 1.26 0.75 
1.24 1.17 1.09 1.06 0.78 0.95 0.97 0.90 0.98 0.88 
0.87 1.00 0.83 1.11 1.15 1.20 1.21 1.00 0.95 0.82 
1.16 0.80 0.99 0.89 1.17 
;

/* round real data to the nearest 0.1 unit */
data Round;
set Channel1;
   Round     = round(Length, 0.1);    /* traditional: round to nearest tenth */
   RoundEven = rounde(Length, 0.1);   /* use round-to-even method to round to nearest thenth */
   /* create a binary indicator variable: Was x rounded up or down? */
   RoundUp     = (Round > Length);    /* 1 if rounded up; 0 if rounded down */
   RoundEvenUp = (RoundEven > Length);
run;

/* compare estimates that use full data and rounded data */
proc means data=Round sum mean ndec=3;
   label Length="True values" Round ="Rounded values" RoundEven="Round-to-even values";
   var Length Round RoundEven;
run;

/* how many obserations were rounded up or down? */
proc means data=Round mean ndec=3;
   label RoundUp    = "Proportion rounded up for ROUND"
         RoundEvenUp= "Proportion rounded up for ROUNDE";
   var RoundUp RoundEvenUp;
run;
