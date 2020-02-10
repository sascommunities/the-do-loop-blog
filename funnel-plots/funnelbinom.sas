/* SAS program to accompany the article 
   "Funnel plots for proportions"
   by Rick Wicklin, published 23NOV2011 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2011/11/23/funnel-plots-for-proportions.html
*/

/* Edit this step to point to the location of the data */
proc import out=Adoptions 
datafile="C:\Users\<userid>\Documents\My SAS Files\Adoptions.csv" 
     dbms=csv replace;
     getnames=yes;
     datarow=2; 
run;

proc iml;
use Adoptions;
read all var {Events Trials};
close Adoptions;

/* 1. Compute observed proportions */
Proportion = Events/ Trials;

/* 2. compute overall proportion */
theta = sum(Events) / sum(Trials); 
Expected = theta * Trials; 

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
/* plot limits at equally spaced points between min & max */
minN = min(Trials); maxN = max(Trials);
n = T( do(minN, maxN, (maxN-minN)/20) );
n = round(n); /* binomial parameter must be integer */
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

/* merge data with proportions */
data Adoptions; 
merge Adoptions Stats; 
label Proportion="Proportion Adopted";
label Trials="Number of Cases";
run;
/* append control limits */
data Funnel; set Adoptions Limits; run;


/* 4. Plot proportions versus sample size. Overlay control limits */
proc sgplot data=Funnel;
title "Proportion of Children Adopted within 12 Months";
band x=N lower=L3sd upper=U3sd / 
     nofill lineattrs=(color=lipk)
/*   fill fillattrs=GraphConfidence2 */
/*   outline lineattrs=GraphPrediction */
     legendlabel="99.8% limits" name="band99";
band x=N lower=L2sd upper=U2sd /
     nofill lineattrs=(color=gray)  
/*   fill fillattrs=GraphConfidence*/
/*   outline lineattrs=GraphPrediction */
     legendlabel="95% limits" name="band95";
refline &AvgProp / axis=y legendlabel="Overall proportion" name="ref";
scatter x=Trials y=Proportion / legendlabel="Local authority" name="scat";
keylegend "scat" "band95" "band99" / location=inside
     across=1 position=bottomright;
yaxis grid values=(0.3 to 1 by 0.1); xaxis grid;
run;

