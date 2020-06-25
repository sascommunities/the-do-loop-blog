/* SAS program to accompany the article 
   "What is a pooled variance?"
   by Rick Wicklin, published 29JUN2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/06/29/pooled-variance.html ‎

   This program shows how to compute the pooled variance
   between groups. If you assume is that the variance is 
   constant across different groups in the data,
   the pooled variance is an estimate of the common variance.
   It is a weighted average of the sample variances,
   where larger groups are weighted more heavily than smaller groups.

   Data modified from the 'AirPoll' example in the PROC UNIVARIATE doc:
   "Example 4.5 Creating Basic Summary Plots"
   Group = Site, Response = 'Ozone level (in ppb)'
*/
data Have(keep = Group Y);
label Y = 'Response' n = 'Sample Size';
length Group $6;
do i = 1 to 3;
   input Group n @@;
   do j = 1 to n;
      input Y @@; output;
   end;
end;
datalines;
A 22   4 6 3 4 7 8 2 3 4 1 3 8 9 5 6 4 6 3 4 7 3 7
B 15   5 3 6 2 1 2 4 3 2 4 6 4 6 3 6
C 18   8 9 7 8 6 7 6 7 9 8 9 8 7 8 5 8 9 7
;

/* visualize the response for each group */
title "Response by Group";
proc sgplot data=Have;
   vbox Y / category=Group;
   yaxis grid;
run;
/*
proc sgpanel data=Have;
   panelby group / columns=1 onepanel;
   histogram Y / binwidth=1;
run;
*/

/* estiamtes the variance and 95% CIs for each group */
proc univariate data=Have CIBASIC;
   class Group;
   var Y;
   ods select BasicIntervals;
   ods output BasicIntervals=Simple1;
run;

/* merge the sample size into the BasicIntervals table */
proc means data=Have noprint;
   class Group;
   var Y;
   output out=OutMeans N=n Var=Var Std=Std;
run;

data Simple;
keep Group n Parameter Estimate LowerCL UpperCL;
merge Simple1 OutMeans(where=(_TYPE_=1));
by Group;
run;

proc print data=Simple noobs; 
where Parameter="Variance";
var Parameter Group n Estimate LowerCL UpperCL;
run;


/* Compute the pooled variance as a weighted average.
   The weight of the i_th sample is (n_i - 1), where
   n_i is the size of the i_th sample.

   "To pool" means to combine information from several sources.
   For example: "If we pool our resources, we can buy Mom a nice gift."
*/

data PooledVar;
Labl = "Pooled Variance";
retain varPooled Wt 0;
set OutMeans(where=(_TYPE_=1)) end=EOF;
varPooled = varPooled + (n-1)*Var;    /* numerator */
Wt = Wt + (n-1);                      /* denominator */
if EOF then do;
   VarPooled = VarPooled / Wt;        /* weighted average */
   StdPooled = sqrt(varPooled);       /* pooled std dev */
   output;
end;
keep Labl VarPooled StdPooled;
run;

/*
proc print data=PooledVar noobs;
run;
*/

/* Although PROC UNIVARIATE computes the confidence intervals, 
   you can also compute them manaually. 
   The formula for confidence interval for Variance in the UNIVARIATE 
   doc, section 
   "Confidence Limits for Parameters of the Normal Distribution"
*/
/*
data OutMeans2;
alpha = 0.05;
length GroupName $ 6.;
set OutMeans(where=(_TYPE_=1));
GroupName = put(Group, 6.);
qChiL = quantile("chisq", 1-alpha/2, n-1);
qChiU = quantile("chisq", alpha/2, n-1);
LowerVar = Var * (n-1)/qChiL;
UpperVar = Var * (n-1)/qChiU;
LowerStd = sqrt(LowerVar);
UpperStd = sqrt(UpperVar);
run;

proc print data=OutMeans2 noobs; 
var Group n Var LowerVar UpperVar Std LowerStd UpperStd;
run;
*/

/* Visualize the variances and the pooled variance.
   Merge the group estimates and the pooled estimate. */
data AllVar;
set Simple(where=(Parameter="Variance") rename=(Estimate=Variance))
    PooledVar;
run;

ods graphics / width=480px height=200px;
title "Group Estimates of Variance for Three Groups";
title2 "Pooled Variance Shown as Horizontal Line";
proc sgplot data=AllVar;
   refline VarPooled / axis=y label=Labl splitchar=" ";
   scatter y=Variance x=Group / yerrorLower=LowerCL yerrorUpper=UpperCL
           markerattrs=(symbol=CircleFilled) name="v" 
           legendlabel="Group Variance";
   legenditem type=line name="ci" / label="95% Confidence Interval" lineattrs=GraphData2;
   keylegend "v" "ci";
run;

/*
title "Group Estimates of Standard Deviation for Three Groups";
title "Pooled Standard Deviation Shown as Vertical Line";
proc sgplot data=AllVar;
   refline StdPooled / axis=x label=Labl;
   scatter x=Std y=GroupName / xerrorLower=LowerStd xerrorUpper=UpperStd;
run;
*/

/* The graph does not answer the question, "is the variance
   of the groups constant?" To answer that question,
   you can perform a homogeneity of variance" test
   by using the HOVFTEST option on the MEANS statement in 
   PROC GLM.

   Notice you can use wildcards to match the ODS table:
   https://blogs.sas.com/content/iml/2018/11/19/select-ods-tables-wildcards-regular-expressions-sas.html
*/
proc GLM data=Have;
   class Group;
   model Y = Group;
   means Group / HOVTEST;
   ods select where=(_label_ ? 'Box Plot') HOVFTest Means;
quit;

title;
