/* SAS program to accompany the article 
   "Bimodal and unimodal beta distributions"
   by Rick Wicklin, published 29APR2024 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2024/04/29/bimodal-unimodal-beta.html

   This program shows how to use the moment-ratio diagram to 
   choose beta distributions that have a specified skewness and 
   kurtosis. The beta distribution has two major shapes:
   - For LOW kurtosis, the beta distribution is bimodal and U-shaped
   - For HIGH kurtosis, the beta distribution in unimodal. It can be
     + L-shaped: The mode is at x=0
     + J-shaped: The mode is at x=1
     + A-shaped: The mode is in (0, 1)
*/

/*****************************/
/* 1. Visualize the moment-ratio diagram and the line that 
   divides the Beta region into unimodal vs bimodal regions */
data BetaRegion;
do sx = 0 to 2.4 by 0.025;
   kLower = 1 + sx**2;     output;  /* boundary of impossible region */
   kUpper = 3 + 1.5*sx**2; output;  /* gamma line = boundry of Beta region */
end;
run;

/* The formulas for the skewness and kurtosis of the Beta(a,b) distribution
   are available at 
   https://en.wikipedia.org/wiki/Beta_distribution
*/
%let maxSkew = 2.4;
%let maxKurt = 12;
data Bimodal;
b = 1;
do a =0.01 to 0.99 by 0.01;  
   skew = ( 2*(b-a)*sqrt(a+b+1) ) /
             ( (a+b+2)*sqrt(a*b) );
   kurt = 3 + 6* ( (a-b)**2 * (a+b+1) - a*b*(a+b+2) ) /
                    ( a*b*(a+b+2)*(a+b+3) );
   if skew >= 0 & skew<= &maxSkew & kurt <= &maxKurt then
      output;
end;
run;

data Text;
text = "Bimodal "; tx = 0.3; ty = 1.5; output;
text = "Unimodal"; tx = 1.0; ty = 3.6; output;
run;

data PlotData2;
set BetaRegion Bimodal Text;
run;

title "Skewness and Kurtosis for the Beta Distribution";
title2 "Curve That Divides Unimodal and Bimodal Distributions";
proc sgplot data=PlotData2 noautolegend;
   band x=sx lower=kLower upper=kUpper / legendlabel="Beta Region" transparency=0.5;
   series x=skew y=kurt;
   text x=tx y=ty text=text / textattrs=(size=12);
   yaxis grid reverse label="Full Kurtosis" max=&maxKurt values=(1 to 12) valueshint;
   xaxis grid label="Skewness" max=&maxSkew min=0 offsetmin=0 offsetmax=0 values=(0 to 2.4 by 0.2) valueshint;
run;

/*****************************/
/* 2. Specify values in the moment-ratio diagram for which
   the Beta distribution has a variety of (s,k) values */
data SKBetaParms;
do s = 0 to 2.4 by 0.1;
   s = round(s,0.1);
   kLow = 1 + s**2;                /* boundary of impossible region */
   kHigh= 3 + 1.5*s**2;            /* gamma line = boundry of Beta region */
   kUnimod = 0.5*kLow + 0.5*kHigh; /* middle of the Beta region */
   kBimod  = 0.9*kLow + 0.1*kHigh; /* close to the impossible region */
   output;
end;
keep s kUnimod kBimod;
run;

data PlotData3;
set BetaRegion Bimodal Text SKBetaParms;
run;

title "Skewness and Kurtosis for the Beta Distribution";
title2 "Curve That Divides Unimodal and Bimodal Distributions";
proc sgplot data=PlotData3 noautolegend;
   band x=sx lower=kLower upper=kUpper / legendlabel="Beta Region" transparency=0.5;
   series x=skew y=kurt;
   scatter x=s y=kUnimod;
   scatter x=s y=kBimod;
   text x=tx y=ty text=text / textattrs=(size=12);
   yaxis grid reverse label="Full Kurtosis" max=&maxKurt values=(1 to 12) valueshint;
   xaxis grid label="Skewness" max=&maxSkew min=0 offsetmin=0 offsetmax=0 values=(0 to 2.4 by 0.2) valueshint;
run;

proc print data=SKBetaParms;
where s=1.0;
run;

/*****************************/
/* 3. You can solve the inverse problem to find the (a,b) values that 
   correspond to 
   (s,k)=(1, 2.25) and (s,k)=(1,3.25)
   For the details, see 
   https://blogs.sas.com/content/iml/2024/04/15/beta-skewness-kurtosis.html

   The results are:
   s   k    a    b 
   1  3.25 0.71 2.29 
   1  2.25 0.09 0.24 
*/
data PDF;
array _a[2] (0.71, 0.09);
array _b[2] (2.29, 0.24);
do i = 1 to dim(_a);
   a = _a[i];
   b = _b[i];
   /* https://blogs.sas.com/content/iml/2018/08/08/plot-curves-two-categorical-variables-sas.html */
   Group = catt("(a,b) = (", putn(a,5.2)) || "," || catt(putn(b,5.2)) || ")";
   do x = 0.0001, 0.001, 0.005, 0.01 to 0.99 by 0.01, 0.999;
      PDF = pdf("Beta", x, a,b);
      output;
   end;
end;
keep a b x PDF Group;
run;

title "Beta Distributions with Specified Skewness and Kurtosis";
title2 "Skewness=1; Kurtosis in {2.25, 3.25}";
proc sgplot data=PDF;
  series x=x y=PDF / group=Group;
  xaxis grid;
  yaxis grid min=0 max=5 label="Density of Beta(a,b)";
run;
