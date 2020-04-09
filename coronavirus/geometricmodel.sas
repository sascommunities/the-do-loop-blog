/* SAS program to accompany the article 
   "On reducing the spread of coronavirus"
   by Rick Wicklin, published 08Apr2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/04/08/reducing-spread-of-coronavirus.html

   This program shows that the public-health mandates to reduce the spread of the 
   coronovirus (COVID-19) are supported by a very simple mathematical 
   model, the cumulative geometric distribution. The cumulative distribution
   shows the probability that an event occurs on the n_th trial when the 
   probability in each trial is the same value, p. 

   The public health mandates do two things:
   1. Reduce the number of interactions: stay at home, social distancing, prohibit congregating,...
   2. Reduce the probability that an interaction results in a new infection: wash hands,
      sanitize surfaces, wear a mask,....
   In addition, the model shows that even if the probability of an infdection is low, it can
   still happen in only a few interactions. Lastly, the model quantifies the how the 
   cumulative probability increase as a function of the number of trials when p << 1.
*/

/* COVID=19 application */
data Geometric;
do p = 0.0001, 0.001, 0.01;
   do n = 0 to 180;    /* event occurs on the n_th trial */
      /* there are n-1 trials for which the event did NOT occur */
      pdf = pdf("Geometric", n-1, p);  /* number of trials before first event */
      cdf = cdf("Geometric", n-1, p);  /* probability that event happens on or before n_th Trial*/
      cdf2 = 1 - (1-p)**n;    /* explicit formula */
      cdfApprox = n*p;
      output;
   end;
end;
run;

data GeoLabel;
set Geometric;
by p;
if last.p then do;
   px = n;
   py = cdf;
   pLabel = " p = " || putn(p, "best6.");
end;
run;

ods graphics / width=640px height=480px;
title "Probability That the Event Occurs on or before n_th Trial";
proc sgplot data=GeoLabel noautolegend;
   series x=n y=cdf / group=p;
   text x=px y=py text=pLabel / group=p position=right textattrs=(size=10);
   label n="Number of Trials (n)" cdf="Cumulative Probability";
   xaxis grid values=(0 to 175 by 25) offsetmax=0.3;
   yaxis grid values=(0 to 0.9 by 0.1);
   refline 0 / axis=y;
run;

proc print data=Geometric;
   where p=0.01 and n in (10, 15, 25, 50, 75, 100);
   var p n cdf;
run;

%let numSamples = 20;           *sample of 20 people;
data RandGeom;
call streaminit(19);            * random seed for covid-19;
p = 0.01;
do ID = 1 to &numSamples;
   x = rand("Geometric", p);    *number of trials before event occurs;
   Zero = 0;
   output;
end;
run;

proc sort data=RandGeom;
   by x;
run;

ods graphics / width=400px height=400px;
title "&numSamples Random Draws from Geom(0.01)";
proc sgplot data=RandGeom;
   highlow y=ID low=Zero high=x / type=bar barwidth=1 CLIPCAP highlabel=x; 
   xaxis grid MAX=90 label="Number of Trials until Event" minor;  
   yaxis type=discrete discreteorder=data reverse valueattrs=(size=7);
run;


proc means data=RandGeom;
   var x;
run;

/*************************************************/
/* Estimates of the probbility of catching the coronavirus through normal activities.

In an interview,
https://www.heraldscotland.com/news/18327916.coronavirus-maths-expert-reveals-chances-catching-covid-19/
Fergus Simpson, an expert in Bayesian statistics, said on 23Mar2020,
"If you were to carry on life as normal" (that is, no social distancing,
the odds are "about 1,000 to 1 each day" of catching the coronavirus
from everyday activities. He then added that "those odds are getting progressively worse each day"
as more people become infected. 

Some aoccupations (for example healthcare workers) might be at an elevated risk.
Public activities in "outbreak areas" can carry an increased risk. 

"By following public health guidelines such as social distancing, then this risk can be greatly reduced."  

Simpson wrote a blog post that has more information about the chances of contracting coronavirus:
https://medium.com/@fergus2/coronavirus-a-personal-risk-assessment-50003c952c03

Also, Keith Devlin
https://www.mathvalues.org/masterblog/making-sense-of-the-covid-19-data

Also Rick McAnulty:
https://inside.uncc.edu/news-features/2020-03-26/psychology-noncompliance-health-care-advisories
"People estimate the probability of risk—and decide the odds outweigh probable danger.:
*/
/*************************************************/
