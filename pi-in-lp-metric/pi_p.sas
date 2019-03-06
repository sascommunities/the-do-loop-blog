/* SAS program to accompany the article 
   "The value of pi depends on how you measure distance"
   by Rick Wicklin, published 13MAR2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/03/13/pi-in-lp-metric.html 

   This program visualizes the unit circle for various Lp metrics.
   For many values of p in 1, 11], the program uses numerical integration
   (the QUAD function in SAS/IML) to compute the length of the unit 
   semicircle, which is the value of pi in the Lp metric. 
   This value is plotted against p to show that 
   1. Different Lp metrics lead to different values of pi
   2. For all p>-1, the value 3.14159... is the unique minimum value of pi

   REFERENCES:
   C. L. Adler and J. Tanton (2000) "&pi; is the Minimum Value for Pi", 
      The College Mathematics Journal, 31:2, 102-106.
      https://www.tandfonline.com/doi/abs/10.1080/07468342.2000.11974122

   J. B. Keller and R. Vakil (2009) "pi_p, the Value of pi in l_p"
      The American Mathematical Monthly, 116:10, 931-935.
      https://www.jstor.org/stable/pdf/40391254.pdf
*/

/* 1. Visualize the unit semicircle for the Lp metric */
data Circles(drop=t);
do p = 1, 1.5, 2, 3, 10;
   do t = -1 to 0 by 0.01;
      x =  t;  y=(1 - abs(x)**p)**(1/p); output; /* quarter circle, second quadrant */
      x = -t;                            output; /* first quadrant */
   end;
end;
run;
proc sort data=Circles; by p x; run;   /* for plotting, sort by p and x */

data Text;               /* add labels along line y=x */
length text $5;
do p = 1, 1.5, 2, 3, 10;
   t = 2**(-1/p);        /* 45-degree line intersects curves at (t, t) */
   text = cats("p=", p); /* label */
   output;
end;
run;

data AllCircles;  set Circles text;  run;

ods graphics / width=640px height=320px;
title "The Unit Circle for Lp Distance";
proc sgplot data=AllCircles aspect=0.5 noautolegend;
   series x=x y=y / group=p;
   text x=t y=t text=Text / group=p textattrs=(size=11) position=SW;
   refline 0 / axis=x;
   refline 0 / axis=y;
   xaxis grid values=(-1 to 1 by 0.2) offsetmin=0.01;
   yaxis grid;
run;


/* 2. Compute the value of pi in the Lp norm as the 
      length of the upper semicircle circle in the Lp metric 
*/
proc iml;
/* arclength of unit circle in Lp metric for each p > 1 */
start ArcLength(x) global(p);
   u = abs(x**(-p) - 1) ** (1-p); 
   return (1 + u)**(1/p);
finish;

s = do(1.1, 2, 0.1);       /* regular grid on [1,2] */
w = 3:10;                  /* regular grid on [3,10] */
t = (s/(s-1));             /* dual values such that 1/s + 1/t = 1 */
parm = union(1, s, t, w);  /* form union of these sets */

pi = j(ncol(parm), 1, 0);
do i = 1 to ncol(parm);
   p = parm[i];            /* set global parameter */
   /* integrate arclength on [0, 2^(-1/p)], which is 1/8th of unit circle */
   call quad(len, "ArcLength", 0||2**(-1/p)); 
   pi[i] = 4*len;          /* pi(p) for each p */
end;
create C var {parm, pi}; append; close;
quit;

data Label; t=2; C=constant('pi'); run;
data PiLP; set C Label; run;

ods graphics / width=640px height=320px;
title "The Value of pi for the Lp Metric";
proc sgplot data=PiLp noautolegend;
   series x=parm y=pi;
   refline 3.14159/axis=y label="3.14159" labelloc=inside;
   scatter x=t y=C / markerattrs=GraphData2(symbol=StarFilled size=10);
   xaxis grid values=(1 to 11) label="p (Exponent for Lp Distance)";
   yaxis grid values=(3 to 4 by 0.2) label="pi(p)";
run;
