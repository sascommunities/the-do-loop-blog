/*************************/
/* pi_th roots of unity */
/*************************/
/* SAS program to accompany the article 
   "The pi_th roots of unity"
   by Rick Wicklin, published 14MAR2022 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2022/03/14/pi_th-roots-of-unity.html

   This program shows how to visualize the pi_th roots of unity.
   There are infinitely many pi_th roots of unity. 
   They are complex numbers of the form exp(2 i n) for any integer n.
   Geometrically, you can visualize these points as the rotations 
   by +/-2 radians of the point (1,0) in the Argand plane. The 
   pi_th roots of unity are dense in the unit circle. 

   If you restrict your attention to only positive rotations 
   (the images of the positive even integers), you can classify 
   each root according to whether it is the image of an even 
   integer in the Beatty sequence for pi or in the complementary 
   sequence. The points that are in the Beatty sequence are 
   mapped to the arcs for angles (pi-1, pi) and (2 pi-1, 2 pi). 
   These arcs are (1/pi)th of the arclength of the unit circle. 
   For more about the Beatty sequence for pi, see
   https://blogs.sas.com/content/iml/2022/03/09/beatty-sequence-pi.html
*/

/* If you wrap the real line around the unit circle and mark the 
   image of the EVEN integers, you get the pi_th roots of unity, 
   as described in 
   https://sxpmaths.wordpress.com/2015/03/14/%CF%80-th-roots-of-unity/
   The pi_th roots of unity are exp(2 i n).
*/
%let maxN = 355;
data Roots;
do n = 1 to &maxN;
   s = 2*n;
   xn = cos(s);   /* image of even integers under mapping to unit circle */
   yn = sin(s);
   output;
end;
run;

%macro PlotCircle();
   ellipseparm semimajor=1 semiminor=1 / xorigin=0 yorigin=0 outline
               lineattrs=(color=cxCCCCCC);
%mend;

ods graphics / width=400px height=400px;
title "Wrap Positive Even Integers onto a Circle";
title2 "First 25 Even Integers";
proc sgplot data=Roots(obs=25) aspect=1 noautolegend;
   %PlotCircle;
   scatter x=xn y=yn / datalabel=s;
   refline 0 / axis=x;
   refline 0 / axis=y;
   xaxis display=none;
   yaxis display=none;
run;

title2 "First &maxN Even Integers";
proc sgplot data=Roots aspect=1 noautolegend;
   %PlotCircle;
   scatter x=xn y=yn;
   refline 0 / axis=x;
   refline 0 / axis=y;
   xaxis display=none;
   yaxis display=none;
run;

/* add information about Beatty sequence */
data Beatty2;
/* find the lessor of the maximum values of the sequences */
pi = constant('pi');
MB = floor(pi* 2*&maxN);
MC = floor(pi/(pi-1)* 2*&maxN);         
maxS = min(MB, MC);
call symputx('maxS', maxS);
drop MB MC maxS pi;

Sequence = 'Beatty       '; 
do n = 1 to 2*&maxN;
   s = floor(pi*n);         /* the Beatty sequence */
   if s <= maxS then output;
end;
Sequence = 'Complementary'; 
do n = 1 to 2*&maxN;
   s = floor(pi/(pi-1)*n);  /* the complementary sequence */
   if s <= maxS then output;
end;
run;

%put &=maxS;   /* display value in SAS log */
proc sort data=Beatty2;
   by s;
run;

data Combine;
merge Beatty2(drop=n) Roots;
by s;
run;

ods graphics / width=400px height=400px;
title "Wrap Even Integers onto a Circle";
title2 "Group by Beatty Sequence";
proc sgplot data=Combine aspect=1;
   %PlotCircle;
   scatter x=xn y=yn / group=Sequence;
   refline 0 / axis=x;
   refline 0 / axis=y;
   xaxis display=none;
   yaxis display=none;
run;

/* add segment at 0 and pi-1 radians */
data Diameters;
do Sequence='Beatty       ',
            'Complementary'; 
   do theta = 0, constant('pi') - 1;
      x0 = cos(theta); y0 = sin(theta); x1 = -x0; y1 = -y0; 
      output;
   end;
end;
run;

data Combine2;
set Combine Diameters;
run;

title "Wrap Even Integers onto a Circle";
title2 "Group by Beatty Sequence";
proc sgplot data=Combine2 aspect=1;
   %PlotCircle;
   vector x=x1 y=y1 / xorigin=x0 yorigin=y0 noarrowheads lineattrs=(color=cxCCCCCC);
   scatter x=xn y=yn / group=Sequence nomissinggroup;
   xaxis display=(nolabel);
   yaxis display=(nolabel);
run;

