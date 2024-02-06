/* SAS program to accompany the article 
   "The elliptical heart"
   by Rick Wicklin, published 07FEB2024 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2024/02/07/elliptical-heart.html

   This program shows how to contruct a heart-shaped region from two ellipses.
   The idea is from 
   Ted Conway, 2023, "Total Ellipse of the heart"
   https://communities.sas.com/t5/Graphics-Programming/Fun-w-SAS-ODS-Graphics-Total-Ellipse-of-the-Heart/m-p/858764
   who adapted the idea from Amanda Newton, "The Heart of Mathematics"
   http://jwilson.coe.uga.edu/EMAT6680Su10/Newton/emat6690/Hearts/Hearts.html

   Rick Wicklin shows how to parameteriz the heart-shaped region and plot it by using the 
   POLYGON statement in PROC SGPLOT in SAS.
*/

/* GOAL: shade ONLY the heart-shaped region by using polar coordinates */
ods graphics / reset;

/* 1. Show Ted Conway's version. He used a different scale than Amanda.
   If you want to duplicate Amanda Newton's ellipse, use maxR and minR */
data textOnly;
retain tx ty  0   msg  "HAPPY/VALENTINE'S/DAY"; 
run;

%let maxR = 1.4142;      *semimajor=sqrt(2);
%let minR = 0.8165;      *semiminir=sqrt(2/3);
%let ellipseAttrs = lineattrs=(color=red) lineattrs=(thickness=2pt) nofill fillattrs=(color=cxC3540C);
title "Ted's Heart from Two Ellipses";
proc sgplot data=textOnly aspect=1 noautolegend;
   ellipseparm semimajor=&maxR semiminor=&minR / slope= 1 &ellipseAttrs; *Ellipse 1, angled at +45 degrees;
   ellipseparm semimajor=&maxR semiminor=&minR / slope=-1 &ellipseAttrs; *Ellipse 2, angled at -45 degrees;
   * "Happy Valentine's Day!" message;
   text x=tx y=ty text=msg / textattrs=(size=22pt color=red) splitchar='/' splitpolicy=splitalways contributeoffsets=none; 
   xaxis display=none;
   yaxis display=none;
run;


/* 2. Parameterize Ellipse #1.
      Draw a slanted ellipse x^2 + y^2 - x*y = 1 in polar coordinates, 
      which is r^2 = 1 / (1 - (1/2)sin(2*theta)) */
data Ellipse;
length label $4;
ID = 1;                                     /* ID var for POLYGON plot */
pi = constant('pi');
do theta = 0 to 2*pi by pi/100;
   r = sqrt(1 / (1 - 0.5*sin(2*theta)));
   x = r*cos(theta);  y = r*sin(theta);     /* plot in Cartesian coords */
   output;
end;
/* add labels for a few values of theta */
do theta = 0 to 2*pi-0.1 by pi/4;
   r = sqrt(1 / (1 - 0.5*sin(2*theta)));
   xx = r*cos(theta); yy = r*sin(theta);
   label = putn(theta, 4.2);
   output;
end;
drop pi;
run;

title "Ellipse in Polar Coordinates";
title2 "x^2 + y^2 - x*y = 1";
proc sgplot data=Ellipse aspect=1 noautolegend;
   polygon x=x y=y ID=ID / fill outline transparency=0.5;
   text x=xx y=yy text=label;
   xaxis grid;
   yaxis grid;
run;


/* 3. Draw only the right half of the ellipse. Note that the polar angle for 
      this region is -pi/2 <= theta <= pi/2 */
%let pi = 3.1415926536;
title "Half-Ellipse";
title2 "-pi/2 <= theta <= pi/2";
proc sgplot data=Ellipse aspect=2 noautolegend;
   where theta <= &pi/2 OR theta >= 3*&pi/2-0.001;
   polygon x=x y=y ID=ID / fill;
   text x=xx y=yy text=label;
   xaxis grid;
   yaxis grid;
run;


/* 4. To get the entire heart-shaped region, reflect (x,y) to (-x,y).
      You also need to adjust the polar angle. If theta is the polar
      angle for (x,y), then pi-theta is the polar angle for (-x,y).
      You can then sort by the polar angle and use the POLYGON statement
      to visualize the reqion. */
data Heart;
retain ID 1;
pi = constant('pi');
do t = -pi/2 to pi/2 by pi/100;
   r = sqrt(1 / (1 - 0.5*sin(2*t)));
   y = r*sin(t);
   /* for this value of y, save both x and -x */
   theta = t;      x = r*cos(theta);   output;  /* right half of heart */
   theta = pi - t; x = -x;             output;  /* left half of heart */
end;
run;

proc sort data=Heart;
   by theta;  /* reorder to get parameterized heart */
run;

data textOnly;  /* message from tc's post */
retain tx ty  0   msg  "HAPPY/VALENTINE'S/DAY"; 
run;

data All;
   set textOnly Heart;
run;

title "Total Ellipse of the Heart";
proc sgplot data=All noborder nowall noautolegend aspect=1;
   polygon x=x y=y ID=ID / fill fillattrs=(color=pink);
   text x=tx y=ty text=msg / textattrs=(size=22pt color=red) splitchar='/' splitpolicy=splitalways contributeoffsets=none;
   xaxis display=none;
   yaxis display=none;
run;


/* 5. In the previous image, the ellipses are gone, so it doesn't look like
      a "total ellipse of the heart" any more. Overlay the ellipses
      to emphasize the method of construction. */
%let maxR = 1.4142;      *semimajor=sqrt(2);
%let minR = 0.8165;      *semiminir=sqrt(2/3);
title "Total Ellipse of the Heart";
title2 "Showing Underlying Ellipses";
proc sgplot data=All aspect=1 noautolegend;
   polygon x=x y=y ID=ID / fill fillattrs=(color=pink);
   ellipseparm semimajor=&maxR semiminor=&minR / slope= 1 lineattrs=(color=lightred) nofill transparency=0.8;
   ellipseparm semimajor=&maxR semiminor=&minR / slope=-1 lineattrs=(color=lightred) nofill transparency=0.8;
   text x=tx y=ty text=msg / textattrs=(size=22pt color=red) splitchar='/' splitpolicy=splitalways contributeoffsets=none;
   xaxis grid;
   yaxis grid;
run;
