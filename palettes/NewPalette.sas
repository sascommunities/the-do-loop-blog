
/* SAS program to accompany the article 
   "A statistical palette of Christmas colors"
   by Rick Wicklin, published 15DEC2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/12/15/statistical-palette-christmas.html

   This program shows how to create a new palette of colors
   by using principal component analysis applied to colors
   from existing (harmonious) palettes.
*/

/* Convert hex colors to RGB. See
   https://blogs.sas.com/content/iml/2014/10/06/hexadecimal-to-rgb.html
*/
data RGB;
length Name $30 color $8;
length palette $450; /* must be big enough to hold all colors */
retain ID 1 palette;
input Name 1-22 n @;
do col = 1 to n;
   input color @;
   R = inputn(substr(color, 3, 2), "HEX2."); /* get RGB colors from hex */
   G = inputn(substr(color, 5, 2), "HEX2.");
   B = inputn(substr(color, 7, 2), "HEX2.");
   palette = catx(' ', palette, color);
   output;
   ID + 1;    /* Assign a unique ID to each color. */
end;
call symput('AllColors', palette);  /* 3. Concatenate colors; store in macro. */
drop n col palette;
/* Palette Name      |n| color1 | color2 | color2 | color4 | ... */
datalines;
Jingle Bell           6 CX44690D CX698308 CXF3C543 CXFFEDC7 CXCA2405 CX9E1007
Holiday Red and Gold  6 CXCF1931 CXAD132D CXD9991A CXEAA61E CXF2BC13 CX216A1B
Green and Gold        6 CX03744B CX15885A CX1E9965 CXFBE646 CXFBC34D CXF69F44
Unwrapping My Gifts   5 CX2A5B53 CX5EB69D CXECEBF1 CXD34250 CX5F1326
Christmas Wrapping    6 CX237249 CX3B885C CXE5D6B5 CXE3CD8E CXDA111E CXC00214
Christmas Wedding     6 CX325C39 CX9C9549 CXDBAA46 CXFFE9D9 CXFF4A4A CXDB2121
Real Christmas Tree   6 CX779645 CX497542 CX274530 CX6E3C3B CXBF403B CXEDB361
;

ods trace off;
ods graphics / reset;
/* use N=3 to show that the 3rd PC separates blue from green */
proc princomp data=RGB N=2 out=PCOut plots(only)=(Score Pattern);
   var R G B;
   ID color;
run;

title "Principal Component Scores for Christmas Colors";
ods graphics / width=640px height=400px ;
proc sgplot data=PCOut noautolegend;
   label Prin1="Component 1 (60.37%)" Prin2="Component 2 (28.62%)";
   styleattrs datacontrastcolors=(&AllColors);
   scatter x=Prin1 y=Prin2 / markerattrs=(symbol=CircleFilled size=16)
        group=ID;
   refline 0 / axis=x;
   refline 0 / axis=y;
run;

/* To look at the paths of specific palettes, 
   we can't use the ALlColors macro. Instead, 
   we need to create a 'SomeColors' macro based on the subset of colors
   WHERE name in ("Jingle Bell" "Holiday Red and Gold" "Unwrapping My Gifts");
*/
%let WhereClause = name in ("Jingle Bell" "Unwrapping My Gifts");
data _null_;
set RGB(where=(&WhereClause)) end=eof;
length pal2 $450; /* must be big enough to hold all colors */
retain pal2 ' ';
pal2 = catx(' ', pal2, color);
if eof then call symput('SomeColors', pal2);    
run;
%put &=SomeColors;

/* visualize two paths */
title "Paths of Select Christmas Palettes";
ods graphics / width=640px height=400px;
proc sgplot data=PCOut noautolegend;
   where &WhereClause;
   label Prin1="Component 1 (60.37%)" Prin2="Component 2 (28.62%)";
   styleattrs datacontrastcolors=(&SomeColors);
   series x=Prin1 y=Prin2 / group=name lineattrs=(color=gray thickness=3) name="path" 
          curvelabel curvelabelpos=start;
   scatter x=Prin1 y=Prin2 / markerattrs=(symbol=CircleFilled size=16)
        group=ID;
   refline 0 / axis=x;
   refline 0 / axis=y;
run;

/* add a circle and a hexagon */
data Poly;
retain cnt 0;
pi = constant('pi');
r = 1.1;                               /* distance from origin in PC space */
/* points on a circle */
PID = 1;
do deg = 0 to 360 by 6;
   angle = deg/180*pi;
   cx = r*cos(angle); cy = r*sin(angle);
   output;
end;
/* points on a regular hexagon */
PID = 2;
do deg = 210, 270, 330, 30, 90, 150;    /* angles in degrees for a regular hexagon */
   cnt + 1;
   angle = deg/180*pi;
   hx = r*cos(angle); hy = r*sin(angle);
   output;
end;
run;

data PC2;
set PCOut Poly;
run;

/* visualize the path of the new Christmas palette */
ods graphics / width=640px height=400px;
title "Path of a New Christmas Palette";
proc sgplot data=PC2 noautolegend;
   label Prin1="Component 1 (60.37%)" Prin2="Component 2 (28.62%)";
   styleattrs datacontrastcolors=(&AllColors);
   scatter x=Prin1 y=Prin2 / markerattrs=(symbol=CircleFilled size=16) group=ID;
   polygon x=cx y=cy ID=PID / lineattrs=(color=gray thickness=1 pattern=solid);
   scatter x=hx y=hy / markerattrs=(symbol=StarFilled size=16 color=darkgray) 
           datalabel=cnt datalabelattrs=(size=12 color=CX444444);
   refline 0 / axis=x;
   refline 0 / axis=y;
run;


/* use PCA and the inverse PCA transformation to 
   create a new Christmas palette */
proc iml;
/* read in all RGB colors */
use PCOut;
read all var {R G B} into X;
close;

/* for correlation matrix; form eigenvector decomposition (spectral decomposition) */
R = corr(X);                    /* correlation matrix (from PROC PRINCOMP) */
call eigen(D, V, R);            /* D=eigenvalues; columns of V are PCs */
*print D V;
/* CHECK */
decomp = V*diag(D)*V`;          /* should equal corr matrix R */
*print R decomp;

/* compute all PC scores: transformation from RGB to PC */
scale = std(X);
c = mean(X);
B = (X - c) / scale;            /* center and scale data about the mean */
score = B * V;                  /* project standardized data onto the PCs */
/* test transformation from RGB to PC */
/* if I did it correctly, the mean score in each component is 0 */
*print (mean(score))[L="mean(score)=(0,0,0)"];

/* Since V is orthogonal, V` = inv(V)
   Let's compute the inverse transformation as
      NewRGB = scale # (score*V`) + c;
   for any set of scores  */

/* Check: do we recover the original data from the scores? */
*orig = scale # (score*V`) + c;
*print (max(abs(X-orig)));    /* should be 0 */

/* let's visualize points on the circle on the plane (PC1, PC2, 0) */
startDeg = 210;
deg = T(do(startDeg, startDeg+360, 5));              /* angles in degrees */
angles = deg * constant('pi')/180;   /* angles in radians */
k = nrow(angles);
r = 1.1;
CircleScores = r*cos(angles) || r*sin(angles) || j(k,1,0);  /* in span(PC1, PC2) */
CircleColors = round( scale # (CircleScores*V`) + c );
/* The range of a color is [0, 255], so make sure you do not get too far from 
   the origin or you could start generating invalid colors. For example, 
   r = 1.4 give invalid colors such as
   267 103 47 
   ^^^
   Too big!
*/
if any(CircleColors > 255) then print "ERROR";
Hex = cats('CX', rowcat(putn(CircleColors,'HEX2.'))); /* convert from RGB to Hex */
ods graphics / width=640px height=160px;
run HeatmapDisc(1:k, Hex) title="Colors on a Circle" ShowLegend=0 displayoutlines=0;

/* choose 6 points on regular hexagon on the circle */
deg = {210, 270, 330, 30, 90, 150};    /* angles in degrees for a regular hexagon */
angles = deg * constant('pi')/180;     /* angles in radians */
k = nrow(angles);
r = 1.1;                               /* distance from origin in PC space */
CircleScores = r*cos(angles) || r*sin(angles) || j(k,1,0);  /* in span(PC1, PC2) */
CircleColors = round( scale # (CircleScores*V`) + c );
print deg CircleColors;

Hex = cats('CX', rowcat(putn(CircleColors,'HEX2.'))); /* convert from RGB to Hex */
ods graphics / width=640px height=160px;
run HeatmapDisc(1:k, Hex) title="New Christmas Palette" ShowLegend=0 
               xvalues=Hex yvalues="Rick's Palette";

QUIT;
