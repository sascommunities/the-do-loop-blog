/* SAS program to accompany the article 
   "Use SAS to create mathematical art"
   by Rick Wicklin, published 04AUG2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/08/04/create-mathematical-art-sas.html

   This program shows how to use orthogonal matrices to induce rotations and 
   reflections on planar polygons. Each matrix corresponds to an action 
   of the dihedral group on the square (D4), which is a finite group of order 8.

   Wicklin generates the data by using SAS/IML software and applying 
   rotations and reflections. If you do not have SAS/IML software, you can use a
   DATA step, which provides the coordinates of each polygon. See the end
   of this program for the DATA step.
*/
/* define the coordinates of one pinwheel shape */
data Shape;
input ID x y @@;
datalines;
0  0 0  0  0.667  0.333  0  1  1  0  0  1  0  0  0  
1  0 0  1 -0.333  0.667  1 -1  1  1 -1  0  1  0  0 
2  0 0  2 -0.667 -0.333  2 -1 -1  2  0 -1  2  0  0 
3  0 0  3  0.333 -0.667  3  1 -1  3  1  0  3  0  0 
;


proc iml;
/* actions of the D4 dihedral group */
start D4Action(v, act);
   /* the subroup of rotations by 90 degrees */
   if      act=0 then M = { 1  0,  0  1};
   else if act=1 then M = { 0 -1,  1  0};
   else if act=2 then M = {-1  0,  0 -1};
   else if act=3 then M = { 0  1, -1  0};
   /* the subgroup of reflections across horz, vert, or diagonals */
   else if act=4 then M = {-1  0,  0  1};
   else if act=6 then M = { 0 -1, -1  0};
   else if act=7 then M = { 1  0,  0 -1};
   else if act=5 then M = { 0  1,  1  0};
   return( v*M` );  /* = (M*z`)` */
finish;

use Shape; read all var {x y} into P; 
           read all var "ID"; close;

OpNames = {"Identity" "R1" "R2" "R3" "S0" "S1" "S2" "S3"};
Name = OpNames[1];
Q = {. . .};
create Panel from Name Q[c={'Name' 'ID' 'x' 'y'}];
do i = 0 to 7;
   R = D4Action(P, i);
   Q = ID || R;
   Name = j(nrow(Q), 1, OpNames[i+1]);
   append from Name Q;
end;
close;
QUIT;

%macro colOpts; colaxis offsetmin=0 offsetmax=0 display=(nolabel noticks novalues); %mend;
%macro rowOpts;  rowaxis offsetmin=0 offsetmax=0 display=(nolabel noticks novalues); %mend;

/* Arrange the 8 images in a lattice. This visualizes the D4 actions
   on the original pinwheel shape. */
ods graphics / width=720px height=480px;
title "Dihedral Group (D4) Actions on Pinwheel Figure";
title2 "Colors=(Teal, Orange, Blue, Salmon)";
proc sgpanel data=Panel noautolegend;
   styleattrs wallcolor=&gray datacolors=(&teal &orange &blue &salmon);
   panelby Name / columns=4 onepanel;
   polygon x=x y=y ID=ID / group=ID fill;
   %colOpts; %rowOpts;
run;

/* switch the colors of the 2nd and 4th vanes */
data Panel2;
set Panel;
if      ID=1 then ID=3;
else if ID=3 then ID=1;
run;

/* arrange the new 8 images in a lattice. These are the images 
   of the second pinwheel under D4. */
title "Dihedral Group (D4) Actions on Pinwheel Figure";
title2 "Colors=(Teal, Salmon, Blue, Orange)";
proc sgpanel data=Panel2 noautolegend;
   styleattrs wallcolor=&gray datacolors=(&teal &orange &blue &salmon);
   panelby Name / columns=4 onepanel;
   polygon x=x y=y ID=ID / group=ID fill;
   %colOpts; %rowOpts;
run;

/* concatenate the 16 pinwheels and plot them in one image */
data Pinwheels;
length Name $20;
set Panel Panel2;
Group = floor((_N_-1)/20);
run;

ods graphics / width=640px height=640px;
title "Dihedral Shadow";
proc sgpanel data=Pinwheels noautolegend;
   styleattrs wallcolor=&gray datacolors=(&teal &orange &blue &salmon);
   panelby Group / columns=4 onepanel noheader;
   polygon x=x y=y ID=ID / group=ID fill;
   %colOpts; %rowOpts;
run;


/*********************************************************************/
/* Data and SAS code for challenge posted on the SAS Support Community
/*********************************************************************/


/* If you do not have SAS/IML software, you can still create the 
   lattice by using PROC SGPANEL to plot the following data, which 
   contains the vertices of the 16 pinwheels */
data Pinwheels;
input Group ID @;
do i = 1 to 5;
   input x y @@;
   output;
end;
drop i;
datalines;
0 0 0 0   0.667  0.333   1  1   0  1  0 0 
0 1 0 0  -0.333  0.667  -1  1  -1  0  0 0 
0 2 0 0  -0.667 -0.333  -1 -1   0 -1  0 0 
0 3 0 0   0.333 -0.667   1 -1   1  0  0 0 
1 0 0 0  -0.333  0.667  -1  1  -1  0  0 0 
1 1 0 0  -0.667 -0.333  -1 -1   0 -1  0 0 
1 2 0 0   0.333 -0.667   1 -1   1  0  0 0 
1 3 0 0   0.667  0.333   1  1   0  1  0 0 
2 0 0 0  -0.667 -0.333  -1 -1   0 -1  0 0 
2 1 0 0   0.333 -0.667   1 -1   1  0  0 0 
2 2 0 0   0.667  0.333   1  1   0  1  0 0 
2 3 0 0  -0.333  0.667  -1  1  -1  0  0 0 
3 0 0 0   0.333 -0.667   1 -1   1  0  0 0 
3 1 0 0   0.667  0.333   1  1   0  1  0 0 
3 2 0 0  -0.333  0.667  -1  1  -1  0  0 0 
3 3 0 0  -0.667 -0.333  -1 -1   0 -1  0 0 
4 0 0 0  -0.667  0.333  -1  1   0  1  0 0 
4 1 0 0   0.333  0.667   1  1   1  0  0 0 
4 2 0 0   0.667 -0.333   1 -1   0 -1  0 0 
4 3 0 0  -0.333 -0.667  -1 -1  -1  0  0 0 
5 0 0 0   0.333  0.667   1  1   1  0  0 0 
5 1 0 0   0.667 -0.333   1 -1   0 -1  0 0 
5 2 0 0  -0.333 -0.667  -1 -1  -1  0  0 0 
5 3 0 0  -0.667  0.333  -1  1   0  1  0 0 
6 0 0 0  -0.333 -0.667  -1 -1  -1  0  0 0 
6 1 0 0  -0.667  0.333  -1  1   0  1  0 0 
6 2 0 0   0.333  0.667   1  1   1  0  0 0 
6 3 0 0   0.667 -0.333   1 -1   0 -1  0 0 
7 0 0 0   0.667 -0.333   1 -1   0 -1  0 0 
7 1 0 0  -0.333 -0.667  -1 -1  -1  0  0 0 
7 2 0 0  -0.667  0.333  -1  1   0  1  0 0 
7 3 0 0   0.333  0.667   1  1   1  0  0 0 
8 0 0 0   0.667  0.333   1  1   0  1  0 0 
8 3 0 0  -0.333  0.667  -1  1  -1  0  0 0 
8 2 0 0  -0.667 -0.333  -1 -1   0 -1  0 0 
8 1 0 0   0.333 -0.667   1 -1   1  0  0 0 
9 0 0 0  -0.333  0.667  -1  1  -1  0  0 0 
9 3 0 0  -0.667 -0.333  -1 -1   0 -1  0 0 
9 2 0 0   0.333 -0.667   1 -1   1  0  0 0 
9 1 0 0   0.667  0.333   1  1   0  1  0 0 
10 0 0 0  -0.667 -0.333  -1 -1  0 -1  0 0 
10 3 0 0   0.333 -0.667   1 -1   1  0  0 0 
10 2 0 0   0.667  0.333   1  1   0  1  0 0 
10 1 0 0  -0.333  0.667  -1  1  -1  0  0 0 
11 0 0 0   0.333 -0.667   1 -1   1  0  0 0 
11 3 0 0   0.667  0.333   1  1   0  1  0 0 
11 2 0 0  -0.333  0.667  -1  1  -1  0  0 0 
11 1 0 0  -0.667 -0.333  -1 -1   0 -1  0 0 
12 0 0 0  -0.667  0.333  -1  1   0  1  0 0 
12 3 0 0   0.333  0.667   1  1   1  0  0 0 
12 2 0 0   0.667 -0.333   1 -1   0 -1  0 0 
12 1 0 0  -0.333 -0.667  -1 -1  -1  0  0 0 
13 0 0 0   0.333  0.667   1  1   1  0  0 0 
13 3 0 0   0.667 -0.333   1 -1   0 -1  0 0 
13 2 0 0  -0.333 -0.667  -1 -1  -1  0  0 0 
13 1 0 0  -0.667  0.333  -1  1   0  1  0 0 
14 0 0 0  -0.333 -0.667  -1 -1  -1  0  0 0 
14 3 0 0  -0.667  0.333  -1  1   0  1  0 0 
14 2 0 0   0.333  0.667   1  1   1  0  0 0 
14 1 0 0   0.667 -0.333   1 -1   0 -1  0 0 
15 0 0 0   0.667 -0.333   1 -1   0 -1  0 0 
15 3 0 0  -0.333 -0.667  -1 -1  -1  0  0 0 
15 2 0 0  -0.667  0.333  -1  1   0  1  0 0 
15 1 0 0   0.333  0.667   1  1   1  0  0 0 
;

%let teal   = CX288c95;
%let orange = CXeba411;
%let blue   = CX0f5098;
%let salmon = CXd5856e;
%let gray   = CX929386;

ods graphics / width=640px height=640px;
%macro colOpts; colaxis offsetmin=0 offsetmax=0 display=(nolabel noticks novalues); %mend;
%macro rowOpts;  rowaxis offsetmin=0 offsetmax=0 display=(nolabel noticks novalues); %mend;

title "Dihedral Shadow";
proc sgpanel data=Pinwheels noautolegend;
   styleattrs wallcolor=&gray datacolors=(&teal &orange &blue &salmon);
   panelby Group / columns=4 onepanel noheader;
   polygon x=x y=y ID=ID / group=ID fill;
   %colOpts; %rowOpts;
run;
