/* SAS program to accompany the article 
   "Visualize palettes of colors in SAS"
   by Rick Wicklin, published 08DEC2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/12/08/visualize-palette-colors.html

   This program shows how to use SAS to visualize
   multiple palettes of colors.
*/

ods graphics/reset;

/* Acknowledgements:
   Jingle Bell palette:  https://www.schemecolor.com/jingle-bells.php
   Holiday Red and Gold: https://www.schemecolor.com/holiday-red-gold.php
   Green and Gold:       https://www.schemecolor.com/green-and-gold.php
   Unwrapping My Gifts:  https://www.schemecolor.com/unwrapping-my-gifts.php
   Christmas Wrapping:   https://www.schemecolor.com/christmas-wrapping.php
   Christmas Wedding:    https://www.schemecolor.com/christmas-wedding.php
   Real Christmas Tree:  https://www.schemecolor.com/real-christmas-tree.php
*/

/* 1. Use the DATA step to read the colors in each palette. */
data Palettes;
length Name $30 color $8;
length palette $450; /* must be big enough to hold all colors */
retain ID 1 row 1 palette;
input Name 1-22 n @;
do col = 1 to n;
   input color @;
   palette = catx(' ', palette, color);
   output;
   ID + 1;    /* Assign a unique ID to each color. */
end;
row + 1;
call symput('AllColors', palette);  /* 3. Concatenate colors; store in macro. */
drop palette;
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

%put &AllColors;   /* check that all colors are here */


%let gray   = CXF8F8F8;
title "First Visualization of Christmas Palettes";
title2 "A Scatter Plot";
proc sgplot data=palettes noautolegend noborder;
   styleattrs backcolor=&gray wallcolor=&gray datacontrastcolors=(&AllColors);
   scatter x=col y=row / group=ID markerattrs=(symbol=SquareFilled size=45);
   text x=col y=row text=color;
run;


/* 4. Create a polygon for each color. 
      The i_th color for the j_th palette is represented by a rectangle 
      that is centered at the coordinate (i,j). */
data Poly;
set palettes;
w = 0.5; h = 0.48;                /* width and height of rectangle */
x = col - w; y = row - h; output; /* four vertices of square */
x = col + w; y = row - h; output;
x = col + w; y = row + h; output;
x = col - w; y = row + h; output;
drop w h;
run;

%let gray   = CXF8F8F8;
title "Second Visualization of Christmas Palettes";
title2 "A Polygon Plot";
proc sgplot data=Poly noautolegend noborder;
   styleattrs backcolor=&gray wallcolor=&gray datacolors=(&AllColors);
   polygon x=x y=y ID=ID / group=ID fill;
   text x=col y=row text=color;
run;

/* 5. Optionally, use the TEXT statement to display the name of each palette. */
data Temp;
set Poly;
by ID;
if ^first.ID then
   color=" ";
run;

/* 6. Optionally, use the TEXT statement to overlay the hex values of the colors */
data PlotPalettes;
set Temp;
by Name notsorted;
label = Name;
if first.Name then do;
   lx = col - 0.6;    ly = row;
end;
else do;
   label=" ";  lx = .;  ly = .;
end;
run;

/* 7. Use the POLYGON and TEXT statements to display the palettes of colors. */
%let gray   = CXF8F8F8;
title "Visualization of Christmas Palettes";
proc sgplot data=PlotPalettes noautolegend noborder;
   styleattrs backcolor=&gray wallcolor=&gray datacolors=(&AllColors);
   polygon x=x y=y ID=ID / group=ID fill;
   text x=col y=row text=color;
   text x=lx y=ly text=label / position=left textattrs=(size=12);
   xaxis display=none offsetmax=0;
   yaxis display=none thresholdmin=0 thresholdmax=0;
run;
