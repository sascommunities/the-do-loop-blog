
/* SAS program to accompany the article 
   "A principal component analysis of color palettes"
   by Rick Wicklin, published 13DEC2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/12/13/principal-components-colors.html
   
   This program shows how to use principal component analysis to visualize
   a 2-D plane that aligns with the maximal  variation in a set of RGB colors. 
   Often this plane also shows similar colors next to each others.
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
