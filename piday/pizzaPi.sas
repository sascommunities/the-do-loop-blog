/* SAS program to accompany the article 
   "Pizza pi"
   by Rick Wicklin, published 11MAR2024 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2024/03/11/pizza-pi.html ‎

   This program demonstrates the geometric "method of exhaustion"
   that Archimedes used to prove that the area of a circle is pi*r^2

   This program uses SAS IML software. In PROC IML, you can 
   compute the coordinates for a sector that is 1/N of the circle.
   You can think of the sector as a polygon.
   You can easily rotate and translate the sector, thus creating
   an alternating up-down sequence of pizza slices.
   For a better effect, you can plot only half of the first slice
   on the left and put the other half on the right, so that the 
   sequence of polygons looks more rectangular.
*/

/* 1. Show a circle cut into 8 equal sectors (or slices). 
      Alternate colors for odd and even slices.
*/
data Circle;
pi = constant('pi');
N=8; r=1;
theta = 2*pi/N;
do j = 1 to N;
   side = mod(j,2);
   t=.; x=0; y=0; output;
   do angle = pi/2-theta/2 to pi/2+theta/2 by theta/20;
      t = (j-1)*theta + angle;
      x = r*cos(t);
      y = r*sin(t);
      output;
   end;
end;
keep j side t x y;
run;

ods graphics / width=400px height=400px;
title "Cut a Circle into Slices";
proc sgplot data=Circle aspect=1 noautolegend;
   polygon x=x y=y ID=j / group=Side fill outline transparency=0.5;
run;

/* 2. Use PROC IML to arrange the slices into a row, alternating the
      orientation up, down, up, down, ...
*/
%let N = 8;
proc iml;
/* Enumerate slices as i=1..N. Draw odd slices with circular arc (the "crust")
   facing up and draw the even slices with arc on the bottom */
pi = constant('pi');
r = 1;
N = &N;          /* number of pieces to slice pizza */
theta = 2*pi/N;
angle = T( do(pi/2-theta/2, pi/2+theta/2, theta/20) );
sector = {0 0} // (r*cos(angle) || r*sin(angle));

/* open SAS data set 'pizza0' for writing */
slice = {. .}; ID = .; Side = .;
create pizza0 from slice ID Side[c={'px' 'py' 'ID' 'Side'}];

a = r/2*cos(theta/2);    /* reflection line */
b = r*sin(theta/2);      /* translation amount */
do j = 1 to N;           /* for each slice ... */
   /* translate to the right by a multiple of r*sin(theta/2) */
   slice = sector + ((j-1)*b || 0);
   if mod(j,2)=0 then   /* j is an even number */
      /* reflect across horizontal line y=a where a = r/2*cos(theta/2). 
         Recall reflection across y=a means that (x,y) --> (x,2*a-y) */
      slice[,2] = 2*a - slice[,2];
   ID = j(nrow(slice), 1, j);
   Side = j(nrow(slice), 1, mod(j,2));
   append from slice ID Side;
end;
close;
QUIT;

ods graphics / width=600px height=250px;
title "Rearranging Pizza Slices";
proc sgplot data=pizza0 aspect=0.333 noautolegend;
   polygon x=px y=py ID=ID / group=Side fill outline transparency=0.5;
   refline 3.1416 / axis=x noclip;
   xaxis grid label='C/2 (=pi*r)'
         values = (0 to 3.1416 by %sysevalf(3.1416/8)) valueshint
         valuesdisplay = ('0' 'pi/8' 'pi/4' '3pi/8' 'pi/2' '5pi/8' '3pi/4' '7pi/8' 'pi');
   yaxis label='r';
run;


/* 3. Generalize the program by cutting the first slice in half.
      Plot half on the left and move the other half to the right.
      If N is odd, you need to rotate the half slice on the right.
*/
%macro PizzaPi(N);
proc iml;
N = &N;          /* number of pieces to slice pizza */
/* Enumerate slices as i=1..N. Draw odd slices with circular arc (the "crust")
   facing up and draw the even slices with arc on the bottom */
pi = constant('pi');
r = 1;
theta = 2*pi/N;
angle = T( do(pi/2-theta/2, pi/2+theta/2, theta/20) );
sector = {0 0} // (r*cos(angle) || r*sin(angle));
halfangle = T( do(pi/2-theta/2, pi/2, theta/10) );   /* half slice */
halfsector = {0 0} // (r*cos(halfangle) || r*sin(halfangle));

/* plot a half slice on the left side */
j = 1;
slice = halfsector; 
ID = j(nrow(slice), 1, j); 
Side = j(nrow(slice), 1, mod(j,2));
create pizza from slice ID Side[c={'px' 'py' 'ID' 'Side'}];
append from slice ID Side;

/* plot the other slices, alternating up and down */
a = r/2*cos(theta/2);    /* reflection line */
b = r*sin(theta/2);      /* translation amount */
do j = 2 to N;
   /* translate to the right by multiple of r*sin(theta/2) */
   slice = sector + ((j-1)*b || 0);
   if mod(j,2)=0 then   /* j is an even number */
      /* reflect across horizontal line y=a where a = r/2*cos(theta/2). 
         Recall reflection across y=a means that (x,y) --> (x,2*a-y) */
      slice[,2] = 2*a - slice[,2];
   ID = j(nrow(slice), 1, j);
   Side = j(nrow(slice), 1, mod(j,2));
   append from slice ID Side;
end;

/* Put the other half slice on the far right side.
   Translate if N is even; flip and translate if N is odd */
isEven = (mod(N,2)=0);
slice = halfsector; 
slice[,1] = -slice[,1];          /* reflect across x=0 */
if isEven then                   /* j is used to assign the color */
   j = N+1;
else do;
   j = N+2;
   slice[,2] = 2*a - slice[,2];  /* reflect across y=a */
end;
slice = slice + (N*b || 0);      /* translate to the right */

ID = j(nrow(slice), 1, j); 
Side = j(nrow(slice), 1, mod(j,2));
append from slice ID Side;
close;
QUIT;

title2 "&N Slices";
proc sgplot data=pizza aspect=0.333 noautolegend;
   polygon x=px y=py ID=ID / group=Side fill outline transparency=0.5;
   refline 3.14 / axis=x noclip;
   xaxis grid max=3.1416  label='C/2 (=pi*r)'
         values = (0 to 3.1416 by %sysevalf(3.1416/8)) valueshint
         valuesdisplay = ('0' 'pi/8' 'pi/4' '3pi/8' 'pi/2' '5pi/8' '3pi/4' '7pi/8' 'pi');
   yaxis label='r';
run;
%mend;

/* 4. Create some pictures of pizza pi */
ods graphics / width=600px height=250px;
title "Rearranging a Sliced Pizza";
%PizzaPi(8);
%PizzaPi(16);
*%PizzaPi(17);

%PizzaPi(56);

/******************************/
/* BONUS: How quickly does the extent (width) of the base approach pi? */
proc iml;
pi = constant('pi');
N = T( 8:72 );
Range = N#sin(pi/N);
title "Width of Pizza Quasi-Rectangle";
call series (N, Range) grid={x y} other="refline 3.14159/axis=y noclip;";

Diff = pi - Range;
ne = nrow(N);
print (N[ne-5:ne])[L='N'] (Diff[ne-5:ne])[L='Diff'];
QUIT;
