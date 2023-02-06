/* SAS program to accompany the article 
   "A random walk inside a heart"
   by Rick Wicklin, published 08FEB2022 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/02/08/random-walk-heart.html
   
   This program visualizes the path of a random walk that stays inside 
   a heart-shaped region. From a poin, x, inside the region, attempt
   to step to x+v, where v is a random bivariate normal vector.
   If x+v is outside the region, truncate the step and step onto the boundary
   of the region.  Repeat many times.

   This program builds on the program in 
   https://blogs.sas.com/content/iml/2022/02/09/billiards-heart-table.html
   which displays the path of a ball on a heart-shaped billiards table.
*/

/*******************************************/
/* Show a random 2-D walk where each step is a vector v, 
   which is bivariate normal v ~ MVN(0, 0.2*I(2)) */
data RandomWalk;
N = 2000;
call streaminit(12345);
x = 0; y = 0;
do i = 1 to N;
   dx = rand("Normal", 0, 0.2);
   dy = rand("Normal", 0, 0.2);
   x = x + dx;
   y = y + dy;
   output;
end;
keep x y;
run;
ods graphics / reset width=480px height=480px;
title "Random Walk with Bivariate Normal Steps";
proc sgplot data=RandomWalk;
   series x=x y=y / lineattrs=(color=lightred);
   xaxis grid;
   yaxis grid;
run;

/* or you can do the same computation in SAS IML */
ods graphics / reset width=480px height=480px;
proc iml;
N = 2000;
call randseed(12345);
vel = randfun(N//2, "Normal", 0, 0.2);  /* N two-dimensional vectors */

x = j(N+1, 2, 0);
do i = 1 to N;
   x[i+1,] = x[i,] + vel[i,];
end;

title "Random Walk with Bivariate Normal Steps";
y = x[,2]; x= x[,1]; 
call series(x, y) option="lineattrs=(color=lightred)" grid={x y};
QUIT;

/*******************************************/

/* Instead of an unconstrained walk, at each step check whether the 
   new location is outside of a region. If so, truncate the step to 
   land on the boundary of the region. */

proc iml;
/* define the heart-shaped region */
start H(xy);
  x = xy[,1]; y = xy[,2];
  return ( (x**2 + y**2 - 1)**3 - x**2 * y**3 );
finish;

/* return 1 if (x,y) is inside the region */
start InRegion(xy, tol=1e-14);
  return (H(xy) <= tol);
finish;

/* given a point, x, and a vector, v, this function returns the function
   f(t) = H(x + t*v), t in [0,1]. If f has a root, then the line segment from 
   x to x+v intersects the boundary of the reqion. 
   This function is used by FROOT to find points on the boundary of the region. */
start OnBoundary(t) global(G_x, G_v);
  return ( H( G_x + t*G_v ) );
finish;

/* Start at x. Try to step to x+v.
   If x+v is in the region, take the step.
   Otherwise, see if we can take part of the step to land on the boundary. */
start StepInRegion(x, v) global(G_x, G_v);
   if InRegion(x + v) then 
      return (x + v);
   if InRegion(x) then do;
      G_x = x; G_v = v;
      /* does the line from x to x+v cross the boundary? */
      t = froot("OnBoundary", {0 1});  /* find intersection of line and region */
      if t > 0 then 
         return (x + t*v);       /* step onto the boundary */
      else
         return (x);             /* cannot step in this direction */
   end;
   /* something is wrong: x is not in the region */
   return ( {. .} );
finish;

N = 2000;
call randseed(12345);
vel = randfun(N//2, "Normal", 0, 0.2);

x = j(1, 2, 0);
create walk from x [c={x y}];
do i = 1 to N;
   x = StepInRegion(x, vel[i,]);
   append from x;
end;
close;
QUIT;

/* visualize only the constrained random walk */
ods graphics / width=480px height=480px;
title "Random Walk Inside a Heart-Shaped Region";
proc sgplot data=Walk;
   series x=x y=y / lineattrs=(color=lightred);
run;


/* add a pink heart in the background to show the constraint region */
data Heart(rename=(x=hx y=hy));
do y = -1.0 to 1.3 by 0.01;
   do x = -1.2 to 1.2 by 0.01;
      z = (x**2 + y**2 - 1)**3 - x**2 * y**3;
      if z < 0 then output;
   end;
end;
run;

data All;
set Walk Heart;
run;


ods graphics / width=480px height=480px;
title "Random Walk Inside a Heart-Shaped Region";
proc sgplot data=All noautolegend nowall;
   scatter x=hx y=hy / markerattrs=(symbol=SquareFilled color=lightpink size=5);
   series x=x y=y  / lineattrs=(color=lightred);
   xaxis display=none min=-1.2 max=1.2;
   yaxis display=none min=-1 max=1.3;
run;


/* Show the evolution of the random walk for N=100, 300, ... */
%macro plotIters();
%do N = 100 %to 700 %by 200;

data All;
set Walk(obs=&N) Heart;
run;
title "Trajectory After N = &N Iterations";
proc sgplot data=All noautolegend nowall;
   scatter x=hx y=hy / markerattrs=(symbol=SquareFilled color=lightpink size=5);
   series x=x y=y / lineattrs=(color=lightred);
   xaxis display=none min=-1.2 max=1.2;
   yaxis display=none min=-1 max=1.3;
run;
%end;
%mend;
ods graphics / width=240px height=240px;
ods layout gridded columns=2 advance=table
      column_gutter=10px row_gutter=0px;
%plotIters;
ods layout end;
