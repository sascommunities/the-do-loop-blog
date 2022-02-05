/* SAS program to accompany the article 
   "Billiards on a heart-shaped table"
   by Rick Wicklin, published 09FEB2022 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2022/02/05/billiards-heart-table.html

   This program builds on the program in 
   https://blogs.sas.com/content/iml/2022/02/07/billiards-elliptical-table.html
   which displays the path of a ball on a circular or elliptical billiard table.

   This program visualizes the path of a rebounding ball on a heart-shaped 
   billiard table.
*/

/* The program is in three parts:
   1. Define the helper functions that compute the path of a ball
   2. Define the function H that defines the heart-shaped table. Create a 
      scatter plot that visualizes the table.
   3. Compute the path of a ball that rebounds across the table. Overlay
      the path on the scatter plot.

   The program ignores the fact that the 
   contour H(x)=0 is not smooth at (0,+1) and (0,-1). Note that 
   grad(H)=(0 0) at those locations, so there is no normal vector and 
   the reflection is undefined.  The gradient also valishes at (+/-1, 0),
   although that is a mathematical artifact that can be handled.
*/


/* 1. Define helper functions for billiards on a table defined by H(x)=0.
      See https://blogs.sas.com/content/iml/2022/02/07/billiards-elliptical-table.html
*/
proc iml;
/* given x such that H(x)=0 and a unit vector, v, find the 
   normal and tangent components of v such that v = v_N + v_T, 
   where v_N is normal to H=0 and v_T is tangent to H=0 */
/* project v onto the normal space at x */
start NormalH(x, v);
   grad = gradH(x);      /* if H(x)=0, grad(H) is normal to H=0 at x */
   N = grad / norm(grad);/* unit normal vector to curve at x */
   v_N = (v*N`) * N;     /* projection of v onto normal vetor */
   return( v_N ); 
finish;
/* project v onto the tangent space at x */
start TangentH(x, v);     
   v_N = NormalH(x, v);
   v_T = v - v_N;
   return( v_T ); 
finish;
/* reflect v across tangent space at x */
start Reflect(x, v);
   v_N = NormalH(x, v);
   v_refl = v - 2*v_N;
   return( v_refl );
finish;
/* To find the intersection of H=0 with 
   the line through G_p in the direction of G_v, see
   https://blogs.sas.com/content/iml/2022/02/02/line-search-sas.html
*/
/* Restriction of H to the line through G_p in the direction of G_v */
start F_Line(t) global(G_p, G_v);
  return( H( G_p + t*G_v ) );
finish;
/* return the intersection of H=0 and the line */
start FindIntersection(p, v, step=20) global(G_p, G_v);
   G_p = p; G_v = v;
   t = froot("F_Line", 1E-3 || step);  /* find F(t) = 0 for t in [1E-3, step] */
   q = p + t*v;                        /* H(q) = 0 */
   return( q );
finish;
/* given a point, p, and a direction vector, v, find the point
   q on the line through p in the direction of v such that H(q)=0.
   Reflect the velocity vector to get a new direction vector, v_refl */
/* redefine the HitTable function to take an optional step length */
start HitTable(q, v_refl, p, v, step=20);
   q = FindIntersection(p, v, step);
   v_refl = Reflect(q, v);
finish;
store module=(NormalH TangentH Reflect F_Line FindIntersection HitTable);
QUIT;

/* 2. Create points in the interior H < 0. */
data Heart(rename=(x=hx y=hy));
do y = -1.0 to 1.3 by 0.01;
   do x = -1.2 to 1.2 by 0.01;
      z = (x**2 + y**2 - 1)**3 - x**2 * y**3;
      if z < 0 then output;
   end;
end;
run;
/*
title "Region where H(x,y) <= 0";
proc sgplot data=Heart;
   scatter x=hx y=hy / markerattrs=(symbol=SquareFilled color=lightred) ;
run;
*/


/* 3. Define the heart-shaped region and gradient vector. Use it to 
      compute the path of a ball as it rebounds across the table. */
proc iml;
load module=(NormalH TangentH Reflect F_Line FindIntersection HitTable);

start H(_x); 
   x = _x[,1]; y = _x[,2];
   z = (x##2 + y##2 - 1)##3 - x##2 # y##3;
   return( z );
finish;
start gradH(_x);
   x = _x[,1]; y = _x[,2];
   dHdx = 6*x #(x##2 + y##2 - 1)##2 - 2*x # y##3;
   dHdy = 6*y #(x##2 + y##2 - 1)##2 - 3*x##2 # y##2;
   return( dHdx || dHdy ); 
finish;

/* computes the trajectory of the ball for the first k bounces,
   starting from position p and moving in the direction of v */
start Billiards(p, _v, NSeg);  
   x = j(NSeg, ncol(p));
   x[1,] = p;
   v = _v;
   FirstTime=1;
   do i = 2 to NSeg;
      tFinal = 3;
      /* NOTE: The root-finding algorithm might jump across the gap near (0, 1).
         To prevent this, find the FIRST value of t for which x_i + t*v is 
         outside of the heart-shaped region. */
      if x[i-1,2] > 0.8 then do;  /* happens mostly when p is near the top */
         tt = T(do(1e-3, 3, 0.1));
         pp = x[i-1,] + tt @ v;
         hh = H( pp );
         posIdx = loc(hh > 0);
         if ncol(posIdx)>0 then 
            tFinal = tt[posIdx[1]];
      end;
      run HitTable( q, v, x[i-1,], v, tFinal ); 
      x[i,] = q;
   end;
   /* write the trajectory to a data set and overlay trajectory on a graph */
   Bounce = T(1:NSeg);
   create Bounce from Bounce x[c={'Bounce' 'x1' 'x2'}];
   append from Bounce x;
   close;
   submit;
      data All;
      set Bounce Heart;
      run;
      proc sgplot data=All noautolegend nowall;
         scatter x=hx y=hy / markerattrs=(symbol=SquareFilled color=lightpink size=5);
         series x=x1 y=x2  / lineattrs=(color=lightred);
         xaxis display=none;
         yaxis display=none;
      run;
   endsubmit; 
finish;

p = {0  0};            /* point p */
v = {1  1}/sqrt(2);    /* vector v */

ods graphics / width=480px height=480px antialias=off;
title "Love Makes My Heart Bounce!";
run Billiards(p, v, 250);

QUIT;
