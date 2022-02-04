/* SAS program to accompany the article 
   "Billiards on an elliptical table"
   by Rick Wicklin, published 07FEB2022 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2022/02/07/billiards-elliptical-table.html

   This program shows how to compute the normal vector to 
   a curve defined by H(x) = 0. Given a point, p, and a 
   velocity vector, v, the program computes the point, q,
   where a ball starting at p will hit the curve H=0.  
   It also computes the reflection of v at the point of impact
   by reflecting v across the normal vector to the curve at q.

   An application of this math is that you can display the 
   path of a ball on a circular or elliptical billiard table.
   This article visualizes the paths of rebounding balls 
   on a nonstandard billiard table.
*/


/* define a circular table */
%let a = 1;
%let b = 1;
data Table;
P=1;
do t = 0 to 2*constant('pi') by 0.06;
   cx = &a*cos(t); cy = &b*sin(t);
   output;
end;
if (&a >= &b) then do; /* output the foci of the ellipse */
   ff = sqrt(&a**2 - &b**2);
   fx =  ff; fy = 0; output;
   fx = -ff; fy = 0; output;
end;
else if (&a < &b) then do;
   ff = sqrt(&b**2 - &a**2);
   fx = 0; fy =  ff; output;
   fx = 0; fy = -ff; output;
end;
run;

/* tangent and normal vectors; reflection of a vector across a boundary */


proc iml;
/* Define the multivariate function H and its gradient.
   Example: H(x)=0 defines a circle and grad(H) is the gradient */
start H(x);               
   return( x[ ,1]##2 + x[ ,2]##2 - 1 );  /* = x##2 + y##2 - 1 */
finish;
start gradH(x);
   return( 2*x );                        /* { dH/dx  dH/dy } */
finish;

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

/* Test the definitions for the following points and direction vector */
p = {0.5 -0.5};           /* initial position of ball */
v = {1 1} / sqrt(2);      /* unit velocity vector: sqrt(2) = 0.7071068 */
q = {1 0};                /* ball hits circle at this point */
N = NormalH(q, v);        /* normal component to curve at q */
T = TangentH(q, v);       /* tangential component to curve at q */
v_refl = Reflect(q, v);   /* reflection of v */
print v, N, T, v_refl;


/* Let's play a game of billiards on a circular table. 
   Start with a point (the ball) and a direction vector.
   See the path of the ball as it reflects around the table. */
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
start HitTable(q, v_refl, p, v);
   q = FindIntersection(p, v);
   v_refl = Reflect(q, v);
finish;

/* Example 1: A periodic path with four segments */
p = {0.5 -0.5};     /* point p */
v = {1 1}/sqrt(2);  /* direction vector v */

NSeg = 6;           /* follow for this many bounces */
x = j(NSeg, ncol(p));
x[1,] = p;
do i = 2 to NSeg;
   run HitTable( q, v, x[i-1,], v ); /* overwrite v on each iteration */
   x[i,] = q;
end;
print x[L="Positions"];

/* computes the trajectory of the ball for the first k bounces,
   starting from position p and moving in the direction of v */
start Billiards(p, _v, NSeg);  
   x = j(NSeg, ncol(p));
   x[1,] = p;
   v = _v;
   do i = 2 to NSeg;
      run HitTable( q, v, x[i-1,], v ); 
      x[i,] = q;
   end;
   /* write the trajectory to a data set and graph */
   Bounce = T(1:NSeg);
   create Bounce from Bounce x[c={'Bounce' 'x1' 'x2'}];
   append from Bounce x;
   close;
   submit;
      data All;
      set Bounce Table;
      run;
      proc sgplot data=All noautolegend;
         polygon x=cx y=cy ID=P;
         series x=x1 y=x2;
         scatter x=fx y=fy / markerattrs=(symbol=CircleFilled size=12 color=red);
         xaxis display=(nolabel) grid;
         yaxis display=(nolabel) grid;
      run;
   endsubmit; 
finish;

ods graphics / width=480px height=480px antialias=off;

/* Example 2. Either periodic orbit or bounces inside an annulus */
p = {0.5 -0.5};     /* point p */
v = {1 1.4};        /* vector v */
v = v / norm(v);    /* standardize to unit vector */

title "50 Bounces on a Round Table";
title2 "p=(0.5, -0.5); v=(0.58, 0.81)";
run Billiards(p, v, 50);  /* NSeg=50 bounces */

/***************************************/

/* if I did everything correctly, all I have to do is change the shape of the 
   table to get new behavior on an elliptical table */
submit;
   /* define an elliptical table */
   %let a = sqrt(5);
   %let b = 2;
   data Table;
   P=1;
   do t = 0 to 2*constant('pi') by 0.06;
      cx = &a*cos(t); cy = &b*sin(t);
      output;
   end;
   if (&a >= &b) then do; /* output the foci of the ellipse */
      ff = sqrt(&a**2 - &b**2);
      fx =  ff; fy = 0; output;
      fx = -ff; fy = 0; output;
   end;
   else if (&a < &b) then do;
      ff = sqrt(&b**2 - &a**2);
      fx = 0; fy =  ff; output;
      fx = 0; fy = -ff; output;
   end;
   run;
endsubmit;

/* define the elliptical region and gradient vector */
start H(x); 
   a = sqrt(5); b = 2;
   return( (x[ ,1]/a)##2 + (x[ ,2]/b)##2 - 1 ); /* ellipse(a,b) */
finish;
start gradH(x);
   a = sqrt(5); b = 2;
   return( 2*x[ ,1]/a##2 || 2*x[ ,2]/b##2 ); 
finish;

/* ellipse has foci at (+/-sqrt(a##2-b##2), 0) 
   For a=sqrt(5) and b=2, foci at (+/-1, 0)
*/
/* There are two types of trajectories: https://maths.ucd.ie/~plynch/Talks/BandB.pdf
   Type 1: Annular region
   Type 2: Hyperbolic region
*/

/* If the ball does not pass between the foci, its path is inside an annular region. */
p = {1.5 0};        /* point p */
v = {0   1};        /* vector v */
title "250 Bounces on an Elliptical Table";
title2 "Ball Does Not Pass Between Foci";
run Billiards(p, v, 250);

/* If the ball passes between the foci, its path is inside a hyperbolic region. */
p = {0  0};         /* point p */
v = {1  1}/sqrt(2); /* vector v */
title "250 Bounces on an Elliptical Table";
title2 "Ball Passes Between Foci";
run Billiards(p, v, 250);

/* If the ball passes over one focal point, it also passes over the second focal point. 
   The path converges to the major axis that contains the foci. */
p = {1 -1};         /* point p */
v = {0  1};         /* vector v */
title "10 Bounces on an Elliptical Table";
title2 "Ball Passes Through a Focal Point";
run Billiards(p, v, 10);

QUIT;

