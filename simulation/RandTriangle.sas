/* SAS program to accompany the article 
   "Generate random points in a triangle"
   by Rick Wicklin, published 19Oct2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/10/19/random-points-in-triangle.html

   This program shows how to generate uniform random points
   inside an arbitrary triangular region of the plane.
*/


/* Optional: You can use the Polygon package to improve the 
   visualization.

   If you want to use the Polygon package, install it once. 
   For information about how to install and use packages, see
   https://blogs.sas.com/content/iml/2016/04/27/packages-share-sas-iml-programs-html.html   
   You can download the Polygon package from
   https://communities.sas.com/t5/SAS-IML-File-Exchange/The-polygon-Package/ta-p/263210

*INSTALL the package;
proc iml;
package install "C:\Packages\polygon.zip";  *put location here;
quit;
*/

proc iml;
package load polygon;

/* How to generate random points in a parallelogram */
n = 1000;
call randseed(1234,1);

a = {3  2};                       /* vector along one side */
b = {1 -2};                       /* vector along adjacent side */
u = randfun(n // 2,  "Uniform");  /* u[,1], u[,2] ~ U(0,1) */
w = u[,1]@a + u[,2]@b;            /* linear combination of a and b */ 
title "Random Points in Parallelogram";
call scatter(w[,1], w[,2]) grid={x,y};

R = {0 0} // a // (a+b) // b;
call PolyDraw(R, w) grid={x y} transparency=0.8;


/* Generate random uniform sample in triangle with vertices
   P1 = (x0,y0), P2 = (x1,y1), and P3 = (x2,y2)
   The triangle is specified as a 3x2 matrix, where each row is a vertex.
*/
start randUnifTriangle(P, n);
   a = P[2,] - P[1,];         /* translate triangle to origin */
   b = P[3,] - P[1,];         /* a and b are vectors at the origin */
   u = randfun(n // 2,  "Uniform");
   idx = loc(u[,+] >= 1);
   if ncol(idx)>0 then 
      u[idx,] = 1 - u[idx,];  /* transform variates */
   w = u[,1]@a + u[,2]@b;     /* linear combination of a and b vectors */ 
   return( P[1,] + w );       /* translate triangle back to original position */
finish;
store module=(randUnifTriangle);

/* triangle contains three vertices */
call randseed(1234,1);
P = {1 2,      /* P1 */
     4 4,      /* P2 */
     2 0};     /* P3 */
n = 1000;
w = randUnifTriangle(P, n);

title "Random Points in Triangle";
ods graphics / width=480px height=480px;
call scatter(w[,1], w[,2]) grid={x,y};

call PolyDraw(P, w) grid={x y} transparency=0.8;


*-------------------;
/* Optional: Illustrate the reflection process that 
   turns points in a parallelogram into points in a triangle.
*/
call randseed(1234,1);
P = {1 2 1,      /* A */
     4 4 1,      /* B */
     2 0 1,      /* C */
     4 4 2,      /* reflection */
     2 0 2,
     5 2 2};

   a = P[2,] - P[1,];         /* translate triangle to origin */
   b = P[3,] - P[1,];         /* a and b are vectors at the origin */
   u = randfun(n // 2,  "Uniform");
   z = 1 + (u[,+] >= 1);
   w = u[,1]@a + u[,2]@b;     /* linear combination of a and b vectors */ 
   v = (P[1,] + w)[,1:2]  || z; 
   call sort(v, 3);

ods graphics / width=640px height=480px;
title "Method of Reflection";
title2 "Generate Random Points in Parallelogram";
call PolyDraw(P, v) grid={x y} transparency=0.8;

