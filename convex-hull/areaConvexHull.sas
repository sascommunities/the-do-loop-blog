/* SAS program to accompany the article 
   "The area and perimeter of a convex hull"
   by Rick Wicklin, published 02NOV2022 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2022/11/02/area-perimeter-convex-hull.html

   This program shows how to compute the area and perimeter of
   a polygon. The formulas are applied to the convex hull polygon
   for a set of points. An n-sided polygon is specified in 
   counterclockwise order as an (n x 2) matrix whose rows are
   the vertices of the polygon.
*/

proc iml;
/* The following functions are modified from 
   Rick Wicklin (2016). The Polygon package. SAS Global Forum 2016 */
/* PolyArea: Return the area of a simple polygon.
   P is an N x 2 matrix of (x,y) values for vertices. N > 2.
   This function uses Gauss' "shoelace formula," which is also 
   known as the "triangle formula."  
   See https://en.wikipedia.org/wiki/Shoelace_formula  */
start PolyArea(P);
   lagIdx = 2:nrow(P) || 1;   /* vertices of the lag; appends first vertex to the end of P */
   xi   = P[ ,1];       yi = P[ ,2];
   xip1 = xi[lagIdx]; yip1 = yi[lagIdx]; 
   area = 0.5 * sum(xi#yip1 - xip1#yi);
   return( area );
finish;

/* PolyPerimeter: Return the perimeter of a polygon.
   P is an N x 2 matrix of (x,y) values for vertices. N > 2.  */
start PolyPeri(_P);
   P = _P // _P[1,];     /* append first vertex */
   /* (x_{i+1} - x_i)##2 + (y_{i+1} - y_i)##2 */
   v = dif(P[,1])##2 + dif(P[,2])##2;  /* squared distance from i_th to (i+1)st vertex */
   peri = sum(sqrt(v));                /* sum of edge lengths */
   return( peri );
finish;
store module=(PolyArea PolyPeri);

/* demonstrate functions by using a 3:4:5 right triangle */
Pts = {1 1,
       4 1,
       4 5};
Area = PolyArea(Pts);
Peri = PolyPeri(Pts);
print Area Peri;

/* You can use the CVEXHULL function in SAS to compute 2-D
   convex hulls. See
   https://blogs.sas.com/content/iml/2021/11/15/2d-convex-hull-sas.html
*/
points = {0  2, 0.5 2, 1 2, 0.5 1, 0 0, 0.5 0, 1  0, 
          2 -1,   2 0, 2 1,   3 0, 4 1,   4 0, 4 -1, 
          5  2,   5 1, 5 0,   6 0 }; 
 
/* Find the convex hull: indices on the convex hull are positive */
Indices = cvexhull( points ); 
hullIdx = indices[loc(indices>0)];   /* the positive indices */
convexPoly = points[hullIdx, ];      /* extract rows to form a polygon */
print hullIdx convexPoly[c={'cx' 'cy'} L=""];

/* what is the area of the convex hull? */
AreaCH = PolyArea(convexPoly);
PeriCH = PolyPeri(convexPoly);
print AreaCH PeriCH;

quit;

