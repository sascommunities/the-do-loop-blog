/* SAS program to accompany the article 
   "Generate random points in a polygon"
   by Rick Wicklin, published 21Oct2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/10/21/random-points-in-polygon.html

   This program shows how to generate uniform random points
   inside an arbitrary polygon in the plane.

   Every polygon can be decomposed into a set of triangles.
   You can then generate random points in each triangle by using
   the RandUnifTriangle function, which is defined at
   https://blogs.sas.com/content/iml/2020/10/19/random-points-in-triangle.html

   The number of points in the i_th triangle 
   should be proportional to
   Area(i_th_triangle) / Area(Polygon).

   It is easy to decompose a convex polygon into triangles:
   Choose a base vertex (P1) and connect it to consecutive
   pairs of other vertices:
   P1 P2 P3
   P1 P3 P4
   P1 P4 P5, etc
*/


/* You can use the Polygon package to improve the visualization.

   To use the Polygon package, install it once. 
   For information about how to install and use packages, see
   https://blogs.sas.com/content/iml/2016/04/27/packages-share-sas-iml-programs-html.html   
   You can download the Polygon package from
   https://communities.sas.com/t5/SAS-IML-File-Exchange/The-polygon-Package/ta-p/263210

*INSTALL the package;
proc iml;
package install "C:\Packages\polygon.zip";  *put location here;
quit;
    
   Alternately, you can just %INCLUDE the source code 
   for the modules in the package.
*/

/* assume the polygon package is installed */
proc iml;
package load polygon;     /* load the polygon package */

/* Decompose a convex polygon into triangles. Return a list
   that contains the vertices for the triangles.
   This function uses a function in the Polygon package, which must be loaded.
*/
start TriangulateConvex(P);            /* input parameter(N x 2): vertices of polygon */
   isConvex = PolyIsConvex(P);
   if ^isConvex then
      return ( [] );                 /* The polygon is not convex */
   numTri = nrow(P) - 2;             /* number of triangles in convex polygon */
   L = ListCreate(numTri);           /* create list to store triangles */
   idx = 2:3;
   do i = 1 to ListLen(L);
      L$i = P[1,] // P[idx,];
      idx = idx + 1;
   end;
   return (L);
finish;

/* Specify a convex polygon and visualize the triangulation.  */
P = { 2 1 ,
      3 1 ,
      4 2 ,
      5 4 ,
      3 6 ,
      1 4 ,
      1 2 };
L = TriangulateConvex(P);

/* visualize the decomposition */
Len = ListLen(L);
Poly = j(3*Len, 3);
j = 1:3;
do i = 1 to Len;
   Poly[j,] = (L$i) || j(3,1,i);
   j = j + 3;
end;
ods graphics / width=400px height=500px;
title "Triangulation of a Convex Polygon";
run PolyDraw(Poly);

/* Given a list of triangles (L), generate N random points in the union,
   where the number of points is proportional to 
   area(triangle) / area(all triangles)
   
   This function uses functions in the Polygon package, which must be loaded.
*/
start RandUnifManyTriangles(N, L);
   numTri = ListLen(L);
   /* compute areas of each triangle in the list */
   AreaTri = j(1, numTri,.);         /* create vector to store areas */
   do i = 1 to numTri;
      AreaTri[i] = PolyArea(L$i);    /* PolyArea is in the Polygon package */
   end;
   /* Numbers of points in the triangles are multinomial with
      probability proportional to Area(triangle)/Area(polygon)
   */
   NTri = RandMultinomial(1, N, AreaTri/sum(AreaTri));
   cumulN = 0 || cusum(NTri);        /* cumulative counts; use as indices */
   z = j(N, 3, .);                   /* columns are (x,y,TriangleID) */
   do i = 1 to numTri;
      k = (cumulN[i]+1):cumulN[i+1]; /* the next NTri[i] elements */
      z[k, 1:2] = RandUnifTriangle(L$i, NTri[i]);
      z[k, 3] = i;                   /* store the triangle ID */
   end;
   return z;
finish;

/* The RandUnifTriangle function is defined at
   https://blogs.sas.com/content/iml/2020/10/19/random-points-in-triangle.html
*/
load module=(RandUnifTriangle);
call randseed(12345);
N = 2000;
z = RandUnifManyTriangles(N, L);
title "Random Points in a Polygon";
title2 "Colors Assigned Based on Triangulation";
call PolyDraw(P, z);
