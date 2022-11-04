/* Rick Wicklin (2016). The Polygon package. SAS Global Forum 2016 */
/* PolyArea: Return the area of a simple polygon.
   Input: 
   P     an N x 2 matrix of (x,y) values for vertices. N > 2.
         Or an N x 3 matrix where the third column is an ID variable
         that identifies different polygons.
*/
/* Area of Nx2 matrix of vertices */
start _PolyArea(P, useAbsVal=1);
   lagIdx = 2:nrow(P) || 1;
   xi   = P[ ,1];       yi = P[ ,2];
   xip1 = xi[lagIdx]; yip1 = yi[lagIdx]; 
   area = 0.5 * sum(xi#yip1 - xip1#yi);
   if useAbsVal then 
      return( abs(area) );
   else 
      return( area );
finish;

start PolyArea(P);
   if ncol(P)=2 then
      return( _PolyArea(P) );

   ID = P[,3];
   u = uniqueby(ID);         /* starting index for each group */
   result = j(nrow(u), 1);   /* allocate vector to hold results */
   u = u // (nrow(ID)+1);    /* append (N+1) to end of indices */
   do i = 1 to nrow(u)-1;    /* for each group... */
      idx = u[i]:(u[i+1]-1); /* get rows in group */
      result[i] = _PolyArea( P[idx, 1:2] );
   end;
   return( result );
finish;
store module=(_PolyArea PolyArea);
