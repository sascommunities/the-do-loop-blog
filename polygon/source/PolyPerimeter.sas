/* Rick Wicklin (2016). The Polygon package. SAS Global Forum 2016 */
/* PolyPerimeter: Return the perimeter of a polygon.
   Input: 
   P        an N x 2 matrix of (x,y) values for vertices. N > 2.
            Or an N x 3 matrix where the third column is an ID variable
            that identifies different polygons.
*/
/* Perimeter of Nx2 matrix of vertices */
start _PolyPerimeter(_P);
   P = _P // _P[1,];     /* append first vertex */
   /* (x_{i+1} - x_i)##2 + (y_{i+1} - y_i)##2 */
   v = dif(P[,1])##2 + dif(P[,2])##2;  
   peri = sum(sqrt(v));
   return( peri );
finish;

start PolyPerimeter(P);
   if ncol(P)=2 then
      return( _PolyPerimeter(P) );
      
   ID = P[,3];
   u = uniqueby(ID);         /* starting index for each group */
   result = j(nrow(u), 1);   /* allocate vector to hold results */
   u = u // (nrow(ID)+1);    /* append (N+1) to end of indices */
   do i = 1 to nrow(u)-1;    /* for each group... */
      idx = u[i]:(u[i+1]-1); /* get rows in group */
      result[i] = _PolyPerimeter( P[idx, 1:2] );
   end;
   return( result );
finish;
store module=(_PolyPerimeter PolyPerimeter);
