/* SAS program to accompany the article 
   "Bilinear interpolation in SAS"
   by Rick Wicklin, published 20MAY2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/05/20/bilinear-interpolation-sas.html

   This program shows how to perform bilinear interpolation. 
   The main function is the bilinearInterp functio in SAS/IML.
   The function takes the following input parameters:

   xgrd : vector of Nx points along x axis
   ygrd : vector of Ny points along y axis
   z    : Ny x Nx matrix. The value Z[i,j] is the value at xg[j] and yg[i]
   t    : k x 2 matrix of points at which to bilinearly interpolate

   The function returns a k x 1 vector, which is the bilinear 
   interpolation at each row of t. The function uses the ideas in 
   the article 
   "What is bilinear interpolation?" by Rick Wicklin;
   https://blogs.sas.com/content/iml/2020/05/14/what-is-bilinear-interpolation.html

*/
title;
ods graphics / width=480px height=400px;

/* BEFORE PROCEDING: DEFINE AND STORE THE blinearInterp FUNCTION */
proc iml;
start bilinearInterp(_xgrd, _ygrd, z, _t);
   xgrd = rowvec(_xgrd);
   ygrd = colvec(_ygrd);
   /* 1. Check for valid parameters */
   /* verify that z is Nx x Ny matrix */
   if any(dimension(z)^= nrow(ygrd)||ncol(xgrd)) then 
      stop "ERROR: Dimension of z is not congruent with dimension of XGRD and YGRD.";
   /* verify that t has two columns */
   if ncol(_t) ^= 2 then 
      stop "ERROR: the T matrix must have two columns.";

   G = j(nrow(_t), 1, .);   /* allocate the result vector */
   /* 2. Exclude any points not in the data range */
   /* BIN returns missing if a value is not in range of grid points */
   ktx = bin(_t[,1], xgrd);            /* which bin each point is in */
   kty = bin(_t[,2], ygrd);
   validIdx = loc(ktx^=. & kty^=.); 
   if ncol(validIdx)= 0 then 
      return( G );  /* no points are in the data range */
   /* restrict to valid points; don't need to use validIdx again until the end */
   t = _t[validIdx,]; 
   kx = ktx[validIdx]; ky = kty[validIdx];  

   /* 3. Each point is in a cell. 
         Get the (x,y) values at the four corners of the cell it is in */
   xL = xgrd[kx];   yL = ygrd[ky];  
   xR = xgrd[kx+1]; yR = ygrd[ky+1];

   /* 4. get the function value at the four corners of the cell for each point */
   zT = T(z);  /* easier indexing if x corresponds to rows and y corresponds to cols */
   idx00 = sub2ndx( dimension(zt), kx    ||ky );
   idx01 = sub2ndx( dimension(zt), (kx+1)||ky );
   idx10 = sub2ndx( dimension(zt), kx    ||(ky+1) );
   idx11 = sub2ndx( dimension(zt), (kx+1)||(ky+1) );
   z00 = zt[idx00];
   z01 = zt[idx01];
   z10 = zt[idx10];
   z11 = zt[idx11];

   /* 5. Bilinear interpolation */
   /* Map each t into fraction of cell coordinates (transform to unit square) */
   s = ((t[,1] - xL) / (xR - xL)) || ((t[,2] - yL) / (yR - yL));
   sx = s[,2];  sy = s[,1];  cx = 1 - sx;
   /* interpolate as if in the unit square */
   F = (z00#cx + z10#sx)#(1-sy) + (z01#cx + z11#sx)#sy;
   G[validIdx] = F;
   return( G );
finish;
store module=bilinearInterp;       /* store the module so it can be loaded later */
QUIT;


/* The fitting data. The response Z is observed on a regular grid.
   These data must be sorted byY then by X
*/
data Have;
input x y z;
datalines;
1 0   6 
2 0   7 
5 0   5 
1 1   9 
2 1   9 
5 1   7 
1 1.5 10 
2 1.5 9 
5 1.5 6 
1 3   12 
2 3   9 
5 3   8 
;

/* If the data are not already sorted, sort by Y and X */
proc sort data=Have; by Y X; run; 

/* MAKE SURE you have ALREADY STORED the bilinearInterp function! */
proc iml;
/* read the fitting data */
use Have; read all var {x y z}; close;
xGrd = unique(x);                      /* find grid points */                     
yGrd = unique(y);
z = shape(z, ncol(yGrd), ncol(xGrd)); /* data must be sorted by y, then x */
print z[r=(char(yGrd)) c=(char(xGrd))];

/* read or otherwise define the scoring data */
t = {0   0,     /* not in data range */
     1   1,
     1   2,
     2.5 2,
     4   0.5,
     4   2,
     5   3,
     6   3};    /* not in data range */

/* LOAD the previously stored function */
load module=bilinearInterp;  
F = bilinearInterp(xGrd, yGrd, z, t);
print t[c={'x' 'y'}] F[format=Best4.];

/* END of bilinear interpolation example */

/* The following is for visualization: 
   you can interpolate on a fine grid of scoring locations */
tx = do(1, 5, 0.1);
ty = do(0, 3, 0.1);
yx = ExpandGrid( ty, tx );  /* vizualization trick: let x change fastest */
t = yx[,{2 1}];
F = bilinearInterp(xGrd, yGrd, z, t);

/* Visualization Options: 
   1. Heat map in SAs/IML
   2. Write data set and use GTL or PROC SGPLOT 
*/

/* 1. Visualize by using HEATMAPCONT function in SAS/IML */
/* Reshape from long to wide */
Q = shape(F, ncol(ty), ncol(tx));   /* reshape in
/* currently, the Y axis is pointing down. Flip it and the labels. */
H = Q[ nrow(Q):1, ];
reverseY = ty[ ,ncol(ty):1];
/* visualize by using a heat map */
call heatmapcont(H) displayoutlines=0
     xValues=tx yValues=reverseY
     colorramp=palette("Spectral", 7)[,7:1]
     title="Bilinear Interpolation";

/* 2. Visualize by using GTL and PROC SGPLOT */
/* 2a. Write to SAS data set */
Out = t || F;
create BLInterp from Out[c={'x' 'y' 'z'}];
append from Out;
close;
QUIT;

/* 2b. Define surface plot by using GTL */
proc template;
  define statgraph surfparm;
   begingraph;
     entrytitle "Bilinear Interpolation";
     layout overlay3d;
       surfaceplotparm x=x y=y z=z / colorresponse=z
         colormodel=(CX3288BD CX99D594 CXE6F598 CXFFFFBF CXFEE08B CXFC8D59 CXD53E4F) /* palette("Spectral", 7)[,7:1] */
         surfacetype=fillgrid;
     endlayout;
   endgraph;
  end;
run;

proc sgrender data=BLInterp template=surfparm;
run;

/* 2c. Use HEATMAPPARM statement in PROC SGPLOT. 
       Optionally, you can overlay the location of the grid 
       points that define the interpolant.
*/
data Combine;
set Have(rename=(x=xx y=yy)) BLInterp;
run;

title "Bilinear Interpolation";
title2 "Interpolation from 12 Points";
proc sgplot data=Combine;
   heatmapparm x=x y=y colorresponse=z / name="HM"
       colormodel=(CX3288BD CX99D594 CXE6F598 CXFFFFBF CXFEE08B CXFC8D59 CXD53E4F);
   refline 1 2 5 / axis=x;
   refline 0 1 1.5 3 / axis=y;
   scatter x=xx y=yy / colorresponse=z markerattrs=(symbol=CircleFilled size=15)
                  filledoutlinedmarkers
                  colormodel=(CX3288BD CX99D594 CXE6F598 CXFFFFBF CXFEE08B CXFC8D59 CXD53E4F);
   xaxis label="x";
   yaxis label="y";
   gradlegend "HM" / title="z";
run;



