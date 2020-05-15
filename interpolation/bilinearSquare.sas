/* SAS program to accompany the article 
   "What is bilinear interpolation?"
   by Rick Wicklin, published 18MAY2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/05/14/what-is-bilinear-interpolation.html

   This program shows how to perform bilinear interpolation on the 
   unit square.

   Assume you know the values 
   f00 = f(0,0), f01=f(0,1), f10=f(1,0), and f11(1,1)
   at the corners of the square.
   Then if (x,y) is in the square,
   f(x,y) ~ f00*(1-x)*(1-y) +
            f10*x*(1-y) + 
            f01*(1-x)*y +
            f11*x*y
*/

/* 1. Visualize the bilinear interpolant for 
      z00 = 0; z10 = 4; z01 = 2; z11 = 1; 
*/
title;
ods graphics / width=480px height=400px;

/* Use GTL to create a template for a contour plot:
   https://blogs.sas.com/content/iml/2012/07/02/create-a-contour-plot-in-sas.html
*/
proc template;
define statgraph ContourPlotParm;
dynamic _X _Y _Z _TITLE;
begingraph;
   entrytitle _TITLE;
   layout overlay;
      contourplotparm x=_X y=_Y z=_Z /
        contourtype=fill nhint=21 name="Contour"
            colormodel=(CX3288BD CX99D594 CXE6F598 CXFFFFBF CXFEE08B CXFC8D59 CXD53E4F) /* palette("Spectral", 7)[,7:1] */
            contourtype=LABELEDLINEGRADIENT;
      continuouslegend "Contour" / title=_Z;
   endlayout;
endgraph;
end;
run;

data Func;
z00 = 0; z10 = 4; z01 = 2; z11 = 1; 
do y = 0 to 1 by 0.01;
   do x = 0 to 1 by 0.01;
      z = z00*(1-x)*(1-y) +
          z10*    x*(1-y) + 
          z01*(1-x)*y +
          z11*    x*y;
      output;
   end;
end;
run;

proc sgrender data=Func template=ContourPlotParm;
dynamic _TITLE="Graph of Interpolant on Unit Square"
        _X="x" _Y="y" _Z="z";
run;

/* 2. Create a SAS/IML function for bilinear interpolation on the unit square */
/* Corners: f00 = f(0,0), f01=f(0,1), f10=f(1,0), and f11(1,1)
   If (x,y) is in the square,
   f(x,y) ~ f00*(1-x)*(1-y) +
            f10*x*(1-y) + 
            f01*(1-x)*y +
            f11*x*y
*/
proc iml;
start bilinearInterpSquare(z, t);
   z00 = z[1]; z10 = z[2]; z01 = z[3]; z11 = z[4]; 
   x = t[,1];  y = t[,2];
   cx = 1-x;                    /* cx = complement of x */
   return( (z00*cx + z10*x)#(1-y) + 
           (z01*cx + z11*x)#y );
finish;

/* specify the fitting data */
z = {0 4,                       /* bottom: values at (0,0) and (1,0) */
     2 1};                      /* top: values at (0,1) and (1,1) */
print z[c={'x=0' 'x=1'} r={'y=0' 'y=1'}];

/* test on some specific points in the square */
t = {0.5 0,                     /* specify the scoring data */
     0    0.25,
     1    0.66666666,
     0.5  0.75,
     0.75 0.5     };
F = bilinearInterpSquare(z, t); /* test the bilinear interpolation function */
print t[c={'x' 'y'} format=FRACT.] F;

/* Test on a regular grid of points */
/* The ExpandGrid function generates a grid where the 
   t[,2] values change the fastest and the t[,1] values change the slowest.
   However, if you view the interpolation on a grid, you want 
   the X variable to change the fastest along each row. Therefore, 
   for visualization purposes, reverse the first and second columns of t. */

/* for visualization, reverse the arguments to ExpandGrid and then swap 1st and 2nd columns of t */
xGrd = do(0, 1, 0.2);
yGrd = do(0, 1, 0.2);
t = ExpandGrid(yGrd, xGrd);  /* want X changing fastest; put in 2nd col */
t = t[ ,{2 1}];              /* reverse the columns from (y,x) to (x,y) */
F = bilinearInterpSquare(z, t);
Q = shape(F, ncol(yGrd), ncol(xGrd));
print Q[r=(char(YGrd,3)) c=(char(xGrd,3)) label="Bilinear Interpolant"];

/* increase the grid resolution and visualize by using a heat map */
xGrd = do(0, 1, 0.1);
yGrd = do(0, 1, 0.1);
t = ExpandGrid(yGrd, xGrd);  /* want X changing fastest; put in 2nd col */
t = t[ ,{2 1}];              /* reverse the columns from (y,x) to (x,y) */
F = bilinearInterpSquare(z, t);
Q = shape(F, ncol(yGrd), ncol(xGrd));

/* currently, the Y axis is pointing down. Flip it and the labels. */
H = Q[ nrow(Q):1, ];
reverseY = yGrd[ ,ncol(yGrd):1];
call heatmapcont(H) range={0 4} displayoutlines=0
     xValues=xGrd yValues=reverseY
     colorramp=palette("Spectral", 7)[,7:1]
     title="Interpolation at Multiple Points in the Unit Square";
QUIT;
