/* SAS program to accompany the article 
   "The distribution of the differences between two beta random variables"
   by Rick Wicklin, published 01MAR2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/03/01/distribution-difference-beta.html

   This program shows how to evaluate Appell's hypergeometric function 
   F1(a,b1,b2,c; x,y)
   for the parameter values c > a > 0
   Appell's F1 function can be evaluated by computing a definite integral.
   The function is defined when |x| < 1 and |y| < 1. 

   References
   Weisstein, Eric W. "Appell Hypergeometric Function." From MathWorld--A Wolfram Web Resource
        https://mathworld.wolfram.com/AppellHypergeometricFunction.html
   Wikipedia, https://en.wikipedia.org/wiki/Appell_series       
*/


/* Appell hypergeometric function of 2 vars 
   F1(a,b1,b2; c; x,y) is a function of (x,y) with parms = a // b1 // b2 // c;
   F1 is defined on the domain {(x,y) | |x|<1 and |y|<1}.
   You can evaluate F1 by using an integral for c > a > 0, as shown at
   https://en.wikipedia.org/wiki/Appell_series#Integral_representations
*/
proc iml;
/* Integrand for the QUAD function */
start F1Integ(u)  global(g_parm, g_x, g_y);
   a=g_parm[1]; b1=g_parm[2]; b2=g_parm[3]; c=g_parm[4];
   x=g_x; y=g_y;
   numer = u##(a-1) # (1-u)##(c-a-1);
   denom = (1-u#x)##b1 # (1-u#y)##b2;
   return( numer / denom );
finish;
 
/* Evaluate the Appell F1 hypergeometric function when c > a > 0
   and |x|<1 and |y|<1
   Pass in parm = {a, b1, b2, c} and
   z = (x1 y1,
        x2 y2,
        ...
        xn yn}; */
start AppellF1(z, parm)  global(g_parm, g_x, g_y);
   /* transfer parameters to global symbols */
   g_parm = parm;
   a = parm[1]; b1 = parm[2]; b2 = parm[3]; c = parm[4];
   n = nrow(z);
   F = j(n, 1, .);
   if a <= 0 | c <= a then do;
      /* print error message or use PrintToLOg function:
         https://blogs.sas.com/content/iml/2023/01/25/printtolog-iml.html */
      print "This implementation of the F1 function requires c > a > 0.",
            parm[r={'a','b1','b2','c'}];
      return( F );
   end;
   /* loop over the (x,y) values */
   do i = 1 to n;
      /* defined only when |x| < 1 */
      g_x = z[i,1]; g_y = z[i,2];
      if g_x^=. & g_y^=. & 
         abs(g_x)<1 & abs(g_y)<1 then do;  
         call quad(val, "F1Integ", {0 1}) MSG="NO";
         F[i] = val;
      end;
   end;
   return( F / Beta(a, c-a) );    /* scale the integral */
finish;
 
/* Example of calling AppellF1 */
parm = {2, 1, 1, 3};
parm = {2, 1, 1, 3};
xy = { 0.9  0,
       0.7  0.3,
      -0.5  0.2,
      -0.9 -0.5,
       0    0};
F1 = AppellF1(xy, parm);
print xy[L="" c={'x' 'y'}] F1;

/* compare with Mathematica's AppellF1 function:
   https://reference.wolfram.com/language/ref/AppellF1.html
*/

/* visualize the 2-D graph of Appell's function */
parm = {2, 1, 1, 3};
xy = expandgrid( do(-0.95, 0.9, 0.05), do(-0.95, 0.9, 0.05) );
F1 = AppellF1( xy, parm);
x=xy[,1]; y=xy[,2];
create Appell var {'x' 'y' 'F1'}; append; close;

QUIT;


/* Use a contour plot to visualize Appell's hypergeometric function.
   Template from "Create a contour plot in SAS"
   https://blogs.sas.com/content/iml/2012/07/02/create-a-contour-plot-in-sas.html
*/
proc template;
define statgraph ContourPlotParm;
dynamic _X _Y _Z _TITLE;
begingraph;
   entrytitle _TITLE;
   layout overlay;
      contourplotparm x=_X y=_Y z=_Z /
        contourtype=fill nhint=19  colormodel=twocolorramp name="Contour";
      continuouslegend "Contour" / title=_Z;
   endlayout;
endgraph;
end;
run;

proc sgrender data=Appell template=ContourPlotParm;
dynamic _TITLE="Graph of Appell's Hypergeometric Function F1(2,1,1,3; x,y)"
        _X="x" _Y="y" _Z="F1";
run;
