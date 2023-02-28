/* SAS program to accompany the article 
   "The distribution of the differences between two beta random variables"
   by Rick Wicklin, published 01MAR2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/03/01/distribution-difference-beta.html

   This program shows how to compute the PDF of the distribution of the 
   random variable d = X-Y where X~beta(a1,b1) and Y~beta(a2,b2)
   The formula for the PDF is given by Pham-gia and Turkkan (1993).
   It requires evaluating Appell's hypergeometric function 
   F1(a,b1,b2,c; x,y)
   The function is defined when |x| < 1 and |y| < 1. 
   For the parameter values c > a > 0, Appell's F1 function can be evaluated 
   by computing a definite integral.

   References
   T. Pham-Gia and N. Turkkan (1993), "Bayesian analysis of the difference of two proportions,"
        Communications in Statistics - Theory and Methods, 22(6), 1755-1771.
        https://doi.org/10.1080/03610929308831114

   Weisstein, Eric W. "Appell Hypergeometric Function." From MathWorld--A Wolfram Web Resource
        https://mathworld.wolfram.com/AppellHypergeometricFunction.html

   Wikipedia, https://en.wikipedia.org/wiki/Appell_series       
*/


/* MOTIVATION: Use simulation to visualize the problem. What is the distribution of 
   d = X-Y where X~beta(a1,b1) and Y~beta(a2,b2)?
*/
/* Case 2 from Pham-Gia and Turkkan, 1993, p. 1765 */
%let a1 = 0.5;
%let b1 = 0.5;
%let a2 = 1;
%let b2 = 1;
%let N = 1E5;
data BetaDiff;
call streaminit(123);
do i = 1 to &N;
   x = rand("beta", &a1, &b1);
   y = rand("beta", &a2, &b2);
   d = x - y;
   output;
end;
drop i;
run;

title "Distribution of X-Y";
title2 "X~beta(&a1,&b1); Y~beta(&a2,&b2)";
proc sgplot data=BetaDiff;
   histogram d;
run;
/*
title "Distribution of X~beta(&a1,&b1)";
proc sgplot data=BetaDiff;
   histogram x;
run;
title "Distribution of Y~beta(&a2,&b2)";
proc sgplot data=BetaDiff;
   histogram y;
run;
*/




/* Appell hypergeometric function of 2 vars 
   F1(a,b1,b2; c; x,y) is a function of (x,y) with 
   parms = a // b1 // b2 // c;
   Can evaluate F1 by using an integral for c > a > 0.
   F1 is defined on the domain {(x,y) | |x|<1 and |y|<1}

   https://en.wikipedia.org/wiki/Appell_series#Integral_representations
   Author: Rick Wicklin, SAS  20FEB2023
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

/* compare with Mathematica's AppellF1 function:
   https://reference.wolfram.com/language/ref/AppellF1.html
*/
parm = {2, 1, 1, 3};
xy = { 0.9  0,
       0.7  0.3,
      -0.5  0.2,
      -0.9 -0.5,
       0    0};
F1 = AppellF1(xy, parm);
print xy[L="" c={'x' 'y'}] F1;

/* visualize the function in (-1, 1) x (-1, 1), but avoid the singularity at x=1 and y=1 */
parm = {2, 1, 1, 3};
xy = expandgrid( do(-0.95, 0.9, 0.05), do(-0.95, 0.9, 0.05) );
F1 = AppellF1( xy, parm);
x=xy[,1]; y=xy[,2];
create Appell var {'x' 'y' 'F1'}; append; close;

/* Use Appell's hypergeometric function to evaluate the PDF
   of the distribution of the difference X-Y between 
   X ~ Beta(a1,b1) and Y ~ Beta(a2,b2)
   Formula from Pham-Gia and Turkkan, 1993
*/
start PDFBetaDiff(_d, _parm);
   a1=_parm[1]; b1=_parm[2]; a2=_parm[3]; b2=_parm[4];
   n = nrow(_d)*ncol(_d);
   PDF = j(n, 1, .);
   if a1<=0 | b1<=0 | a2<=0 | b2<=0 then do;
      print "All parameters must be positive", a1 b1 a2 b2;
      return(PDF);
   end;
   A = Beta(a1,b1) * Beta(a2,b2);
   do i = 1 to n;
      d = _d[i];
      *print d a1 b1 a2 b2;
      if d = . then 
         PDF[i] = .;
      if abs(d) < 1 then do;
         if abs(d)<1e-14 then do;
            if a1+a2 > 1 & b1+b2 > 1 then 
               PDF[i] = Beta(a1+a2-1, b1+b2-1) / A;
            *print "d=0" (a1+a2-1)[L='a1+a2-1'] (b1+b2-1)[L='b1+b2-1'] (PDF[i])[L='PDF'];
         end;
         else if d > 0 then do;   /* (0,1] */
            parm = b1            //
                   a1+b1+a2+b2-2 //
                   1-a1          //
                   a2+b1;
            xy = 1-d || 1-d##2;
            F1 = AppellF1(xy, parm);
            PDF[i] = Beta(a2,b1) # d##(b1+b2-1) # (1-d)##(a2+b1-1) # F1 / A;
            *print "d>0" xy (PDF[i])[L='PDF'];
         end;
         else do;            /* [-1,0) */
            parm = b2            //
                   1-a2          //
                   a1+b1+a2+b2-2 //
                   a1+b2;
            xy = 1-d##2 || 1+d;
            F1 = AppellF1(xy, parm);
            PDF[i] = Beta(a1,b2) # (-d)##(b1+b2-1) # (1+d)##(a1+b2-1) # F1 / A;
            *print "d<0" xy (PDF[i])[L='PDF'];
         end;
      end;
   end;
   return PDF;
finish;

/* Case 2 from Pham-Gia and Turkkan, 1993, p. 1765 */
a1 = 0.5;
b1 = 0.5;
a2 = 1;
b2 = 1;
parm = a1//b1//a2//b2;

print "*** Case 2 in Pham-Gia and Turkkan, p. 1767 ***";
d = {-0.9, -0.5, -0.25, 0, 0.25, 0.5, 0.9};
PDF = PDFBetaDiff(d, parm);
print d PDF;

/* The PDF simplifies for these parameters */
check = j(nrow(d), 1, 1);
do i = 1 to nrow(d);
   di = d[i];
   if di < 0 then 
     check[i] = 2*sqrt(-di)#sqrt(1+di)#AppellF1(1-di##2||1+di, {1,0,1,1.5}) / constant('pi');
   else if di > 0 then 
     check[i] = 2*sqrt(di#(1-di))#AppellF1(1-di||1-di##2, {0.5,1,0.5,1.5}) / constant('pi');
end;
print PDF check (PDF-check)[L='diff'];


/* graph the distribution of the difference */
d = T(do(-0.975, 0.975, 0.025));
PDF = PDFBetaDiff(d, parm);

title "X-Y for X ~ Beta(0.5,0.5) and Y ~ Beta(1,1)";
title2 "Fig 1, p. 1766";
call series(d, PDF) grid={x y}
      other="yaxis grid values=(0 to 1 by 0.2) min=0;";

/* write to data set for overlay with histogram */
create PDF var {'d' 'PDF'};
append;
close;


/***************************************************/
/* REPEAT THE ANALYSIS FOR OTHER PARAMETERS        */
/* Case 5 from Pham-Gia and Turkkan, 1993, p. 1767 */
a1 = 3;
b1 = 5;
a2 = 2;
b2 = 8;
parm = a1//b1//a2//b2;

/* again, the hypergeometric function simplifies */
print "*** CHECK CASE 5 in PHAM-GIA/TURKKAN p. 1767 ***";
d = {-0.9, -0.5, -0.25, 0, 0.25, 0.5, 0.9};
do i = 1 to nrow(d);
   di = d[i];
   f = PDFBetaDiff(di, parm);
   if di < 0 then 
      check = 21*(-di)**12*(1+di)**10*AppellF1(1-di##2||1+di, {8,-1,16,11});
   else if di > 0 then 
      check = 252*di**12*(1-di)**6*AppellF1(1-di||1-di##2, {5,16,-2,7});
   else
      check = 18/13;
   print f check (f-check)[L='diff'];
end;

/* graph the distribution of the difference */
d = T(do(-0.975, 0.975, 0.025));
pdf = PDFBetaDiff(d, parm);

title "X-Y for X ~ Beta(3,5) and Y ~ Beta(2,8)";
title2 "Fig 2, p. 1767";
call series(d, pdf) grid={x y}
      other="yaxis grid values=(0 to 2 by 0.2) min=0;";

/* write to data set and overlay with histogram of differences */
create PDF2 var {'d' 'PDF'};
append;
close;

QUIT;


/* overlay histogram of simulated data and the PDF 
   https://blogs.sas.com/content/iml/2013/04/24/overlay-density-curve-on-a-histogram.html
*/
proc template;
define statgraph ContPDF;
dynamic _X _T _Y _Title _HAlign _binstart _binstop _binwidth;
begingraph;
   entrytitle halign=center _Title;
   layout overlay /xaxisopts=(linearopts=(viewmax=_binstop));
      histogram _X / name='hist' SCALE=DENSITY binaxis=true 
         endlabels=true xvalues=leftpoints binstart=_binstart binwidth=_binwidth;
      seriesplot x=_T y=_Y / name='PDF' legendlabel="PDF" lineattrs=(thickness=2);
      discretelegend 'PDF' / opaque=true border=true halign=_HAlign valign=top 
            across=1 location=inside;
   endlayout;
endgraph;
end;
run;


/* CASE 1 */
data All;
set BetaDiff(keep=d) PDF(rename=(d=t));
label d="Difference X-Y";
run;
 
proc sgrender data=All template=ContPDF;
dynamic _X="d" _T="t" _Y="PDF" _HAlign="right"
        _binstart=-1.05 _binstop=1.05 _binwidth=0.1
        _Title="PDF of X-Y for X ~ Beta(0.5,0.5) and Y ~ Beta(1,1)";
run;


/***************************************************/
/* REPEAT THE ANALYSIS FOR OTHER PARAMETERS        */
/* CASE 2 */
%let a1 = 3;
%let b1 = 5;
%let a2 = 2;
%let b2 = 8;
%let N = 1E5;
data BetaDiff2;
call streaminit(123);
do i = 1 to &N;
   x = rand("beta", &a1, &b1);
   y = rand("beta", &a2, &b2);
   d = x - y;
   output;
end;
drop i;
run;

data All;
set BetaDiff2(keep=d) PDF2(rename=(d=t));
label d="Difference X-Y";
run;
 
proc sgrender data=All template=ContPDF;
dynamic _X="d" _T="t" _Y="PDF" _HAlign="right"
        _binstart=-1.05 _binstop=1.05 _binwidth=0.1
        _Title="PDF of X-Y for X ~ Beta(3,5) and Y ~ Beta(2,8)";
run;

/***********************************************/
/***********************************************/
/***********************************************/


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
