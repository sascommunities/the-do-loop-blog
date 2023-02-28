/* SAS program to accompany the article 
   "Compute hypergeometric functions in SAS"
   by Rick Wicklin, published 27FEB2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/02/27/hypergeometric-function-sas.html

   This program shows how to evaluate Gauss's hypergeometric function 
   F(a,b;c;x)
   for c > b > 0 and a > 0. For those parameters, you can evaluate the 
   hypergeometric function by solving a definite integral. Gauss's 
   hypergeometric function is defined when |x| < 1. 

   References
   Beukers, F. (2014) "Hypergeometric Functions, How Special Are They?" Notices of the AMS, 61(1), p. 48-56.
         https://www.ams.org/notices/201401/rnoti-p48.pdf
   Weisstein, Eric W. "Hypergeometric Function." From MathWorld--A Wolfram Web Resource.
         https://mathworld.wolfram.com/HypergeometricFunction.html
*/
/* 1-D hypergeometric function 1F2(a,b;c;x) */
proc iml;
/* Integrand for the QUAD function. Parameters are in the GLOBAL clause. */
start HG1(u)  global(g_parm, g_x);
   a=g_parm[1]; b=g_parm[2]; c=g_parm[3]; x=g_x;
   /* Note: integrand unbounded as x->1 */
   F = u##(b-1) # (1-u)##(c-b-1) / (1-x#u)##a;
   return( F );
finish;

/* Evaluate the 2F1 hypergeometric function when c > b > 0 and |x|<1 */
start Hypergeometric2F1(x, parm)  global(g_parm, g_x);
   /* transfer the parameters to global symbols for QUAD */
   g_parm = parm;
   a = parm[1]; b = parm[2]; c = parm[3]; 
   n = nrow(x)* ncol(x);
   F = j(n, 1, .);
   if b <= 0 | c <= b then do;
      /* print error message or use PrintToLOg function:
         https://blogs.sas.com/content/iml/2023/01/25/printtolog-iml.html */
      print "This implementation of the F1 function requires c > b > 0.",
            parm[r={'a','b','c'}];
      return( F );
   end;
   /* loop over the x values */
   do i = 1 to n;
      if x[i]^=. & abs(x[i])<1 then do;  /* defined only when |x| < 1 */
         g_x = x[i];
         call quad(val, "HG1", {0 1}) MSG="NO";
         F[i] = val;
      end;
   end;
   return( F / Beta(b, c-b) );    /* scale the integral */
finish;

/*
In WolframAlpha:
hypergeometric 2F1[ 1/2,1/2, 1, [-0.5,-0.25, 0, 0.25, 0.5, 0.9]]
{0.901286, 0.945006, 1, 1.07318, 1.18034, 1.64126}
*/

/* TEST */
a = 0.5;
b = 0.5;
c = 1;
parm = a//b//c;

x = {-0.5, -0.25, 0, 0.25, 0.5, 0.9};
F = Hypergeometric2F1(x, parm);
print x F[L='2F1(0.5,0.5;1;x)'];
/*
 x    F 
-0.5  0.9012863 
-0.25 0.9450063 
 0    1 
 0.25 1.073182 
 0.5  1.1803406 
 0.9  1.6412644 
*/

/* draw a graph for x in (-1, 1) */
x = T( do(-0.99, 0.99, 0.01) );
F = Hypergeometric2F1(x, parm);

title "2F1(0.5,0.5; 1; x)";
call series(x, F) grid={x y} label={'x' '2F1(x)'};

/* Note: from 
   https://mathworld.wolfram.com/HypergeometricFunction.html
   2F1(1,1; 2, x) = -log(1-x)/x
*/
x = round( T( do(-0.9, 0.9, 0.1) ), 0.01);
parm = {1, 1, 2};
F = Hypergeometric2F1(x, parm);
F2 = -log(1-x)/x;
/* when x=0, log(1-x)/x is undefined, but IML evaluates 
   log(1-eps)/eps = 0 for small eps */
print (max(F-F2))[L="MaxDiff"];
print x F F2[L='formula'] (F-F2)[L='Diff'];
title "Hypergeometric Function 2F1(1,1, 2; x)";
call series(x, F);
QUIT;
