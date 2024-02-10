/* SAS program to accompany the article 
   "An exact formula for the sampling distribution of the correlation coefficient"
   by Rick Wicklin, published 12FEB2024 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2024/02/12/distribution-correlation.html

   This program shows how to evaluate an exact formula for the 
   distribution of the Pearson correlation in random samples of 
   size n drawn from a bivariate normal distribution with correlation rho.

   The program evaluates two formulas:
   1. Fisher's formula requires evaluating the integral of 
      1/(cosh(w) - r*rho)^(n-1)
      on the domain w in (0,Infinity)
   2. Hotteling's formula requires evaluating the hypergeometric function
*/

/* Fisher's exact density function */
proc iml;

/* set parameters for an example */
n = 20;
rho = 0.6;

/*  Integrand is 1 / (cosh(w) - r*rho)##(n-1) = (cosh(w) - r*rho)##-(n-1)  
    Specify constants as GLOBAL variables prefixed by 'g_'. */
start Integrand(w) global(g_r, g_rho, g_n);
   return (cosh(w) - g_r*g_rho)##(-(g_n-1));
finish;

/* define the global variables for  r = 0.5 */
g_r =  0.5;  g_rho = rho;  g_n = n;
domain = {0  .P};
call quad(integral, "Integrand", domain);
print integral;

/* Fisher's exact PDF */
start CorrPDFFisher(r, rho, n) global(g_r, g_rho, g_n);
   /* copy values to global vars for the integration */
   g_r = r;   g_rho = rho;   g_n = n;
   domain = {0  .P};
   call quad(integral, "Integrand", domain); /* compute integral */

   pi = constant('pi');
   C = (n-2) * (1-rho**2)##((n-1)/2) * (1-r**2)##((n-4)/2) / pi;
   return C*integral;
finish;

/* note that it is more probably to get r=0.65 and r=0.7, even though rho=0.6 */
/* evaluate the PDF for a range of r values and graph the PDF */
r = do(-0.2, 1, 0.01);
PDF = j(1, ncol(r), .);
do i = 1 to ncol(r);
   PDF[i] = CorrPDFFisher(r[i], rho, n );
end;

title "Distribution of Sample Correlation Coefficient";
title2 "rho = 0.6; n = 20; Fisher's Formula";
call series(r, PDF) grid={x y};
/* save for later plotting */
create Fisher var {"r" "PDF"}; append; close;

/* where is the mode? */
idx = PDF[<:>];
print (r[idx])[L='mode'];

/* for kicks, estimate the expected value */
start momentInteg(r) global(g_r, g_rho, g_n);
   rho = g_rho; n = g_n;
   return r*CorrPDFFisher(r, rho, n);
finish;
start ExpectedR(rho, n) global(g_rho, g_n);
   /* copy values to global vars for the integration */
   g_rho = rho;   g_n = n;
   domain = {-1  1};
   call quad(mu, "momentInteg", domain) peak=rho; /* compute expected value */
   return mu;
finish;
/* expected value does not equal rho, so the sample correlation coefficient 
   is a biased estimator */
ER = ExpectedR(rho, n);
print ER;

/* let's compute the PDF for rho = 0, 0.5, and 0.8 and save to a data set */
rhoVec = {0 0.5 0.8};
n = 20;
r = do(-0.5, 1, 0.01);
PDF = j(1, ncol(r), .);
results = {. . .};
create CorrPDF from results[c={'rho' 'r' 'PDF'}];
do j = 1 to ncol(rhoVec);
   rho = rhoVec[j];
   do i = 1 to ncol(r);
      PDF[i] = CorrPDFFisher(r[i], rho, n );
   end;
   results = j(ncol(r), 1, rho) || r` || PDF`;
   append from results;
end;
close;
QUIT;

title "Distribution of Sample Correlation Coefficient";
title2 "Random Sample (n=20) from Bivariate Normal Distribution (rho)";
proc sgplot data=CorrPDF;
   series x=r y=PDF / group=rho lineattrs=(thickness=2);
   xaxis grid values=(-1 to 1 by 0.1) valueshint label="Sample Correlation (r)";
   yaxis grid label="PDF (Density)";
run;


/* Let's make sure we compute this correctly. Use simulation to generate
   10,000 correlation estimates from random samples of size n=20 
   drawn from bivariate normal distribution.

   You can use the Wishart distribution to simulate the sample correlation:
   https://blogs.sas.com/content/iml/2017/10/11/simulate-correlations-wishart-distribution.html
*/
proc iml;
call randseed(1234);
mu = {0 0};
Sigma = {1    0.6, 
         0.6  1  };
numSamples = 10000;
n = 20;
/* Wishart */
DF = n - 1;               /* X ~ N obs from MVN(0, Sigma) */
S = RandWishart(NumSamples, DF, Sigma); /* each row is 2x2 scatter matrix. Cov is 1/DF * S */
r = j(numSamples, 1, .);  /* allocate vector for correlation estimates */
do i = 1 to nrow(S);      /* convert to correlation; extract off-diagonal */
   r[i] = cov2corr( shape(S[i,], 2, 2) )[1,2];
end;
title "Distribution of Sample Correlation Coefficient";
title2 "rho = 0.6; n = 20";
call histogram(r);

/* write simulated distribution to data set */
x = r;
create Sim var "x"; append; close;
QUIT;

/* Overlay the histogram of the simulated correlations and Fisher's
   exact formula for the PDF */
data All;
   set Sim Fisher;
run;

/* GTL for overlaying density curve on histogram 
   https://blogs.sas.com/content/iml/2013/04/24/overlay-density-curve-on-a-histogram.html
*/
proc template;
define statgraph ContPDF;
dynamic _X _T _Y _Title _HAlign _binstart _binstop _binwidth;
begingraph;
   entrytitle halign=center _Title;
   layout overlay / xaxisopts=(linearopts=(viewmax=_binstop
                    /* I've added a specific tick list */
                    tickvaluelist= (-0.2 -0.1 0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1)
                    ));
      histogram _X / name='hist' SCALE=DENSITY binaxis=true 
         endlabels=true xvalues=leftpoints binstart=_binstart binwidth=_binwidth;
      seriesplot x=_T y=_Y / name='PDF' legendlabel="PDF" lineattrs=(thickness=2);
      discretelegend 'PDF' / opaque=true border=true halign=_HAlign valign=top 
            across=1 location=inside;
   endlayout;
endgraph;
end;
run;

proc sgrender data=All template=ContPDF;
dynamic _X="X" _T="r" _Y="PDF" _HAlign="right"
        _binstart=-0.2 _binstop=1 _binwidth=0.02
        _Title="Distribution of Correlation Coefficient (n=20; rho=0.6)";
run;


/* NOTE: There is also a formula that uses the hypergeometric distribution */
/* HOTELLING (1953) */
/* to evaluate the hypergeometri 1F2 function in SAS, see
   https://blogs.sas.com/content/iml/2023/02/27/hypergeometric-function-sas.html
*/

proc iml;
/******************************************************
   Copy functions from 
   https://blogs.sas.com/content/iml/2023/02/27/hypergeometric-function-sas.html
 ******************************************************/

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

/* Now use the previous definitions to implement Hotteling's formula */
start CorrPDFHotel(r, rho, n);
   numer = (n-2) * (1-rho**2)##((n-1)/2) * (1-r**2)##((n-4)/2);
   denom = sqrt(2) * (n-1) * Beta(1/2, n-1/2) * (1-rho*r)##(n-3/2);
   param = 1/2 //  1/2 // n-1/2;
   hyper = Hypergeometric2F1((1+rho*r)/2, param);
   return numer / denom * hyper;
finish;

/* test it out by graphing the PDF for an example */
n = 20;
rho = 0.6;
r = do(-0.2, 1, 0.01);

PDFH = j(1, ncol(r), .);
do i = 1 to ncol(r);
   PDFH[i] = CorrPDFHotel(r[i], rho, n);
end;

title "Distribution of Sample Correlation Coefficient";
title2 "rho = 0.6; n = 20; Hotelling's Formula";
call series(r, PDFH) grid={x y};

QUIT;

