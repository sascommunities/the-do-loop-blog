/* SAS program to accompany the article 
   "Estimate polychoric correlation by maximum likelihood estimation"
   by Rick Wicklin, published 13DEC2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/12/13/polychoric-correlation-mle.html

   This program shows how to use MLE to manually estimate the 
   polychoric correlation for two ordinal variables.
   NOTE: You can use the PLCORR in PROC FREQ to obtain the
   polychoric correlation, so this method reproduces
   a statistic that is easily available in SAS.
*/

/****************************************************/
/* DEFINE AND STORE THE CDFBN and PROBBVN FUNCTIONS */
/* SAS IML definitions of functions that compute the bivariate normal probability
   for rectangular regions in the plane.
   
   "Bivariate normal probability in SAS"
   by Rick Wicklin, published 04DEC2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/12/04/bivariate-normal-probability-sas.html
*/

proc iml;
/* Extend the bivariate CDF, which is PROBBNRM(a,b, rho) to support an 
   upper limit of infinity (.P) for either argument:
   Pr(u,.P; rho) = Phi(u) = cdf("Normal", u) is probability over left half plane
   Pr(.P,v; rho) = Phi(v) = cdf("Normal", v) is probability over lower half plane
*/
start CDFBN(u,v,rho);
   if missing(u) & missing(v) then 
      return 1;
   if missing(u) then 
      return cdf("Normal", v);
   if missing(v) then 
      return cdf("Normal", u);
   return probbnrm(u, v, rho);
finish;

start ProbBVN(a,b,c,d,rho);
   ma = missing(a);   mb = missing(b);
   mc = missing(c);   md = missing(d);
   /* 1. complete plane */
   if ma & mb & mc & md then
      return 1;
   /* 2. lower half plane */
   if ma & mb & mc & ^md then
      return CDFBN(.P,d,rho);
   /* 3. upper half plane */
   if ma & mb & ^mc & md then
      return 1 - CDFBN(.P,c,rho);
   /* 4. horiz strip */
   if ma & mb & ^mc & ^md then
      return CDFBN(.P,d,rho) - CDFBN(.P,c,rho);
   /* 5. left half plane */
   if ma & ^mb & mc & md then 
      return CDFBN(b,.P,rho); 
   /* 6. SW quadrant */
   if ma & ^mb & mc & ^md then
      return CDFBN(b,d,rho);
   /* 7. NW quadrant */
   if ma & ^mb & ^mc & md then 
      return CDFBN(b,.P,rho) - CDFBN(b,c,rho); 
   /* 8. left strip (W) */
   if ma & ^mb & ^mc & ^md then
      return CDFBN(b,d,rho) - CDFBN(b,c,rho);
   /* 9. right half plane */
   if ^ma & mb & mc & md then
      return 1 - CDFBN(.P,a,rho);
   /* 10. SE quadrant = lower half - (SW quad)*/
   if ^ma & mb & mc & ^md then
      return CDFBN(.P,d,rho) - CDFBN(a,d,rho);
   /* 11. NE quadrant = right half - (SE quad) */
   if ^ma & mb & ^mc & md then 
      return 1 - CDFBN(.P,a,rho)
               - (CDFBN(.P,c,rho) - CDFBN(a,c,rho));
   /* 12. right strip (E) */
   if ^ma & mb & ^mc & ^md then
      return CDFBN(.P,d,rho) - CDFBN(.P,c,rho) 
           - CDFBN(a,d,rho) + CDFBN(a,c,rho);
   /* 13. vert strip */
   if ^ma & ^mb & mc & md then
      return CDFBN(.P,b,rho) - CDFBN(.P,a,rho);
   /* 14. lower strip (S) */
   if ^ma & ^mb & mc & ^md then
      return CDFBN(b,d,rho) - CDFBN(a,d,rho);
   /* 15. upper strip (N) */
   if ^ma & ^mb & ^mc & md then
      return CDFBN(b,.P,rho) - CDFBN(a,.P,rho) /*upper strip _|_  */
           - CDFBN(b,c,rho) + CDFBN(a,c,rho);
   /* 16. rectangle */
   if ^ma & ^mb & ^mc & ^md then
      return CDFBN(b,d,rho) - CDFBN(a,d,rho) 
           - CDFBN(b,c,rho) + CDFBN(a,c,rho); 
   return( . );  /* should never execute this statement */
finish;
store module=(CDFBN ProbBVN);
quit;
/****************************************************/



/* Encode the ordinal variables as 1,2,...,k.
   Use a SAS format to display the class levels.
*/
proc format;
value QualFmt
      1="Strongly Disagree"   2="Disagree"
      3="Agree"               4="Strongly Agree";
value PerfFmt
      1="Low"   2="Medium"    3="High";
run;
data Agree;
format Quality QualFmt. Performance PerfFmt.;
input Quality Performance Count @@;
label Quality="Quality Course?" Performance="Performance on Tests";
datalines;
1 1 131   1 2  71   1 3  20 
2 1 217   2 2 207   2 3 112 
3 1 213   3 2 337   3 3 257 
4 1  52   4 2 139   4 3 244 
;

/* two easy ways to obtain the polychoric correlation: PROC FREQ or PROC CORR */
proc freq data=Agree;
   tables Quality * Performance / PLCORR CL norow nocol nopercent out=FreqOut;
   weight Count;
run;

proc corr data=Agree polychoric;
   var Quality Performance;
   freq Count;
run;



/* How to understanding polychoric correlation.
   This program computes the polychoric correlation according to 
   Olsson (1979)
   MAXIMUM LIKELIHOOD ESTIMATION OF THE POLYCHORIC CORRELATION COEFFICIENT
   PSYCHOMETRIKA--VOL. 44, NO. 4,
   https://www.researchgate.net/profile/Ulf-Olsson-4/publication/24062390_Maximum_likelihood_estimation_of_the_polychoric_correlation_coefficient/links/53f30b5e0cf2dd48950c8372/Maximum-likelihood-estimation-of-the-polychoric-correlation-coefficient.pdf

   1: Simulate bivariate normal data X ~ N(0, Sigma) where
         Sigma = {1 0.4, 0.4 1};
   2: Create cut points to construct ordinal variable. Visualize.
         Cut points for x1 are {-0.5, 0.5}
         Cut points for x2 are {-1.2, -0.3, 0.8}
   3: Use PROC FREQ to obtain estimate of polychoric correlation.
         Save frequencies and marginal frequencies to data sets.
   4. Use PROC IML to estimate cut points from marginal cumulative frequencies.
   5. Define LL as sum( n[i,j]*pi[i,j] ), where 
         pi[i,j] is the BVN probability integrated over the region defined by the cut points.
         Plot LL as a function of rho.
   6. Use one-parameter MLE to estimate the polychoric correlation (simple estimate)
   7. Use multi-parameter MLE to estimate polychoric correlation (more complicated)
*/

/* Use PROC FREQ to get the frequencies and marginal cumulative percentages */
proc freq data=Agree(rename=(Performance=Y1 Quality=Y2));
   tables Y1      / out=Y1Out outcum;
   tables Y2      / out=Y2Out outcum;
   tables Y2 * Y1 / out=FreqOut norow nocol nopercent;
   weight Count;
run;

/* perform MLE for polychoric correlation */
proc iml;
load module=(CDFBN ProbBVN);

/* The cutoff points a[i] and b[i] are estimated from the 
   empirical cumulative proportions. */
use Y1Out; read all var {Y1 'cum_pct'}; close;
u = cum_pct / 100;                  /* convert to proportion */
k1 = nrow(u);
a = .M // quantile("Normal", u[1:k1-1]) // .P;
use Y2Out; read all var {Y2 'cum_pct'}; close;
v = cum_pct / 100;                  /* convert to proportion */
k2 = nrow(v);
b = .M // quantile("Normal", v[1:k2-1]) // .P;

ods layout gridded advance=table columns=2 column_gutter=10px;
   print a[F=7.4], b[F=7.4];
ods layout end;


/* evaluate the probability that a BVN(0,rho) variate
   is in each region. For example, choose rho=0.4 */
rho = 0.4;
pi = j(nrow(b)-1, nrow(a)-1, .);
do i = 1 to nrow(b)-1;
   do j = 1 to nrow(a)-1;
      pi[i,j] = ProbBVN(a[j], a[j+1], b[i], b[i+1], rho);
   end;
end;
*print pi[r=a c=b F=4.2];
/* reverse the rows so Y1 points up */
P = pi[nrow(pi):1, ];
print P[r={0.7807 -0.3081 -1.2212 .M} 
         c={.M -0.5058 0.4775} F=4.2 L="Bivariate Normal Probability"];

/* compute the standard BVN probability for the rectangular
   regions defined by the vectors a and b. 
   The a vector defines the k1-1 intervals
      (a1,a2], (a2,a3], ... (a[k1-1], a[k1]),  k1 >= 2
   and the b vector defines the k2-1 intervals
      (b1,b2], (b2,b3], ... (b[k2-1], a[k2]),  k2 >= 2
   The (k1-1)*(k2-1) rectangular regions are the product of these intervals.
   The BVN probability is returned in a (k2-1) x (k1-1) matrix
   where rows correspond to b and columns correspond to a.
   It is common for a1    = b1    = .M = -Infinity 
                and a[k1] = b[k2] = .I = +Infinity
*/
start ProbBVNMatrix(rho, _a, _b);
   a = colvec(_a); b = colvec(_b);
   pi = j(nrow(b)-1, nrow(a)-1, .);
   do i = 1 to nrow(b)-1;
      do j = 1 to nrow(a)-1;
         pi[i,j] = ProbBVN(a[j], a[j+1], b[i], b[i+1], rho);
      end;
   end;
   return pi;
finish;
/* test calling the ProbBVNMatrix function */
pi2 = ProbBVNMatrix(rho, a, b);
print (pi-pi2);


/* Olsson (1979) "Maximum Likelihood Estimation of the Polychoric Correlation Coefficient", Psychometrika
   Example: Find LL for the value of rho
*/
/* read the observed frequencies */
use FreqOut; read all var {'Count'}; close;
n = shape(Count, nrow(Y2), nrow(Y1));
*print n[r=Y2 c=Y1];

/* evaluate LL for this value of rho */
LL = sum( n#log(pi) );
print rho LL;


/* plot LL vs rho */
pi = j(nrow(n), ncol(n), .);
rhos = T( do(0, 0.6, 0.01) );  /* choose some values of rho */
LL = j(nrow(rhos), 1);
do k = 1 to nrow(rhos);        
   pi = ProbBVNMatrix(rhos[k], a, b); /* probs for each rho value */
   LL[k] = sum( n#log(pi) );          /* evaluate LL(rho) */
end;

ods graphics/reset;
title "Loglikelihood for polychoric correlation";
call series(rhos, LL) grid={x y};

maxIdx = LL[<:>];
maxLLrho = rhos[maxIdx];
print maxLLrho;

print "=================================";
/* MLE optimization of the LL objective function. See
   https://blogs.sas.com/content/iml/2011/10/12/maximum-likelihood-estimation-in-sasiml.html */
start LLPolychorRHO(rho) global( g_a, g_b, g_Freqs );
   a = g_a; b = g_b;
   pi = ProbBVNMatrix(rho, a, b); 
   LL = sum( g_Freqs#log(pi) );
   return LL;
finish;

/* define global variables used in the optimization */
g_a = a;
g_b = b;
g_Freqs = n;

/* set-up for Newton-Raphson optimization NLPNRA */
con = { -0.99,  /* lower bounds: -1 < rho     */
         0.99}; /* upper bounds:      rho < 1 */
rho0 = 0.1;       /* initial guess for solution */
opt = {1};        /* find maximum of function   */

/* test */
LL0 = LLPolychorRHO(rho0);
*print LL0;

/* find optimal value of rho, given the data and cut points */
call nlpnra(rc, polycorrEst1, "LLPolychorRHO", rho0) opt=opt; * blc=con;
print polycorrEst1;


print "=================================";

/* the full objective function for LL of polychoric correlation  */
start LLPolychor(param) global( g_Freqs );
   /* param = rho // a // b */
   k1 = ncol(g_Freqs)-1;  /* num threshold params for X1 */
   k2 = nrow(g_Freqs)-1;  /* num threshold params for X2 */
   rho = param[1];
   a = .M // param[2:k1+1] // .P;
   b = .M // param[k1+2:1+k1+k2] // .P;
   pi = ProbBVNMatrix(rho, a, b); 
   if any(pi<0) then return( . );
   LL = sum( g_Freqs#log(pi) );
   return LL;
finish;

reset wide;

/* Maximize the full log-likelihood function */
/* First, create initial guess from empirical thresholds and best rho value */
call nlpnra(rc, polycorr0, "LLPolychorRHO", rho0) opt=opt blc=con;
param0 = polycorr0      //
        a[2:nrow(a)-1] //
        b[2:nrow(b)-1];

/* define global variables used in the optimization */
g_Freqs = n;

/* set-up for Newton-Raphson optimization NLPNRA */
/* Second, set up contraint matrix 
   -0.99 < rho < 0.99
   0 < a[1] < a[2] < 1
   0 < b[1] < b[2] < b[3] < 1
*/   
/*      rho  a1 a2 b1 b2 b3 op RHS*/
con = {-0.99 .  .  .  .  .   . .,   /* lower bounds */
        0.99 .  .  .  .  .   . .,   /* upper bounds */
        0    1 -1  0  0  0  -1 0,   /* a1 < a2 */
        0    0  0  1 -1  0  -1 0,   /* b1 < b2 */
        0    0  0  0  1 -1  -1 0};  /* b2 < b3 */
opt = {1,    /* find maximum of function */
       1};   /* print some intermediate optimization steps */

/* Third, optimize the LL to estimate polychoric correlation */
call nlpnra(rc, polycorrEst2, "LLPolychor", param0) opt=opt blc=con;
print polycorrEst2[c={'rho' 'a1' 'a2' 'b1' 'b2' 'b3'}];

LLopt = LLPolychor(polycorrEst2);
print LLopt;

QUIT;
