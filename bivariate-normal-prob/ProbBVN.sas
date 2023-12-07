/* SAS IML definitions of functions that compute the bivariate normal probability
   for rectangular regions in the plane.
   
   "Bivariate normal probability in SAS"
   by Rick Wicklin, published 04DEC2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/12/04/bivariate-normal-probability-sas.html
*/

/* DEFINE AND STORE THE CDFBN and PROBBVN FUNCTIONS */
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
