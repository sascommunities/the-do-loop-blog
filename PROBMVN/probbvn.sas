proc iml;

/* Rectangular probability for 1-D standard normal distribution.
   For the interval (L,U), use the CDF function (left-tailed probability)
   to compute the probability P( L < X < U | X ~ N(0,1) )
   L and U are scalars. A missing value indicates infinity.
   Return the probability P( L < X < U | X ~ N(0,1) ) */
start probuvn_mod(a, b, sigma, mu=0);
   L = (a - mu)/sigma;
   U = (b - mu)/sigma;
   /* Determine the type of limits and compute probability accordingly */
   if ^missing(L) & ^missing(U) then 
     return( CDF("Normal", U) - CDF("Normal", L) );
   if  missing(L) & ^missing(U) then
     return( CDF("Normal", U) );
   if ^missing(L) &  missing(U) then 
     return( SDF("Normal", L) );
   return( 1 );
finish;


/* bivariate normal probabilities on rectangular domains for 
   X~BVN(mu, Sigma) with upper integration limits U=(U1,U2). 
   The function standardizes the upper limits U and the covariance 
   matrix Sigma to get the corresponding standardized values for PROBBNRM.
   The function uses missing values in U to indicate infinity.
   .M indicates negative infinity
   .I indicates positive infinity.
*/
start probbvn_mod(L, U, Sigma, mu=j(1,2,0));
   R = cov2corr(Sigma);
   D = rowvec(sqrt(vecdiag(Sigma)));
   L_std = (L - mu)/ D;
   U_std = (U - mu)/ D;
   return probbvn_std(L_std, U_std, R[1,2]);
finish;

/* Extend the standard bivariate CDF, which is PROBBNRM(a,b, rho), to support an 
   upper limit of infinity (.P) for either argument:
   Pr(u,.P; rho) = Phi(u) = cdf("Normal", u) is probability over left half plane
   Pr(.P,v; rho) = Phi(v) = cdf("Normal", v) is probability over lower half plane
   See https://blogs.sas.com/content/iml/2023/12/04/bivariate-normal-probability-sas.html
*/
start cdfbvn_std(a,b, rho);
   if missing(a) & missing(b) then 
      return 1;
   if missing(a) then 
      return cdf("Normal", b);
   if missing(b) then 
      return cdf("Normal", a);
   return probbnrm(a, b, rho);
finish;

/* probabilities on rectangular regions for BVN(0, rho) */
start probbvn_std(L, U, rho);
   a = L[1]; c = L[2];
   b = U[1]; d = U[2];
   ma = missing(a);   mb = missing(b);
   mc = missing(c);   md = missing(d);
   if ^ma & ^mb & ^mc & ^md then       * 16. rectangle;
      return cdfbvn_std(b,d,rho) - cdfbvn_std(a,d,rho) 
           - cdfbvn_std(b,c,rho) + cdfbvn_std(a,c,rho); 
   if ^ma & ^mb & ^mc & md then        * 15. upper strip (N);
      return cdfbvn_std(b,.P,rho) - cdfbvn_std(a,.P,rho)  
           - cdfbvn_std(b,c,rho) + cdfbvn_std(a,c,rho);
   if ^ma & ^mb & mc & ^md then        * 14. lower strip (S);
      return cdfbvn_std(b,d,rho) - cdfbvn_std(a,d,rho);
   if ^ma & ^mb & mc & md then         * 13. vert strip;
      return cdfbvn_std(.P,b,rho) - cdfbvn_std(.P,a,rho);
   if ^ma & mb & ^mc & ^md then        * 12. right strip (E);
      return cdfbvn_std(.P,d,rho) - cdfbvn_std(.P,c,rho) 
           - cdfbvn_std(a,d,rho) + cdfbvn_std(a,c,rho);
   if ^ma & mb & ^mc & md then         * 11. NE quadrant;
      return 1 - cdfbvn_std(.P,a,rho)
               - (cdfbvn_std(.P,c,rho) - cdfbvn_std(a,c,rho));
   if ^ma & mb & mc & ^md then         * 10. SE quadrant;
      return cdfbvn_std(.P,d,rho) - cdfbvn_std(a,d,rho);
   if ^ma & mb & mc & md then          * 9. right half plane;
      return 1 - cdfbvn_std(.P,a,rho);
   if ma & ^mb & ^mc & ^md then        * 8. left strip (W);
      return cdfbvn_std(b,d,rho) - cdfbvn_std(b,c,rho);
   if ma & ^mb & ^mc & md then         * 7. NW quadrant;
      return cdfbvn_std(b,.P,rho) - cdfbvn_std(b,c,rho); 
   if ma & ^mb & mc & ^md then         * 6. SW quadrant;
      return cdfbvn_std(b,d,rho);
   if ma & ^mb & mc & md then          * 5. right half plane;
      return cdfbvn_std(b,.P,rho); 
   if ma & mb & ^mc & ^md then         * 4. horiz strip;
      return cdfbvn_std(.P,d,rho) - cdfbvn_std(.P,c,rho);
   if ma & mb & ^mc & md then          * 3. upper half plane;
      return 1 - cdfbvn_std(.P,c,rho);
   if ma & mb & mc & ^md then          * 2. lower half plane;
      return cdfbvn_std(.P,d,rho);
   if ma & mb & mc & md then           * 1. complete plane;
      return 1;
   return( . );          * should never reach this statement;
finish;

store module=(
probuvn_mod
probbvn_mod 
cdfbvn_std
probbvn_std
);
QUIT;
