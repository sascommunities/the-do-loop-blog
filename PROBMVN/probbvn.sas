proc iml;

/* Rectangular probability for 1-D standard normal distribution.
   For the interval (L,U), use the CDF function (left-tailed probability)
   to compute the probability P( L < X < U | X ~ N(0,1) )
   L and U are scalars. A missing value indicates infinity.
   Return the probability P( L < X < U | X ~ N(0,1) ) */
start probuvn_mod(L, U, sigma=1, mu=0);
   L_std = (L - mu)/sigma;
   U_std = (U - mu)/sigma;
   /* Determine the type of limits and compute probability accordingly */
   prob = j(1, nrow(L_std), .);
   do i = 1 to nrow(L_std);
      a = L_std[i];
      b = U_std[i];
      if ^missing(a) & ^missing(b) then
        prob[i] = CDF("Normal", b) - CDF("Normal", a);
      else if  missing(a) & ^missing(b) then
        prob[i] = CDF("Normal", b);
      else if ^missing(a) &  missing(b) then
        prob[i] = SDF("Normal", a);
      else
        prob[i] = 1;
   end;
   return prob;
finish;


/* bivariate normal probabilities on rectangular domains for 
   X~BVN(mu, Sigma) with upper integration limits U=(U1,U2). 
   The function standardizes the upper limits U and the covariance 
   matrix Sigma to get the corresponding standardized values for PROBBNRM.
   The function uses missing values in U to indicate infinity.
   .M indicates negative infinity
   .I indicates positive infinity.
   L, U and mu are row vectors of length 2. Sigma is a 2x2 covariance matrix.
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
   Pr(a<.P; rho) = Phi(a) = cdf("Normal", a) is probability over left half plane
   Pr(.P,b; rho) = Phi(b) = cdf("Normal", b) is probability over lower half plane
   See https://blogs.sas.com/content/iml/2023/12/04/bivariate-normal-probability-sas.html
   This function is for SCALAR arguments a, b and rho.
*/
start cdfbvn_std(a,b, rho);
   if ^missing(a) & ^missing(b) then
      return probbnrm(a, b, rho);
   if missing(a) then 
      return cdf("Normal", b);
   if missing(b) then 
      return cdf("Normal", a);
   return 1;
finish;

/* probabilities on rectangular regions for BVN(0, rho)
   L and U are kx2 matrices, where each rows is the lower or upper limit of integration.
   rho is a scalar correlation coefficient */
start probbvn_std(L, U, rho);
   prob = j(nrow(L), 1, .);
   do i = 1 to nrow(L);
      a = L[i,1]; c = L[i,2];
      b = U[i,1]; d = U[i,2];
      ma = missing(a);   mb = missing(b);
      mc = missing(c);   md = missing(d);
      if ^ma & ^mb & ^mc & ^md then            * 16. rectangle;
         prob[i] = cdfbvn_std(b,d,rho) - cdfbvn_std(a,d,rho)
                 - cdfbvn_std(b,c,rho) + cdfbvn_std(a,c,rho);
      else if ^ma & ^mb & ^mc & md then        * 15. upper strip (N);
         prob[i] = cdfbvn_std(b,.P,rho) - cdfbvn_std(a,.P,rho)
                 - cdfbvn_std(b,c,rho) + cdfbvn_std(a,c,rho);
      else if ^ma & ^mb & mc & ^md then        * 14. lower strip (S);
         prob[i] = cdfbvn_std(b,d,rho) - cdfbvn_std(a,d,rho);
      else if ^ma & ^mb & mc & md then         * 13. vert strip;
         prob[i] = cdfbvn_std(.P,b,rho) - cdfbvn_std(.P,a,rho);
      else if ^ma & mb & ^mc & ^md then        * 12. right strip (E);
         prob[i] = cdfbvn_std(.P,d,rho) - cdfbvn_std(.P,c,rho)
                  - cdfbvn_std(a,d,rho) + cdfbvn_std(a,c,rho);
      else if ^ma & mb & ^mc & md then         * 11. NE quadrant;
         prob[i] = 1 - cdfbvn_std(.P,a,rho)
                  - (cdfbvn_std(.P,c,rho) - cdfbvn_std(a,c,rho));
      else if ^ma & mb & mc & ^md then         * 10. SE quadrant;
         prob[i] = cdfbvn_std(.P,d,rho) - cdfbvn_std(a,d,rho);
      else if ^ma & mb & mc & md then          * 9. right half plane;
         prob[i] = 1 - cdfbvn_std(.P,a,rho);
      else if ma & ^mb & ^mc & ^md then        * 8. left strip (W);
         prob[i] = cdfbvn_std(b,d,rho) - cdfbvn_std(b,c,rho);
      else if ma & ^mb & ^mc & md then         * 7. NW quadrant;
         prob[i] = cdfbvn_std(b,.P,rho) - cdfbvn_std(b,c,rho);
      else if ma & ^mb & mc & ^md then         * 6. SW quadrant;
         prob[i] = cdfbvn_std(b,d,rho);
      else if ma & ^mb & mc & md then          * 5. right half plane;
         prob[i] = cdfbvn_std(b,.P,rho);
      else if ma & mb & ^mc & ^md then         * 4. horiz strip;
         prob[i] = cdfbvn_std(.P,d,rho) - cdfbvn_std(.P,c,rho);
      else if ma & mb & ^mc & md then          * 3. upper half plane;
         prob[i] = 1 - cdfbvn_std(.P,c,rho);
      else if ma & mb & mc & ^md then          * 2. lower half plane;
         prob[i] = cdfbvn_std(.P,d,rho);
      else if ma & mb & mc & md then           * 1. complete plane;
         prob[i] = 1;
   end;
   return( prob );
finish;

store module=(
probuvn_mod
probbvn_mod 
cdfbvn_std
probbvn_std
);
QUIT;
