/* Translation of the Alan Genz MVN probability algorithm from MATLAB to SAS IML.
   See the original file multinorm_CDF.m for comments and references.

   Implements equation (14) in Section 3.2 of Genz (2004), integrating each
   term in (14) separately in terms of theta between 0 and arsin(rho_j1), using
   adaptive quadrature.
   **************************************/

proc iml;

/* Truncate a value x to be within the bounds ab=[a,b] */
start Clip(x, ab);
   a = ab[1]; b = ab[2];
   return( choose(x=., ., a <> (x >< b)) );
finish;

/* Convert from centered Covariance scale to Correlation scale.
   This function OVERWRITES the arguments, so be sure to 
   send in a COPY if you need to preserve the original data.
   Assume that we have already verified that Sigma is SPD.
   Also assumes that the problem is centered.
   For the uncentered problem, call as 
      b_new = b - mu;
      R = Sigma;
      run Convert_Cov_Limits(b_new, R);  
   */
start Convert_Cov_Limits(b, Sigma);
   n = nrow(Sigma);
   D = sqrt(vecdiag(Sigma));
   Sigma = cov2corr(Sigma);
   /* assume b is a row vector or a matrix of row vectors */
   b = b / rowvec(D);
finish;

/* Integrand function using vectorized globals */
start tvnIntegrand(theta) 
      GLOBAL(g_b, g_rho_vec);
   
   /* rho_vec = {rho_j1, rho_k1, rho_32} */
   rho_j1 = g_rho_vec[1];
   rho_k1 = g_rho_vec[2];
   rho_32 = g_rho_vec[3];
   if abs(rho_j1)=0 then return(0);
   asin_rho = arsin(rho_j1);
   Left = min(0, asin_rho);
   Right = max(0, asin_rho);
   if theta < Left | theta > Right then 
      return( 0 );   /* or should it return 0? */

   /* g_b = {b1, bj, bk} */
   c1 = ((g_b[1] - g_b[2])**2) / 2;
   c2 = g_b[1] * g_b[2];

   sintheta = sin(theta);
   cossqtheta = cos(theta)**2;

   if cossqtheta < 1e-12 then 
      return (0);

   expon = (c1 + c2 * (1 - sintheta)) / cossqtheta;
   if -expon < constant('logsmall') then 
      term1 = 0; 
   else 
      term1 = exp(-expon);

   sinphi = sintheta * rho_k1 / rho_j1;
   numeru = g_b[3] * cossqtheta 
            - g_b[1] * (sinphi - rho_32 * sintheta) 
            - g_b[2] * (rho_32 - sintheta * sinphi);
   
   denom_sq = cossqtheta * (cossqtheta - sinphi*sinphi - 
              rho_32 * (rho_32 - 2 * sintheta * sinphi));
   
   if denom_sq > 0 then 
      denomu = sqrt(denom_sq);
   else 
      denomu = 1e-12;

   term2 = cdf("Normal", numeru / denomu);
   return (term1 * term2);
finish;

/* Helper function to compute the integral for a specific rho component.
   NOTE: If you call the QUAD subroutine on [c,d] and c > d,
   IML silently sorts the limits of integration. For TVN, the limits are always specified as 
   [0, arsin(rho)], and the arcsine function can be negative. The fix is to add a negative sign 
   to the result if c > d. See
   https://blogs.sas.com/content/iml/2014/08/04/reversing-the-limits-of-integration.html
*/
start tvnComputeTerm(b, rho_vec, tol)
      GLOBAL(g_b, g_rho_vec);
   
   n = nrow(b);
   p_term = j(n, 1, 0);
   rho_j1 = rho_vec[1];

   if abs(rho_j1) > 0 then do;
      loLimit = 0;
      hiLimit = arsin(rho_j1);
      
      /* Set correlation global vector for this integration term */
      g_rho_vec = rowvec(rho_vec); 
      
      do i = 1 to n;
         /* Ensure g_b is a column vector as per preferences */
         g_b = colvec(b[i,]);
         
         if ^any(missing(g_b)) then do;
            call quad(val, "tvnIntegrand", loLimit || hiLimit) eps=(tol/3);
            p_term[i] = choose(loLimit<=hiLimit, val, -val);
         end;
      end;
   end;
   return(p_term);
finish;


/* ------------------------------------------------------------------
   Main Function: cdftvn
   ------------------------------------------------------------------ */
/* Return the trivariate CDF for MVN(Sigma, mu) at each row of b.
   Validate the parameters, standardize to correlation scale, and call cdftvn_impl */
start cdftvn(b, Sigma, mu={0 0 0});
   IsValid = IsValidParmsTVN(b, Sigma, mu);
   if ^IsValid then
      return(.);
   b_center = b - mu;
   R = Sigma;
   run Convert_Cov_Limits(b_center, R);
   tol = 1e-6; /* tolerance for numerical integration */
   cdf = cdftvn_impl(b_center, R, tol);
   return ( cdf );
finish;

/* This subroutine permutes the 3-element b and rho 
   vectors so that the largest magnitude is R[3,2].
   This routine overwrites the arguments in place, so 
   if you want to preserve the original values, send in a COPY of b and rho. 
   EXAMPLE:
   b_new = b;
   rho_new = rho;
   run PermuteTVNCorr(b_new, rho_new);
*/
start PermuteTVNCorr(b, rho);
   /* Permutation Logic: Maximize abs(rho_32) */
   imax = abs(rho)[<:>];    
   
   if imax = 1 then do; 
      rho_idx = {3 2 1}; 
      b_perm = b[, {3 2 1}];
   end;
   else if imax = 2 then do; 
      rho_idx = {1 3 2};
      b_perm = b[, {2 1 3}];
   end;
   else do; 
      rho_idx = {1 2 3};
      b_perm = b;
   end;

   /* Overwrite the arguments */
   rho = rho[rho_idx];
   b = b_perm;
finish;


/* the argumente have been validated and scaled. Return CDF for TVN(b; R, mu=0) */
start cdftvn_impl(b, R, tol);
   rho_vec = R[{2 3 6}]; /* Initial {rho_21, rho_31, rho_32} */
   
   /* Use a copy of b to avoid modifying the caller's matrix */
   b_perm = b; 
   run PermuteTVNCorr(b_perm, rho_vec);
   
   /* Clip limits after permutation for numerical stability */
   b_perm = Clip(b_perm, {-10, 10});

   /* Term 1: Standard bivariate normal calculation */
   p1 = cdf("Normal", b_perm[,1]) # probbnrm(b_perm[,2], b_perm[,3], rho_vec[3]);

   /* Term 2: Integration for rho_21 in the permuted set */
   p2 = tvnComputeTerm(b_perm, rho_vec, tol);

    /* Term 3: Integration for rho_31 in the permuted set */
   p3 = tvnComputeTerm(b_perm[,{1 3 2}], rho_vec[{2 1 3}], tol);

   pi = constant("pi");
   p = p1 + (p2 + p3) / (2 * pi);
   return (p);
finish;

store module=(
    Clip 
    Convert_Cov_Limits 
    tvnIntegrand 
    tvnComputeTerm 
    PermuteTVNCorr
    cdftvn_impl 
    cdftvn);

start VerifyPD(R_vec);
    if abs(R_vec[3]) < max(abs(R_vec)) then print "R_vec does not have R[2,3] as largest element", R_vec;
    R = I(3);
    R[ {2 3 6}] = R_vec;
    R[ {4 7 8}] = R_vec;
    G = root(R, "NoError");
    notPD = any(G=.);
    if notPD then print "rho_vec NOT PD!", R_vec;
    return ^notPD;
finish;

start VerifyAndPrint(tol = 1E-6) global(g_b, g_rho_vec);
   isPD = VerifyPD(g_rho_vec);
   L = 0;
   R = arsin(g_rho_vec[1]);
   theta = L // (L+R)/2 // R;
   z = j(nrow(theta), 1, 0);
   do i = 1 to nrow(theta);
      z[i] = tvnIntegrand(theta[i]);
   end;
   p_term = 0;
   rho_j1 = g_rho_vec[1];
   if abs(rho_j1) > 0 then do;
   loLimit = 0;
      hiLimit = arsin(rho_j1);
      call quad(val, "tvnIntegrand", loLimit || hiLimit) eps=(tol/3);
      p_term = choose(loLimit<=hiLimit, val, -val);
   end;
   print theta z, p_term;
finish;

/* unit tests for tvn_integrand. Make sure rho[3] is the largest in magnitude */

cnt = 1;
print cnt;
g_b = {0 0 0};
g_rho_vec = {0.1 0.2 0.5};
run VerifyAndPrint();

cnt = cnt+1;
print cnt;
g_b = {0 1 -1};
g_rho_vec = {-0.1 0.2 -0.5};
run VerifyAndPrint();

cnt = cnt+1;
print cnt;
g_b = {1 -1 0};
g_rho_vec = {0.3 0.3 0.5};
run VerifyAndPrint();

cnt = cnt+1;
print cnt;
g_b = {1 1 1};
g_rho_vec = {0.5 -0.5 0.5};
run VerifyAndPrint();

cnt = cnt+1;
print cnt;
g_b = {1 1 -1};
g_rho_vec = {-0.7 0.2 -0.7};
run VerifyAndPrint();

cnt = cnt+1;
print cnt;
g_b = {1 2 3};
g_rho_vec = {0.9 0.7 0.9};
run VerifyAndPrint();

cnt = cnt+1;
print cnt;
g_b = {0 0 0};
g_rho_vec = {0 0 0};
theta = 0;
z = tvnIntegrand(theta);
print theta z;
run VerifyAndPrint();



*quit;
