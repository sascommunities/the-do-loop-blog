/* downloaded 25AUG2017 from 
   https://www.biostat.uni-hannover.de/fileadmin/institut/probmvn.sas

   For an overview, see https://www.biostat.uni-hannover.de/89.html?&L=1
   For orthant probabilities, see https://www.biostat.uni-hannover.de/91.html?&L=1#c140
*/

/* SAS/IML program for the calculation of multivariate normal probabilities. 
   The code uses the RANDGEN function for the generation of uniform random variables. 
   The program evaluates the multivariate normal integral by applying randomised 
   lattice rule on the transformed integral as described by Genz (1992, 1993). 
   For the evaluation of singular integrals the method follow the representation 
   of Genz and Kwong (2000). Furthermore, variable priorization and anthitec 
   sampling is used. The program is based on the Fortran code of Alan Genz.
   The mathematics are described in the monograph
   Genz and Bretz (2009), "Computation of Multivariate Normal and t Probabilities",
   Springer, New York.

   The program computes multivariate normal probabilities for positive 
   semi-definite covariance matrices until dimension 100. 
   Original Author : ALAN GENZ and FRANK BRETZ 
   Date : 21.06.00 - start version 
   Input : LOWER : lower integration limits (N-rowvector) 
           UPPER : upper integration limits (N-rowvector) 
           COVAR : positive semi-definie covariance matrix (N*N-matrix)  

   Output : ERROR : estimated absolute error, with 99% confidence level 
                  VALUE : estimated integral value

   Completely rewritten by Rick Wicklin, April 2026, to modernize, modularize, 
   and vectorize the code.
*/
/* Example call: 
   ************* 
   The statements at the end of the program 
   N = 5; 
   LOWER = J(1,N,-2); UPPER = J(1,N,2); 
   COVAR = I(5);
   RUN MVN_DIST( LOWER, UPPER, COVAR, ERROR, VALUE ); 
   PRINT ERROR VALUE; 

   lead to the following output: 
   ERROR     VALUE
   0.0000756 0.9030463
*/
options nodate ps=32000;
/* define and STORE univariate, bivariate, and other helper functions */
/*
%include "probbvn.sas";
%include "mvn_validate.sas";
%include "probmvn_validate.sas";
*/
/* assume that we have stored the matrix validation helper functions. When you 
   run 
   mvn_validate.sas
   and
   probmvn_validate.sas
   the functions are stored in a library where they can be LOADed by the LOAD MODULE= statement.
*/
/* Define the PROBMVN_MOD function for rectangular bivariate normal probabilities.
   Let X~MVN(mu, Sigma) be a multivariate normal random vector, where
   Sigma is an nxn covariance matrix and mu is a 1xn row vector.
   Then 
   prob = PROBMVN_MOD(L, U, Sigma<, mu>) 
   returns the probability that X falls within the rectangular region defined by L and U:
   P(L1 < X1 < U1 & ... & Ln < Xn < Un) 
   where a missing element for L means -Infinity and a missing element for U means +Infinity. */

proc iml;
load module=_all_;

/* probmvn_mod: Main top-level routine for rectangular multivariate normal probabilities. 
   L and U are 1xn row vectors of lower and upper limits, respectively. 
   Sigma is an nxn covariance matrix. mu is an optional 1xn mean vector. 
   The function standardizes the limits and covariance matrix
   then calls PROBMVN_STD, which computes the probability for the 
   standardized problem X~MVN(0,R). 
   The function uses missing values in L and U to indicate infinity. 
   .M indicates negative infinity.
   .I indicates positive infinity.
*/
start probmvn_mod(L, U, Sigma, mu=j(1,ncol(Sigma),0));
   D = rowvec(sqrt(vecdiag(Sigma)));
   L_std = (L - mu)/ D;
   U_std = (U - mu)/ D;
   R = cov2corr(Sigma);
   /* ensure that diagonal is EXACTLY 1 */
   diagIdx = do(1,nrow(R)*ncol(R), ncol(R)+1);
   R[diagIdx] = 1;             /* set diagonal elements */
   return probmvn_std(L_std, U_std, R); /* Note: From here on, we deal only with correlation matrices */
finish;

/* Define some constants and call mvn_dist for the standardized problem X~MVN(0,R). */
start probmvn_std(L0, U0, R0);
   /* Validate standardized arguments once so downstream routines can assume validity. */
   isValid = IsValidParmsPROBMVN(L0, U0, R0);
   if ^isValid then 
      return( j(nrow(L0),1,.) );

   /* If all limits are infinite, the probability is 1 by definition. */
   new_params = RemoveInfiniteLimits(L0, U0, R0);
   L = new_params$1;
   U = new_params$2;
   R = new_params$3;
   /* if the effective dimensions are 0, 1, or 2, solve the problem exactly */
   if ncol(L) = 0 then
      return(1);     /* If all limits are infinite, the probability is 1 by definition. */
   if ncol(L) = 1 then 
      value = probuvn_std(L, U);
   else if ncol(L)=2 then 
      value = probbvn_std(L, U, R[1,2]);
   else 
      run mvn_dist(L, U, R,
                   error, value );
   return(value);
finish;

start RemoveInfiniteLimits(L, U, R);
   /* Keep all dimensions except those with both bounds infinite/missing. */
   idx = loc( ^(L= . & U= .) );
   if ncol(idx)>0 then do;
      L_new = L[,idx];
      U_new = U[,idx];
      R_new = R[idx, idx];
   end;
   else do;
      L_new = {};
      U_new = {};
      R_new = {};
   end;
   return ( [L_new, U_new, R_new] );
finish;

/* FUNCTION: mvn_dist
High-level driver. 
1. Initializes evaluation counters.
2. Calls 'mvndnt' to sort variables and compute Cholesky decomposition.
3. Only called when effective dimension > 2 after handling infinities, so 
   calls the integration routine 'dkbvrc'.

Input Arguments:
- lower: 1xn row vector of lower limits (missing values indicate -Infinity)
- upper: 1xn row vector of upper limits (missing values indicate +Infinity)
- covar: nxn covariance matrix (in practice, is a correlation matrix)
Output Arguments:
- error: output scalar for error estimate
- value: output scalar for probability value estimate
*/
start mvn_dist( lower, upper, covar, 
                error, value );
   /* Phase 1: Setup, Pivoting, and Cholesky Factorization.
               If all of the checks pass, proceed to numerical integration
               by running dkbvrc. */
   run mvndnt( lower, upper, covar, value, error );

   n = ncol(covar);
   maxpts = 2000*n**3;
   if n < 10 then abseps = 1E-4;
   else abseps = 1E-3;
   releps = 0;
   /* perform numerical integration */
   run dkbvrc( n-1, 0, maxpts, abseps, releps, error, value );
finish mvn_dist;


/* FUNCTION: mvn_dfn
Computes the integrand value for a specific point 'w' in the unit hypercube.
This implements the Genz transformation:
1. Maps uniform point 'w' to the normal domain using Inverse Normal CDF (Quantile).
2. Computes the integration limits (ai, bi) for the current dimension conditional on previous dimensions.
3. Accumulates the probability mass (ei - di).
*/
start mvn_dfn( n, w ) 
      global( g_corr, g_done, g_eone, g_a, g_b, g_y );
   value = 1;
   infa = 0;
   infb = 0;
   ik = 1;

   /* Loop through dimensions to compute conditional probabilities */
   do i = 1 to n+1;
      vsum = 0;
      /* Compute inner product of Cholesky row 'i' with previous realizations 'y' */
      if ( ik > 1 ) then
         /* vsum = vsum + sum( g_corr[i,1:ik-1]#t(g_y[1:ik-1]));*/
         vsum = vsum + g_corr[i,1:ik-1] * g_y[1:ik-1];

      /* Update lower limit (ai) based on conditional mean */
      if ( g_a[i] ^= . ) then do;
         if ( infa = 1 ) then
            ai = max( ai, g_a[i] - vsum );
         else do;
            ai = g_a[i] - vsum;
            infa = 1;
         end;
      end;

      /* Update upper limit (bi) based on conditional mean */
      if ( g_b[i] ^= . ) then do;
         if ( infb = 1 ) then
            bi = min( bi, g_b[i] - vsum );
         else do;
            bi = g_b[i] - vsum;
            infb = 1;
         end;
      end;

      /* Check for end of integration or singular dimension */
      if ( i < n + 1 & ik < n + 1 ) then
         aaa = g_corr[i+1,ik+1];
      else aaa = 0;

      /* Calculate Probability for current dimension */
      if ( i = n+1 | aaa > 0 ) then do;
         if ( i = 1 ) then do;
            di = g_done;  /* this seems to be the only place where g_done and g_eone are used. They are set in the covsrt routines by ProbIntegralTransform. */
            ei = g_eone;
         end;
         else do;
            di = 0;
            ei = 1;
            if ( infa = 1 ) then
               di = cdf("Normal", ai);

            if ( infb = 1 ) then
               ei = cdf("Normal", bi);

            ei = max( ei, di );
         end;

         /* Accumulate the volume */
         if ( di >= ei ) then do;
            value = 0;
            i = n+1; /* Break loop */
         end;
         else do;
            value = value*( ei - di );

            /* Transform uniform w[ik] to Normal g_y[ik] for next iteration */
            if ( i <= n ) then
               g_y[ik] = quantile("Normal", di + w[ik]*( ei - di ) );
            ik = ik + 1;
            infa = 0;
            infb = 0;
         end;
      end;
   end;

   mvndfn = value;
   return( mvndfn );
finish mvn_dfn;


/* FUNCTION: mvndnt
Initialization and Analytic Handling.
1. Defines GLOBAL variables.
2. Calls 'covsrt' to pivot variables and compute Cholesky.
*/
start mvndnt( lower, upper, covar, value, error ) 
      global( g_corr, g_a, g_b, g_y, g_done, g_eone, g_nn, g_hisum, g_olds );
   n = ncol(covar);
   if n<2 then do;
      print "ERROR: 1-D and 2-D problems are handled directly before calling this routine.";
      value=.; error=.;
      return;
   end;
   /* Initialize global work arrays. 
      Fully infinite dimensions were already removed.
      After the reduction, 1-D and 2-D were handled directly, so n > 2 from here forward. */
   g_y = j( 1, n, . );
   g_a = lower;
   g_b = upper;
   g_corr = covar;

   /* Reset quasi-random sequence state for each top-level probability call. */
   g_nn = j(1, 51, 0);
   g_hisum = .;
   g_olds = 0;

   /* The parameters have been validated by IsValidParmsPROBMVN, which 
      performs all necessary checks, including that COVAR is SPD. */

   /* Perform Cholesky Decomposition with Variable Reordering */
   /* FEATURE FLAG: 0 = original covsrt (Bretz: greedy pivot + singularity handling)
                    1 = covsrt_naive (Rick: ROOT + ridge, no permuting) 
                    2 = covsrt_static (Rick: permute vars according to 1-D marginal probabilities) */
   use_new_covsrt = 2;  /* 0=Greedy, 1=None, 2=Static Sort */
   if use_new_covsrt = 2 then
      run covsrt_static( lower, upper, covar );
   else if use_new_covsrt = 1 then
      run covsrt_naive( lower, upper, covar );
   else
      run covsrt( lower, upper, covar );
finish mvndnt;


/* FUNCTION: ProbIntegralTransform
Convert one-dimensional normal limits on the z-scale into
probability-scale limits for the Genz transformation.
Given scalar bounds a and b, return lower = Phi(a) and upper = Phi(b),
with lower=0 for -Infinity and upper=1 for +Infinity.
Missing values in a or b are treated as infinities.
The routine is used to form the interval mass upper-lower for pivoting,
initialization, and low-dimensional special cases. In the original
Genz/Bretz code, this routine was named "mvnlms".
*/
start ProbIntegralTransform( a, b, lower, upper );
   lower = 0;
   upper = 1;
   if ( a ^= . ) then
      lower = cdf("Normal", a);
   if ( b ^= . ) then
      upper = cdf("Normal", b);
finish ProbIntegralTransform;

/* CopyLowerTriToUpper */
start CopyLowerTriToUpper(M);
   v = vech(M);          /* extract lower triangular elements (including diagonal) */
   return( sqrvech(v) ); /* create dense symmetric matrix */
finish;

/* FUNCTION: covsrt_static
   Cholesky decomposition with a STATIC Pre-Sort on 1-D marginal probabilities.
   Strategy:
   1. Move fully infinite dimensions to the end (same as original).
   2. STATIC SORT: Compute 1-D marginal probability for each active dimension, 
      sort them ascending, and symmetrically permute the active limits and covariance block.
   3. Compute Cholesky via ROOT with ridge policy (same as nopivot).
   4. Scale factor and compute conditional expectations (same as nopivot).
*/
start covsrt_static( lower, upper, covar )
      global( g_corr, g_a, g_b, g_y, g_done, g_eone );
   n = ncol(covar);
v = vecdiag(covar);
IF ANY(v^=1) THEN PRINT "WARNING: covsrt_static is designed for correlation matrices. Results may be inaccurate.", v;
   /* STATIC PRE-SORT on 1-D marginals before Cholesky factorization */
   p = j(n, 1, .);
   do i = 1 to n;
      run ProbIntegralTransform( g_a[i], g_b[i], dd, ee );
      p[i] = ee - dd;   /* marginal probability CDF(b)-CDF(a) */
   end;

   call sortndx(idx, p, 1);
   idx = t(idx);

   g_a   = g_a[idx];
   g_b   = g_b[idx];
   g_corr = g_corr[idx, idx];

   run CholeskyAndStd( n );
   run ComputeConditionalExpectations( n );
   run ProbIntegralTransform( g_a[1], g_b[1], g_done, g_eone );/* sets these globals...necessary? */
finish covsrt_static;

/* FUNCTION: covsrt_naive
Cholesky decomposition without pivoting, using ROOT with ridge policy. 
I tried this method, but in increased the variance of the estimates.
Strategy:
   1. Mirror covsrt exactly: copy inputs into global work arrays g_corr,g_a,g_b,g_y.
   2. Permute infinite dimensions to the end using the same rcswp calls as covsrt.
   3. Compute Cholesky of the active submatrix via ROOT (two-pass ridge policy).
   4. Scale the factor and limits so each diagonal entry of C equals 1.
   5. Compute conditional expectations g_y[i] for the Genz transform.
Input arguments lower, upper, and covar are not modified (same as covsrt).
*/
start covsrt_naive( lower, upper, covar )
      global( g_corr, g_a, g_b, g_y, g_done, g_eone );
   n = ncol(covar);

   run CholeskyAndStd( n );
   run ComputeConditionalExpectations( n );
   run ProbIntegralTransform( g_a[1], g_b[1], g_done, g_eone );/* sets these globals...necessary? */
finish covsrt_naive;


/* FUNCTION: covsrt
Computes Cholesky Decomposition with Greedy Pivoting.
This is critical for numerical stability and efficiency.
It reorders variables so that the outermost integration variables 
have the smallest expected integration intervals.
*/
start covsrt( lower, upper, covar )
      global( g_corr, g_a, g_b, g_y, g_done, g_eone );
   n = ncol(covar);
v = vecdiag(covar);
IF ANY(v^=1) THEN PRINT "WARNING: covsrt_static is designed for correlation matrices. Results may be inaccurate.", v;
   eps = 1e-10;       /* numerical tolerance in pivot/singularity checks */
   sqtwpi = sqrt( 2*constant("PI") );

   *print "DBG covsrt: after permutation, before Cholesky loop";
   *print n;
   *print g_corr, g_a, g_b;
   /* MAIN CHOLESKY LOOP WITH PIVOTING */
   do i = 1 to n;
      demin = 1;
      jmin = i;
      cvdiag = 0;
      epsi = i*i*eps;

      /* SEARCH LOOP: Find variable 'j' that minimizes expected interval width */
      do j = i to n;
         if ( g_corr[j,j] > epsi ) then do;
            sumsq = sqrt( g_corr[j,j] );
            vsum = 0;

            if ( i > 1 ) then
               vsum = g_corr[j,1:i-1] * g_y[1:i-1];
            aj = ( g_a[j] - vsum )/sumsq;
            bj = ( g_b[j] - vsum )/sumsq;
            run ProbIntegralTransform( aj, bj, dd, ee );

            /* Pivot Criterion: Choose j with smallest probability mass (ee-dd) */
            if ( demin >= ee - dd ) then do;
               jmin = j;
               amin = aj;
               bmin = bj;
               demin = ee - dd;
               cvdiag = sumsq;
            end;
         end;
      end;

      /* Swap rows/cols if a better variable was found */
      if ( jmin > i ) then
         run rcswp( i, jmin, g_a, g_b, n, g_corr );

      /* Perform Cholesky Update for current row */
      if ( cvdiag > 0 ) then do;
         g_corr[i,i] = cvdiag;

         do l = i+1 to n;
            g_corr[l,i] = g_corr[l,i]/cvdiag;
            g_corr[l,i+1:l] = g_corr[l,i+1:l] - g_corr[l,i] # t(g_corr[i+1:l,i]);
         end;

         /* Calculate expected value 'g_y[i]' to center next iteration */
         if ( demin > epsi ) then do;
            yl = 0;
            yu = 0;

            if ( g_a[i] ^= . ) then
               yl = -exp( -amin**2/2 )/sqtwpi;

            if ( g_b[i] ^= . ) then
               yu = -exp( -bmin**2/2 )/sqtwpi;
            g_y[i] = ( yu - yl )/demin;
         end;
         else do;
            if ( g_a[i] = . & g_b[i] ^= . ) then          /* (-inf, b) */
               g_y[i] = bmin;

            if ( g_a[i] ^= . & g_b[i] = . ) then          /* (a, inf) */
               g_y[i] = amin;

            if ( g_a[i] ^= . & g_b[i] ^= . ) then         /* (a, b) */
               g_y[i] = ( amin + bmin )/2;
         end;

         /* Normalize current row */
         g_corr[i,1:i] = g_corr[i,1:i]/cvdiag;
         g_a[i] = g_a[i]/cvdiag;
         g_b[i] = g_b[i]/cvdiag;
      end;
      else do;
         /* Handling Singularities / Semi-definite cases */
         if ( g_corr[i,i] > -epsi ) then do;
            g_corr[i:n,i] = 0;
            aaa = 0;

            /* Back-substitute to resolve dependency */
            do j = i-1 to 1 by -1;
               if ( abs( g_corr[i,j] ) > epsi ) then do;
                  g_a[i] = g_a[i]/g_corr[i,j];
                  g_b[i] = g_b[i]/g_corr[i,j];

                  if ( g_corr[i,j] < 0 ) then do;
                     tmp = g_a[i];
                     g_a[i] = g_b[i];
                     g_b[i] = tmp;
                  end;

                  g_corr[i,1:j] = g_corr[i,1:j]/g_corr[i,j];

                  /* Re-sort rows to maintain triangular structure after singularity fix */
                  do l = j+1 to i-1;
                     if( g_corr[l,j+1] > 0 ) then do;
                        do k = i-1 to l by -1;
                           /* Large block swapping logic for singularity handling */
                           tmp = g_corr[k,1:k];
                           g_corr[k,1:k] = g_corr[k+1,1:k];
                           g_corr[k+1,1:k] = tmp;
                           tmp = g_a[k];
                           g_a[k] = g_a[k+1];
                           g_a[k+1] = tmp;
                           tmp = g_b[k];
                           g_b[k] = g_b[k+1];
                           g_b[k+1] = tmp;
                        end;
                        l = i-1;
                     end;
                  end;
                  j = 1;
                  aaa = 1;
               end;
               
               /* Zero out remaining elements if singularity resolved */
               if aaa = 1 then;
               else g_corr[i,j] = 0;
            end;
            g_y[i] = 0;
         end;
         else do;
            /* Error: Shouldn't happen b/c COVAR is pre-checked to be positive definite */
            i = n;
         end;
      end;
   end;

   *print "DBG covsrt: final state before ProbIntegralTransform";
   *print g_corr, g_a, g_b, g_y;

   run ProbIntegralTransform( g_a[1], g_b[1], g_done, g_eone ); /* sets these globals...necessary? */
finish covsrt;


/* FUNCTION: rcswp
Utility to swap Row 'p' and Row 'q' in the covariance matrix 
and associated limit vectors.
** NOTE: This section contains the user-patched recovery logic. **
*/
start rcswp( p, q, a, b, n, c );
   /* Swap Limits */
   tmp = a[p];   a[p] = a[q];         a[q] = tmp;
   tmp = b[p];   b[p] = b[q];         b[q] = tmp;

   /* Swap Diagonal Elements */
   tmp = c[p,p]; c[p,p] = c[q,q];     c[q,q] = tmp;

   /* Swap Columns to the left of p */
   if (p>1) then do;
      tmp = c[q,1:p-1];
      c[q,1:p-1] = c[p,1:p-1];
      c[p,1:p-1] = tmp;
   end;

   /* Swap Elements between p and q (the rectangular block) */
   do i = p+1 to q-1;
      tmp = c[i,p];
      c[i,p] = c[q,i];
      c[q,i] = tmp;
   end;

   /* Swap columns below q */
   if (q<n) then do;
      tmp = c[q+1:n,p];
      c[q+1:n,p] = c[q+1:n,q];
      c[q+1:n,q] = tmp;
   end;
finish rcswp;


/* FUNCTION: dkbvrc
Main Adaptive Integration Loop using Randomized Korobov Rules.
1. Estimates integral value and error.
2. Increases number of points or changes prime base if error > tolerance.
*/
start dkbvrc( ndim, minvls, maxvls, abseps, releps,   abserr, finest );

   minsmp = 8;           /* minimum number of random shifts per stage */

      /* Korobov lattice primes and generator coefficients (local to QMC driver). */
      p_vector = { 31 47 73 113 173 263 397 593 907 1361 2053 3079 4621 6947 10427 15641 23473 35221 52837 79259 118891 178349 267523 401287 601942};
      mat = { 12 9 9 13 12 12 12 12 12 12 12 12 3 3 3 12 7 7 12,
         13 11 17 10 15 15 15 15 15 15 22 15 15 6 6 6 15 15 9,
         27 28 10 11 11 20 11 11 28 13 13 28 13 13 13 14 14 14 14,
         35 27 27 36 22 29 29 20 45 5 5 5 21 21 21 21 21 21 21,
         64 66 28 28 44 44 55 67 10 10 10 10 10 10 38 38 10 10 10,
         111 42 54 118 20 31 31 72 17 94 14 14 11 14 14 14 94 10 10,
         163 154 83 43 82 92 150 59 76 76 47 11 11 100 131 116 116 116 116,
         246 189 242 102 250 250 102 250 280 118 196 118 191 215 121 121 49 49 49,
         347 402 322 418 215 220 339 339 339 337 218 315 315 315 315 167 167 167 167,
         505 220 601 644 612 160 206 206 206 422 134 518 134 134 518 652 382 206 158,
         794 325 960 528 247 247 338 366 847 753 753 236 334 334 461 711 652 381 381,
         1189 888 259 1082 725 811 636 965 497 497 1490 1490 392 1291 508 508 1291 1291 508,
         1763 1018 1500 432 1332 2203 126 2240 1719 1284 878 1983 266 266 266 266 747 747 127,
         2872 3233 1534 2941 2910 393 1796 919 446 919 919 1117 103 103 103 103 103 103 103,
         4309 3758 4034 1963 730 642 1502 2246 3834 1511 1102 1102 1522 1522 3427 3427 3928 915 915,
         6610 6977 1686 3819 2314 5647 3953 3614 5115 423 423 5408 7426 423 423 487 6227 2660 6227,
         9861 3647 4073 2535 3430 9865 2830 9328 4320 5913 10365 8272 3706 6186 7806 7806 7806 8610 2563,
         10327 7582 7124 8214 9600 10271 10193 10800 9086 2365 4409 13812 5661 9344 9344 10362 9344 9344 8585,
         19540 19926 11582 11113 24585 8726 17218 419 4918 4918 4918 15701 17710 4037 4037 15808 11401 19398 25950,
         34566 9579 12654 26856 37873 38806 29501 17271 3663 10763 18955 1298 26560 17132 17132 4753 4753 8713 18624,
         31929 49367 10982 3527 27066 13226 56010 18911 40574 20767 20767 9686 47603 47603 11736 11736 41601 12888 32948,
         40701 69087 77576 64590 39397 33179 10858 38935 43129 35468 35468 2196 61518 61518 27945 70975 70975 86478 86478,
         103650 125480 59978 46875 77172 83021 126904 14541 56299 43636 11655 52680 88549 29804 101894 113675 48040 113675 34987,
         165843 90647 59925 189541 67647 74795 68365 167485 143918 74912 167289 75517 8148 172106 126159 35867 35867 35867 121694,
         130365 236711 110235 125699 56483 93735 234469 60549 1291 93937 245291 196061 258647 162489 176631 204895 73353 172319 28881};

      plim = ncol(p_vector);
      klim = ncol(mat) + 1;

   vk   = j( 1, klim, . );
   intvls = 0;
   klimi  = klim;
   
   /* RECOVERED SECTION: Entry check */
   if ( minvls >= 0 ) then do;      
      finest = 0;
      varest = 0;
      sampls = minsmp;

      /* Determine starting prime index based on input 'minvls' */
      do i = 1 to plim;
         np = i;
         if ( minvls < 2*sampls*p_vector[i] ) then
            i = plim;
      end;

      if ( minvls >= 2*sampls*p_vector[plim] ) then
         sampls = minvls/( 2*p_vector[plim] );
   end;

   value = j( 1, sampls, . );
   exit = 0;

   /* ADAPTIVE LOOP */
   do until( exit = 1);
      /* Initialize Lattice Rule Generator */
      vk[1] = 1/p_vector[np];
      do i = 2 to min( ndim, klim );
         vk[i] = mod( mat[np, min(ndim-1,klim-1)]*vk[i-1], 1 );
      end;

      finval = 0;
      varsqr = 0;

      /* Compute Integral Estimate (value) */
      do i = 1 to sampls;
         value[i] = dksmrc( ndim, klimi, p_vector[np], vk );
      end;

      /* Compute Variance and Error Estimates */
      finval = value[:];
      varsqr = (value[##] - value[+]##2/sampls) / (sampls # (sampls-1));
      intvls = intvls + 2*sampls*p_vector[np];
      varprd = varest*varsqr;
      
      /* Weighted update of final estimate */
      finest = finest + ( finval - finest )/( 1 + varprd );

      if ( varsqr > 0 ) then
         varest = ( 1 + varprd )/varsqr;
      if varsqr <= 0 then
         abserr = 0;
      else 
         abserr = 3*sqrt( varsqr/( 1 + varprd ) );

      /* Check Convergence Criteria */
      if ( abserr > max( abseps, abs(finest)*releps ) ) then do;
         /* If not converged, increase Prime index (np) or Sample size */
         if ( np < plim ) then
            np = np + 1;
         else do;
            sampls = min( 3*sampls/2, ( maxvls - intvls )/( 2*p_vector[np] ) );
            sampls = max( minsmp, sampls );
         end;

         if ( intvls + 2*sampls*p_vector[np] > maxvls ) then
            exit = 1;
      end;
      else do;
         exit = 1;
      end;
   end;
finish dkbvrc;


/* FUNCTION: dksmrc
Computes the sum for the Randomized Lattice Rule.
Uses the Periodizing Transformation (Baker's/Tent Transform).
Includes Antithetic Variates (x and 1-x).
*/
start dksmrc( ndim, klim, prime, vk );
   x = j( 1, ndim, . );
   sumkro = 0;
   nk = min( ndim, klim );
   
   /* Random Shift (Cranley-Patterson Rotation) */
   u = randfun(nk-1, "Uniform");
   do j = 1 to nk-1;
      jp = j + u[j]*( nk + 1 - j );
      xt = vk[j];
      vk[j] = vk[int(jp)];
      vk[int(jp)] = xt;
   end;

   xp = randfun(ndim, "Uniform");
   diff = ndim - klim;

   /* Loop over lattice points */
   do k = 1 to prime;
      /* Handle dimensions beyond the lattice limit using Richtmyer generators */
      if ( diff>0 ) then do;
         run dkrcht( diff, x );
         do j = 1 to ndim-klim;
            x[nk+j] = x[j];
         end;
      end;

      x[1:nk] = mod( k*vk[1:nk], 1 );

      /* Apply Shift and Baker (Tent) Transformation */
      do j = 1 to ndim;
         xt = x[j] + xp[j];

         if ( xt > 1 ) then
            xt = xt - 1;
         x[j] = abs( 2*xt - 1 ); /* Tent Transform */
      end;

      /* Evaluate Integrand + Antithetic Variate (1-x) */
      mvndfn = mvn_dfn( ndim, x );
      sumkro = sumkro + ( mvndfn - sumkro )/( 2*k - 1 );
      x = 1 - x;
      mvndfn = mvn_dfn( ndim, x );
      sumkro = sumkro + ( mvndfn - sumkro )/( 2*k );
   end;

   return (sumkro);
finish dksmrc;


/* FUNCTION: dkrcht
Richtmyer Generator.
Generates quasi-random numbers using square roots of primes.
Used for dimensions exceeding the optimal lattice rule coefficients.
*/
start dkrcht( s, quasi ) global( g_nn, g_hisum, g_olds );
   /* Local read-only constants for quasi-random generation */
   mxdim = 80;        /* max number of dimensions supported */
   mxhsum = 50;       /* max depth of digit expansion */
   bb = 2;            /* base for Richtmyer counter */
   primes = {2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71
      73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173
      179 181 191 193 197 199 211 223 227 229 233 239 241 251 257 263 269 271 277 281
      283 293 307 311 313 317 331 337 347 349 353 359 367 373 379 383 389 397 401 409};
   psqt = sqrt(primes);
   
   if ( s ^= g_olds | s < 1 ) then do;
      g_olds = s;
      g_nn[1] = 0;
      g_hisum = 0;
   end;

   i = 0;
   crit = 0;

   /* Update counters for the sequence */
   do until( crit = 1 | i = g_hisum + 1 );
      g_nn[i + 1] = g_nn[i + 1] + 1;

      if ( g_nn[i + 1] < bb ) then do;
         crit = 1;
         i = i - 1;
      end;
      else g_nn[i + 1] = 0;
      i = i + 1;
   end;

   if ( i > g_hisum ) then do;
      g_hisum = g_hisum + 1;

      if ( g_hisum > mxhsum ) then
         g_hisum = 0;
      g_nn[g_hisum + 1] = 1;
   end;

   rn = 0;
   do i = g_hisum to 0 by -1;
      rn = g_nn[i + 1] + bb * rn;
   end;

   quasi[1:s] = mod( rn # psqt[1:s], 1 );
finish dkrcht;


/* define helper functions used in the various covsrt* modules */
/* Helper: Cholesky factorization via ROOT with ridge policy, followed by a standardization */
start CholeskyAndStd( n )
   global( g_corr, g_a, g_b );
   
   R = g_corr;
   
   /* Pass 1: try ROOT without ridge */
   G = root(R, "NoError");

   /* Pass 2: Retry with small ridge if not PD; prevent division by zero */
   ridge_eps = 1E-12;
   if any(vecdiag(G) < ridge_eps) then do;
      G = root(R + ridge_eps*I(n), "NoError");
      /* TO DO: Gracefully return without calling ABORT */
      if any(G=.) then
         ABORT "CholeskyAndStd: ROOT fails even with ridge=1E-12";
   end;

   /* C is the lower-triangular Cholesky factor.
      Scale C and limits so diagonal of C is 1 */
   C = t(G); 
   v = vecdiag(C);
   do i = 1 to n;
      C[i,] = C[i,] / v[i];
      g_a[i]  = g_a[i]  / v[i];
      g_b[i]  = g_b[i]  / v[i];
   end;
   g_corr = C;
   /* this looks wrong; the Cholesky should be triangular */
   g_corr = CopyLowerTriToUpper(g_corr); /* from this point on, g_corr is the Cholesky factor */
finish CholeskyAndStd;

/* Helper: Compute conditional expectations for the Genz transformation
           These are put into the g_y array */
start ComputeConditionalExpectations( n )
   global( g_corr, g_a, g_b, g_y, g_done, g_eone );
   
   sqtwpi = sqrt( 2*constant("PI") );
   do i = 1 to n;
      if (i > 1) then vsum = g_corr[i,1:i-1] * g_y[1:i-1];
      else vsum = 0;
      
      ai = g_a[i] - vsum;
      bi = g_b[i] - vsum;

      /* this is the same calc as ProbIntegralTransform */
      di = 0; ei = 1;
      if ( g_a[i] ^= . ) then di = cdf("Normal", ai);
      if ( g_b[i] ^= . ) then ei = cdf("Normal", bi);
      ei = max(ei, di);

      if ( ei - di > 1e-10 ) then do;
         yl = 0; yu = 0;
         /* be careful here: better to compute the pdf values directly from ai and bi.
            The following calculations might overflow. Consider using logs. */
         if ( g_a[i] ^= . ) then yl = -exp(-ai**2/2)/sqtwpi;
         if ( g_b[i] ^= . ) then yu = -exp(-bi**2/2)/sqtwpi;
         g_y[i] = (yu - yl) / (ei - di);
      end;
      else do;
         if      ( g_a[i] = .  & g_b[i] ^= . ) then g_y[i] = bi;  /* flag 0: (-inf, b) */
         else if ( g_a[i] ^= . & g_b[i] = .  ) then g_y[i] = ai;  /* flag 1: (a, inf) */
         else if ( g_a[i] ^= . & g_b[i] ^= . ) then g_y[i] = (ai + bi)/2;  /* flag 2: (a, b) */
      end;
   end;
   
   run ProbIntegralTransform( g_a[1], g_b[1], g_done, g_eone );
finish ComputeConditionalExpectations;

store module=(
probmvn_mod
probmvn_std
RemoveInfiniteLimits
mvn_dist
mvn_dfn
mvndnt
ProbIntegralTransform
covsrt
covsrt_naive
covsrt_static
rcswp
dkbvrc
dksmrc
dkrcht
CopyLowerTriToUpper
CholeskyAndStd
ComputeConditionalExpectations
);
