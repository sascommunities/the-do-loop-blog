/* This program defines validation functions for the PROBMVN_MOD function, which
   computes P(L < X < U) for a multivariate normal random vector X ~ MVN(mu, Sigma).

   Key difference from CDFMVN: L and U are allowed to contain missing values.
   A missing value in L[i] represents -Infinity, and a missing value in U[i]
   represents +Infinity. All other parameters (Sigma, mu) cannot be missing.

   Functions defined here:
   IsValidRectLimits(L, U)   -- returns 1 if L and U are compatible limit vectors;
                                missing values in L or U are allowed
   IsValidParmsPROBMVN(L, U, Sigma, mu) -- returns 1 if all parameters are valid
*/
proc iml;
/* IsValidRectLimits: validate the lower (L) and upper (U) limit vectors for PROBMVN_MOD.
   Missing values are explicitly allowed: a missing L[i] represents -Infinity
   and a missing U[i] represents +Infinity.
   Returns 1 if valid, 0 otherwise.
   Checks:
   1. L and U have the same number of elements.
   2. For every dimension i where both L[i] and U[i] are non-missing, L[i] <= U[i].
*/
start IsValidRectLimits(L, U);
   if ncol(L) ^= ncol(U) then do;
      run ErrorToLog("The L and U parameters must have the same number of columns.");
      return( 0 );
   end;
   if nrow(L) ^= nrow(U) then do;
      run ErrorToLog("The L and U parameters must have the same number of rows.");
      return( 0 );
   end;
   n = ncol(L);
   do j = 1 to nrow(L);
       do i = 1 to n;
          if L[j,i] ^= . & U[j,i] ^= . then do;
             if L[j,i] > U[j,i] then do;
                run ErrorToLog("L[j,i] must be <= U[j,i] for every non-missing pair of limits.");
                return( 0 );
             end;
          end;
       end;
   end;
   return( 1 );
finish;

/* IsValidParmsPROBMVN: validate all input arguments to PROBMVN_MOD.
   Parameters:
     L     -- 1xn row vector of lower integration limits (missing = -Infinity)
     U     -- 1xn row vector of upper integration limits (missing = +Infinity)
     Sigma -- nxn covariance matrix; must be symmetric and positive definite
     mu    -- 1xn mean vector; no missing values allowed
   Returns 1 if all checks pass, 0 otherwise.
*/
start IsValidParmsPROBMVN(L, U, Sigma, mu=j(1,ncol(Sigma),0));
   n = ncol(Sigma);

   /* 1. Sigma must be symmetric */
   if ^IsSym(Sigma) then do;
      run ErrorToLog("The Sigma parameter must be symmetric.");
      return( 0 );
   end;

   /* 2. Sigma cannot contain missing values */
   if any(Sigma=.) then do;
      run ErrorToLog("The Sigma parameter cannot contain missing values.");
      return( 0 );
   end;

   /* 3. mu cannot contain missing values */
   if any(mu=.) then do;
      run ErrorToLog("The mu parameter cannot contain missing values.");
      return( 0 );
   end;

   /* 4. Dimensional compatibility: L, U, mu vs Sigma */
   if ncol(L) ^= n then do;
      run ErrorToLog("The L and Sigma parameters have incompatible dimensions.");
      return( 0 );
   end;
   if ncol(U) ^= n then do;
      run ErrorToLog("The U and Sigma parameters have incompatible dimensions.");
      return( 0 );
   end;
   if ncol(mu) ^= n then do;
      run ErrorToLog("The mu and Sigma parameters have incompatible dimensions.");
      return( 0 );
   end;

   /* 5. L and U must be internally consistent (missing values allowed) */
   if ^IsValidRectLimits(L, U) then do;
      run ErrorToLog("The L and U parameters must satisfy L[i] < U[i].");
      return( 0 );
   end;

   /* 6. Sigma must be positive definite */
   if ^IsSPD(Sigma) then do;
      run ErrorToLog("The Sigma parameter must be positive definite.");
      return( 0 );
   end;

   /* 7. Dimension must be between 1 and 100 */
   if n < 2 | n > 100 then do;
      run ErrorToLog("PROBMVN supports problems between 2 and 100 dimensions.");
      return( 0 );
   end;

   return( 1 );
finish;

store module=(IsValidRectLimits IsValidParmsPROBMVN);
QUIT;
