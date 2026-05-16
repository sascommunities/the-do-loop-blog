/* Common matrix validation functions (IsSym, IsSPD, IsCorr) and PrintToLog/ErrorToLog
   are defined and stored by mvn_validate.sas. */
*%include "mvn_validate.sas";

/* CDF-specific validation functions for CDFMVN and CDFTVN. */
proc iml;
start IsValidParmsCDF(b, Sigma, mu);
   if ^IsSym(Sigma) then do;
      run ErrorToLog( "The Sigma parameter must be symmetric.");
      return( 0 );
   end;
   if ncol(b) ^= ncol(Sigma) then do;
      run ErrorToLog( "The b and Sigma parameters are not compatable dimensions.");
      return( 0 );
   end;
   if ncol(mu) ^= ncol(Sigma) then do;
      run ErrorToLog( "The mu and Sigma parameters are not compatable dimensions.");
      return( 0 );
   end;
   if any(b=.) | any(Sigma=.) | any(mu=.) then do;
      run ErrorToLog( "The parameters cannot contain missing values.");
      return( 0 );
   end;
   if ^IsSPD(Sigma) then do;
      run ErrorToLog( "The Sigma parameter must be positive definite.");
      return( 0 );
   end;
   return( 1 );
finish;

start IsValidParmsTVN(b, Sigma, mu);
   isValid = IsValidParmsCDF(b, Sigma, mu);
   if ^isValid then return( 0 );
   if ncol(b) ^= 3 then do;
      run ErrorToLog( "cdftvn only supports 3-dimensional input.");
      return( 0 );
   end;
   return( 1 );
finish;

start IsValidParmsCDFMVN(b, Sigma, mu);
   isValid = IsValidParmsCDF(b, Sigma, mu);
   if ^isValid then return( 0 );
   if ncol(b)<2 | ncol(b) > 32 then do;
      run ErrorToLog( "CDFMVN supports problems between 2 and 32 dimensions.");
      return( 0 );
   end;
   return( 1 );
finish;

store module=(IsValidParmsCDF IsValidParmsTVN IsValidParmsCDFMVN);
QUIT;
