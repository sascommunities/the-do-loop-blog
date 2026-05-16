/* COMMON utility functions for validating parameters to PROBMVN and CDFMVN.
   These functions are loaded by tests for PROBMVN and CDFMVN.
   Functions defined here:
   1. PrintToLog(msg, errCode) -- prints a message to the SAS log with an optional error code (0=note, 1=warning, 2=error)
   2. ErrorToLog(msg)          -- prints an error message to the SAS log
   3. IsSym(A)                 -- returns 1 if A is a symmetric square matrix
   4. IsSPD(M)                 -- returns 1 if M is symmetric and positive definite
   5. IsCorr(M)                -- returns 1 if M is a correlation matrix (SPD with unit diagonal)
*/
/* Check the SYSVER macro to see if SAS 9.4 is running.
   In SAS Viya, the macro is empty and does nothing.
   In SAS 9.4, the macro defines a function that emulates the PrintToLog call.
   The syntax is as follows:
   call PrintToLog("This is a log message.");
   call PrintToLog("This is a note.", 0);
   call PrintToLog("This is a warning.", 1);
   call PrintToLog("This is an error.", 2);
*/
%macro DefinePrintToLog;
%if %sysevalf(&sysver = 9.4) %then %do;
start PrintToLog(msg,errCode=-1);
   if      errCode=0 then prefix = "NOTE: ";
   else if errCode=1 then prefix = "WARNING: ";
   else if errCode=2 then prefix = "ERROR: ";
   else prefix = "";
   stmt = '%put ' + prefix + msg + ';';
   call execute(stmt);
finish;
store module=(PrintToLog);
%end;
start ErrorToLog(msg);
   run PrintToLog(msg, 2);
finish;
store module=(ErrorToLog);
%mend;

/* Common validation functions shared by CDFMVN and PROBMVN.
   Functions defined here:
   IsSym(A)   -- returns 1 if A is a symmetric square matrix
   IsSPD(M)   -- returns 1 if M is symmetric and positive definite
   IsCorr(M)  -- returns 1 if M is a correlation matrix (SPD with unit diagonal)
*/
proc iml;
%DefinePrintToLog;

start IsSym(A);
   if nrow(A) ^= ncol(A) then return(0);    /* A is not square */
   c = max(abs(A));
   sqrteps = constant('SqrtMacEps');
   return( all( abs(A-A`) < c*sqrteps ) );
finish;

start IsSPD(M);
   if ^IsSym(M) then return( 0 );
   U = root(M, "NoError");
   if any(U=.) then return( 0 );
   return( 1 );
finish;

start IsCorr(M);
   if ^IsSPD(M) then return( 0 );
   if any(vecdiag(M) ^= 1) then return ( 0 );
   return( 1 );
finish;

/* Sigma is a kxk covariance matrix and b is a row vector with k elements.
   Convert b to inv(D)*(b-mu), where D is the diagonal matrix of
   standard deviations: sqrt(vecdiag(Sigma))
   Note: b can have multiple rows. b can also have missing values.
*/
start Xform_Limits_Cov2Corr( b, Sigma, mu=j(1,ncol(Sigma),0) );
   D = sqrt(vecdiag(Sigma));
   c = (b - mu) /rowvec(D);
   return c;
finish;

store module=(IsSym IsSPD IsCorr Xform_Limits_Cov2Corr);
QUIT;
