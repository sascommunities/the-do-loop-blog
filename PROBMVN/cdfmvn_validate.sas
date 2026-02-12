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

/* this program runs in SAS 9.4 or in SAS Viya */
proc iml;
%DefinePrintToLog; 

/* validate the arguments for CDFMVN:
   b does not support missing values
   Sigma is SPD
*/
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

store module=(IsSym IsSPD IsCorr IsValidParmsCDF IsValidParmsTVN);
QUIT;
