/* SAS program to accompany the article 
   "Write to the log from SAS IML programs"
   by Rick Wicklin, published 25JAN2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/01/25/printtolog-iml.html

   The SAS IML product in SAS Viya supports the PrintToLog function, which writes
   messages to the log from SAS IML programs (PROC IML or the iml action).
   This program defines a PrintToLog function that runs in PROC IML in SAS 9.4. The 
   function has the same syntax and functionality.
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
%mend;
%DefinePrintToLog;
