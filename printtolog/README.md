# The PrintToLog Subroutine

## Overview

The PrintToLog subroutine is discussed at ["Write to the log from SAS IML programs"](https://blogs.sas.com/content/iml/2023/01/25/printtolog-iml.html).
The subroutine is supported on SAS Viya in both PROC IML and the iml action. In SAS 9.4, you can use the %INCLUDE statement to include this file in a PROC IML program.

Example"

```sas
proc iml;
/* define the PrintToLog subroutine */
%include "&path/printtolog.sas";
/* call the PrintToLog subroutine */
run PrintToLog("This is a test message.", 0);
quit;
```
