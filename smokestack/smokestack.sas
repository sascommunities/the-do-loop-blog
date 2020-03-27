/* SAS program to accompany the article 
   "Smokestack plots: A visualization technique for comparing cumulative curves"
   by Rick Wicklin, published 30MAR2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/03/30/smokestack-plots-cumulative-curves.html

   This program shows how to create a "smokestack plot," which is a plot of
   multiple cumulative curves.

   To generate the data, run the file smokestackDATA.sas.
*/

proc contents data=Have short varnum; run;

proc freq data=Have order=freq;
table ID;
run;

/* Compute cumulative quantities */
data Cumul;
set Have;                  /* assume sorted by Date for each ID variable */
by ID notsorted;           /* use NOTSORTED option if not sorted by ID */
if first.ID then do;
   Cumul=0;                /* initialize the cumulative count */
   output;                 /* Optional: to display 0, include the OUTPUT statements */
end;
Cumul + Value;             /* increment the cumulative count */
output;
run;

ods graphics / reset width=640px height=200px;
title "Events versus Time";
title2 "Views for Blog Post ID=25559";
proc sgplot data=Cumul;
   where ID=25559;
   needle x=Date y=Value / lineattrs=(thickness=5);;
   xaxis display=(nolabel);
   yaxis grid label="New Views";
run;


/* with < 10 groups, you can panel these. At 150 pixels
   per plot, you need 150*numGroups pixels */
ods graphics / width=640px height=800px;
title "Events versus Time";
title2 "Views for Blog Posts, excluding ID=25559";
proc sgpanel data=cumul;
   where ID^=25559;
   panelby ID / layout=rowlattice onepanel sort=data;
   needle x=Date y=Value/ lineattrs=(thickness=5);
   colaxis display=(nolabel);
   rowaxis grid label="New Views";
run;

ods graphics / width=640px height=480px;
title "Smokestack Plot";
title2 "Overlay Multiple Cumulative Curves";
proc sgplot data=Cumul;
*   where ID^=25559;
   series x=Date y=Cumul / group=ID lineattrs=(pattern=solid) curvelabel; 
   xaxis grid display=(nolabel);
   yaxis grid label="Cumulative Views" values=(0 to 1000 by 100) valueshint;
run;

data Origin;
retain BaseDate;
set Cumul;                 /* start with the cumulative data */
by ID notsorted;           /* use the NOTSORTED option if not sorted by ID */
if first.ID then           /* remember the first day for this ID */
   BaseDate = Date;
Day = Date - BaseDate;     /* days since first day */
drop BaseDate;
run;

title2 "Cumulative Curves with a Common Origin";
proc sgplot data=Origin;
   series x=Day y=Cumul / group=ID lineattrs=(pattern=solid) curvelabel; 
   xaxis grid label="Days Since Publication" values=(0 to 100 by 10) valueshint;
   yaxis grid label="Cumulative Views" values=(0 to 1000 by 100) valueshint;
run;

title2 "Cumulative Curves (log scale)";
proc sgplot data=Origin;
   where Cumul > 0;
   series x=Day y=Cumul / group=ID lineattrs=(pattern=solid) curvelabel; 
   xaxis grid label="Days Since Publication" values=(0 to 100 by 10) valueshint;
   yaxis grid type=log logbase=10     /* LOG10 transformation of axis */
         label="Cumulative Views" values=(100 to 1000 by 100) valueshint;
run;

