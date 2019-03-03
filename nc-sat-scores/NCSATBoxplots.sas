/* SAS program to accompany the article 
   "Use PROC BOXPLOT to display hundreds of box plots"
   by Rick Wicklin, published 06MAR2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/03/06/proc-boxplot-hundreds.html

   This program uses PROC BOXPLOT to visualize the distribution of 
   SAT scores for public high schools in NC in 2018. The graphs include
   1. A box plot that displays a table of descriptive statitics
   2. A panel of box plots that displays 115 boxes across five graphs

   Data from NC Department of Public Instruction. The web site
   http://www.ncpublicschools.org/accountability/reporting/sat/2018
   contains an Excel spread sheet "2018 SAT Performance by District and School"
   You should download that fileas a CSV file to
   NCSAT2018.csv
   in the MYPATH direcory
*/
option PS=999;                             /* make pages extra long to support tall graphs */
%let MYPATH = C:\Users\frwick\Downloads;   /* specify directory for CSV file */
/* Import CSV to SAS data set */
data NCSAT;
infile "&MYPATH\NCSAT2018.csv" delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
informat District School $40. ;
informat _Tested _PctTested _Total $4. ;
informat _ReadWrite _Math $3. ;
length Charter $ 3;
input District $ School $ _Tested _PctTested _Total _ReadWrite _Math;
/* missing values are represented by strings '<10' and '*' */
if _Total='*' then do;
  Tested = .F; PctTested=.; Total=.; ReadWrite=.; Math=.;
end;
else do;
  Tested = _Tested; PctTested=_PctTested; Total=_Total; 
  ReadWrite=_ReadWrite; Math=_Math;
end;
drop _:;
District = tranwrd(District, "Public Schools Of R", "R");
District = tranwrd(District, "Chapel Hill", "Chapel-Hill"); /* handle Chapel Hill */
DistrictAbbr = scan(District, 1, ' ');   /* grab first word as abbreviation */
DistrictAbbr = tranwrd(DistrictAbbr, "-", " ");   /* handle Chapel Hill */
if School ^= "District Average" AND Total ^= .;
Charter = ifc(District="Charter", "Yes", "No ");  /* create indicator variable */
label Total="Average Total SAT Score"
      Math="Average SAT Math Score"
      ReadWrite="Average SAT Reading/Writing Score"
      DistrictAbbr = "District";
run;

/* If you want, you can sort by a statistic, such as median score */
proc sort data=NCSAT; by DistrictAbbr; run;
proc means data=NCSAT noprint;
   class DistrictAbbr;
   var Total;
   output out=MeansOut(where=(_Type_=1)) Median=MedianSAT;
run;
data SATSortMerge;
   merge NCSAT MeansOut;
   by DistrictAbbr;
run;
proc sort data=SATSortMerge; by DESCENDING MedianSAT DistrictAbbr; run;

/***************************************************/
/* Add a statistical table to a graph of box plots */
/***************************************************/
proc freq data=NCSAT order=freq;
   table DistrictAbbr / maxlevels=20;
run;

ods graphics / width=700px height=480px;
title "Average SAT Scores for Large NC School Districts";
proc boxplot data=SATSortMerge;
   where _FREQ_ >= 7;
   plot Total*DistrictAbbr / grid odstitle=title nohlabel
      boxstyle=schematicID vaxis=800 to 1450 by 50;
   insetgroup Q2 N;
run;

ods graphics / width=640px height=400px;
title "Average SAT Scores for NC School Districts";
proc boxplot data=SATSortMerge;
   plot Total*DistrictAbbr / grid odstitle=title nohlabel
      boxstyle=schematicID vaxis=800 to 1450 by 50;
run;
