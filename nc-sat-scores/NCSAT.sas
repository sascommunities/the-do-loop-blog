/* SAS program to accompany the article 
   "Visualize SAT scores in North Carolina"
   by Rick Wicklin, published 04MAR2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/03/04/visualize-sat-scores-nc.html 

   This program visualizes the distribution of SAT scores for 
   public high schools in NC in 2018. The graphs include
   1. A histogram of total SAT scores
   2. A comparative histogram of math and English SAT scores
   3. A scatter plot of math and English SAT scores
   4. A scatter plot of total SAT score by school district

   Data from NC Department of Public Instruction, The web site
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

/* transpose data for the comparative histogram */
data SATTrans;
set NCSAT;
length Test $ 7;
Test = "Math";    Score = Math;      output;
Test = "English"; Score = ReadWrite; output;
run;

/* If you want, you can sort by a statistic, such as median score */
proc sort data=NCSAT; by DistrictAbbr; run;
proc means data=NCSAT noprint;
   class DistrictAbbr;
   var Total;
   output out=MeansOut(where=(_Type_=1)) Median=MedianSAT;
run;
data All;
merge NCSAT MeansOut;
by DistrictAbbr;
run;
proc sort data=All; by DESCENDING MedianSAT DistrictAbbr; run;

/************************************************/
/* Visualize histogram of Total SAT scores */
/************************************************/
ods graphics/reset;
title "Average SAT Scores at NC High Schools";
title2 "Class of 2018";
footnote J=L 'Source: NC Department of Public Instruction';
proc sgplot data=NCSAT;
   histogram Total / binwidth=25 binstart=800; 
   xaxis values=(800 to 1500 by 100);
run;

footnote;
proc means data=NCSAT;
var Total Math ReadWrite;
run;

title "High Performing Schools";
proc print data=NCSAT;
where Total > 1275;
var District School Tested Total;
run;

proc freq data=NCSAT;
tables Charter;
run;


/************************************************/
/* Visualize comparative histogram of Total SAT scores */
/************************************************/
ods graphics / width=640px height=480px Attrpriority=NONE IMAGEMAP;
title "Average SAT Scores at NC High Schools";
footnote J=L 'Source: NC Department of Public Instruction';

/* https://blogs.sas.com/content/iml/2016/03/09/comparative-panel-overlay-histograms-sas.html */
/* compare math vs English scores */
proc sgpanel data=SATTrans;
   panelby Test / columns=1;
   histogram Score / binwidth=25; 
   colaxis grid values=(400 to 750 by 25) valueshint min=375;
run;

/* compare charter vs non-charter schools */
proc sgpanel data=SATTrans;
   panelby Charter / columns=1;
   histogram Score / binwidth=25; 
   colaxis grid values=(400 to 750 by 25) valueshint min=375;
run;

/* visualize joint distribution of math and English scores */
title "Average SAT Scores at NC High Schools";
proc sgplot data=NCSAT;
   styleattrs datasymbols=(Circle TriangleFilled);
   scatter x=Math y=ReadWrite / group=Charter
           tip=(School Total Math ReadWrite District Tested); 
   xaxis grid values=(400 to 750 by 50);
   yaxis grid values=(400 to 750 by 50);
   keylegend / location=inside;
run;
footnote;

/* various counts: number of school districts, range, ... */
proc iml;
use NCSAT;
read all var {"DistrictAbbr" "Total"};
close;
nSchools=nrow(Total);
print nSchools;

u = unique(DistrictAbbr);
numDistricts = ncol(u);
print numDistricts;

idx = loc(Total >= 1000 & Total <= 1200);
nSchoolsInRange = ncol(idx);
print nSchoolsInRange;
pct = nSchoolsInRange / nSchools;
print pct;

use All where (MedianSAT >= 1035);
read all var {"DistrictAbbr" "Total"};
close;
nSchools = nrow(Total);
nDistricts = ncol(unique(DistrictAbbr));
print nSchools nDistricts;
quit;

/************************************************/
/* Visualize Total SAT scores by school district*/
/************************************************/
  
title "Average SAT Scores at NC High Schools";
title2 "By School District (where Median SAT > 1035)";
footnote;
ods graphics / width=600px height=1000px;

/* viz Y axis with many discrete categories:
   https://blogs.sas.com/content/iml/2017/08/16/pairwise-correlations-bar-chart.html
*/
proc sgplot data=All;
where MedianSAT >= 1035;
scatter x=Total y=DistrictAbbr / markerattrs=(symbol=CircleFilled) transparency=0.4;
xaxis grid values=(800 to 1500 by 100) valueattrs=(size=8pt) labelattrs=(size=8pt);
yaxis discreteorder=data reverse display=(nolabel) 
      valueattrs=(size=6pt) fitpolicy=none 
      offsetmin=0.012 offsetmax=0.012  /* half of 1/k, where k=number of catgories */
      colorbands=even colorbandsattrs=(color=gray transparency=0.9);
run;
