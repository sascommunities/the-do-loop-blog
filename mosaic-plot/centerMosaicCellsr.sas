/*********************************************************/
/* SAS program to accompany the article 
   "Find the center of each cell in a mosaic plot"
   by Rick Wicklin, published 10JUL2019 on The DO Loop blog:
   https:

   This program shows how to find the center of cells
   in a mosaic plot in SAS. You can us the centers to add annotations.
*/

/* Get the summary statistics. Use FreqOut for mosaic plot
   but use FreqList (same info, different format) 
   to obtain frequencies, percentages, and labels for annotation */

%let hVar = Type;             /* name of horizontal variable */
%let vVar = Origin;           /* name of vertical variable */

proc freq data=Sashelp.Cars;
   where Type ^= 'Hybrid';
   tables &vVar * &hVar / crosslist sparse 
                          out=FreqOut(where=(Percent^=.));
   ods output CrossList=FreqList;
run;

proc contents data=FreqList varnum ;run;

title "Horizontal Percentages: &hVar";
proc print data=FreqList;
where F_&vVar='Total' & F_&hVar^='Total';
var F_&vVar F_&hVar Frequency Percent;
run;

title "Vertical Percentages: &vVar";
proc print data=FreqList;
by F_&vVar notsorted;
var F_&vVar F_&hVar Frequency Percent ColPercent;
where F_&vVar ^= 'Total' & F_&hVar ^= 'Total';
run;

/* Read CROSSLIST data set and write data set that contains centers of cells in mosaic plot.
  Main idea: If a categorical variable has three levels and observed proportions are
  20, 30, and 50, then the midpoints of the bars are
       10 = 20/2
       35 = 20 + 30/2
       75 = 20 + 30 + 50/2
*/
proc iml;
/* 1. read percentages for the horizontal variable */
use FreqList;
read all var {Percent F_&hVar F_&vVar}
         where (F_&vVar='Total' & F_&hVar^='Total');
hPos = Percent;            
nH = nrow(hPos);
hLabel = F_&hVar;

/* 2. horizontal centers
   h1/2,  h1 + h2/2, h1 + h2 + h3/2, h1 + h2 + h3 + h4/2, ... */
halfW = hPos / 2;
hCenter = halfW + cusum(0 // hPos[1:(nH-1)]);
print (hCenter`)[c=hLabel L="hCenter" F=5.2];

/* 3. For each column, read cell percentages and frequencies */
read all var {Frequency Percent ColPercent F_&vVar F_&hVar}
     where (F_&vVar^='Total' &  F_&hVar ^= 'Total');
close;
FhVar = F_&hVar;              /* levels of horiz var */
FvVar = F_&vVar;              /* lavels of vert var */
*print FhVar fVVar Frequency Percent ColPercent;

/* 4. Get the counts and percentages for each cell.
   Vertical centers are 
   v1/2,  v1 + v2/2, v1 + v2 + v3/2, ... */
vLabel = shape( FvVar, 0, nH )[,1];
vPos = shape( ColPercent, 0, nH );
nV = nrow(vPos);
halfW = vPos / 2;
vCenters = j(nrow(vPos), ncol(vPos), .);
do i = 1 to nH;
   vCenters[,i] = halfW[,i] + cusum(0 // vPos[1:(nV-1),i]);
end;
print vCenters[r=vLabel c=hLabel F=5.2];

/* 5. convert to a long format: (hPos, vPos, Freq, Pct) 
   and write to SAS data set */
hCenters = repeat(hCenter, nV);
CellFreq = shape( Frequency, 0, nH );
CellPct = shape( Percent, 0, nH );
result = colvec(hCenters) || colvec(vCenters) || 
         colvec(CellFreq) || colvec(CellPct);
*print FhVar FvVar result[c={"hCenter" "vCenter" "Freq" "Pct"}];

/* Optional: You might want to get rid of any labels that account to fewer 
   than 2% (or 1%) of the total. The criterion is up to you. For example: 
   idx = loc( result[,4] > 2 );  * keep if greater than 2% of total;
   result = result[idx, ];
   print result[c={"hCenter" "vCenter" "Freq" "Pct"}];
*/

/* write character and numeric vars separately, then merge together */
/* Character vars: pairwise labels of horiz and vertical categories */
create labels var {FhVar FvVar}; append; close;
/* Numeric vars: centers of cells, counts, and percentages */
create centers from result[c={"hCenter" "vCenter" "Freq" "Pct"}];
   append from result;
close;
QUIT;

data annoData;
merge labels centers;
run;


/* If we could use the WALLPERCENT drawing space, we could
   use (hCenter, vCenter) as the drawing coordinates:
   x1 = hCenter;
   y1 = vCenter;

   Unfortunately, we have to use LAYOUTPERCENT, so perform a
   linear transformation from "wall coordinates" to layout
   coordinates.

   For information about variable names, see
   https://blogs.sas.com/content/iml/2019/07/08/add-annotation-mosaic-plot-sas.html
*/
data anno;
set annoData;
length label $12;
retain function 'text' 
       y1space 'layoutpercent' x1space 'layoutpercent'
       width 4 
       anchor 'center';
/* Guess a linear transform to LAYOUTPERCENT coordinates.
   Need to move inside graph area, so shrink and 
   translate to correct for left and bottom axes areas */
x1 = 0.9*hCenter + 10;
y1 = 0.9*vCenter + 10;   
label = put(Freq, 5.);
run;

/* Note: PROC FREQ reverses Y axis b/c it sorts the 
   FreqOut data in descending order. To get the 
   reversed data set, call PROC SORT or use
      ods output MosaicPlot=MP;
   in the PROC FREQ call.
*/
proc template;
  define statgraph mosaicPlotParm;
  dynamic _VERTVAR _HORZVAR _FREQ _TITLE;
    begingraph;
      entrytitle _TITLE;
      layout region;    /* REGION layout, so can't overlay text! */
      MosaicPlotParm category=(_HORZVAR _VERTVAR) count=_FREQ / 
             datatransparency=0.5
             colorgroup=_VERTVAR name="mosaic";
      endlayout;
      annotate;         /* required for annotation */
    endgraph;
  end;
run;

proc sgrender data=FreqOut template=mosaicPlotParm sganno=anno;
dynamic _VERTVAR="Origin" _HORZVAR="Type" _FREQ="Count"
        _TITLE="Basic Mosaic Plot with Counts";
run;

