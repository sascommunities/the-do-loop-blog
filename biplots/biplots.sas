/* SAS program to accompany the article
   "How to create biplots in SAS"
   by Rick Wicklin, published 13NOV2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/11/30/create-biplots-in-sas.html

   This program shows how to create biplots in SAS by using modern ODS
   statistical graphics (PROC SGPLOT). Biplots are discussed 
   in the article https://blogs.sas.com/content/iml/2019/11/06/what-are-biplots.html

   This program is in two parts:
   1. Use Michael Friendly's %BIPLOT macro. The OUT= option writes the 
      coordinates for the markers and vectors to a data set. 
      Transpose the data from long to wide form and 
      use PROC SGPLOT to create the biplot.
   2. Use SAS/IML modules directly. Plot by using PROC SGPLOT.
      
   You can download the %BIPLOT macro from 
         http://www.datavis.ca/books/vcd/biplot.html
   You can download the %EQUATE macro from 
         http://www.datavis.ca/books/vcd/equate.html
*/
title;
data iris;
set Sashelp.iris;
id = " ";        /* create empty label variable */
/* or use obs numbers:*/
*id = put(_N_, 3.);     
run;

/***************************************************/
proc prinqual data=iris plots=(MDPref) 
              n=2       /* project onto Prin1 and Prin2 */
              mdpref=1; /* use COV scaling */
   transform identity(SepalLength SepalWidth PetalLength PetalWidth);  /* identity transform */
   id ID;
   ods select MDPrefPlot;
run;

/***************************************************/
/* A. Use %BIPLOT macro, which uses SAS/IML to compute the biplot coordinates. 
      Use the OUT= option to get the coordinates for the markers and vectors.
   B. Transpose the data from long to wide form.
   C. Use PROC SGPLOT to create the biplot
*/
%let FACTYPE = SYM;   /* options are GH, COV, JK, SYM */
title "Biplot: &FACTYPE, STD";
%biplot(data=iris, 
        var=SepalLength SepalWidth PetalLength PetalWidth, 
        id=id,
        factype=&FACTYPE,  /* GH, COV, JK, SYM */
        std=std,           /* NONE, MEAN, STD  */
        scale=1,           /* if you do not specify SCALE=1, vectors are auto-scaled */
        out=biplotFriendly,/* write SAS data set with results */ 
        symbols=circle dot, inc=1);
/* looks like the p1 and p2 macro variables are not availble */
/* %put &=_user_; */

/* transpose from long to wide */
data Biplot;
set biplotFriendly(where=(_TYPE_='OBS') rename=(dim1=Prin1 dim2=Prin2 _Name_=_ID_))
    biplotFriendly(where=(_TYPE_='VAR') rename=(dim1=vx dim2=vy _Name_=_Variable_));
run;

proc sgplot data=Biplot aspect=1 noautolegend;
   refline 0 / axis=x; refline 0 / axis=y;
   scatter x=Prin1 y=Prin2 / datalabel=_ID_;
   vector  x=vx    y=vy    / datalabel=_Variable_
                             lineattrs=GraphData2 datalabelattrs=GraphData2;
   xaxis grid offsetmin=0.1 offsetmax=0.2;
   yaxis grid;
run;   

/***************************************************/

/* NOTE: Before you can use the SAS/IML modeuls, first run the 
   program that defines the modeuls.
   Download the module definitions from 
   https://github.com/sascommunities/the-do-loop-blog/blob/master/biplots/biplot-iml-modules.sas
*/

ods graphics / width=480px height=480px;

proc iml;
/* assumes the modules have been previously stored */
load module=(CalcBiplot WriteBiplot Biplot);
use sashelp.iris;
read all var _NUM_ into X[rowname=Species colname=varNames];
close;

/* Generate all four biplots */

/* GH biplot: the principal component variables have unit variances. 
   The vectors for the original variables might be relatively large.
   The COV Biplot is a rescaled version of the GH biplot, 
   which is often more useful.
*/
title "GH Biplot: Relationships between Variables";
run Biplot(X, Species, varNames, "GH", "Std");
title "COV Biplot: Rescaled GH Biplot";
run Biplot(X, Species, varNames, "COV", "Std");

/* JK biplot: the variances of the PC variables are equal to the 
   corresponding eigenvalues. */
title "JK Biplot: Relationships between Observations";
run Biplot(X, Species, varNames, "JK", "Std");

/* SYM biplot: displays approximate relationships between observations
   and variables. */
title "SYM Biplot: Relationships between Observations and Variables";
run Biplot(X, Species, varNames, "SYM", "Std");


/* Demonstate using StdMethod = "Mean" and "None" */
title "GH Biplot: No Scaling";
run Biplot(X, Species, varNames, "GH", "Mean");
title "GH Biplot: Automatic Scaling";
run Biplot(X, Species, varNames, "GH", "None") scale=0;

/* label obs by Species */
title "COV Biplot with Labels";
run Biplot(X, Species, varNames, "COV", "Std") labelPoints=1;   
/* Label obs by observation number */
title "JK Biplot: Relationships between Observations";
run Biplot(X, NULL, varNames, "JK", "Std") labelPoints=1;

/* autoscale; empty ID var */
title "JK Biplot: Automatic Scaling of Vectors";
run Biplot(X, NULL, varNames, "JK", "Std") scale=0;            

/* scale vectors by 0.25 */
title "SYM Biplot: Vectors Scaled by 0.25";
run Biplot(X, NULL, varNames, "SYM", "Std") scale=0.25;

/* Optional: Iterate over all biplot types */
/*
FacType = {JK SYM COV};
StdMethod = "Std";
title "Biplot";
do i = 1 to ncol(FacType);
   run Biplot(X, Name, varNames, FacType[i], StdMethod);
end;
*/
QUIT;

/************************************************/

/* Show how to use the WriteBiplot module to create a data set
   and then use PROC SGPLOT to create the biplot. 

   You can easily create a macro from this set of statements. */
%let FACTYPE=JK;
%let STD=Std;
title "&FACTYPE Biplot: Automatic Scaling of Vectors";
title2 "FacType=&FACTYPE; Std=&STD";
proc iml;
/* assumes the modules have been previously stored */
load module=(CalcBiplot WriteBiplot);
use sashelp.iris;
read all var _NUM_ into X[rowname=Species colname=varNames];
close;

run WriteBiplot(X, NULL, varNames, "&FACTYPE", "&STD") scale=0;
QUIT;

data Biplot;
   set _Scores _Vectors;
run;

proc sgplot data=Biplot aspect=1 noautolegend;
   refline 0 / axis=x; refline 0 / axis=y;
   scatter x=Prin1 y=Prin2 / ; 
   vector  x=vx    y=vy    / datalabel=_Variable_
                             lineattrs=GraphData2 datalabelattrs=GraphData2;
   xaxis grid offsetmin=0.1 offsetmax=0.1 min=&minAxis max=&maxAxis;
   yaxis grid min=&minAxis max=&maxAxis;
run;

/***********************************/

/* Additional example in which the biplot data is merged with the 
   original data so that you can color-code the observations by 
   using the Species variable.
*/
data _Obs;
   merge iris _Scores;
run;
data Biplot;
   set _Obs _Vectors;
run;
proc sgplot data=Biplot aspect=1 noautolegend;
   refline 0 / axis=x; refline 0 / axis=y;
   scatter x=Prin1 y=Prin2 / GROUP=Species; 
   vector  x=vx    y=vy    / datalabel=_Variable_
                             lineattrs=(color=black) datalabelattrs=(color=black);
   xaxis grid offsetmin=0.1 offsetmax=0.1 min=&minAxis max=&maxAxis;
   yaxis grid min=&minAxis max=&maxAxis;
run;

/***********************************/
