/* SAS program to accompany the article 
   "Compute the silhouette statistic in SAS"
   by Rick Wicklin, published 15MAY2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/05/15/silhouette-statistic-cluster.html

   This program shows how to compute the silhouette statistic in SAS.
*/


/* Implement the silhouette statistics for n points in k clusters.
   See p. 5 of
   Abbey, R. (2019) "How to Evaluate Different Clustering Results"
   https://support.sas.com/resources/papers/proceedings19/3409-2019.pdf
   The silhouette statistic was defined by 
   Peter Rousseeuw (1987) "Silhouettes: A graphical aid to the interpretation 
       and validation of cluster analysis", https://doi.org/10.1016/0377-0427(87)90125-7
*/
proc iml;

/* Define two helper functions */
/* AvgDistSelf: For each point, p, in a cluster this function returns the 
   average distance between p and other points in the same cluster.
*/
start AvgDistSelf(Y);
   n = nrow(Y);
   if n=1 then return(0);    /* convention */
   D = distance(Y);          /* D[i,j] is distance from Y[i,] to Y[j,] */
   AvgDist = D[ ,+] / (n-1); /* exclude D[i,i] */
   return AvgDist;
finish;
/* AvgDistOther: Let Y contains points in one cluster and Z contain
   points in a different cluster. For each point p in Y, this function 
   returns the average distance between p and the points in Z.
*/
start AvgDistOther(Y, Z);
   D = distance(Y,Z);        /* D[i,j] is distance from Y[i,] to Z[j,] */
   return( D[ ,:] );         /* average of the D[i,j] */
finish;

/* Silhouette: The main computational routine.
   For each point, p:
   1. Define AvgDistIn = average distance between p and other points in the same cluster.
      Note: The divisor for this cluster is N-1, not N.
   2. Find average distance between p and points in each of the other clusters.
      Define AvgDistOut to be the MINIMUM of the average distance.
   3. Define s(p) = (AvgDistOut(p) - AvgDistIn(p)) / max(AvgDistOut(p), AvgDistIn(p))

   Input: X         is an (n x d) matrix of observations 
          ClusterID is an (n x 1) vector that associates each point to a cluster. Typically,
                    the values are in the range {1,2,...,k}, where k is the total number of clusters.
   Output: S        is an (n x 1) vector that contains the silhouette statistic for each point in X.
*/
start Silhouette(X, ClusterID);
   n = nrow(X);
   AvgDistIn  = j(n, 1, .);
   AvgDistOut = j(n, 1, .);

   /* Use the Unique-Loc technique to iterate over clusterID:
   https://blogs.sas.com/content/iml/2011/11/01/the-unique-loc-trick-a-real-treat.html */
   uID = unique(ClusterID);
   k = ncol(uID);
   do i = 1 to k;
      idxIn = loc(ClusterID = uID[i]);
      Y = X[idxIn, ];            /* obs in the i_th cluster */
      v = AvgDistSelf(Y);        /* average inter-cluster distance */
      AvgDistIn[idxIn,] = v;
      dMin = j(ncol(idxIn), 1, constant('BIG'));   /* retain the min distances */
      h = setdif(1:k, i);        /* indices for the other clusters */
      do j = 1 to ncol(h);       /* loop over the other clusters */
         jdx = loc(ClusterID = uID[h[j]]); 
         Z = X[jdx, ];           /* obs in the j_th cluster, j ^= i */
         v = AvgDistOther(Y, Z); /* average between cluster distances */
         dMin = (dMin >< v);     /* retain the elementwise minimum */
      end;
      AvgDistOut[idxIn,] = dMin; /* min average distance to another cluster */
   end;
   s = (AvgDistOut - AvgDistIn) / (AvgDistOut <> AvgDistIn);
   return s;
finish;
store module=(AvgDistSelf AvgDistOther Silhouette);
QUIT;

/***************************************************/
/* These macros are not necessary, but they reduce 
   the amount of code in the blog post */
/***************************************************/

/* The SilhouetteMergeSort macro
   1. Merges the original data and the silhouette statistics
   2. Computes the overall average silhouette measure and stores it in &SILH macro variable
   3. Sorts data by cluster and (descending) by silhouette value
   4. Adds an _ObsNum variable for ease of plotting the silhouette values as a bar chart

   Args: DSname      = Original data
         SilhDSName  = Output from PROC IML, containing silhouette values
         ClusName    = Name of indicator variable for cluster membership
         SilhVarName = Name of variable that contains silhouette values
         OutName     = Name of merged data set
*/
%macro SilhouetteMergeSort(DSname, SilhDSName, ClusName, SilhVarName, OutName);
   %GLOBAL Silh;
   data _temp; 
      merge &DSName &SilhDSName end=EOF; 
      _mean + &SilhVarName;
      if EOF then
         call symputx( "Silh", round(_mean/_N_, 0.001));
   run;
   proc sort data=_temp; by &ClusName  descending &SilhVarName; run;
   data &OutName / view=&OutName; set _temp; _ObsNum=_N_; label ObsNum="Observation"; run;
%mend;

/* The SilhouetteScatter macro:
   To run this macro, you must first run the SilhouetteMergeSort macro!
   This macro creates a scatter plot of two variables in the cluster data
   and colors the markers by the silhouette statistic.
   The marker shapes indicate the cluster.

   Note: The macro hardcodes a spectral color ramp and uses only three shapes.
   Feel free to modify these attributes.

   Args: OutName     = Name of merged data set
         xVarName, yVarName = Names of X and Y variables for scatter plot
         ClusName    = Name of indicator variable for cluster membership
         SilhVarName = Name of variable that contains silhouette values
*/
%macro SilhouetteScatter(OutName, xVarName, yVarName, ClusName, SilhVarName);
   %let Spectral7 = (CX3288BD CX99D594 CXE6F598 CXFFFFBF CXFEE08B CXFC8D59 CXD53E4F );
   ods graphics / PUSH attrpriority=none;
   proc sgplot data=&OutName;
      /* add more symbols if you have more than 3 clusters */
      styleattrs datasymbols=(CircleFilled SquareFilled TriangleFilled);
      scatter x=&xVarname y=&yVarName / group=&ClusName name="scat"
                 colormodel=&Spectral7 colorresponse=&SilhVarName 
                 filledoutlinedmarkers markerattrs=(size=14);
      xaxis grid; yaxis grid;
      keylegend "scat" / type=Marker;  gradlegend "scat"; 
   run;
   ods graphics / POP;
%mend;


/* The SilhouettePlot macro:
   To run this macro, you must first run the SilhouetteMergeSort macro!
   This macro creates a silhouette plot, which is a panel of bar charts
   of the silhouette values.

   Args: OutName     = Name of merged data set
         ClusName    = Name of indicator variable for cluster membership
         SilhVarName = Name of variable that contains silhouette values
*/
%macro SilhouettePlot(OutName, ClusName, SilhVarName);
   proc sgpanel data=&OutName;
      panelby &ClusName / layout=rowlattice onepanel spacing=0 rowheaderpos=right uniscale=column;
      hbar _ObsNum / response=&SilhVarName group=&ClusName;
      rowaxis display=(noticks novalues);
      colaxis grid label="Silhouette Statistic";
   run;
%mend;

/* The SilhouetteHistogram macro:
   To run this macro, you must first run the SilhouetteMergeSort macro!
   This macro creates a silhouette plot as a panel of histograms.

   Args: OutName     = Name of merged data set
         ClusName    = Name of indicator variable for cluster membership
         SilhVarName = Name of variable that contains silhouette values
         SilhRef (optional) = Value of X variable at which to place a vertical reference line
                              This is often the overall silhouette value for the data set.
*/
%macro SilhouetteHistogram(OutName, ClusName, SilhVarName, SilhRef=);
   proc sgpanel data=&OutName;
      panelby &ClusName / layout=rowlattice onepanel spacing=0 rowheaderpos=right uniscale=column;
      histogram &SilhVarName;
      %if %length(&SilhRef) %then %do;
          refline &SilhRef / axis=x;
      %end;
      rowaxis grid; 
   run;
%mend;

/* SAS program to accompany the article 
   "What is the silhouette statistic in cluster analysis?"
   by Rick Wicklin, published 15MAY2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/05/15/silhouette-statistic-cluster.html

   This program demonstrates the geometry and visualization of the silhouette statistic.
*/

/*********************************************************/
/* 1. Visualize the geometry of the silhouette statistic */
/*********************************************************/
data Have;
input ID $ x y;
datalines;
A  0  0
A  0  1
A  1  0
A -1  0
A  0 -1
B 12 -1
B  6  0
B 12  1
C  1 11
C -1 11
C -1  9
C  1  9
;


/* run some preliminary tests */
/*
proc iml;
p = {0 0};
A = {0  1, 1  0, -1  0, 0 -1};
B = {12 -1, 6  0, 12  1};
C = {1 11, -1 11, -1  9, 1  9};

DA = distance(p, A);
meanDA = DA[,:];
print meanDA, DA;
DB = distance(p, B);
meanDB = DB[,:];
print meanDB, DB;
DC = distance(p, C);
meanDC = DC[,:];
print meanDC, DC;
QUIT;
*/

/* k-means cluster analysis */
proc fastclus data=Have maxclusters=3 out=ClusOut noprint;
  var x y;
run;

ods graphics / reset attrpriority=color;
ods graphics / push width=480px height=480px;

/* Visualize data and clusters */
title "Example Data in Three Clusters";
proc sgplot data=ClusOut aspect=1;
  scatter x=x y=y / group=Cluster markerattrs=(symbol=CircleFilled size=14);
  xaxis grid values=(0 to 12 by 2) valueshint;
  yaxis grid values=(0 to 12 by 2) valueshint;
run; 

/* show vectors from 1st point to others in 1st cluster */
data Within1;
set ClusOut;
xf = .; yf = .;
if Cluster=1 then do;
   xf = x; yf = y;
end;
run;

title "Within-Cluster Distances for One Point in Cluster 1";
proc sgplot data=Within1 aspect=1;
  scatter x=x y=y / group=Cluster markerattrs=(symbol=CircleFilled size=14);
  vector x=xf y=yf / group=Cluster noarrowheads lineattrs=(thickness=3) transparency=0.5;
  xaxis grid values=(0 to 10 by 2) valueshint;
  yaxis grid values=(0 to 10 by 2) valueshint;
run; 

/* show vectors from 1st point to points in other clusters */
data Between1;
set ClusOut;
xf = .; yf = .;
if Cluster^=1 then do;
   xf = x; yf = y;
end;
run;

title "Within-Cluster Distances for One Point in Cluster 1";
proc sgplot data=Between1 aspect=1;
  scatter x=x y=y / group=Cluster markerattrs=(symbol=CircleFilled size=14);
  vector x=xf y=yf / group=Cluster noarrowheads lineattrs=(thickness=3) transparency=0.5;
  xaxis grid values=(0 to 12 by 2) valueshint;
  yaxis grid values=(0 to 12 by 2) valueshint;
run; 


/***************************************************************/
/* 2. Create scatter and silhouette plots for the example data */
/***************************************************************/
data Have;
input ID $ x y;
datalines;
A  0  0
A  0  1
A  1  0
A -1  0
A  0 -1
B 12 -1
B  6  0
B 12  1
C  1 11
C -1 11
C -1  9
C  1  9
;

/* k-means cluster analysis */
proc fastclus data=Have maxclusters=3 out=ClusOut noprint;
  var x y;
run;

/* compute vector of silhouette statistics */
proc iml;
load module=(Silhouette);
varNames = {'x' 'y'};
use ClusOut;
   read all var varNames into X;
   read all var "Cluster" into ClusterID;
close;
silhouette = Silhouette(X, ClusterID);

print ClusterID X[c=varNames L=""] silhouette;
m = mean(silhouette);
print m[L="Silhouette Measure"];

/* write silhouette values to data set for graphing */
create S var "Silhouette"; append; close;
QUIT;

/* merge data and create macro variable */
%SilhouetteMergeSort(ClusOut, S, Cluster, Silhouette, _All);
%put &=Silh;

/* scatter plot; markers colored by silhouette values */
title "Clusters and Silhouette Values";
title2 "Average Silhouette Value = &Silh";
%SilhouetteScatter(_All, x, y, Cluster, Silhouette);

/* for this example, LABEL the markers with silhouette values */
%let Spectral7 = (CX3288BD CX99D594 CXE6F598 CXFFFFBF CXFEE08B CXFC8D59 CXD53E4F);
ods graphics / PUSH attrpriority=none;
proc sgplot data=_All;
   /* add more symbols if you have more than 3 clusters */
   format Silhouette 4.2;
   styleattrs datasymbols=(CircleFilled SquareFilled TriangleFilled);
   scatter x=x y=y / group=Cluster name="scat" DATALABEL=Silhouette 
              colormodel=&Spectral7 colorresponse=Silhouette 
              filledoutlinedmarkers markerattrs=(size=14);
   xaxis grid; yaxis grid;
   keylegend "scat" / type=Marker;  gradlegend "scat"; 
run;
ods graphics / POP;

/* silhouette plot (bar chart) */
ods graphics /push width=480px height=480px;
title "Silhouette Plot";
title2 "Average Silhouette Value = &Silh";
%SilhouettePlot(_All, Cluster, Silhouette);
ods graphics /pop;

/* summary statistics by cluster for the silhouette values */
proc means data=_All Min Q1 Median Q3 Max ndec=3;
   class Cluster;
   var Silhouette;
run;


/***************************************************************/
/* 2. Create scatter and silhouette plots for the iris data    */
/***************************************************************/
/* k-means cluster analysis */
proc fastclus data=Sashelp.Iris maxclusters=3 out=ClusOut noprint;
  var SepalLength SepalWidth PetalLength PetalWidth;
run;

/* compute vector of silhouette statistics */
proc iml;
load module=(Silhouette);
varNames = {'SepalLength' 'SepalWidth' 'PetalLength' 'PetalWidth'};
use ClusOut;
   read all var varNames into X;
   read all var "Cluster" into ClusterID;
close;
silhouette = Silhouette(X, ClusterID);
create S var "Silhouette"; append; close; /* write silhouette values to data set */
QUIT;

/* merge data and create macro variable */
%SilhouetteMergeSort(ClusOut, S, Cluster, Silhouette, _All);
%put &=Silh;

/* scatter plot; markers colored by silhouette values */
title "Clusters and Silhouette Values";
title2 "Average Silhouette Value = &Silh";
%SilhouetteScatter(_All, SepalLength, PetalLength, Cluster, Silhouette);

/* 4-D data, so use a scatter plot matrix to visualize other projections */
ods graphics / PUSH attrpriority=none width=720px height=640px;
footnote J=L "Marker shapes indicate clusters";
proc sgscatter data=_All datasymbols=(CircleFilled SquareFilled TriangleFilled);
compare x=(SepalLength SepalWidth)
        y=(PetalLength PetalWidth) / grid group=Cluster
        colormodel=&Spectral7 colorresponse=Silhouette filledoutlinedmarkers
        markerattrs=(size=14);
run;
footnote;
ods graphics / POP;

/* silhouette plot (bar chart) */
title "Silhouette Plot";
title2 "Average Silhouette Value = &Silh";
%SilhouettePlot(_All, Cluster, Silhouette);

/* summary statistics by cluster for the silhouette values */
proc means data=_All Min Q1 Median Q3 Max;
   class Cluster;
   var Silhouette;
run;

/* silhouette plot (histogram) */
ods graphics / PUSH width=480px height=480px;
title "Distribution of Silhouette Values by Cluster";
%SilhouetteHistogram(_All, Cluster, Silhouette, SilhRef=&Silh);
ods graphics / pop;


/*********************************************/
/* Additional Examples */
/*********************************************/


/*-----------------------------------------------------------------*/
/* Example: Assign cluster ID alternating on a circle */
/*-----------------------------------------------------------------*/
data AltCircle;
pi = constant('pi');
n = 8;
do i = 0 to n-1;
   theta = i * 2*pi/n;
   x = cos(theta);
   y = sin(theta);
   /* alternate assignment to cluster 0 and cluster 1 */
   Cluster = mod(i,2);
   output;
end;
drop pi n i;
run;


/*-----------------------------------------------------------------*/
/* Example: reflect point through origin and change the cluster ID */
/*-----------------------------------------------------------------*/
data ReflectOut;
call streaminit(1234);
do i = 1 to 100;
   x = rand("Normal", 0, 1);
   y = rand("Normal", 0, 1);
   Cluster = 1;
   output;
   x = -x; y = -y; Cluster = 2;
   output;
end;
drop i;
run;


/*-----------------------------------------------------------------*/
/* Example: Assign cluster ID randomly                             */
/* We expect an overall silhouette measure close to 0              */
/*-----------------------------------------------------------------*/
data ClusOut;
call streaminit(1234);
do i = 1 to 100;
   x = rand("Normal", 0, 1);
   y = rand("Normal", 0, 1);
   /* random assignment to either cluster 0 or cluster 1 */
   Cluster = rand("Bernoulli", 0.5);
   output;
end;
drop i;
run;

proc means data=ClusOut N Mean;
   class Cluster;
   var x y;
run;

proc iml;
load module=(Silhouette);
varNames = {'x' 'y'};
use ClusOut;
   read all var varNames into X;
   read all var "Cluster" into ClusterID;
close;
silhouette = Silhouette(X, ClusterID);
create S var "Silhouette"; append; close;
QUIT;


%SilhouetteMergeSort(ClusOut, S, Cluster, Silhouette, _All);
%put &=Silh;

title "Silhouette Values for Random Assignment to Clusters";
title2 "Average Silhouette Value = &Silh";
%SilhouetteScatter(_All, x, y, Cluster, Silhouette);

title "Silhouette Plot";
title2 "Average Silhouette Value = &Silh";
%SilhouettePlot(_All, Cluster, Silhouette);

title "Distribution of Silhouette Values by Cluster";
%SilhouetteHistogram(_All, Cluster, Silhouette, SilhRef=&Silh);
