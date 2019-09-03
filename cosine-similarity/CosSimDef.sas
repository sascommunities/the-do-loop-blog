/* SAS program to accompany the article 
   "Cosine similarity of vectors"
   by Rick Wicklin, published 03SEP2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/09/03/cosine-similarity.html

   This program shows how to compute the cosine similarity in SAS.
   You can 
   1. Use PROC DISTANCE to compute the cosine similarities of rows
   2. Use PROC IML to compute the cosine similarity of rows or columns.

   PROC IML also provides an easy way to create a heat map that visualizes
   the cosine similarity matrix.
*/

/********************************************/
/* Very simple 2-D example                  */
/********************************************/
data Vectors;
length Name $1;
input Name x y;
datalines;
A  0.5 1
B  3   5
C  3   2.8
D  5   1
;

/* plot the vectors */
ods graphics / width=400px height=400px;
title "Four Row Vectors";
proc sgplot data=Vectors aspect=1;
   vector x=x y=y / datalabel=Name datalabelattrs=(size=14);
   xaxis grid;
   yaxis grid;
run;

/*
proc distance data=Vectors out=Dist method=EUCLID shape=square;
   var ratio(_NUMERIC_);
   id Name;
run;
proc print data=Dist noobs; run;
*/
proc distance data=Vectors out=Cos method=COSINE shape=square;
   var ratio(_NUMERIC_);
   id Name;
run;
proc print data=Cos noobs; run;

/* you can also compute the cosine similarity by using SAS/IML */
proc iml;
/* complete cases:
   https://blogs.sas.com/content/iml/2015/02/23/complete-cases.html
   exclude any row with a missing value */
start ExtractCompleteCases(X);
   if all(X ^= .) then return X;
   idx = loc(countmiss(X, "row")=0);
   if ncol(idx)>0 then return( X[idx, ] );
   else                return( {} ); 
finish;

/* cosine similarity of variables */
start CosSimCols(X, checkForMissing=1);
   if checkForMissing then do;    /* by default, check for missing and exclude */
      Z = ExtractCompleteCases(X);
      Y = Z / sqrt(Z[##,]);       /* stdize each column */
   end;
      else Y = X / sqrt(X[##,]);  /* skip check if you know all values are valid */
   cosY = Y` * Y;                 /* pairwise inner products */
   /* because of finite precision, elements could be 1+eps or -1-eps */
   idx = loc(cosY> 1); if ncol(idx)>0 then cosY[idx]= 1; 
   idx = loc(cosY<-1); if ncol(idx)>0 then cosY[idx]=-1; 
   return cosY;
finish;

/* cosine similarity of observations */
start CosSimRows(X);
   Z = ExtractCompleteCases(X);  /* check for missing and exclude */
   return T(CosSimCols(Z`, 0));  /* transpose and call CosSimCols */
finish;
store module=(ExtractCompleteCases CosSimCols CosSimRows);

/* test the computations on the simple data */
use Vectors; read all var _NUM_ into X[r=Name]; close;
cosY = CosSimRows(X);
print cosY[r=Name c=Name format=7.5];

/* If you want the actual angles, apply the inverse cosine function */
/*
deg = arcos(cosY)*180/constant('pi');
reset fuzz; print deg[f=5.2];
*/

/* plot the standarized vectors */
Y = X / sqrt(X[,##]);   /* stdize each row */
create StdVectors from Y[r=Name c={'x' 'y'}];
append from Y[r=Name];
close;
QUIT;

title "Standardized Vectors";
proc sgplot data=StdVectors aspect=1;
   vector x=x y=y / datalabel=Name datalabelattrs=(size=14);
   xaxis grid min=0 max=1;
   yaxis grid min=0 max=1;
run;

/***************************************************/
/* second example, vehicles with very large/small horsepower */
data Vehicles;
   set sashelp.cars(where=(Origin='USA'));
   if Horsepower < 140 OR Horsepower >= 310;
run;
proc sort data=Vehicles; by Type; run;
proc means data=Vehicles;run;

ods graphics / reset;
proc iml;
load module=(CosSimCols CosSimRows);
use Vehicles;
read all var _NUM_ into X[c=varNames r=Model]; 
read all var {Model Type}; close;

labl = compress(Type + ":" + substr(Model, 1, 10));
cosRY = CosSimRows(X);
*call heatmapcont(cosRY) xvalues=labl yvalues=labl
          title="Cosine Similarity between Vehicles";

cosCY = CosSimCols(X);
call heatmapcont(cosCY) xvalues=varNames yvalues=varNames
          title="Cosine of Angle Between Variables";

/* If desired, you can compare cosine similarity for correlation.
   For these data, the cosine similarity is often in the range [0.8, 1],
   The real advantage of cosine similarity is when the data have a lot
   of zeros, because then the inner product between two vectors 
   only contains the product of the mutually nonzero elements.
*/
corr=corr(X);
call heatmapcont(corr) xvalues=varNames yvalues=varNames
          title="Correlation Between Variables";

QUIT;

/****************************/
