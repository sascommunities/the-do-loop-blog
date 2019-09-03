/* SAS program to accompany the article 
   "Use cosine similarity to make recommendations"
   by Rick Wicklin, published 05SEP2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/09/05/cosine-similarity-recommendations.html

   This program shows computes the cosine similarity of rows for six recipes that 
   have a total of 35 ingredients among them. The cosine similarity tells you 
   which recipes are similar to each other. You can use a heat map to visualize
   the cosine similarity matrix.

   The cosine similarity is one tool that is available for recommender engines.
   If a user likes one item, the user might be interested in a similar item.
*/
data recipes;
   input Recipe $ 1-20
      (Tomato Garlic Salt Onion TomatoPaste OliveOil Celery Broth 
       GreenPepper Cumin Flour BrownSugar BayLeaf GroundBeef 
       BlackPepper ChiliPowder Cilantro Carrot CayennePepper Oregano 
       Oil Parsley PorkSausage RedPepper Paprika Thyme Tomatillo 
       JalapenoPepper WorcestershireSauce Lime
       Eggplant GreenOlives Capers Sugar
       ) (1.);
datalines;
Spag Sauce          1111110000000000000101000000000000
Spag Meat Sauce     1111111010001100010000110000000000
Eggplant Relish     0111110000000000000000000000001111
Creole Sauce        1011111110000010000000001100100000
Salsa               1111000000000000100000000011010000
Enchilada Sauce     1000000101110001001010000000000000
;


/* you can compute the cosine similarity by using SAS/IML */
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
QUIT;


ods graphics/reset;
proc iml;
load module=(CosSimRows CosSimCols);
use Recipes; read all var _NUM_ into X[c=varNames r=Recipe]; close;

cosSim = CosSimRows(X);
Colors = palette('BRBG', 7);
call heatmapcont(cosSim) xvalues=Recipe yvalues=Recipe range={0 1} 
            colorramp=Colors title="Cosine Similarity of Recipes";
print CosSim[F=4.2 r=Recipe c=Recipe];

/* compare to the correlation matrix */
/*
corr = Corr(X`);
call heatmapcont(corr) xvalues=Recipe yvalues=Recipe 
            colorramp=Colors title="Correlation of Recipes";
print Corr[F=5.2 r=Recipe c=Recipe];

cosC = CosSimCols(X);
call heatmapcont(cosC)
            colorramp=Colors title="Cosine Similarity of Ingredients";
*/

/* pairwise similarities in a bar chart 
   https://blogs.sas.com/content/iml/2017/08/16/pairwise-correlations-bar-chart.html
*/
ColNames = Recipe;
numCols = ncol(CosSim);                /* number of variables */
numPairs = numCols*(numCols-1) / 2;
length = 2*nleng(ColNames) + 5;       /* max length of new ID variable */
PairNames = j(NumPairs, 1, BlankStr(length));
i = 1;
do row= 2 to numCols;                 /* construct the pairwise names */
   do col = 1 to row-1;
      PairNames[i] = strip(ColNames[col]) + " vs. " + strip(ColNames[row]);
      i = i + 1;
   end;
end;
 
lowerIdx = loc(row(CosSim) > col(CosSim));  /* indices of lower-triangular elements */
CosSim = CosSim[ lowerIdx ];
create CosPairs var {"PairNames" "CosSim"};  append;  close;
QUIT;

proc sort data=CosPairs;  by CosSim;  run;
 
%macro HalfWidth(nCat);
   %sysevalf(0.5/&nCat)
%mend;

ods graphics / width=500px height=400px;
title "Pairwise Cosine Similarities";
proc sgplot data=CosPairs;
format CosSim 4.2;
hbar PairNames / response=CosSim datalabel;
refline 0 / axis=x;
yaxis discreteorder=data display=(nolabel) 
      labelattrs=(size=6pt) fitpolicy=none 
      offsetmin=%HalfWidth(15) offsetmax=%HalfWidth(15) /* half of 1/k, where k=number of catgories */
      colorbands=even colorbandsattrs=(color=gray transparency=0.9);
xaxis grid display=(nolabel);
run;

/*************************************/
