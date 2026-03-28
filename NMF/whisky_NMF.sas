
/* SAS program to accompany the article 
   "An NMF analysis of Scotch whiskies"
   by Rick Wicklin, published 30MAR2026 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2026/03/30/nmf-whisky.html

   This program performs an NMF factoriztions on whisky data.
   This data analysis was first presented in 
   Young, Fogel, and Hawkins (2006),
   "Clustering scotch whiskies using non-negative matrix factorization",
   SPES Newsletter.
   https://www.niss.org/sites/default/files/ScotchWhisky.pdf
*/

/********************************/

/* Perform a principal component analysis (PCA) of the 
Scotch Whisky data set in Young, Fogel, and Hawkins 
(SPES, 2006), who analyzed the data by using a nonnegative 
matrix factorization (NMF). The article is available at 
https://www.niss.org/sites/default/files/ScotchWhisky.pdf

I downloaded the data from  
https://www.niss.org/sites/default/files/ScotchWhisky01.txt
and corrected a few typos:
- ID=12: Replaced 'Belvenie' with 'Balvenie' 
- ID=56: Replaced 'Laphroig' with 'Laphroaig' 
- ID-85: Replaced 'Tomore' with 'Tormore'
*/
data Whisky_Orig;
infile datalines dsd dlm=','; 
length Distillery $20;
input
RowID Distillery Body Sweetness Smoky Medicinal Tobacco Honey Spicy Winey Nutty Malty Fruity Floral;
datalines;
01,Aberfeldy,2,2,2,0,0,2,1,2,2,2,2,2
02,Aberlour,3,3,1,0,0,4,3,2,2,3,3,2
03,AnCnoc,1,3,2,0,0,2,0,0,2,2,3,2
04,Ardbeg,4,1,4,4,0,0,2,0,1,2,1,0
05,Ardmore,2,2,2,0,0,1,1,1,2,3,1,1
06,ArranIsleOf,2,3,1,1,0,1,1,1,0,1,1,2
07,Auchentoshan,0,2,0,0,0,1,1,0,2,2,3,3
08,Auchroisk,2,3,1,0,0,2,1,2,2,2,2,1
09,Aultmore,2,2,1,0,0,1,0,0,2,2,2,2
10,Balblair,2,3,2,1,0,0,2,0,2,1,2,1
11,Balmenach,4,3,2,0,0,2,1,3,3,0,1,2
12,Balvenie,3,2,1,0,0,3,2,1,0,2,2,2
13,BenNevis,4,2,2,0,0,2,2,0,2,2,2,2
14,Benriach,2,2,1,0,0,2,2,0,0,2,3,2
15,Benrinnes,3,2,2,0,0,3,1,1,2,3,2,2
16,Benromach,2,2,2,0,0,2,2,1,2,2,2,2
17,Bladnoch,1,2,1,0,0,0,1,1,0,2,2,3
18,BlairAthol,2,2,2,0,0,1,2,2,2,2,2,2
19,Bowmore,2,2,3,1,0,2,2,1,1,1,1,2
20,Bruichladdich,1,1,2,2,0,2,2,1,2,2,2,2
21,Bunnahabhain,1,2,1,1,0,1,1,1,1,2,2,3
22,Caol Ila,3,1,4,2,1,0,2,0,2,1,1,1
23,Cardhu,1,3,1,0,0,1,1,0,2,2,2,2
24,Clynelish,3,2,3,3,1,0,2,0,1,1,2,0
25,Craigallechie,2,2,2,0,1,2,2,1,2,2,1,4
26,Craigganmore,2,3,2,1,0,0,1,0,2,2,2,2
27,Dailuaine,4,2,2,0,0,1,2,2,2,2,2,1
28,Dalmore,3,2,2,1,0,1,2,2,1,2,3,1
29,Dalwhinnie,2,2,2,0,0,2,1,0,1,2,2,2
30,Deanston,2,2,1,0,0,2,1,1,1,3,2,1
31,Dufftown,2,3,1,1,0,0,0,0,1,2,2,2
32,Edradour,2,3,1,0,0,2,1,1,4,2,2,2
33,GlenDeveronMacduff,2,3,1,1,1,1,1,2,0,2,0,1
34,GlenElgin,2,3,1,0,0,2,1,1,1,1,2,3
35,GlenGarioch,2,1,3,0,0,0,3,1,0,2,2,2
36,GlenGrant,1,2,0,0,0,1,0,1,2,1,2,1
37,GlenKeith,2,3,1,0,0,1,2,1,2,1,2,1
38,GlenMoray,1,2,1,0,0,1,2,1,2,2,2,4
39,GlenOrd,3,2,1,0,0,1,2,1,1,2,2,2
40,GlenScotia,2,2,2,2,0,1,0,1,2,2,1,1
41,GlenSpey,1,3,1,0,0,0,1,1,1,2,0,2
42,Glenallachie,1,3,1,0,0,1,1,0,1,2,2,2
43,Glendronach,4,2,2,0,0,2,1,4,2,2,2,0
44,Glendullan,3,2,1,0,0,2,1,2,1,2,3,2
45,Glenfarclas,2,4,1,0,0,1,2,3,2,3,2,2
46,Glenfiddich,1,3,1,0,0,0,0,0,0,2,2,2
47,Glengoyne,1,2,0,0,0,1,1,1,2,2,3,2
48,Glenkinchie,1,2,1,0,0,1,2,0,0,2,2,2
49,Glenlivet,2,3,1,0,0,2,2,2,1,2,2,3
50,Glenlossie,1,2,1,0,0,1,2,0,1,2,2,2
51,Glenmorangie,2,2,1,1,0,1,2,0,2,1,2,2
52,Glenrothes,2,3,1,0,0,1,1,2,1,2,2,0
53,Glenturret,2,3,1,0,0,2,2,2,2,2,1,2
54,Highland Park,2,2,3,1,0,2,1,1,1,2,1,1
55,Inchgower,1,3,1,1,0,2,2,0,1,2,1,2
56,Isle of Jura,2,1,2,2,0,1,1,0,2,1,1,1
57,Knochando,2,3,1,0,0,2,2,1,2,1,2,2
58,Lagavulin,4,1,4,4,1,0,1,2,1,1,1,0
59,Laphroaig,4,2,4,4,1,0,0,1,1,1,0,0   
60,Linkwood,2,3,1,0,0,1,1,2,0,1,3,2
61,Loch Lomond,1,1,1,1,0,1,1,0,1,2,1,2
62,Longmorn,3,2,1,0,0,1,1,1,3,3,2,3
63,Macallan,4,3,1,0,0,2,1,4,2,2,3,1
64,Mannochmore,2,1,1,0,0,1,1,1,2,1,2,2
65,Miltonduff,2,4,1,0,0,1,0,0,2,1,1,2
66,Mortlach,3,2,2,0,0,2,3,3,2,1,2,2
67,Oban,2,2,2,2,0,0,2,0,2,2,2,0
68,OldFettercairn,1,2,2,0,1,2,2,1,2,3,1,1
69,OldPulteney,2,1,2,2,1,0,1,1,2,2,2,2
70,RoyalBrackla,2,3,2,1,1,1,2,1,0,2,3,2
71,RoyalLochnagar,3,2,2,0,0,2,2,2,2,2,3,1
72,Scapa,2,2,1,1,0,2,1,1,2,2,2,2
73,Speyburn,2,4,1,0,0,2,1,0,0,2,1,2
74,Speyside,2,2,1,0,0,1,0,1,2,2,2,2
75,Springbank,2,2,2,2,0,2,2,1,2,1,0,1
76,Strathisla,2,2,1,0,0,2,2,2,3,3,3,2
77,Strathmill,2,3,1,0,0,0,2,0,2,1,3,2
78,Talisker,4,2,3,3,0,1,3,0,1,2,2,0
79,Tamdhu,1,2,1,0,0,2,0,1,1,2,2,2
80,Tamnavulin,1,3,2,0,0,0,2,0,2,1,2,3
81,Teaninich,2,2,2,1,0,0,2,0,0,0,2,2
82,Tobermory,1,1,1,0,0,1,0,0,1,2,2,2
83,Tomatin,2,3,2,0,0,2,2,1,1,2,0,1
84,Tomintoul,0,3,1,0,0,2,2,1,1,2,1,2
85,Tormore,2,2,1,0,0,1,0,1,2,1,0,0
86,Tullibardine,2,3,0,0,1,0,2,1,1,2,2,1
;

/* For later analysis, identify some of the most popular single-malt 
   Scotch whiskies (eg, Glenlivet, Glenfiddich, Macallan,...) and 
   some that have extreme scores in a PCA of the data.
*/
data Whisky / view=Whisky;
length ID $20;
set Whisky_Orig;
selected = 0;
BestSeller = 0;
if Distillery in (
   'Glenlivet' 'Glenfiddich' 'Macallan' 'Glenmorangie' 'Balvenie' 
   'Laphroaig' 'Aberlour' 'Lagavulin' 'Ardbeg' 'Talisker' 
   ) then BestSeller=1;
/* also select a few whiskies that have extreme values in a PC */
if BestSeller | RowID in (07 11 32 33 35 43 60 62 82 85) then 
   selected = 1;
if selected then 
   ID = Distillery;
else 
   ID = put(RowID, Z2.);
run;


/*
Here are 10 of the most popular single-malt Scotch whiskies:
The Glenlivet (12 Year Old): Often cited as the best-selling single malt in the US, known for a smooth, fruity style.
Glenfiddich (12 Year Old): A global leader in volume, famous for its accessible, pear-forward, and green apple flavor profile.
The Macallan (12 or 18 Year Old): Highly coveted for its rich Sherry Oak maturation, considered a premium staple.
Glenmorangie (The Original 10 Year Old): A popular Highland malt known for its light, delicate, and citrusy character.
Balvenie [MISPELLED?] (DoubleWood 12 Year Old): Known for aging in two different wood types, offering honeyed, nutty, and vanilla notes.
Laphroaig (10 Year Old): The most popular heavily-peated Islay malt, recognized for its medicinal, smoky, and maritime character.
Aberlour (12 or 16 Year Old): Highly regarded for its Sherried Speyside style, with rich dried fruit and spice notes.
Lagavulin (16 Year Old): A cult-favorite, intensely smoky, and complex Islay malt.
Ardbeg (10 Year Old): Another major Islay player, celebrated for a smoky, peppery, and balanced flavor.
Talisker (10 Year Old): Known for its maritime, peppery, and coastal character from the Isle of Skye.
*/
%let varNames = Tobacco Medicinal Smoky Body Spicy Winey Nutty Honey Malty Fruity Sweetness Floral;

/* DEFINE AN IML IMPLEMNTATION OF NMF */
/* Use the Lee and Seung "multiplicative update" algorithm for NMF.
   This isn't the best algorithm, but it is simple to implement in IML
   by folloeing the pseudocode at Wikipedia:
   https://en.wikipedia.org/wiki/Non-negative_matrix_factorization

   This implementation is a work in progress.
   Stay tuned for other NMF algorithms in IML.
*/
proc iml;
start Initialize_LR_SVD(L, R, X, k);
   call SVD(A, D, B, X); 
   L =   A[ ,1:k] # sqrt(D[1:k])`  <> 0;  /* clip negative values */
   R = T(B[ ,1:k] # sqrt(D[1:k])`) <> 0;
finish;
start Initialize_LR_Random(W, H, X, k);
   /* Initialize W and H to the mean values of X, then add random noise */
   mu = X[:];          /* mean of all elements in X */
   sd = stddev( colvec(X) );
   c = sqrt(mu / k);
   L = randfun(n//k, "Normal", c, sd/3 );  /* most values are within 1 SD of c */
   R = randfun(k//m, "Normal", c, sd/3 );
   L = L <> 0;  /* clip negative values */
   R = R <> 0;
finish;

/* COMPARE SVD+CLIP to SVD+MEAN */
start Initialize_LR(L, R, X, k, method="SVD");
   if upcase(method)="SVD" then 
      run Initialize_LR_SVD(L, R, X, k);
   else 
      run Initialize_LR_Random(L, R, X, k);
finish;

/* standardize the rows of R, then scale L accordingly */
start Stdize_LR(L, R, method="L2");
   if upcase(method)="L2" then
      wts = R[,##];                   /* sum of squares of rows of R */
   else if upcase(method)="L1" then
      wts = abs(R)[,+];               /* sum of abs values of rows of R */
   else return;

   wts = choose(wts=0, 1, wts);       /* avoid division by 0 */
   R = R / wts;
   L = L # wts`;
finish;

start Permute_LR(L, R, method="L2");
   if upcase(method)="L2" then
      stat = L[##, ];                  /* sum of squares of cols of L */
   else if upcase(method)="L1" then
      stat = abs(L)[+,];               /* sum of abs values of rows of R */
   else return;
   /* sort columns of L in descending order by statistic */
   call sortndx(idx, stat`, 1, 1);   
   L = L[ ,idx];
   R = R[idx, ];
finish;
 
/* Lee and Seung "multiplicative update" algorithm for NMF */
start nmf_mult(W, H,  /* OUTPUT args */
                X, k, convergenceCrit = 1E-6, maxIter=10000);
   eps = constant("MACEPS");    /* factor to prevent division by 0 */
   minIter = 100;  /* do at least this many iterations to avoid premature convergence */

   run Initialize_LR_SVD(W, H, X, k);     /* initialize W and H */
   run Stdize_LR(W, H);
   Xr = W*H;
   prev_error = norm(X-Xr, "Frob");   /* Frobenius Norm */
   delta = convergenceCrit + 1;       /* ensure delta > convergenceCrit */
   do iter = 1 to maxIter  while(delta >= convergenceCrit);
      /* H Update */
      F1 = (W` * X) / ((W` * W) * H + eps); 
      H = H # F1;
       
      /* W Update */
      F2 = (X * H`) / (W * (H * H`) + eps);
      W = W # F2;

      run Stdize_LR(W, H);
      /* Check convergence  */
      Xr = W * H;
      curr_error = norm(X-Xr, "Frob"); /* Frobenius Norm */           
      /* delta = relative change, but ensure a minimum number of iterations */
      if iter > minIter then do;
         delta = abs(prev_error - curr_error) / prev_error;
      end;
      prev_error = curr_error;
   end;
   *print "NMF algorithm exits after iteration = " iter;
   run Permute_LR(W, H);
finish;

store module = (
Initialize_LR_SVD
Initialize_LR_Random
Initialize_LR
Stdize_LR Permute_LR
nmf_mult
);

/*******/

/* Explore the whisky analysis by following
   Young, S., Fogel, P., and Hawkins, D. (2006)
*/

%let varNames = Tobacco Medicinal Smoky Body Spicy Winey Nutty Honey Malty Fruity Sweetness Floral;
proc iml;
/* Apply NMF to the Scotch whisky data */
varNames = propcase({&varNames});
use Whisky;
   read all var varNames into X;
   read all var {"Distillery" "Selected"};
close;

/* The PCA analysis of the whisky data kept four PCs.
   Load the nmf_mult function and compute a four-factor NMF */
load module = _ALL_;
names = 'NMF 1':'NMF 4';
run nmf_mult(W1, H1, X, ncol(names));

/* scale the H matrix to use the same scale used in the PCA article */
maxValPCA = 0.8506;             /* max value in PCA eigenbasis */
scale = maxValPCA / H1[,<>];    /* scale each row by maximum, then multiply by maxVal */
W = W1 / scale`;
H = H1 # scale;
*print  (W[<>,])[L="Max W Cols"],   (H[,<>])[L="Max H Rows"];

%let Reds = CXF7F7F7 CXFDDBC7 CXF4A582 CXD6604D CXB2182B ;
ods graphics / width=360px height=480px;
call heatmapcont(H`) title="NMF Factors (H`) for Whiskies"
                     colorramp={&Reds} range=(0||max(H))
                     xvalues=names yvalues=varNames;

/* show that the profile vectors are not orthogonal */
R = corr(H`);
print R[c=names r=names F=Bestd5.];

ods graphics / width=480px height=640px;
title "Comparison of W` (NMF)";
/* choose a subset of whiskies */
selIdx = loc(Selected);
W_sel = W[selIdx, ];
Whisky_name = Distillery[selIdx];
call heatmapcont(W_sel) title="NMF Characteristics of Selected Whiskies (W)"
         colorramp={&Reds}
         range=(0||max(W))
         xvalues=names yvalues=Whisky_name;

/* To compute the variance explained by an arbitrary linear subspace, 
   you can project the covariance matrix of the data onto that subspace.
   Because the rows of H are not orthogonal, you 
   must first construct an orthonormal basis for the rowspace 
   of H.
*/
/* The total variance is the sum of diag of S = cov(X) */
S = cov(X);
totalVar = trace(S);
print totalVar;

/* Let S be nxn covariance matrix. Let A be nxk matrix with linearly indep columns.
   What proportion of the variance in S is explained by column space of A?
   Method: Let Q be an orthonomal basis for the column space of A.
   The projection matrix onto this subspace is P = Q*Q`.
   The covariance of the projected data turns out to be P*S*P.
   So, explained variance is trace(P*S*P).
   The trace of a product of square matrices is invariant under a cyclic
   shift. See https://en.wikipedia.org/wiki/Trace_(linear_algebra)#Cyclic_property
   Thus
*/
start VarExplained(A, S);
   /* Find orthonormal basis, Q, for the span of A */
   call svd(U, D, V, A);
   k = ncol(A);
   Q = U[, 1:k]; 

   var_explained = trace(Q`*S*Q);   /* variance explained by colspace(W) */
   total_var = trace(S);            /* total variance in S */
   prop_var = var_explained / total_var; /* proportion of var explained */
   return prop_var;
finish;

/* For the PCA, the variance explained by the first 4 PCs is <70%.
   What proportion of variance is explained by the rowspace of H? */
call eigen(D, U, S);                 /* PCA decomp */
U = U[, 1:4];
prop1 = trace(U`*S*U) / trace(S); print prop1;
prop_var_PCA = VarExplained(U, S);
prop_var_NMF = VarExplained(H`, S);
print prop_var_PCA[F=PERCENT7.1] prop_var_NMF[F=PERCENT7.1];

QUIT;
