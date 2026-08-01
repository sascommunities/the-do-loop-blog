/* The first part of this files is
   a SAS program to accompany the article 
   "An NMF analysis of Scotch whiskies"
   by Rick Wicklin, published 30MAR2026 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2026/03/30/nmf-whisky.html

   The second part of this file examines 
   Hoyer's method for NMF with sparseness constraints 
   Hoyer (2004) https://www.jmlr.org/papers/volume5/hoyer04a/hoyer04a.pdf
   GitHub repo: https://github.com/aludnam/MATLAB/tree/master/nmfpack
*/

/********************************/
/* Part 1: Unconstrained NMF    */
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
      wts = sqrt(R[,##]);             /* L2 norm rows of R */
   else if upcase(method)="L1" then
      wts = abs(R)[,+];               /* L1 norm rows of R */
   else return;

   wts = choose(wts=0, 1, wts);       /* avoid division by 0 */
   R = R / wts;
   L = L # wts`;
finish;

start Permute_LR(L, R, method="L2");
   if upcase(method)="L2" then
      stat = sqrt(L[##, ]);            /* L2 norm of cols of L */
   else if upcase(method)="L1" then
      stat = abs(L)[+,];               /* L1 norm of rows of R */
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

/*****************************************/
/* Part 2: Sparseness-constrained NMF    */
/*****************************************/

/**********************************************************
   Add Hoyer's method for NMF with sparseness constraints 
   Hoyer (2004) https://www.jmlr.org/papers/volume5/hoyer04a/hoyer04a.pdf
   GitHub repo: https://github.com/aludnam/MATLAB/tree/master/nmfpack
 **********************************************************/
proc iml;
/* The ProjFunc function translates the logic in Hoyer's MATLAB function projfunc.m. 
   It forces a vector to satisfy a specific L1 norm (k1) and 
   squared L2 norm (k2), while projecting the coefficients to be strictly 
   nonnegative. 

   In other words, given a vector, v, find the vector v having 
   ||v||_1=k1 and ||v||_2^2=k2 which is closest to v in the Euclidian sense.
*/
start ProjFunc(s, k1, k2);
   v = colvec(s);
   N = nrow(v);
   v = v + (k1 - v[+]) / N;   /* project v to the L1 hyperplane */
   
   maxIter = 1000;
   converged = 0;
   /* move the vector v radially outward within the L1 hyperplane until 
      it intersects the L2=k2 sphere. */
   do j = 1 to maxIter while(^converged);
      /* Identify zeroed coefficients */
      zerocoeff = loc(v <= 0);
      num_zeros = ncol(zerocoeff);
      
      /* Compute the proposed projection operator onto L2 constraint */
      midpoint = j(N, 1, k1 / (N - num_zeros));
      if num_zeros > 0 then midpoint[zerocoeff] = 0;
      
      /* solve quadratic equation */
      w = v - midpoint;
      a = w[##]; 
      b = 2 * (w` * v);
      c = v[##] - k2;
      discrim = b**2 - 4*a*c;
      if discrim < 0 then discrim = 0; /* Equivalent to real() in MATLAB */
      alphap = (-b + sqrt(discrim)) / (2*a);
      
      v = alphap*w + v;   /* now v satisfies the L2 constraint */
      
      /* if solution is nonnegative, leave. Vector satisfies L1 and L2 constraints. */
      if min(v) >= 0 then
         converged = 1;
      else do;
         /* Set negs to zero, subtract appropriate amount from rest */
         zerocoeff = loc(v <= 0);
         v[zerocoeff] = 0;
         num_zeros = ncol(zerocoeff);
         v = v + (k1 - v[+]) / (N - num_zeros);
         v[zerocoeff] = 0;
      end;
   end;
   return(v);     /* always return a column vector */
finish;

/* Implement NMF with sparseness constraints (Hoyer, 2004).
   Translated from Hoyer's MATLAB code in nmfsc.m.
   Initialization Note: 
   Hoyer projects the random matrices to have the correct initial sparseness. 
   I apply this same logic after Initialize_LR_SVD and Stdize_LR routines. 
   However, I scale W's projection by its initial column norms to avoid 
   accidentally destroying the SVD-derived scale.  

   The 'sparseness' parameter is a 2-element vector where
   sparseness[1] = s_W and sparseness[2] = s_H in Hoyer's paper.
   You can pass in a missing value for one or both elements to skip
   the constraint.  In the Scotch whisky example, I constrain 
   only H, not W,
*/
start nmf_mult_sc(W, H, X, k, sparseness, convergenceCrit = 1E-6, maxIter=10000);
   eps = constant("MACEPS");    
   minIter = 100;  
   
   vdim = nrow(X);
   samples = ncol(X);
   
   sW = sparseness[1];
   sH = sparseness[2];
   
   has_sW = (sW ^= .);
   has_sH = (sH ^= .);
   
   /* Initialize W and H using existing methodology */
   run Initialize_LR_SVD(W, H, X, k);     
   run Stdize_LR(W, H);  
   
   /* Apply initial sparseness projections if constraints exist */
   if has_sW then do;
      L1a = sqrt(vdim) - (sqrt(vdim)-1)*sW;
      norms = sqrt(W[##, ]); /* L2 norm of columns */
      do i = 1 to k;
         W[,i] = ProjFunc(W[,i], L1a*norms[i], norms[i]**2);
      end;
   end;
   if has_sH then do;
      L1s = sqrt(samples) - (sqrt(samples)-1)*sH;
      do i = 1 to k;
         /* ProjFunc returns a col vector; transpose to place in row */
         H[i,] = ProjFunc(H[i,], L1s, 1)`; 
      end;
   end;
   
   /* Compute initial objective error */
   Xr = W*H;
   prev_error = 0.5 * norm(X-Xr, "Frob")**2;   
   delta = convergenceCrit + 1;       
   
   stepsizeW = 1;
   stepsizeH = 1;
   min_step = 1E-10; /* Threshold for abandoning line search */
   
   do iter = 1 to maxIter while(delta >= convergenceCrit);
      /* --- H Update --- */
      if has_sH then do;
         dH = W` * (W*H - X);
         begobj = 0.5 * norm(X - W*H, "Frob")**2;
         H_flag = 0;
         do while(^H_flag);
            Hnew = H - stepsizeH * dH;
            do i = 1 to k;
               Hnew[i,] = ProjFunc(Hnew[i,], L1s, 1)`;
            end;
            newobj = 0.5 * norm(X - W*Hnew, "Frob")**2;
            
            if newobj <= begobj then
               H_flag = 1;
            else do;
               stepsizeH = stepsizeH / 2;
               if stepsizeH < min_step then H_flag = 1; 
            end;
         end;
         stepsizeH = stepsizeH * 1.2;
         H = Hnew;
      end;
      else do;
         /* Unconstrained multiplicative update for H */
         F1 = (W` * X) / ((W` * W) * H + eps); 
         H = H # F1;
         /* Renormalize so rows of H have constant energy */
         run Stdize_LR(W, H); 
      end;
      
      /* --- W Update --- */
      if has_sW then do;
         dW = (W*H - X) * H`;
         begobj = 0.5 * norm(X - W*H, "Frob")**2;
         W_flag = 0;
         do while(^W_flag);
            Wnew = W - stepsizeW * dW;
            norms = sqrt(Wnew[##, ]); /* L2 norm of columns */
            do i = 1 to k;
               Wnew[,i] = ProjFunc(Wnew[,i], L1a*norms[i], norms[i]**2);
            end;
            newobj = 0.5 * norm(X - Wnew*H, "Frob")**2;
            
            if newobj <= begobj then
               W_flag = 1;
            else do;
               stepsizeW = stepsizeW / 2;
               if stepsizeW < min_step then W_flag = 1;
            end;
         end;
         stepsizeW = stepsizeW * 1.2;
         W = Wnew;
      end;
      else do;
         /* Unconstrained multiplicative update for W */
         F2 = (X * H`) / (W * (H * H`) + eps);
         W = W # F2;
      end;
      
      /* Check convergence */
      Xr = W * H;
      curr_error = 0.5 * norm(X-Xr, "Frob")**2; 
      if iter > minIter then
         delta = abs(prev_error - curr_error) / prev_error;
      prev_error = curr_error;
   end;
   run Permute_LR(W, H);
finish;

/* Hoyer's sparseness measure.
   For a column vector, v, this function returns 
   the Hoyer (2004) measure of sparseness
       s = (sqrt(n) - ||v||_1 / ||v||_2) / (sqrt(n)-1)
   See https://blogs.sas.com/content/iml/2026/07/27/hoyer-sparseness.html
  */
start Sparseness(v);
   u = colvec(v);    /* always treat input as a column vector */
   L2 = norm(u, "L2");
   if L2=0 then return(.);
   L1 = norm(u, "L1");
   sn = sqrt( countn(u) );
   sparseness = (sn - L1/L2) / (sn - 1);
   return( clip01(sparseness) );
finish;
/* Helper function: clip values into [0,1]
   If x is any matrix, return a new matrix whose i_th element has the value 
   0    if x[i] < 0
   x[i] if 0 <= x[i] <= 1
   1    if x[i] > 1
   See https://blogs.sas.com/content/iml/2026/02/04/clip-values.html */
start clip01(x);
   return( 0 <> (x >< 1) );
finish;

store module = (
ProjFunc
nmf_mult_sc
Sparseness
clip01
);

/*******/

/* Example: implement a constrained NMF analysis of the whisky data 
   (first analyzed by Young, S., Fogel, P., and Hawkins, D. (2006))
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
k = ncol(names);
run nmf_mult(W1, H1, X, k);

/* Recall that the NMF is not unique. In particular, you 
   can permute the rows of H and the corresponding columns of W 
   without changing the approximation. That is, if P is any permutation
   matrix, then X = W*H = (W*P)*(P`*H). Let's permute the 
   rows of H to maintain the same columns that were shows in the article 
   https://blogs.sas.com/content/iml/2026/03/30/nmf-whisky.html 
*/
perm = {3,1,2,4};
W1 = W1[,perm];
H1 = H1[perm,];

H = H1;
print (H1`)[c=names r=varNames];

/* Compute baseline sparseness of unconstrained H.
   (Assumes your Sparseness() function from the earlier post is loaded.)
   Note that EVERY row of H has the SAME sparseness! */
s_H1 = j(1, nrow(H), .);
do i = 1 to nrow(H);
   s_H1[i] = Sparseness(H[i,]);
end;
print s_H1[F=5.3 L="Baseline Sparseness (H1 rows)"];

%let Reds = CXF7F7F7 CXFDDBC7 CXF4A582 CXD6604D CXB2182B ;
ods graphics / width=360px height=480px;
call heatmapcont(H`) title="Unconstrained NMF Factors (H`) for Whiskies"
                     colorramp={&Reds} range={0 1}
                     xvalues=names yvalues=varNames;
/* interpretation:
   H (4x12): The right factor contains the flavor profiles. 
             If you want distinct, easily interpretable profiles, 
             you want a higher S_H. For example, a classic Islay 
             flavor profile should load heavily on just a few 
             flavors like "Smoky" and "Medicinal," 
             while leaving "Floral" or "Sweetness" near zero. 
             Maybe aim for a sparseness constraint of 0.6-0.75 to 
             force these distinct profiles.
   W (86x4): These are the weights indicating how much of each 
             flavor profile belongs to a specific whisky. 
             If you constrain S_W, you are forcing the algorithm 
             to describe each whisky using only one or two profiles, 
             rather than a blend of all four. If you leave S_W 
             unconstrained (by using a missing value), a whisky 
             can be an equal blend of multiple profiles.
   So let's force the flavor profiles (H) to be sparse, but leave the 
   weights (W) unconstrained.
*/
/* What is a good value of sparseness to use? Well, there are 
   12 variables and we are using k=4, so a perfectly sparse factorization would
   use 3 or 4 nonzero components for each factor. A 12-D vector that has all
   zeros except for 3 or 4 equal components has sparseness measure
   between 0.6-0.7.
*/
sparseness = {., 0.6};  /* Unconstrained W, highly sparse H */
run nmf_mult_sc(W2, H2, X, k, sparseness);
/* permute the rows of H to maintain the same columns that were shows in the article 
   https://blogs.sas.com/content/iml/2026/03/30/nmf-whisky.html 
*/
perm = {4,1,2,3};
W2 = W2[,perm];
H2 = H2[perm,];
*print (H2`)[c=names r=varNames];
H = H2;

%let Reds = CXF7F7F7 CXFDDBC7 CXF4A582 CXD6604D CXB2182B ;
ods graphics / width=360px height=480px;
call heatmapcont(H`) title="Sparseness-Constrained NMF Factors (H`) for Whiskies"
                     colorramp={&Reds} range={0 1}
                     xvalues=names yvalues=varNames;

s_H2 = j(1, nrow(H), .);
do i = 1 to nrow(H);
   s_H2[i] = Sparseness(H[i,]);
end;
print s_H2[F=5.3 L="New Sparseness (H2 rows)"];

/* let's look at the errros in each approximation. 
   A theorem says that the least possible error comes from the rank-k SVD,
   so let's use the rank-k SVD to standardize the errors. */
call SVD(U, Q, V, X);
U = U[,1:k];
Q = Q[1:k];
V = V[,1:k];
errorSVD = norm(X - U*diag(Q)*V`, "Frob");

errorUnconstrained = norm(X - W1*H1, "Frob");
error2 = norm(X - W2*H2, "Frob");
relError1 = (errorUnconstrained - errorSVD) / errorSVD;
relError2 = (error2 - errorSVD) / errorSVD;
print relError1[L='RelErr Unconstrained NMF'] relError2[L='RelErr Sparse NMF'];

/* let's look at how the relative error changes as we constrain the sparsess in NMF */
sparseness = {., .};  /* Unconstrained W; sparse H */
sparse = do(0.05, 0.95, 0.025);
sparseError = j(1, ncol(sparse), .);
do i = 1 to ncol(sparse);
   sparseness[2] = sparse[i];
   run nmf_mult_sc(W, H, X, k, sparseness);
   sparseError[i] = norm(X - W*H, "Frob");
end;

ods graphics / width=640px height=480px;
title "Approximation Error: Rank-4 Sparse NMF";
refStmt = catx(" ","refline ",char(errorUnconstrained),"/axis=y noclip label='Unconstrained NMF';");
call series(sparse, sparseError) grid={x y} label={'Hoyer Sparseness' 'Relative Error'}
            other=(refStmt +"yaxis grid min=0;") option="lineattrs=(thickness=3)";
            
/* and let's generate a few heat maps of H1 for various sparseness values */
ods graphics / width=360px height=480px;
sparse = {0.2  0.4  0.6  0.8};
do i = 1 to ncol(sparse);
   sparseness[2] = sparse[i];
   run nmf_mult_sc(W, H, X, k, sparseness);
   t = "H' Sparseness=" + strip(char(sparse[i])) + " for Whisky Data";
   call heatmapcont(H`) title=t
                     colorramp={&Reds} range={0 1}
                     xvalues=names yvalues=varNames;
end;

QUIT;
