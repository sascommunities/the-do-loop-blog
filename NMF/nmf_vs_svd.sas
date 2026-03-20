/* SAS program to accompany the article 
   "A comparison of SVD and NMF factorizations"
   by Rick Wicklin, published 23MAR2026 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2026/03/23/nmf-vs-svd.html

   This program compares the SVD and NMF factoriztions on simulated data.
   The example is a small modification of the example in  
   Section 4.2 "Realistic synthetic mixture," p. 212, of 
   Fogel, P., Hawkins, D. M., Beecher, C., Luta, G., & Young, S. S. (2013),
   "A tale of two matrix factorizations", The American Statistician, 67(4), 207-218.
   https://doi.org/10.1080/00031305.2013.845607
*/

/* 
A comparison of SVD and NMF factorizations

INTRODUCTION:
Example from 
Fogel, P., Hawkins, D. M., Beecher, C., Luta, G., & Young, S. S. (2013),
"A tale of two matrix factorizations", The American Statistician, 67(4), 207-218.
https://doi.org/10.1080/00031305.2013.845607
This example is from Section 4.2 "Realistic synthetic mixture," p. 212.

In this article, I use SAS IML to simulate data from a clinical study 
that involves gene expression levels for 30 genes across 40 patient samples.
Because the data are simulated, we know that the data
matrix is rank 2, and we know a set of true rank-2 matrices whose product 
equals the data matrix. The purpose of the article is to compare the 
singular value decomposition (SVD) and the nonnegative matrix factorization
(NMF). Although the SVD can produce factors that exactly generate the rank-2
data matrix, the factors contain both positive and negative values, which 
can make interpreting the factors difficult in some applications. 
In contrast, the NMF produce factors that contain only nonnegative values.
However, the product of the factors might not exactly equal the 
data matrix.

SIMULATING PATIENT-GENE DATA FOR A DISEASE

Assume a disease (such as diabetes) has two distinct genetic etiologies.
The L matrix represents 40 patients.
The first column of L represents a measure of one etiology.
The second column of L represents a measure of a second etiology.
In this example, 10 patients have only the first etiology,
10 have only the second etiology, and 10 have both.

The R matrix represents expression levels for 30 genes. In this example, we assume the first 10 
genes are related to the first etiology. The next 10 genes are related to the 
second etiology. The last 10 genes are not related to the disease.

We simulate the L and R matrix by randomly assigning a value between 1 and 2 to 
specific cells. Notice that these factors are entirely nonnegative,
so the product matrix, X = L*R, is also nonnegative.
The following SAS IML statements simulate the L and R matrices 
and create heat maps to visualize them.
*/
title;
proc iml;
reset fuzz;
/*** Simulation of L and R for clinical example of two etiologies ***/
call randseed(123);
nrow = 40;
ncol = 30;
L = j(nrow, 2, 0);
idx1 = (11:20) || (31:40);
L[idx1, 1] = randfun( ncol(idx1), "Uniform", 1, 2); /* 1st col */
idx2 =  21:40;
L[idx2, 2] = randfun( ncol(idx2), "Uniform", 1, 2); /* 2nd col */

R = j(2, ncol, 0);
idx1 = 1:10;
R[1, idx1] = randfun( 1|| ncol(idx1), "Uniform", 1, 2); /* 1st row */
idx2 = 11:20;
R[2, idx2] = randfun( 1||ncol(idx2), "Uniform", 1, 2);  /* 2nd row */

/* simplify L and R by rounding to 2 decimal places */
L = round(L, 0.01);
R = round(R, 0.01);
/*** end simulation of L and R ***/

/* create blue-white-red color ramp from ColorBrewer */
%let BlueReds = CX2166AC CX4393C3 CX92C5DE CXD1E5F0 CXF7F7F7 CXFDDBC7 CXF4A582 CXD6604D CXB2182B ;
%let Reds = CXF7F7F7 CXFDDBC7 CXF4A582 CXD6604D CXB2182B ;

/* Visualize L (40 rows, transposed for viewing) and R (30 columns) */
ods graphics / width=640px height=84px; 
names = {'Real 1' 'Real 2'};
call heatmapcont(L`) colorramp={&BlueReds} yvalues=names range={-2.2 2.2}
     title="Original L`";

ods graphics / width=480px height=84px;
names = {'Real 1' 'Real 2'};
call heatmapcont(R) colorramp={&BlueReds} yvalues=names range={-2.2 2.2}
     title="Original R";

/* Note that the heat maps are entirely red or white because all values of 
   L and R are nonnegative. The product, L*R, simulates a data matrix
   that is also nonnegative. The X matrix represents data in a 
   clinical study that involves gene expression levels for 
   30 genes across 40 patient samples.
*/
X = L*R;
ods graphics / width=360px height=480px; 
call heatmapcont(X) colorramp={&Reds} range={-0.01 4.04} /* Note: Dark red now represents 4, not 2.2 */
     title="Original (Synthetic) Data Matrix";


/* In practice, the order of the patients and the order of the genese would
   be permuted arbitrarily. Also, there would be random noise that 
   appears in the 40x30 panel of measurements. Regardless, when presented with this data,
   a researcher would like to recover the underlying factors that 
   generate this data.  If the researchers know or suspect that the disease
   has two etiologies, they might try to approximate the matrix by using a 
   rank-2 decomposition of factors. In this example, we know the true factors 
   whose product forms the X matrix, but in the real world those factors are unknown
   and must be estimated from X.


THE SVD FACTORIZATION

   There are two popular ways to decompose X into rank-2 factors. The first
   is the classical singular value decomposition (SVD). The second 
   is a newer algorithm called the nonnegative matrix factorization (NMF).
   Following Fogel et al. (2013), the remainder of this article compares 
   these two factorization.

   I have previous discussed how a low-rank factorization is not unique, but
   can be arbitrarily rescaled and permuted. Fogel et al. did not worry about 
   this fact and visualized the results by using red-blue color ramps that are
   individually scaled to each factorization. That makes it hard to compare the 
   two models. In this article, I carefully scale the factors and set the 
   color ramps so that red colors represent positive values, white represents zero,
   and blue represents negative values.


   Let's begin by looking at an SVD analysis of the simulated data matrix, X.
   For more about the SVD, see 
   https://blogs.sas.com/content/iml/2017/08/30/svd-low-rank-approximation.html
   The following statements call the SVD subroutine in SAS IML:
   X = U*D*V`

   To get the rank-2 approximation of X, the first two colunms of U and V are 
   extracted. For this example, the rank-2 product exactly equals X because we 
   didn't add noise to X and because we constructed X as the product of two rank-2 matrices.
*/
k = 2;
call SVD(A, D, B, X);
i = 1:k;   
L_SVD =   A[, i] # sqrt(D[i])`;   /* we can scale arbitrarily. Assign help the scale to each factor */
R_SVD = T(B[, i] # sqrt(D[i])`);
X2 = L_SVD*R_SVD;

/* show that the SVD recovers X exactly */
maxDiff = max(abs(X-X2));
print maxDiff;

/* INTERPRETING THE SVD FACTORS

   Although the product of the SVD exactly equals X, the SVD will NOT usually recover the 
   underlying factors (not even up to scaling).
   The reason is that the factors in an SVD are orthogonal, which means that the 
   factors must contains both positive and negative elements when k &ge; 2.
   (Technically, if X is block diagonal,
   then the SVD can recover the original features, but that is not the case in this example.)
   Because X contains strictly nonnegative values and isn't block diagonal, 
   it is mathematically 
   impossible for the SVD to return orthogonal factors that are entirely
   nonnegative.

   This is shown in the following heat maps, which 
   show the original factors and the SVD factors on the same color scale.
*/

/* Visualize true L` (40 columns) and the left factor from the SVD */
ods graphics / width=640px height=84px; 
ods layout gridded columns=1 advance=table row_gutter=5px;
call heatmapcont(L`) colorramp={&BlueReds} yvalues={'Real 1' 'Real 2'}
     range={-2.2 2.2} title="Real Factor L`";
call heatmapcont(L_SVD`) colorramp={&BlueReds} yvalues={'SVD 1' 'SVD 2'}
     range={-2.2 2.2} title="SVD Factor U";
ods layout end;

/* As seen previously, the true left factors for X have strictly nonnegative values.
   However, the second factor in the SVD contains negative values
   because the inner product of the factors must equal 0 by orthogonality.

   In a similar way, you can visualize the right factors:
*/
   
/* Visualize true R (30 columns) and the right factor from the SVD */
ods graphics / width=480px height=84px;
ods layout gridded columns=1 advance=table row_gutter=5px;
call heatmapcont(R) colorramp={&BlueReds} yvalues={'Real 1' 'Real 2'}
     range={-2.2 2.2} title="Real Factor R";
call heatmapcont(R_SVD) colorramp={&BlueReds} yvalues={'SVD 1' 'SVD 2'}
     range={-2.2 2.2} title="SVD Factor V";
ods layout end;

/* Again, one of the SVD factors contains negative values, as required by 
   orthogonality.
*/

/*********************************************************/
/* THE NMF FACTORIZATION

   Because the SVD cannot return nonnegative positive factors,
   researchers created the NMF factorization.  The NMF 
   decomposes a nonnegative matrix X into the product W*H
   where both W and H are nonnegative.   

   Let's decompose X by using NMF to produce rank-2 factors.
   In SAS, you can obtain the NMF by using PROC NMF in SAS Viya, or by implementing
   the algorithm in SAS IML. I will show both of these approaches in a 
   subsequent article. For now, I'll show the result of the NMF
   factorization from running an IML algorithm that I wrote. 
   The NMF returns two rank-2 matrices, W and H, such that 
   X = W * H. As with the SVD, the product is a rank-2 approximation
   to X. However, for this example, the product does not exactly equal X
   even though X is rank-2. 
*/

/*******/
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
/*******/

/* Run the Lee and Seung "multiplicative update" algorithm for NMF.
   The IML implementation can be downloaded from GitHub. */
run nmf_mult(W1, H1, X, 2);

/* First, let's see whether the product W1*H1 is a reasonable
   approximation of the data matrix.
*/
X_est = W1*H1;
ods graphics / width=360px height=480px; 
call heatmapcont(X_est) colorramp={&Reds} range={-0.01 4.04}
     title="NMF Estimate of Matrix"; /* Note: Dark red now represents ~4, not 2.2 */

/* The product W1*H1 does a good job approximating X. 
   An obvious difference is seen for Patient 32 and 39.
   Where these patients actually have positive gene expression levels 
   for genese 1-10, the rank-2 NMF does not recover that fact.
*/

/* COMPARING THE NMF FACTORS

   As explained in a previous article, there are two inherent ambiguities
   in the factors of the NMF.
   1. The column of W and the rows of H can be arbitrarily scaled, provided
      that they are scaled together in an appropriate manner.
   2. The column of W and the rows of H can be arbitrarily permuted,
      provided that they are permuted together.

   To compare the factors from the NMF with the original factors, it is
   best that they are represented on the same scale. We could standardize the
   factors L and R, but since we have already visualized them, I will instead
   adjust the scaling of W and H so that the maximum elements of W and H are 
   approximately equal to 2, which was the maximum value in the factors L and R.

   The following statements rescale W and H to try to put all values in [0,2].
*/
maxVal = 2;
scale = maxVal / W1[<>,];    /* scale each col by maximum, then multiply by maxVal */
W = W1 # scale;
H = H1 / scale`;
*print  (W[<>,])[L='max(scaled W)'] (H[,<>])[L='max(scaled H)'];

/* Depending on the data and the initialization of the NFM algorithm,
   you might need to permute the NMF factor, 
   which can come out in the opposite order from the true factors.
W = W[ ,{2 1}];
H = H[{2 1}, ];
*/

/* Now, let's compare the factors L and R to the NMF factors, W and H.
   As before, we use heat maps. We also transpose L and W so that they 
   fit better on the screen.
*/   
ods layout gridded columns=1 advance=table row_gutter=5px;
ods graphics / width=640px height=70px; 
call heatmapcont(L`) colorramp={&BlueReds} yvalues={'Real 1' 'Real 2'}
     range={-2.2 2.2} title="Real Factor L`";
call heatmapcont(W`) colorramp={&BlueReds} yvalues={'NMF 1' 'NMF 2'}
     range={-2.2 2.2} title="NMF Factor W`";
ods layout end;

/* As seen previously, the true left factors for X have strictly nonnegative values.
   However, the second factor in the SVD contains negative values
   because the inner product of the factors must equal 0 by orthogonality.

   In a similar way, you can visualize the right factors:
*/
   
/* Visualize true R (30 columns) and the right factor from the NMF */
ods layout gridded columns=1 advance=table row_gutter=5px;
ods graphics / width=480px height=84px;
call heatmapcont(R) colorramp={&BlueReds} yvalues={'Real 1' 'Real 2'}
     range={-2.2 2.2} title="Real Factor R";
call heatmapcont(H) colorramp={&BlueReds} yvalues={'NMF 1' 'NMF 2'}
     range={-2.2 2.2} title="NMF Factor H";
ods layout end;

/* This output shows two features of the NMF algorithm:
   1. Both factors contain only nonnegative elements. It follows that the 
      product W*H must also be nonnegative. 
   2. The NMF decomposition does a much better job of recovering the true 
      factors. If you compare the true factors (L` and R) with the 
      NMF factors (W and H, when appropriately scaled), you can see that the 
      factors are very similar. The most obvious difference is that 
      the first row of the W factor from the NMF does not correctly 
      estimate the values for Patients 32 and 39. The NMF estimate indicates
      that these patients only have the second etiology for the disease,
      whereas in reality they have both etiologies. 
*/

/* SUMMARY 
   
   This article uses SAS to create an example that is similar to 
   Section 4.2 in Fogel, et al, (2013).

   The example simulates a data matrix that represents
   the expression profiles of 30 target genes for 40 clinical study participants.
   A classical SVD analysis produced factors whose product equals X,
   but the factors contain both positive and negative values.
   That is because the SVD factors are essential contrasts.
   In contrast, the NMF produce factors that contain only nonnegative values.
   When the true underlying features of the data involve positive quantities
   whose effects are additive in nature,
   the NMF factors are more helpful for interpreting the features.  
*/

