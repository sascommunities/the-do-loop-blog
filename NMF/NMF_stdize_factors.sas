/* SAS program to accompany the article 
   "On the nonuniqueness of low-rank matrix factors"
   by Rick Wicklin, published 16MAR2026 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2026/03/16/on-the-nonuniqueness-of-low-rank-matrix-factors.html

   Thesis: You can remove some (but not all!) ambiguities in a rank-k 
   two-factor approximation to a matrix if you choose a way to standardize the factors.

   I want to discuss rank-k products and how they can be represented as a product L*R. 
   Let Y be an nxm matrix that is rank k. 
   Then SVD says that Y = U*Sigma*V` where U is nxk and V is mxk. 
   Immediately, we see that there is nonuniqueness in the factorization of Y, 
   because I can assign L=U*Sigma^{\alpha} and R=Sigma^{1-\alpha}*V` for any 
   choice of \alpha. 

   Nonuniqueness is a problem. It is hard for researchers to communicate. 
   Standardization is a way to eliminate or diminish the problem.  

   This shows a general way that we can use to scale the L and R matrices. 
   If Y=L*R and D is ANY diagonal matrix with positive elements on the diagonal, 
   then Y=(L*D)*(D^{-1}*R), which means that you can scales the rows of R 
   or the columns of L. This helps to standardize the process and remove 
   an ambiguity when representing a rank-k matrix as a product. 
   Some convenient choices are:
   (1) Standardize the rows of R to have unit sums-of-squares (L2 standardization), or 
   (2) Standardize the rows of R to have unit sums-of-absolute values (L2 standardization). 
       If a row of R is entirely zero, use 1 as the scale factor.
*/
proc iml; 
reset fuzz;
/* Define a simple matrix with an underlying rank of 2 */
Y = { 8 11 5,
     11 12 5,
     12  9 3,
     17 14 5 };
 
/* Use SVD to find the factors */
call svd(U, Sigma, V, Y);
 
/* Extract the first k=2 columns to get a rank-2 product */
i = 1:2;
U2 = U[,i];
S  = Sigma[i];
V2 = V[,i];
 
/* You can distribute the scaling of S into the left factor, 
   the right factor, or both factors by varying alpha */
alpha = {0, 0.5, 1};
do j = 1 to nrow(alpha);
   L = U2 * diag(S##alpha[j]);
   R = T( V2 * diag(S##(1-alpha[j])) );
   /* The product L * R is identical for all choices of alpha! */
   prod = L*R;
   print prod;
end;

/* Permute the columns of L and rows of R */
P = {0 1,
     1 0};
L_new = L * P;
R_new = P` * R;
Y_check = L_new * R_new; /* perfectly reconstructs Y */
print Y_check;

/* A can be ANY invertible matrix. */
A = { 0.9  0.1,
      0.1  0.9};
L_new = L * A;
R_new = inv(A) * R;
Y_check = L_new * R_new; /* the product is preserved */
print Y_check;

L1 = {1 2 ,
      3 8 ,
      5 4 };
R1 = {22 21 20 19,
      12 22 11 10};
Y = L1 * R1;
ods layout gridded columns=3 advance=table;
print Y, L1, R1;
ods layout end;

A = {0.9  0.1, 0.1  0.9};
invA = inv(A);
ods layout gridded columns=2 advance=table;
print A, invA;     /* Note: inv(A) is not nonnegative! */
ods layout end;

/* nevertheless, the products (L*A) and (inv(A)*R) can still be nonnegative */
L2 = L1*A;
R2 = inv(A)*R1;
Y_check = L2 * R2;
ods layout gridded columns=3 advance=table;
print Y_check, L2, R2;
ods layout end;

/* What if we use different algorithms and one algorithm returns
   (L1,R1) whereas the other returns (L2,R2)? Can we find an invertible 
   matrix, A, that maps one solution to the other? */
A_from_L = ginv(L1) * L2;
A_from_R = R1 * ginv(R2);
/* Do the two LS solutions give the same matrix? */
max_diff = max(abs(A_from_L - A_from_R));
print max_diff;
ods layout gridded columns=2 advance=table;
print A_from_L, A_from_R;     /* Note: inv(A) is not nonnegative! */
ods layout end;

QUIT;
