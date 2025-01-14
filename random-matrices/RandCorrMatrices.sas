/* RandCorrMatrices.sas */

/* SAS code to accompany "Generate correlation matrices with specified eigenvalues"
by Rick Wicklin. Published on The DO Loop blog on 18DEC2024.
https://blogs.sas.com/content/iml/2024/12/18/correlation-matrix-eigenvalues.html

Use the method of alternating projections (MAP) to find a correlation 
matrix that has a specified spectrum (eigenvalues). The projections are:
1. Choose a spectrum Lambda = {L1, L2, ..., Ln} where L_i>0 and sum(Lambda)=n.
2. Create a covariance matrix, S, that has the specified eigenvalues.
3. Convert the matrix to a correlation matrix, R.
4. If R has the desired eigenvalues, quit. Otherwise,
   use the eigenstructure of R to create a new covariance matrix that has
   the same eigenvectors but the desired eigenvalues. Go to Step 3.

I have previously use MAP to find the nearest correlation matrix:
https://blogs.sas.com/content/iml/2012/11/28/computing-the-nearest-correlation-matrix.html
This implementation of MAP uses similar functions:
Given a spectrum L={L1, L2, ..., LN} start with a random SPD
matrix Q*diag(L)*Q` where Q is an arbitrary orthogonal matrix.
Then use 
ProjU: project a covariance matrix, S, onto U={correlation matrix for S}
ProjS: project a correlation matrix onto 
       S(lambda)={positive semidefinite matrices that has spectrum L}

To use these modules, run or %INCLUDE this file. Then in a PROC IML 
program, use
LOAD module=_all_;
*/
proc iml;
/* To get a random orthogonal matrix, use 
   https://blogs.sas.com/content/iml/2012/03/28/generating-a-random-orthogonal-matrix.html
   To get a random SPD matrix, use 
   https://blogs.sas.com/content/iml/2012/03/30/geneate-a-random-matrix-with-specified-eigenvalues.html
*/
/* Generate random orthogonal matrix G. W. Stewart (1980).    
   Algorithm adapted from QMULT MATLAB routine by Higham (1991) */
start RndOrthog(n);
A = I(n);
d = j(n,1,0);
d[n] = sgn(RndNormal(1,1));
do k = n-1 to 1 by -1; 
   /* generate random Householder transformation */
   x = RndNormal(n-k+1,1);
   s = sqrt(x[##]); /* norm(x) */
   sgn = sgn( x[1] );
   s = sgn*s;
   d[k] = -sgn;
   x[1] = x[1] + s;
   beta = s*x[1];
   /* apply the transformation to A */
   y = x`*A[k:n, ];
   A[k:n, ] = A[k:n, ] - x*(y/beta);
end;
A = d # A; /* change signs of i_th row when d[i]=-1 */
return(A);
finish;
 
/**** Helper functions ****/
/* Compute SGN(A). Return matrix of same size as A with 
   m[i,j]= {  1 if A[i,j]>=0
           { -1 if A[i,j]< 0           */
start sgn(A);
return( choose(A>=0, 1, -1) );
finish;
 
/* return (r x c) matrix with standard normal variates */
start RndNormal(r,c);
   x = j(r,c);
   call randgen(x, "Normal");
   return(x);
finish;

start RndMatWithEigenval(_lambda);
   lambda = rowvec(_lambda);   /* ensure lambda is row vector */
   n = ncol(lambda);
   Q = RndOrthog(n);
   return( Q`*diag(lambda)*Q ); 
finish;

/* efficient computation of A*D*A` where D = diag(v) and v[i]>0 for all i */
start QuadFormDiagPos(A, v);
   D = rowvec(v);     /* ensure D is row vector */
   W = A#sqrt(D);     /* left half of A*diag(v)*A` */
   return( W*W` );    /* W*W` = A*diag(v)*A` */
finish;

/* Project symmetric A onto S = {positive semidefinite matrices} */
start ProjS(A, lambda);
   call eigen(D, Q, A);             /* A = Q*D*Q` */
   S = QuadFormDiagPos(Q, lambda);  /* S = Q*diag(lambda)*Q` */
   return S;
finish;
 
/* project SPD matrix S onto U={correlation matrices}.
   Return correlation matrix for S. */
start ProjU(S);
   return cov2corr(S);
finish;
 
/* Helper function: the L2 norm of the difference between 
   targetLambda and the eigenvalues of a symmetric matrix A.
   Assume targetLambda is sorted from largest to smallest. */
start DistBetweenEigenvalues(A, targetLambda);
   lambda = eigval(A);
   return( sqrt(SSQ(targetLambda - lambda)) );
finish;

/* for column vector v, standardize so sum(v)=nrow(v) */
start StdizeEigenval(v);
    L = colvec(v);
    return ( L * nrow(L)/sum(L) );  /* make sum(v)=dimension */
finish;

/* Encapsulate the MAP algorithm into a function */
start CorrWithEigen(targetLambda, tol=1e-10, maxIters=500);
   lambda = StdizeEigenval(targetLambda);
   /* generate a random orthogonal matrix and use to construct PSD matrix with spectrum */
   S = RndMatWithEigenval(lambda); /* eigenvalues are lambda */
   delta = 1E6;
   do k = 1 to maxIters while(delta > tol);
      /* project onto space of correlation matrices (or matrices with unit diagonal) */
      R = ProjU(S);
      /* project R onto matrices with the spectrum */
      S = ProjS(R, lambda);
      delta = DistBetweenEigenvalues(R, lambda);
   end;
   R = ProjU(S);   /* final projection */
   return R;
finish;

/* MORE helper functions! */
/* extract elements below the diagonal. 
   For a symmetric matrix, this is equivalent to extracting elements above the diagonal */
start StrictLowerTriangular(X);
   v = cusum( 1 || (ncol(X):2) );
   return( remove(vech(X), v) );
finish;
/* extract the name "Rij" for the cells below the diagonal */
start StrictUpperTriLabels(n);
   rows = row(J(n));
   cols = col(J(n));
   v1 = StrictLowerTriangular(rows);
   v2 = StrictLowerTriangular(cols);
   L = compress("R" + char(v2) + char(v1));
   return L;
finish;
/* Concatenate values into a string separated by a delimiter (by default, 
   a blank). Create a macro variable with the specified name.
   https://blogs.sas.com/content/iml/2016/01/18/create-macro-list-values.html
*/
start CreateMacro(values, macroName, delimiter=' ');
   if type(values)='N' then          
      y = rowvec( char(values) );   /* convert numeric to character */
   else 
      y = rowvec(values);
   s = rowcat(y + delimiter);       /* delimit and concatenate */
   s = substr(s, 1, nleng(s)-nleng(delimiter)); /* remove delimiter at end */
   call symputx(macroName, s);      /* create macro variable */
finish;
store module = _all_;
QUIT;
