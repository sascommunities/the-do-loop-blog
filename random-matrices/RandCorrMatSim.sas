/* SAS code to accompany "Generate correlation matrices with specified eigenvalues"
by Rick Wicklin. Published on The DO Loop blog on 18DEC2024.
https://blogs.sas.com/content/iml/2024/12/18/correlation-matrix-eigenvalues.html
*/
/**********************************************/
/* STORE modules by running RandCorrMatrices.sas */
/**********************************************/
%include "RandCorrMatrices.sas";


proc iml;
/* start with a symmetric positive definite matrix (a covariance matrix) */
S = {2.0   0.25  0.16,
     0.25  0.6   0.09,
     0.16  0.09  0.4};
/* The naive attempt to obtain a correlation matrix is to scale the 
   covariance matrix to get unit diagonals, but that destroys the eigenvalues */
R = cov2corr(S);
eigValS = eigval(S);
eigValR = eigval(R);
print R;
print eigValS eigValR;
QUIT;


proc iml;
load module=_all_;   /* modules are defined in the Appendix */
call randseed(1234, 1);
Spectrum = {3, 0.5, 0.3, 0.2};
 
/* call the MAP algorithm multiple times */
R1 = CorrWithEigen(Spectrum);
R2 = CorrWithEigen(Spectrum);
R3 = CorrWithEigen(Spectrum);
 
/* verify that each correlation matrix has the required spectrum */
eigR1 = eigval(R1);
eigR2 = eigval(R2);
eigR3 = eigval(R3);
print Spectrum eigR1 eigR2 eigR3;
/* the spectrum is the same, but the correlations are different */
labl = 'X1':'X4';
print R1[r=labl c=labl], R2[r=labl c=labl], R3[r=labl c=labl];
QUIT;


/* Theorem: If A and B have the same spectrum, then 
   || A || = || B || = || Spectrum ||
   Where || A || is the Frobenius matrix norm and 
         || v || is the Euclidean vector norm sqrt(SSQ(v))
*/
proc iml;
S1 = {
1.6777914558 -0.387406153 -0.443734244 0.0068766222 ,
-0.387406153 0.9137286028 -0.01915772 -0.058629084  ,
-0.443734244 -0.01915772 0.842547134 -0.348769043   ,
0.0068766222 -0.058629084 -0.348769043 0.5659328073 };
 
S2 = {
0.7912845768 0.0216179596 -0.193768339 -0.350178629 ,
0.0216179596 0.7775696813 -0.039066959 -0.295580783 ,
-0.193768339 -0.039066959 1.7418846881 0.4392865998 ,
-0.350178629 -0.295580783 0.4392865998 0.6892610538 };
 
eval1 = eigval(S1);
eval2 = eigval(S2);
print eval1 eval2;

L2Norm1 = norm(S1, "L2");
L2Norm2 = norm(S2, "L2");
print L2Norm1 L2Norm2;
 
FrobNorm1 = norm(S1, "Frob");
FrobNorm2 = norm(S2, "Frob");
print FrobNorm1 FrobNorm2;

FrobNormSq = FrobNorm1**2;
SSQLambda = SSQ(eval1);
print FrobNormSq SSQLambda;
QUIT;
