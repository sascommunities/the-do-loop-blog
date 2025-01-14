/* RandCovMatSim.sas */

/* SAS code to accompany 
   "Generating a random orthogonal matrix"
   by Rick Wicklin. Published on The DO Loop blog on 28MAR2012
   https://blogs.sas.com/content/iml/2012/03/28/generating-a-random-orthogonal-matrix.html
   and 
   "Generate a random matrix with specified eigenvalues"
   by Rick Wicklin. Published on The DO Loop blog on 30MAR2012.
   https://blogs.sas.com/content/iml/2012/03/30/geneate-a-random-matrix-with-specified-eigenvalues.html

   To use these modules, run or %INCLUDE the RandCorrMatrices.sas file. Then use
   LOAD module=_all_;
*/
/**********************************************/
/* STORE modules by running RandCorrMatrices.sas */
/**********************************************/
%include "RandCorrMatrices.sas";


proc iml;
LOAD module=_all_;

call randseed(1);
Q=RndOrthog(4);
reset fuzz;
Ident = Q`*Q;
print Ident, Q;

N = 200;
evals = j(N, 2, 0); /* 2*N complex values */
p = 2;
do i = 1 to N by p;
   e = eigval( RndOrthog(p) );
   evals[i:i+p-1, 1:ncol(e)] = e;
end;
create C2 from evals[c={"x" "y"}]; append from evals; close;

/* eigenvalues of 1 200x200 matrix */
evals = eigval( RndOrthog(N) );
create C200 from evals[c={"x" "y"}]; append from evals; close;

proc kde data=c2; /* Eigenvalues of 100 2x2 Matrices */
bivar x y / plots=ContourScatter bwm=0.5;
run;
 
proc kde data=c200; /* Eigenvalues of 200x200 Matrix */
bivar x y / plots=ContourScatter bwm=0.5;
run;

proc iml;
LOAD module = _all_;
 
lambda = {2 1 0.75 0.25};
call randseed(1234);
S1 = RndMatWithEigenval(lambda); /* eigenvalues are lambda */
eval = eigval(S1);
print eval;
print S1;

S2 = RndMatWithEigenval(lambda); /* eigenvalues are lambda */
print S2;
