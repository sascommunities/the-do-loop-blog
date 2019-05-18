/* SAS program to accompany the article
   "Gershgorin discs and the location of eigenvalues"
   by Rick Wicklin, published 22MAY2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/05/22/gershgorin-discs-location-eigenvalues.html ?

   The Gershgorin Disc Theorem is: 
   If A is square matrix, every eigenvalue lies within one of the disks
      D_i = {z \in C | |z - A[i,i]| < r_i }
   where r_i = \Sigma |A[i,j]| and j \neq i

   The program computes the centers and radii of the Gershgorin discs
   that are associated with a matrix. The discs, along with the
   true eigenvalues of the matrix, are plotted in the complex plane
   to illustrate the Gershgorin Disc Theorem.

   References:
   https://en.wikipedia.org/wiki/Diagonally_dominant_matrix
*/

/*********************************/

%macro MakeGershgorin;
data GershgorinDiscs ;
retain Min 1e308 Max -1e308;
set CR;
Min = min(Min, xC - Radius, yC - Radius);
Max = max(Max, xC + Radius, yC + Radius);
run;

data All;
set GershgorinDiscs Evals;
run;

proc sgplot data=All aspect=1 noautolegend;
ellipseparm semimajor=Radius semiminor=Radius / fill outline
     xorigin=xC yorigin=yC transparency=0.5 name="c";
scatter x=xLambda y=yLambda / markerattrs=(symbol=X size=12) name="ev" legendlabel="Eigenvalues";
scatter x=Min y=Min / markerattrs=(size=0);
scatter x=Max y=Max / markerattrs=(size=0);
refline 0 / axis=x;      refline 0 / axis=y;
xaxis grid label="Real"; yaxis grid label="Imaginary";
keylegend "ev" / opaque location=inside position=topright across=1;
run;
%mend;

title "Gershgorin Discs";
proc iml;
A = { 200  30 -15  5,
       30 100   5  5,
      -15   5  55  0, 
        5   5   0 15};

evals = eigval(A);                 /* compute the eigenvalues */
center = vecdiag(A);               /* centers = diagonal elements */
radius = abs(A)[,+] - abs(center); /* sum of abs values of off-diagonal elements of each row */
discs = center || radius || round(evals,0.01);
print discs[c={"Center" "Radius" "Eigenvalue"} L="Gershgorin Discs"];

/* check to see if all eigenvalues must be positive */
if all(center >= radius) then 
   print "Diagonally dominant";
else 
   print "Not diagonally dominant";

xC = center; yC = 0#xC;          /* matrix is real, so centers are x + i*0 */
create CR var {"xC" "yC" "Radius"};  append;  close;

if ncol(evals)=1 then evals = evals || j(nrow(evals),1,0); /* make complex */
create Evals from evals[c={"xLambda" "yLambda"}];
append from evals; close;

submit;
title2 "Diagonally Dominant Matrix";
%MakeGershgorin
endsubmit;

/***************/

/* return centers and radii of Gershgorin theorem. Input is real square matrix. */
start GershgorinCentersRadii(A);
   c = vecdiag(A);          /* centers = diagonal elements */
   R = abs(A)[,+] - abs(c); /* sum of abs values of off-diagonal elements of each row */
   return c || R;
finish;

use Sashelp.Cars;
read all var _NUM_ into Y[c=varNames];
close;
A = corr(Y);

CR = GershgorinCentersRadii(A);
center = CR[,1]; radius = CR[,2];
xC = center; yC = 0#xC;          /* matrix is real, so centers are x + i*0 */
create CR var {"xC" "yC" "Radius"};  append;  close;

evals = eigval(A);
if ncol(evals)=1 then evals = evals || j(nrow(evals),1,0); /* make complex */
create Evals from evals[c={"xLambda" "yLambda"}];
append from evals; close;

submit;
title2 "Correlation Matrix";
%MakeGershgorin
endsubmit;

/***************/
/* example from IML doc */
A = {-1  2  0       0       0       0       0  0,
     -2 -1  0       0       0       0       0  0,
      0  0  0.2379  0.5145  0.1201  0.1275  0  0,
      0  0  0.1943  0.4954  0.1230  0.1873  0  0,
      0  0  0.1827  0.4955  0.1350  0.1868  0  0,
      0  0  0.1084  0.4218  0.1045  0.3653  0  0,
      0  0  0       0       0       0       2  2,
      0  0  0       0       0       0      -2  0 };
CR = GershgorinCentersRadii(A);
center = CR[,1]; radius = CR[,2];
xC = center; yC = 0#xC;          /* matrix is real, so centers are x + i*0 */
create CR var {"xC" "yC" "Radius"};  append;  close;

evals = eigval(A);
print evals;

if ncol(evals)=1 then evals = evals || j(nrow(evals),1,0); /* make complex */
create Evals from evals[c={"xLambda" "yLambda"}];
append from evals; close;

submit;
title2 "Unsymmetric Matrix with Complex Eigenvalues";
%MakeGershgorin
endsubmit;

/******************/
/* If A is singular, what is a small parameter d that we can add to 
   A so that (A + d*I) is nonsingular? 
   Or if A is not positive definite, what is a small parameter so that
   (A + d*I) is positive definite?
   The parameter d is called a "ridge parameter." Gershgorin's theorem
   provides a bound for d.
*/

