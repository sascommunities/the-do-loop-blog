/* Motivating Example: scatter ellipses and PCA 
   https://blogs.sas.com/content/iml/2014/07/23/prediction-ellipses-from-covariance.html
*/
/**********************************************/
/* STORE modules by running RandCorrMatrices.sas */
/**********************************************/
%include "RandCorrMatrices.sas";

%let rho = 0.6;

proc iml;
call randseed(1234);
/* 1. Generate correlated standardized data */
X = randnormal(1000, {0 0}, {1 &rho, &rho 1});
X = (X-mean(X)) / std(X);
S = cov(X);
call eigen(eval,evec, S); print S, eval evec;
call eigen(Spectrum, U, S);
V = -sqrt(Spectrum)` # U;  /* change sign; evecs not unique */
/* scale to match 95% confidence ellipse */
n = nrow(X);
c = 2 * (n-1)/n * (n+1)/(n-2);          /* adjust for finite sample size */
p = 0.95;                               /* confidence levels */
F = sqrt(c * quantile("F", p, 2, n-2)); /* p = 1-alpha */
V = F*V`; /* transpose for plotting */
labl = {'e1', 'e2'};
create Eigen from X V labl [c={'X1' 'X2' 'v1' 'v2' 'name'}];
append from X V labl;
close;
QUIT;

ods graphics / antialiasmax=50000;
ods graphics / width=480px height=480px;
title "Original Data and Eigenstructure";
title2 "Correlation ~ &rho";
proc sgplot data=Eigen aspect=1 noautolegend;
   lineparm x=0 y=0 slope=1;
   lineparm x=0 y=0 slope=-1;
   scatter x=X1 y=X2;
   ellipse x=X1 y=X2;
   vector x=v1 y=v2 / lineattrs=(thickness=3 color=darkred) 
          datalabel=name datalabelattrs=(size=16 color=darkred);
   xaxis grid min=-4 max=4;
   yaxis grid min=-4 max=4;
run;

/* what happens if we change the sign of X1? Mirror image reflection */
data eigen2;
set eigen;
X1 = -X1;
v1 = -v1;
run;

title "Flipped Data and Eigenstructure";
title2 "Correlation ~ -&rho";
proc sgplot data=Eigen2 aspect=1 noautolegend;
   lineparm x=0 y=0 slope=1;
   lineparm x=0 y=0 slope=-1;
   scatter x=X1 y=X2;
   ellipse x=X1 y=X2;
   vector x=v1 y=v2 / lineattrs=(thickness=3 color=darkred) 
          datalabel=name datalabelattrs=(size=16 color=darkred);
   xaxis grid min=-4 max=4;
   yaxis grid min=-4 max=4;
run;

/**************************/
/* first, visualize the 2-D case */
proc iml;
load module=_all_;   /* modules are defined in the Appendix */
call randseed(12345);

rho = 0.5;
Sigma = I(2); 
Sigma[{2 3}] = rho;
Spectrum = eigval(Sigma);
print Spectrum;
/* simulate many matrices and count how many have +rho vs -rho */
NSim = 10000;
result = j(NSim, 1);
do i = 1 to NSim;
   R = CorrWithEigen(Spectrum);
   result[i] = R[2];
end;
result = round(result, 1E-10);  /* fuzz numerical results like 0.5 vs 0.4999999 */
call tabulate(Values, Freq,  result);
print Freq[c=Values];
QUIT;

/****************************/

/* Visualize the off-diagonal elements for a 3x3 correlation 
   matrix that has the same spectrum as an AR1(rho) matrix */
proc iml;
/* modules are defined in the Appendix at 
   https://blogs.sas.com/content/iml/2024/12/18/correlation-matrix-eigenvalues.html
*/
load module=_all_;   
call randseed(12345);
/* AR1 correlation structure:
   https://blogs.sas.com/content/iml/2012/11/05/constructing-common-covariance-structures.html */
Sigma = {1.0  0.5  0.25,
         0.5  1.0  0.5,
         0.25 0.5  1.0};
Spectrum = eigval(Sigma);
print Spectrum;

NSim = 5000;          /* simulate NSim random matrices */
result = j(NSim, 3);
do i = 1 to NSim;
   R = CorrWithEigen(Spectrum);
   result[i,] = rowvec( R[{2 3 6}] ); /* save the off-diagonal values */
end;
/* write the values to a data set */
create RSim from result[c={R12 R13 R23}];
   append from result;
close;

R2 = ssq(Spectrum);
print R2;
QUIT;

/* visualize in 2-D */
data Orig;
O12 = 0.5;
O13 = 0.25;
O23 = 0.5;
run;

data RSim2;
set RSim Orig;
run;

title "Cell Values for 3x3 Matrices That Have a Common Spectrum";
proc sgplot data=RSim2 aspect=1;
   scatter x=R12 y=R13 / markerattrs=(symbol=CircleFilled size=3) colorresponse=R23 name="S1";
   scatter x=O12 y=O13 / markerattrs=(symbol=StarFilled color=black size=12);
   refline 0 /axis=x; 
   refline 0 /axis=y;
   xaxis grid;
   yaxis grid;
   gradlegend "S1";
run;

/****************************/
/* use PROC G3D to visualize in 3-D. First, generate some lines
   of longitude and latitude on the sphere with radius R.
   Use R^2 as the parameter to the subroutine.
*/

proc iml; 
/* create a sphere SSQ(x,y,z) = Rsq */
start MakeSphere(Rsq);
   pi = constant('pi');
   theta = T( do(0, 2*pi, pi/100) );
   equator = cos(theta) || sin(theta) || j(nrow(theta),1,0);
   latitude = {. . . };
   create sphere from latitude[c={'R12' 'R13' 'R23'}];
   radius = sqrt(Rsq);
   dh =  radius / 5;
   dhLow = dh*1.01;
   h = do( -dhLow, -radius, -dhLow) || do( 0, radius -dh, dh);
   do i = 1 to ncol(h);
      rad = sqrt( Rsq - h[i]##2 );
      *print i rad (h[i])[L="h[i]"];
      z = {0 0 } || h[i];
      latitude = rad*equator + z;
      append from latitude;
      latitude = latitude[ ,{3,2,1}];
      append from latitude;
   end;
   close;
finish;

Sigma = {1.0  0.5  0.25,
         0.5  1.0  0.5,
         0.25 0.5  1.0};
Spectrum = eigval(Sigma);
S2 = SSQ(Spectrum);
SumOffDiag2 = (S2 - ncol(Sigma)) / 2;
print Spectrum, S2 SumOffDiag2;
Rsq = SumOffDiag2; 
print Rsq (sqrt(RSq))[L="R"];
run MakeSphere( Rsq );
quit;


data All;
set RSim2(in=data) sphere;
if data then do;
   group=1; shape="Balloon"; color="blue"; size=1;
end;
else do;
   group=2; shape="points"; color="gray"; size=1;
end;
run;

goptions hsize=5in vsize=5in;
title "Cell Values That Have a Common Spectrum";
proc g3d data=All;
   scatter R13*R12=R23 / shape=shape color=color size=size noneedle grid;
run;
quit;


/********************************/
/* What are some other matrices on the sphere? */

/* what other matrices are on the sphere? */
proc iml;
pt = {0.75  0  0};

R12 = pt[1]; R13 = pt[2]; R23 = pt[3];
Sigma = (1.0 || R12 || R13) //
        (R12 || 1   || R23) //
        (R13 || R23 || 1  );
Spectrum = eigval(Sigma);
S2 = SSQ(Spectrum);
SumOffDiag2 = (S2 - ncol(Sigma)) / 2;
print pt, Spectrum, S2 SumOffDiag2;


pts = {0.75   0      0,
       0      0.75   0,
       0      0      0.75};
R12 = 0.4; R23 = 0.4;
R13 = sqrt(SumOffDiag2 - R12**2 - R23**2); 
pts = pts // 
      (R12 || R13 || R23);

do i = 1 to nrow(pts);
   pt = pts[i,];
   R12 = pt[1]; R13 = pt[2]; R23 = pt[3];
   Sigma = (1.0 || R12 || R13) //
           (R12 || 1   || R23) //
           (R13 || R23 || 1  );
   Spectrum = eigval(Sigma);
   S2 = SSQ(Spectrum);
   SumOffDiag2 = (S2 - ncol(Sigma)) / 2;
   print pt, Spectrum, S2 SumOffDiag2;
end;

S2 = SSQ(Spectrum);
N2 = norm(Sigma)**2;
print S2[L="SSQ(Spectrum)"] N2[L="||Sigma||^2"];
Rsq = (SSQ(Spectrum) - 3) / 2; 
print Rsq (sqrt(RSq))[L="R"];
QUIT;
