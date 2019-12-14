/* SAS program to accompany the article 
   "Math-ing around the Christmas tree: Can the SVD de-noise an image?"
   by Rick Wicklin, published 16DEC2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/12/16/svd-de-noise-image.html

   This program creates a binary matrix where the 1's are in a triangular region 
   that looks like a Christmas tree. It then "adds noise" by 
   switching 10% of the values in the matrix.

   The matrix is factored by using the SVD. A rank-4 approximation is used
   to estimate the noisy matrix. When you use a threashold value, you can convert the 
   low-rank approximation to a binary matrix that estimates the original (pre-noise) 
   image. For this experiment, 98% of the cells are recovered when the threshold
   value is 0.5.
*/

ods graphics / width=300px height=400px;
proc iml;
/* Create the Christmas tree matrix by using 
   https://blogs.sas.com/content/iml/2013/12/18/christmas-tree-matrix.html */

/* parameters for the tree dimensions: 101 x 100 */
h = 100; w = h+1; b = int(h/10);
M = j(w, h, 0);         /* initialize h x w matrix to missing */
x = col(M);             /* column numbers */
y = w + 1 - row(M);     /* reverse Y axis */
 
/* define the leafy portion of the tree */
TreeIdx = loc(y>b & y<=2*x & y<=-2*x+2*h); /* a triangle */
M[TreeIdx] = 1;
 
/* define the trunk */
width = int(b/2);
TrunkIdx = loc( y<= b & abs(x-nrow(M)/2)<width );
M[TrunkIdx] = 1;

/* save the original (true) tree matrix */
Tree = M;

/* add random noise */
call randseed(12345);
PctNoise = 0.1;
idx = sample(1:w*h, 0.2*w*h, "noReplace");
M[idx] = ^M[idx];

ramp = {WHITE GREEN};
s = "Christmas Tree with " + strip(putn(PctNoise, 'PERCENT7.')) + " Noise";
call HeatmapDisc(M) colorramp=ramp displayoutlines=0 showlegend=0 title=s;

/******************************/

/* Factor by using SVD:
   https://blogs.sas.com/content/iml/2017/08/30/svd-low-rank-approximation.html 
   Plot largest singular values */
call svd(U, D, V, M);  /* A = U*diag(D)*V` */
title "Singular Values";
ods graphics / width=400px height=300px;
call series(1:20, D) grid={x y} xvalues=1:nrow(D) label={"Component" "Singular Value"};
ods graphics / width=300px height=400px;

/* based on graph, use rank-4 approximation of matrix */
keep = 4;        /* how many components to keep? */
idx = 1:keep;
Ak = U[ ,idx] * diag(D[idx]) * V[ ,idx]`;

Ak = (Ak - min(Ak)) / range(Ak); /* for plotting, standardize into [0,1] */
s = "Rank-" + strip(char(keep)) + " Approximation of Noisy Christmas Tree Image";
call heatmapcont(Ak) colorramp=ramp displayoutlines=0 showlegend=0 title=s range={0, 1};

/* what is the distribution of cell values in low-rank approximation? */
ods graphics / width=400px height=300px;
call histogram(Ak);

/* use threshold parameter to convert low-rank approx into binary matrix */
ods graphics / width=300px height=400px;
t = 0.5;  /* threshold parameter */
s = "Discrete Filter: Threshold = " + strip(char(t,4,2)) ;
call HeatmapDisc(Ak > t) colorramp=ramp displayoutlines=0 showlegend=0 title=s;

/* explore other cutoff values */
ods layout gridded columns=2 advance=table;
do t = 0.4 to 0.6 by 0.2;
   s = "Discrete Filter: Threshold = " + strip(char(t,4,2)) ;
   call HeatmapDisc(Ak > t) colorramp=ramp displayoutlines=0 showlegend=0 title=s;
end;
ods layout end;

/* what proportion of the original cells are the same as the new binary estimate? */
A4 = (Ak > 0.5);
propSame = 1 - (Tree - A4)[:];
print propDiff;

/* Create ROC curve for the rank-4 approx:
   https://blogs.sas.com/content/iml/2018/11/14/compare-roc-curves-sas.html */

create Score var {'Tree' 'Ak'};
append;
close;
QUIT;

ods graphics / reset;
proc logistic data=Score;
model Tree(event='1') = Ak / nofit;
roc 'Expert Predictions' pred=Ak;
ods select ROCcurve;
run;
