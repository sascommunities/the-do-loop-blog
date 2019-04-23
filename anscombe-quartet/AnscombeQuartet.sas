/* SAS program to accompany the article 
   "Create your own version of Anscombe's quartet: 
          Dissimilar data that have similar statistics"
   by Rick Wicklin, published 17APR2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/04/17/create-version-of-anscombes-quartet.html

   This program shows how to create data that are similar to 
   Anscombe's quartet
   https://en.wikipedia.org/wiki/Anscombe%27s_quartet

   Create 
   1. Linear data with normal errors
   2. Quadratic data, no errors
   The remaining Anscombe data sets are
   3. Linear with one extreme outlier
   4. Constant X, varying Y, except for one outlier
*/

/* THE (original) Anscombe's quartet data:
   For all pairs (x, y1), (x, y2), (x, y3), and (x4, y4)
   the correlation between is 0.816 
   and the linear regression line is y = 3 + 0.5x
*/
data Anscombe;
input x y1 y2 y3 x4 y4;
datalines;
10.0    8.04    9.14    7.46    8.0     6.58
8.0     6.95    8.14    6.77    8.0     5.76
13.0    7.58    8.74    12.74   8.0     7.71
9.0     8.81    8.77    7.11    8.0     8.84
11.0    8.33    9.26    7.81    8.0     8.47
14.0    9.96    8.10    8.84    8.0     7.04
6.0     7.24    6.13    6.08    8.0     5.25
4.0     4.26    3.10    5.39    19.0    12.50
12.0    10.84   9.13    8.15    8.0     5.56
7.0     4.82    7.26    6.42    8.0     7.91
5.0     5.68    4.74    5.73    8.0     6.89
;

ods graphics / width=300px height=200px;
ods layout gridded columns=2 advance=table;
title "Anscombe's Quartet: Linear Data Set";
proc sgplot data=Anscombe noautolegend;
   scatter x=x y=y1;
   lineparm x=4 y=5 slope=0.5;
run;
title "Anscombe's Quartet: Quadratic Data Set";
proc sgplot data=Anscombe noautolegend;
   scatter x=x y=y2;
   lineparm x=4 y=5 slope=0.5;
run;
ods layout end;


proc iml;
/* 1. Create first data set (randomly) */
call randseed(12345);
x = T( do(4, 14, 0.2) );                              /* evenly spaced X */
eps = round( randfun(nrow(x), "Normal", 0, 1), 0.01); /* normal error */
y = 3 + 0.5*x + eps;                                  /* linear Y + error */
 
/* Helper function. Return paremater estimates for linear regression. Args are col vectors */
start LinearReg(Y, tX);
   X = j(nrow(tX), 1, 1) || tX;
   b = solve(X`*X, X`*Y);       /* solve normal equation */
   return b;
finish;
 
/* 2. compute target statistics. Other data sets have to match these */
targetB = LinearReg(y, x);          /* compute regression estimates */
targetCorr = corr(y||x)[2];         /* compute sample correlation */
print (targetB`||targetCorr)[c={'b0' 'b1' 'corr'} F=5.3 L="Target"];


/* Define system of simultaneous equations:
   https://blogs.sas.com/content/iml/2018/02/28/solve-system-nonlinear-equations-sas.html */
/* This function returns linear regression estimates (b0, b1) and correlation for a choice of beta */
start LinearFitCorr(beta) global(x);
   y2 = beta[1] + beta[2]*x + beta[3]*x##2;    /* construct quadratic Y */
   b = LinearReg(y2, x);      /* linear fit */
   corr = corr(y2||x)[2];     /* sample corr */
   return ( b` || corr);      /* return row vector */
finish;
 
/* This function returns the vector quantity (beta - target). 
   Find value that minimizes Sum | F_i(beta)-Target+i |^2 */
start Func(beta) global(targetB, targetCorr);
   target = rowvec(targetB) || targetCorr;
   G = LinearFitCorr(beta) - target;
   return( G );              /* return row vector */
finish;
 
/* 3. Let's solve for quadratic parameters so that same 
   linear fit and correlation occur */
beta0 = {-5 1 -0.1};         /* initial guess */
con = {.  .  .,              /* constraint matrix */
       0  .  0};             /* quadratic term is negative */
optn = ncol(beta0) || 0;     /* LS with 3 components || amount of printing */
/* minimize sum( beta[i] - target[i])**2 */
call nlphqn(rc, beta, "Func", beta0, optn) blc=con;  /* LS solution */
print beta[L="Optimal beta"];
 
/* How nearly does the solution solve the problem? Did we match the target values? */
Y2Stats = LinearFitCorr(beta);
print Y2Stats[c={'b0' 'b1' 'corr'} F=5.3];


/* 4. Visualize the new versions of the first two data set
      in Anscombe's quartet */
y1 = y;
y2 = beta[1] + beta[2]*x + beta[3]*x##2;
create Anscombe2 var {x y1 y2};
append;
close;
QUIT;

proc print data=Anscombe2 noobs;run;

ods layout gridded columns=2 advance=table;
proc sgplot data=Anscombe2 noautolegend;
   scatter x=X y=y1;
   lineparm x=0 y=3.561 slope=0.447 / clip;
run;
proc sgplot data=Anscombe2 noautolegend;
   scatter x=x y=y2;
   lineparm x=0 y=3.561 slope=0.447 / clip;
run;
ods layout end;


/**********************************************************/
/* PROC AUTOREG models a time series with autocorrelation */
/**********************************************************/

/* SAS program to accompany the article 
   "A CUSUM test for autregressive models"
   by Rick Wicklin, published 24APR2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/04/24/cusum-test-autregressive-models.html

   This program shows how to use PROC AUTOREG in SAS to generate CUSUM
   charts for recursive residuals. You can use the CUSUM charts as a
   diagnostic goodness-of-fit tool for a misspeccified model.

   Data from "Create your own version of Anscombe's quartet: 
   Dissimilar data that have similar statistics" by Rick Wicklin
   https://blogs.sas.com/content/iml/2019/04/17/create-version-of-anscombes-quartet.html
*/
data Anscombe2;
input X Y1 Y2 ;
datalines;
4.0 5.26 3.42619 
4.2 6.17 3.74636 
4.4 6.02 4.05710 
4.6 4.75 4.35843 
4.8 6.94 4.65034 
5.0 4.27 4.93282 
5.2 5.46 5.20589 
5.4 6.74 5.46955 
5.6 5.87 5.72378 
5.8 7.13 5.96859 
6.0 5.85 6.20399 
6.2 6.39 6.42996 
6.4 5.05 6.64652 
6.6 5.82 6.85366 
6.8 5.95 7.05138 
7.0 6.77 7.23968 
7.2 6.38 7.41857 
7.4 6.92 7.58803 
7.6 6.38 7.74808 
7.8 7.15 7.89870 
8.0 7.91 8.03991 
8.2 8.58 8.17170 
8.4 8.86 8.29407 
8.6 5.62 8.40702 
8.8 6.64 8.51056 
9.0 8.13 8.60467 
9.2 9.03 8.68937 
9.4 8.06 8.76465 
9.6 7.78 8.83050 
9.8 8.37 8.88694 
10.0 8.74 8.93397 
10.2 7.54 8.97157 
10.4 9.57 8.99975 
10.6 8.38 9.01852 
10.8 7.09 9.02786 
11.0 6.57 9.02779 
11.2 9.20 9.01830 
11.4 9.52 8.99939 
11.6 9.85 8.97106 
11.8 7.90 8.93332 
12.0 10.08 8.88615 
12.2 8.60 8.82957 
12.4 8.04 8.76356 
12.6 8.84 8.68814 
12.8 10.34 8.60330 
13.0 7.53 8.50904 
13.2 9.17 8.40536 
13.4 9.75 8.29227 
13.6 8.76 8.16975 
13.8 10.69 8.03782 
14.0 10.39 7.89647 
;

/**********************************************************/
/* PROC AUTOREG models a time series with autocorrelation */
/**********************************************************/

ods graphics on;
proc autoreg data=Anscombe2;
  Linear: model y1 = x;     /* Y1 is linear. Model is oorrectly specified. */
  output out=CusumLinear cusum=cusum cusumub=upper cusumlb=lower 
                         r=Residual recres=RecursiveResid;
run;
proc autoreg data=Anscombe2;
  Quadratic: model y2 = x;     /* Y1 is quadratic. Model is misspecified. */
  output out=CusumQuad cusum=cusum cusumub=upper cusumlb=lower 
                       r=Residual recres=RecursiveResid;
run;

title "Cumulative Recursive Residuals"
title2 "Linear Data (left) and Quadratic Data (right)";

ods graphics / width=400px height=300px;
ods layout gridded columns=2 advance=table;
 proc sgplot data=CusumLinear noautolegend;
    band x=x lower=lower upper=upper;
    series x=x y=cusum / break markers;
    refline 0  /axis=y noclip;
    xaxis grid; yaxis grid;
 run;
 proc sgplot data=CusumQuad noautolegend;
    band x=x lower=lower upper=upper;
    series x=x y=cusum / break markers;
    refline 0  /axis=y noclip;
    xaxis grid; yaxis grid;
 run;
ods layout end;

/* CUSUM test of the linearity assumption */
proc iml;
start seriesBounds(x,y, YBounds);
    create _temp var {'x' 'y'}; append; close;
    lower=yBounds[1]; upper=YBounds[2];
    submit lower upper; 
       proc sgplot data=_temp noautolegend;
          band x=x lower=&lower upper=&upper;
          series x=x y=y / break markers;
          refline 0 &lower &upper /axis=y noclip;
          xaxis grid; yaxis grid;
       run;
    endsubmit;
finish;
start scatterLineparm(x,y, beta);
    create _temp var {'x' 'y'}; append; close;
    intercept=beta[1]; slope=beta[2];
    submit intercept slope; 
       proc sgplot data=_temp noautolegend;
          series x=x y=y / break markers;
          lineparm x=0 y=&intercept slope=&slope /clip;
          xaxis grid; yaxis grid;
       run;
    endsubmit;
finish;
store module=(seriesBounds scatterLineparm);

proc iml;
load module=(seriesBounds scatterLineparm);
use Anscombe; read all var {'x' 'y1' 'y2'}; close;
beta = {3, 0.5};
use Anscombe2; read all var {'x' 'y1' 'y2'}; close;
beta = {3.561, 0.447};

title "Linear Data";
call scatterLineParm(x, y1, beta);
title "Quadratic Data";
call scatterLineParm(x, y2, beta);

/* Passing-Bablok test for linearity in 
   method-comparison studies: p 712-713 */
start CusumTestForLinearity(x, y, beta);
   resid = y - (beta[1] + beta[2]*x);
   bAbove = (resid > 0); bBelow = (resid < 0);
   NAbove = sum(bAbove);
   NBelow = sum(bBelow);
   r = j(nrow(resid), 1, 0);
   r[loc(bAbove)] =  sqrt(NBelow / NAbove);
   r[loc(bBelow)] = -sqrt(NAbove / NBelow);
 
   D = (y + x/beta[2] - beta[1]) / sqrt(1 + 1/beta[2]**2);
   rank = rank(D);   
   /* order by rank 
      https://blogs.sas.com/content/iml/2011/04/06/how-to-rank-values.html 
   */ 
   tmp = r;        /* copy values */
   r[rank] = tmp; /* re-order */
   c = 0 // cusum(r);

   KS_crit = 1.36; /* critical value of KS distrib for alpha=0.05 */
   Bound = KS_crit * sqrt(NBelow + 1);
   RejectH0 = ( max(abs(c)) > bound ); 
   if RejectH0 then 
        print "Linear relationship rejected at 0.05 signifiance.";
   else print "Linear relationship cannot be rejected at 0.05 signifiance.";
   print (max(abs(c)))[L="Max Cusum"] Bound;
   test = [#'MaxCusum'=max(abs(c)), #'Bound'=Bound, #'Cusum'=c];
   return (test);
finish;

Test = CusumTestForLinearity(x, y1, beta);
c = Test$'Cusum';
bound = Test$'Bound';
call seriesBounds(0:nrow(x), c, bound // -bound);

Test = CusumTestForLinearity(x, y2, beta);
c = Test$'Cusum';
bound = Test$'Bound';
call seriesBounds(0:nrow(x), c, bound // -bound);




/* translate and rotate coordinates */
y2 = y - beta[1];
/* clockwise rotation by NEGATIVE theta */
Rot = ( 1       || beta[2]) //
      (-beta[2] || 1);
Rot = Rot / sqrt(1+beta[2]**2);  /* rotation by negative theta where tan(theta)=b */

uv = (x||y2) * Rot`;
call scatter(uv[,1], uv[,2]) datalabel=rank grid={x y}; 
print x y y2 uv;



/********************/
/* CUSUM test for binary {-1, +1} sequences */
proc iml;
eps = {1 1 0 0 1 0 0 1 0 0 0 0 1 1 1 1 1 1 0 1
       1 0 1 0 1 0 1 0 0 0 1 0 0 0 1 0 0 0 0 1 
       0 1 1 0 1 0 0 0 1 1 0 0 0 0 1 0 0 0 1 1 
       0 1 0 0 1 1 0 0 0 1 0 0 1 1 0 0 0 1 1 0 
       0 1 1 0 0 0 1 0 1 0 0 0 1 0 1 1 1 0 0 0 }; 
x = 2*eps - 1;            /* convert to {-1, +1} sequence */

call tabulate(values, freq, eps);
print freq[c=(char(values))];
p = 2*cdf("Binomial", min(freq), 0.5, 100);
print p;

S = cusum(x);
title "Cumulative Sums of Sequence of {-1, +1} Values";
call series(1:ncol(S), S) option="markers" other="refline 0 / axis=y" label="Observation Number";

/* NIST test for randomness in binary sequence:
   https://nvlpubs.nist.gov/nistpubs/legacy/sp/nistspecialpublication800-22r1a.pdf
   Section 2.13.  Pages 2-31 through 2-33.
   CUSUM test for binary {-1, +1} sequence.
   INPUT: x is sequence of {-1, +1} values.
*/
start BinaryCUSUMTest(x, PrintTable=1, alpha=0.01);
   S = colvec( cusum(x) );     /* cumulative sums */
   n = nrow(S);
   z = max(abs(S));            /* test statistic = maximum deviation */

   /* compute probability of this test statistic for a sequence of this length */
   zn = z/sqrt(n);
   kStart = int( (-n/z +1)/4 );
   kEnd   = int( ( n/z -1)/4 );
   k = kStart:kEnd;
   sum1 = sum( cdf("Normal", (4*k+1)*zn) - cdf("Normal", (4*k-1)*zn) );

   kStart = int( (-n/z -3)/4 );
   k = kStart:kEnd;
   sum2 = sum( cdf("Normal", (4*k+3)*zn) - cdf("Normal", (4*k+1)*zn) );
   pValue = 1 - sum1 + sum2;

   /* optional: print the test results in a nice human-readable format */
   cusumTest = z || pValue;
   if PrintTable then do;
      print cusumTest[L="Result of CUSUM Test" c={"Test Statistic" "p Value"}];
      labl= "H0: Sequence is a random binary sequence";
      if pValue <= alpha then 
         msg = "Reject H0 at alpha = " + char(alpha); /* sequence seems random */
      else
         msg = "Do not reject H0 at alpha = " + char(alpha); /* sequence does not seem random */
      print msg[L=labl];
   end;
   return ( cusumTest );
finish;

/* call the function for the {-1, +1} sequence */
cusumTest = BinaryCUSUMTest(x);


eps = {0 1 0 0 1 0 0 1 0 0 0 0 1 1 0 0 0 1 0 1
       0 1 0 0 1 1 0 0 0 1 0 0 1 1 0 0 0 1 1 0 
       0 1 1 0 1 0 0 0 1 1 0 0 0 0 1 0 0 0 1 1 
       1 0 0 0 1 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 
       0 1 1 0 0 0 1 0 1 0 0 0 1 0 0 1 1 0 0 0 }; 
x = 2*eps - 1;

eps = mod(1:100, 2);

eps = j(100, 1); eps[1:10]=0; eps[31:30]=0; eps[41:50]=0; 
eps[61:70]=0; eps[81:90]=0; 
x = 2*eps - 1;
cusumTest = BinaryCUSUMTest(x);

eps = j(1, 29,1) || j(1, 29,0) || mod(1:42, 2);
x = 2*eps - 1;
cusumTest = BinaryCUSUMTest(x);

/* find critical max deviation for sequence of length N=100 (28 for alpha=0.01) */
zz = 1:50;
pValue = 0*zz;
n = 100;
do i = 1 to ncol(zz);
   z = zz[i];
   /* compute probability of this test statistic for a sequence of this length */
   zn = z/sqrt(n);
   kStart = int( (-n/z +1)/4 );
   kEnd   = int( ( n/z -1)/4 );
   k = kStart:kEnd;
   sum1 = sum( cdf("Normal", (4*k+1)*zn) - cdf("Normal", (4*k-1)*zn) );

   kStart = int( (-n/z -3)/4 );
   k = kStart:kEnd;
   sum2 = sum( cdf("Normal", (4*k+3)*zn) - cdf("Normal", (4*k+1)*zn) );
   pValue[i] = 1 - sum1 + sum2;
end;
print (zz // pValue);

QUIT;

