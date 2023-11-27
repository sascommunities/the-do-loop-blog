/* SAS program to accompany the articles
   "Bivariate normal probability in SAS: Rectangular regions"
    by Rick Wicklin, published 29NOV2023 on The DO Loop blog:
    https://blogs.sas.com/content/iml/2023/11/29/bivariate-normal-rectangle.html  

   "Bivariate normal probability in SAS"
   by Rick Wicklin, published 04DEC2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/12/04/bivariate-normal-probability-sas.html

   This program shows how to compute the probability of random points
   from the bivariate standard normal distribution with correlation rho.
   The distribution is denoted BVN(0, rho) for -1 < rho < 1.
*/


/**************************************************/
/* FIRST ARTICLE: INTRODUCE THE PROBBNRM FUNCTION */
/**************************************************/

/* show examples of using the PROBBNRM function */
/* prob in SW quadrant */
data ProbSW;
rho = 0.4;
a = -0.5; b =  0.5;  /* cut points in X direction */
c = -1.2; d =  0.8;  /* cut points in Y direction */
SWProb = probbnrm(a, c, rho);
run;

proc print data=ProbSW noobs; 
   var rho a c SWProb;
run;

/* prob in finite rectangle (a,b) x (c,d) */
data ProbRect;
rho = 0.4;
a = -0.5; b =  0.5;  /* cut points in X direction */
c = -1.2; d =  0.8;  /* cut points in Y direction */
bdProb = probbnrm(b, d, rho);   /* upper right corner */
adProb = probbnrm(a, d, rho);   /* upper left corner  */
bcProb = probbnrm(b, c, rho);   /* lower right corner */
acProb = probbnrm(a, c, rho);   /* lower left corner  */
rectProb = bdProb - adProb - bcProb + acProb;
run;

proc print data=ProbRect noobs; 
   var rho a b c d rectProb;
run;


/* The probability depends on the correlation, rho.
   Graph probability over rectangle as a function of rho */
data ProbRho;
a = -0.5; b =  0.5;  /* cut points in X direction */
c = -1.2; d =  0.8;  /* cut points in Y direction */
do rho = -0.95 to 0.95 by 0.01;
   bdProb = probbnrm(b, d, rho);   /* upper right corner */
   adProb = probbnrm(a, d, rho);   /* upper left corner  */
   bcProb = probbnrm(b, c, rho);   /* lower right corner */
   acProb = probbnrm(a, c, rho);   /* lower left corner  */
   rectProb = bdProb - adProb - bcProb + acProb;
   output;
end;
keep a b c d rho rectProb;
run;

title "Probability in a Rectangular Region";
title2 "(-0.5, 0.5) x (-1.2, 0.8)";
proc sgplot data=ProbRho;
   series x=rho y=rectProb;
   xaxis grid;
   yaxis grid;
   label rho="Correlation (rho)" rectProb="P(a<X<b,c<Y<d; rho)";
run;

/* Use a different region (0,1) x (0,1) 
   Graph probability over rectangle as a function of rho */
data ProbRho;
a = 0; b =  1;
c = 0; d =  1; 
do rho = -0.9999, -0.999, -0.99 to 0.99 by 0.01, 0.999, 0.9999;
   bdProb = probbnrm(b, d, rho);   /* upper right corner */
   adProb = probbnrm(a, d, rho);   /* upper left corner  */
   bcProb = probbnrm(b, c, rho);   /* lower right corner */
   acProb = probbnrm(a, c, rho);   /* lower left corner  */
   rectProb = bdProb - adProb - bcProb + acProb;
   output;
end;
keep a b c d rho rectProb;
run;

/* You can compute the upper bound on the probability, which is the 
   limit as rho -> 1 */
data _null_;
prob = CDF("normal",1) - CDF("Normal",0);
put prob=;
run;

title "Probability in a Rectangular Region";
title2 "(0, 1) x (0, 1)";
proc sgplot data=ProbRho;
   series x=rho y=rectProb;
   refline 0.341 / axis=y noclip;
   xaxis grid;
   yaxis grid;
   label rho="Correlation (rho)" rectProb="P(a<X<b,c<Y<d; rho)";
run;

/* NOTE: There is a symmetry (x,y) --> (y,x). 
   There is also a symmetry when rho --> -rho */
data ProbTest;
a = -0.5; b =  0.5;  /* cut points in X direction */
c = -1.2; d =  0.8;  /* cut points in Y direction */
do rho = 0 to 0.95 by 0.01;
   bdProb = probbnrm(b, d, rho);   /* upper right corner */
   adProb = probbnrm(a, d, rho);   /* upper left corner  */
   bcProb = probbnrm(b, c, rho);   /* lower right corner */
   acProb = probbnrm(a, c, rho);   /* lower left corner  */
   rhoProb = bdProb - adProb - bcProb + acProb;

   bdProb = probbnrm(b, d, -rho);   /* upper right corner */
   adProb = probbnrm(a, d, -rho);   /* upper left corner  */
   bcProb = probbnrm(b, c, -rho);   /* lower right corner */
   acProb = probbnrm(a, c, -rho);   /* lower left corner  */
   negrhoProb = bdProb - adProb - bcProb + acProb;

   Diff = rhoProb - negrhoProb;
   output;
end;
keep rho rhoProb negrhoProb diff;
run;


/***********************************************/
/* SECOND ARTICLE: COMPUTE BVN PROBABILITY OVER 
   ANY REGION (a,b) x (c,d), where we allow 
   a = -Infinity, b = +Infinity, c = -Infinity, 
   and d = +Infinity 
   If we use the missing values .M and .P to represent 
   -Infinity and +Infinity, respectively, then there 
   are 16 possible regions to consider.
*/
/***********************************************/

/* DEFINE AND STORE THE CDFBN and PROBBVN FUNCTIONS */
proc iml;
/* Extend the bivariate CDF, which is PROBBNRM(a,b, rho) to support an 
   upper limit of infinity (.P) for either argument:
   Pr(u,.P; rho) = Phi(u) = cdf("Normal", u) is probability over left half plane
   Pr(.P,v; rho) = Phi(v) = cdf("Normal", v) is probability over lower half plane
*/
start CDFBN(u,v,rho);
   if missing(u) & missing(v) then 
      return 1;
   if missing(u) then 
      return cdf("Normal", v);
   if missing(v) then 
      return cdf("Normal", u);
   return probbnrm(u, v, rho);
finish;

start ProbBVN(a,b,c,d,rho);
   ma = missing(a);   mb = missing(b);
   mc = missing(c);   md = missing(d);
   /* 1. complete plane */
   if ma & mb & mc & md then
      return 1;
   /* 2. lower half plane */
   if ma & mb & mc & ^md then
      return CDFBN(.P,d,rho);
   /* 3. upper half plane */
   if ma & mb & ^mc & md then
      return 1 - CDFBN(.P,c,rho);
   /* 4. horiz strip */
   if ma & mb & ^mc & ^md then
      return CDFBN(.P,d,rho) - CDFBN(.P,c,rho);
   /* 5. left half plane */
   if ma & ^mb & mc & md then 
      return CDFBN(b,.P,rho); 
   /* 6. SW quadrant */
   if ma & ^mb & mc & ^md then
      return CDFBN(b,d,rho);
   /* 7. NW quadrant */
   if ma & ^mb & ^mc & md then 
      return CDFBN(b,.P,rho) - CDFBN(b,c,rho); 
   /* 8. left strip (W) */
   if ma & ^mb & ^mc & ^md then
      return CDFBN(b,d,rho) - CDFBN(b,c,rho);
   /* 9. right half plane */
   if ^ma & mb & mc & md then
      return 1 - CDFBN(.P,a,rho);
   /* 10. SE quadrant = lower half - (SW quad)*/
   if ^ma & mb & mc & ^md then
      return CDFBN(.P,d,rho) - CDFBN(a,d,rho);
   /* 11. NE quadrant = right half - (SE quad) */
   if ^ma & mb & ^mc & md then 
      return 1 - CDFBN(.P,a,rho)
               - (CDFBN(.P,c,rho) - CDFBN(a,c,rho));
   /* 12. right strip (E) */
   if ^ma & mb & ^mc & ^md then
      return CDFBN(.P,d,rho) - CDFBN(.P,c,rho) 
           - CDFBN(a,d,rho) + CDFBN(a,c,rho);
   /* 13. vert strip */
   if ^ma & ^mb & mc & md then
      return CDFBN(.P,b,rho) - CDFBN(.P,a,rho);
   /* 14. lower strip (S) */
   if ^ma & ^mb & mc & ^md then
      return CDFBN(b,d,rho) - CDFBN(a,d,rho);
   /* 15. upper strip (N) */
   if ^ma & ^mb & ^mc & md then
      return CDFBN(b,.P,rho) - CDFBN(a,.P,rho) /*upper strip _|_  */
           - CDFBN(b,c,rho) + CDFBN(a,c,rho);
   /* 16. rectangle */
   if ^ma & ^mb & ^mc & ^md then
      return CDFBN(b,d,rho) - CDFBN(a,d,rho) 
           - CDFBN(b,c,rho) + CDFBN(a,c,rho); 
   return( . );  /* should never execute this statement */
finish;

/* what is the correct probability for a given rho? 
   To validate the ProbBVN function, use Monte Carlo simulation.
   Assume X ~ BVN(0, Sigma) where Sigma = {1 rho, rho 1} */
start BVNMC(a, b, c, d, XMat);
   ma = missing(a);   mb = missing(b);
   mc = missing(c);   md = missing(d);
   n = nrow(XMat);
   x = XMat[,1]; y = XMat[,2];

   /* 1. Complete plane */
   if ma & mb & mc & md then
      return 1;
   /* 2. lower half plane */
   if ma & mb & mc & ^md then
      return sum(y < d) / n;
   /* 3. upper half plane */
   if ma & mb & ^mc & md then
      return sum(c < y) / n;
   /* 4. horiz strip */
   if ma & mb & ^mc & ^md then
      return sum(c < y & y < d) / n;
   /* 5. left half plane */
   if ma & ^mb & mc & md then 
      return sum(x < b) / n;
   /* 6. SW quadrant */
   if ma & ^mb & mc & ^md then
      return sum(x < b & y < d) / n;
   /* 7. NW quadrant */
   if ma & ^mb & ^mc & md then 
      return sum(x < b & c < y) / n;
   /* 8. Horiz Left Strip */
   if ma & ^mb & ^mc & ^md then
      return sum(x < b & c < y & y < d) / n;
   /* 9. right half plane */
   if ^ma & mb & mc & md then
      return sum(a < x) / n;
   /* 10. SE quadrant = lower half - (SW Quad)*/
   if ^ma & mb & mc & ^md then
      return sum(a < x & y < d) / n;
   /* 11. NE quadrant = right half - (SE Quad) */
   if ^ma & mb & ^mc & md then 
      return sum(a < x & c < y) / n;
   /* 12. Horiz Right Strip */
   if ^ma & mb & ^mc & ^md then
      return sum(a < x & c < y & y < d) / n;
   /* 13. vert strip */
   if ^ma & ^mb & mc & md then
      return sum(a < x & x < b) / n;
   /* 14. vert lower strip */
   if ^ma & ^mb & mc & ^md then
      return sum(a < x & x < b & y < d) / n;
   /* 15. vert upper strip */
   if ^ma & ^mb & ^mc & md then
      return sum(a < x & x < b & c < y) / n;
   /* 16. Rectangle */
   if ^ma & ^mb & ^mc & ^md then
      return sum(a < x & x < b & c < y & y < d) / n;
   return( . );  /* should never execute this statement */
finish;

store module=(CDFBN ProbBVN BVNMC);
quit;


/* ------------------------------------ */
/* TEST THE CDFBN and ProbBVN functions */
proc iml;
load module=(CDFBN);
rho = 0.4;
x0 = 0.5; y0 = 0.8;

/* compute prob on half planes */
probLeft  = CDFBN(x0,.P, rho); /* left half plane {(x,y) | x < x0} */
probRight = 1 - probLeft;      /* right half plane */
probLower = CDFBN(.P,y0, rho); /* lower half plane {(x,y) | y < y0} */
probUpper = 1 - probLower;     /* upper half plane */
print probLeft probRight probLower probUpper;

/* use half planes to compute prob on quadrants */
P1 = CDFBN(x0, y0, rho);        /* probability SW from (x0,y0) */
P2 = CDFBN(.P, y0, rho) - P1;   /* probability SE from (x0,y0) */
P3 = CDFBN(x0,.P, rho) - P1;    /* probability NW from (x0,y0) */
P4 = 1 - (P1 + P2 + P3);        /* probability NE from (x0,y0) */
print P1 P2 P3 P4;

/* compute prob on strips */
a = -0.5; b =  0.5;
c = -1.2; d =  0.8;
probW = CDFBN(a,d,rho) - CDFBN(a,c,rho);  /* left strip  ---] */
probE = CDFBN(.P,d,rho) - CDFBN(.P,c,rho) /* right strip [--- */ 
      - CDFBN(b,d,rho) + CDFBN(b,c,rho); 
probS = CDFBN(b,c,rho) - CDFBN(a,c,rho);  /* lower strip  T  */
probN = CDFBN(b,.P,rho) - CDFBN(a,.P,rho) /*upper strip _|_  */
      - CDFBN(b,c,rho) + CDFBN(a,c,rho); 
print probW probE probS probE;
QUIT;

/* You can use the ProbBVN function to evaluate the BVN probability
   on any of the 9 common shapes */
proc iml;
load module=(CDFBN ProbBVN);   /* the CDFBN and ProbBVN functions must be stored by using the STORE statement */

rho = 0.4;
a = -0.5; b =  0.5;
c = -1.2; d =  0.8;
Prob = j(9, 1, .);
Prob[1] = ProbBVN(.M,  a, .M,  c, rho);   /* SW */
Prob[2] = ProbBVN( a,  b, .M,  c, rho);   /* S  */
Prob[3] = ProbBVN( b, .P, .M,  c, rho);   /* SE */
Prob[4] = ProbBVN(.M,  a,  c,  d, rho);   /* W  */
Prob[5] = ProbBVN( a,  b,  c,  d, rho);   /* C  */
Prob[6] = ProbBVN( b, .P,  c,  d, rho);   /* E  */
Prob[7] = ProbBVN(.M,  a,  d, .P, rho);   /* NW */
Prob[8] = ProbBVN( a,  b,  d, .P, rho);   /* N  */
Prob[9] = ProbBVN( b, .P,  d, .P, rho);   /* NE */

desc = {'SW', 'S', 'SE', 'W ', 'C', 'W ', 'NW', 'N', 'NE'};
labl={'x1'  'x2'  'y1'  'y2'};
Pt = (.M || a  || .M ||  c) //   /* SW */
     ( a || b  || .M ||  c) //   /* S  */
     ( b || .P || .M ||  c) //   /* SE */
     (.M || a  ||  c ||  d) //   /* W  */
     ( a || b  ||  c ||  d) //   /* C  */
     ( b || .P ||  c ||  d) //   /* E  */
     (.M || a  ||  d || .P) //   /* NW */
     ( a || b  ||  d || .P) //   /* N  */
     ( b || .P ||  d || .P);     /* NE */

ods layout gridded advance=table columns=2 column_gutter=1px;
   print Pt[c=labl r=desc], Prob[L=" " c="Prob" F=bestd5.];
ods layout end;
print (sum(Prob))[L='Sum Prob'];
QUIT;

/* additional tests for correctness. Not in the blog posts. */
proc iml;
load module=(CDFBN ProbBVN BVNMC);

rho = 0.4;
a = -0.5; b =  0.5;
c = -1.2; d =  0.8;

/* Test BVN(x,y; rho) = BVN(y,x; rho) for all rho */
rho = 0.4;
p1 = CDFBN(3,2, rho);
p2 = CDFBN(2,3, rho);
print p1 p2 (p1-p2);
p1 = CDFBN(0.5,-2, rho);
p2 = CDFBN(-2,0.5, rho);
print p1 p2 (p1-p2);

/* NOTE: BVN(x,.P; rho) = CDF("Normal",x) and 
         BVN(.P,y; rho) = CDF("Normal",y) for all rho 
   Test by using 6 to approximate +Infinity */

rho = 0.4;
p1 = CDFBN(2, 6, rho);
p2 = CDFBN(2,.P, rho);
print p1 p2 (p1-p2);
p1 = CDFBN(-2, 6, rho);
p2 = CDFBN(-2,.P, rho);
print p1 p2 (p1-p2);
p1 = CDFBN(6 , 1, rho);
p2 = CDFBN(.P, 1, rho);
print p1 p2 (p1-p2);
p1 = CDFBN(6 , -1, rho);
p2 = CDFBN(.P, -1, rho);
print p1 p2 (p1-p2);


/* Test the probability computation for all 16 regions 
   Compare the results with a Monte Carlo estimate.
 */
M = expandgrid( (.M || a), (.P || b), (.M || c), (.P || d) );
desc = {'1. Plane', '2. Lower Half', '3. Upper Half', '4. Horiz Strip',
        '5. Left Half', '6. SW', '7. NW', '8. Horiz Left Strip', 
        '9. Right Half', '10. SE', '11. NE', '12. Horiz Right Strip',
        '13. Vert Strip', '14. Vert Lower Strip', '15. Vert Upper Strip', '16. Rect'};
NMC = 400000;
Mu = {0 0};
Sigma = (1  || rho) // (rho || 1);
XMC = randnormal(NMC, Mu, Sigma);

P  = j(nrow(M), 1, .);
MC = j(nrow(M), 1, .);
do i = 1 to nrow(M);
   a = M[i,1];   b = M[i,2];
   c = M[i,3];   d = M[i,4];
   P[i] = ProbBVN(a,b,c,d, rho);
   MC[i] = BVNMC(a,b,c,d, XMC);
end;
print (M ||P || MC || (P-MC))[c={'a' 'b' 'c' 'd' 'Prob' 'MC' 'Diff'} r=desc];


/* Verify: Some areas must sum to 1 or 0 */
s1 = P[2] + P[3] - P[4];  /* horiz half planes - horiz strip */
print s1;
s2 = P[5] + P[9] - P[13]; /* vert half planes - vert strip */
print s2;
/* left half plane = -Hor L Strip + NW + SW */
s3 = -P[8] + P[7] + P[6] - P[5];
print s3;
/* right half plane = -Hor R Strip + NE + SE */
s3 = -P[12] + P[11] + P[10] - P[9];
print s3;
/* upper half plane = -Vert Up Strip + NW + NE */
s4 = -P[15] + P[7] + P[11] - P[3];
print s4;
/* lower half plane = -Vert Low Strip + SW + SE */
s5 = -P[14] + P[6] + P[10] - P[2];
print s5;
/* Verify: The 9 checkerboard regions must sum to 1 */
a = 0;  b = 1;
c = -1; d = 0.5;

R=j(9,1,.);
R[1] = ProbBVN(.M,a ,d,.P, rho); /* UpLeft */
R[2] = ProbBVN( a,b ,d,.P, rho); /* UpMiddle */
R[3] = ProbBVN( b,.P,d,.P, rho); /* UpRight */
R[4] = ProbBVN(.M,a ,c, d, rho); /* MidLeft */
R[5] = ProbBVN( a,b ,c, d, rho); /* Middle */
R[6] = ProbBVN( b,.P,c, d, rho); /* MidRight */
R[7] = ProbBVN(.M,a ,.M,c, rho); /* BotLeft */
R[8] = ProbBVN( a,b ,.M,c, rho); /* BotMiddle */
R[9] = ProbBVN( b,.P,.M,c, rho); /* BotRight */
*print R;
s6 = sum(R);
print s6;
QUIT;




/****************************************************/
/* VISUALIZATION: Use series plots and shaded tails 
   how univariate CDF computes left-tailed probabilities, but
   you can use the CDF to compute the probability on a finite interval.
*/
/*********************************************/
/* Draw univariate P(a<X<b) = CDF(b) - CDF(a) */
data density;
low  = quantile("Normal", 0.20);
high = quantile("Normal", 0.75);
drop low high;
/* if step size is dx, location of quantile can be off by dx/2 */
do x = -4 to 4 by 0.01;            /* for x in [xMin, xMax] */
   pdf = pdf("Normal", x);         /* call PDF to get density */
   if x<=low  then up1 = pdf; else up1 = .;
   if x<=high then up2 = pdf; else up2 = .;
   output;
end;
run; 
data Arrows;
yLoc = 0.5; low=-4; high = quantile("Normal", 0.75); text="CDF(b)"; xlabl='b'; 
output;
yLoc = 0.45; low=-4; high = quantile("Normal", 0.2); text="CDF(a)"; xlabl='a';
output;
run;
data All;
set density Arrows;
run;

ods graphics / reset;
title "Univariate Normal Probability";
title2 "P(a < X < b) for X ~ N(0,1)";
proc sgplot data=All noautolegend;
   series x=x y=pdf / legendlabel="Normal Density";
   band x=x upper=up2 lower=0 / fillattrs=GraphData2 legendlabel="CDF(b)" transparency=0.5;
   band x=x upper=up1 lower=0 / fillattrs=GraphData1 legendlabel="CDF(a)" transparency=0.5;
   highlow y=yLoc low=low high=high / lowcap=FilledArrow highlabel=text;
   refline high / axis=x label=xlabl labelpos=min;
   xaxis label='x';
   yaxis label='Density' offsetmin=0;
run;


/****************************************************/
/* VISUALIZATION: Use heat maps to show BVN(0, 0.4) and overlay 
   reference lines that define rectangular regions, quadrants, and strips
*/

/* This example uses
rho = 0.4;
a = -0.5; b =  0.5;
c = -1.2; d =  0.8;
*/
/* visualize the case where rho=0 */
data BivarDensity ;
rho = 0.4;
twopi = 2*constant('pi');
step = 0.1;
do y = -3 to 3 by step; 
   do x = -3 to 3 by step;
      Density = exp( -(x**2-2*rho*x*y+y**2)/(2*(1-rho**2)) )
                / (twopi*sqrt(1-rho*rho));   
      output;
   end;
end;
run;

data Rect ;
a = -0.5; b =  0.5;
c = -1.2; d =  0.8;
ID = 1;
Rx = a; Ry = c; output;
Rx = b; Ry = c; output;
Rx = b; Ry = d; output;
Rx = a; Ry = d; output;
drop a b c d;
run;
data Region;
a = -0.5; b =  0.5;
c = -1.2; d =  0.8;
Text="SW"; Tx = a-1; Ty = c-1; output;
Text="S "; Tx = (a+b)/2; Ty = c-1; output;
Text="SE"; Tx = b + 1; Ty = c-1; output;
Text="W "; Tx = a-1; Ty = (c+d)/2; output;
Text="C "; Tx = (a+b)/2; Ty = (c+d)/2; output;
Text="E "; Tx = b + 1; Ty = (c+d)/2; output;
Text="NW"; Tx = a-1; Ty = d + 1; output;
Text="N "; Tx = (a+b)/2; Ty = d + 1; output;
Text="NE"; Tx = b + 1; Ty = d + 1; output;
drop a b c d;
run;

data All;
set BivarDensity Rect Region;
run;
%let REFOPT = lineattrs=(thickness=2) curvelabelposition=min curvelabellocation=outside;

proc template;
define statgraph MyContour;
dynamic _X _Y _DENSITY _RX _RY _TX _TY;
begingraph / designwidth=640 designheight=480;
   entrytitle halign=center 'Bivariate Normal Density: rho=0.4';
   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
      layout overlay /;
         contourplotparm x=_X y=_Y z=_DENSITY / name='contour' contourtype=GRADIENT colormodel=TwoColorRamp gridded=true;
         referenceline x=eval(COLN(-2, -1, 0, 1, 2)) / datatransparency=0.85;
         referenceline y=eval(COLN(-2, -1, 0, 1, 2)) / datatransparency=0.85;
         referenceline x=-0.5 / curvelabel="a" &REFOPT;
         referenceline x= 0.5 / curvelabel="b" &REFOPT;
         referenceline y=-1.2 / curvelabel="c" &REFOPT;
         referenceline y= 0.8 / curvelabel="d" &REFOPT;
         polygonplot x=_RX y=_RY ID=ID / name='poly' display=(fill outline) fillattrs=(transparency=0.9) 
                 outlineattrs=(color=black thickness=2);
         scatterplot x=_RX y=_RY / name='scat';
         textplot text=Text x=_TX y=_TY / name='text' textattrs=(size=16);
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.ALL template=MyContour;
dynamic _X="X" _Y="Y" _DENSITY="DENSITY" _RX="RX" _RY="RY" _TX="TX" _TY="TY";
run;


/****************************************************/
/* VISUALIZATION: Use heat map to visualize four 
   quadrants defined by the point (x0,y0) = (0.5, 1) */
data Region2;
a = 0.5; b = 1;
Text="P1"; Tx = a-1; Ty = b-1; output;
Text="P2"; Tx = a+1; Ty = b-1; output;
Text="P3"; Tx = a-1; Ty = b+1; output;
Text="P4"; Tx = a+1; Ty = b+1; output;
run;

data R;
set BivarDensity Region2;
run;
%let REFOPT = lineattrs=(thickness=2) curvelabelposition=min curvelabellocation=outside;

proc template;
define statgraph MyContourB;
dynamic _X _Y _DENSITY _TX _TY;
begingraph / designwidth=640 designheight=480;
   entrytitle halign=center 'Bivariate Normal Density: rho=0.4';
   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
      layout overlay /;
         contourplotparm x=_X y=_Y z=_DENSITY / name='contour' contourtype=GRADIENT colormodel=TwoColorRamp gridded=true;
         referenceline x=0.5 / curvelabel="x0" &REFOPT;
         referenceline y=1   / curvelabel="y0" &REFOPT;
         scatterplot x=a y=b / name='scat';
         textplot text=Text x=_TX y=_TY / name='text' textattrs=(size=16);
      endlayout;
   endlayout;
endgraph;
end;
run;

proc sgrender data=WORK.R template=MyContourB;
dynamic _X="X" _Y="Y" _DENSITY="DENSITY" _TX="TX" _TY="TY";
run;


