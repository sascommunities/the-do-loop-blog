/* SAS program to accompany the article 
   "The Kolmogorov D distribution and exact critical values"
   by Rick Wicklin, published 24JUN2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/06/24/kolmogorov-d-distribution-exact.html

   This program shows how to compute the CDF of the Kolmogorov D distribution
   for any value of n (sample size) and at any value of D in (0, 1).

   The program uses an algorithm from Facchinetti (2009) 
   "A Procedure to Find Exact Critical Values of Kolmogorov-Smirnov Test"
   Statistica Applicata – Italian Journal of Applied Statistics Vol. 21 n. 3-4 2009
   http://luk.staff.ugm.ac.id/stat/ks/IJAS_3-4_2009_07_Facchinetti.pdf

   For a sample of size n, the method constructs a system of 
   2n+2 linear equations for 2n+2 unknowns. The system is singular, so the system
   is solved by using the Moore-Penrose pseudo-inverse matrix.
   
   The first PROC IML program defines the KolmogorovCDF function and saves it.
   The second program shows how to use the function to 
   1. Find the probability CDF(D; n)
   2. Compute and graph the entire CDF function for various values of n
   3. Display the density (PDF) by using forward-difference approximations
   4. Construct a table of critical values. 
*/

proc iml;
/* Given n and D, this routine returns the probability  Pr(D_n > D)
   Based on MATLAB program in appendix of Facchinetti (2009) 
*/
start KolmogorovCDF(n, D, vizB=0);
   nD = n*D;
   m1 = round(n*D+0.5);
   m2 = round(n-n*D-0.5);
   l1 = round(2*n*D+0.5);
   n1 = n*(1+D);
   n2 = n*(1-D);
   /* B matrices */
   B = j(2*(n+1), 2*(n+1), 0);
   /* B11 matrix */
   do k = m1+1 to n+1;
      B[k,k]=1;
   end;
   do r = m1 to n-1;
      k = T((r+1):n);
      p = (k-r)/(n1-r);   p = (p><1);
      B[k+1,r+1] = pdf("binomial", k-r, p, n-r);
   end;

   /* B22 matrix */
   do k = n+2 to n+2+m2;
      B[k,k]=1;
   end;
   do r = 0 to m2-1;
      k = T((r+1):m2);
      p = (k-r)/(n2-r);    p = (p><1);
      B[k+n+2,r+n+2] = pdf("binomial", k-r, p, n-r);
   end;

   /* B21 matrix */
   do r = m1 to m2;
      k = T(r:m2);
      p = (k-r+2*nD)/(n1-r);     p = (p><1);
      kmr = k-r; 
      nmr = repeat(n-r, nrow(k)); /* make into vector */
      do i = 1 to nrow(k);
         pi = pdf("binomial", k[i]-r, p[i], nmr[i]);
      end;
      B[k+n+2,r+1] = pdf("binomial", k-r, p, nmr);
   end;

   /* B12 matrix */
   do r = 0 to n-l1;
      k = T((l1+r):n);
      kmr = k-r; 
      nmr = repeat(n-r, nrow(k)); /* make into vector */
      p = (k-r-2*nD)/(n2-r);     p = (p><1);
      B[k+1,r+n+2] = pdf("binomial", k-r, p, nmr);
   end;

   /* C vector */
   C = j(2*(n+1), 1, 0);
   /* C1 vector */
   k = T(m1:n);
   nn = repeat(n, nrow(k));
   p = (k-nD)/n;     p = (p><1);   /* min(p,1) b/c p cannot exceed 1 */
   C[k+1] = pdf("binomial", k, p, nn);

   /* C2 vector */
   k = 0:m2;
   nn = repeat(n, nrow(k));
   p = (k+nD)/n;     p = (p><1);   /* min(p,1) b/c p cannot exceed 1 */
   C[n+2+k] = pdf("binomial", k, p, nn);

   /* Optional: visualize the structure of B, as in Fig 6 of Facchinetti (2009) */
   if vizB then do; 
      btitle = "Structure of B Matrix (n=" + strip(char(n)) + ")";
      call heatmapdisc(B^=0) colorramp={white steelblue} title=btitle;
   end;

   /* solve the linear system */
   Binv = ginv(B);
   Z = Binv*C;
   alpha = sum(Z);
   cdf = 1-alpha;
   return cdf;
finish;
store module=KolmogorovCDF;
QUIT;

/* You can compare the results of the KolmogorovCDF function with 
   the critical values in tables, such as 
   http://www.matf.bg.ac.rs/p/files/69-[Michael_J_Panik]_Advanced_Statistics_from_an_Elem(b-ok.xyz)-776-779.pdf
   Kolmogorov–Smirnov Table. from L.H. Miller, 'Tables of Percentage Points of Kolmogorov Statistics,' JASA, 51, 1956, 111–121
*/
proc iml;
load  module=KolmogorovCDF;

/* simple example of calling the function */
n = 20;
D = 0.294;
cdf = KolmogorovCDF(n, D);
print n D cdf;

ods graphics / width=400px height=400px;
cdf = KolmogorovCDF(n, D, 1);  /* Use 1 as third argument to create a heat map */

ods graphics/reset;
/* Interpretations of the calculation:
   1. The critical value of the D statistic for alpha=0.05 is D=0.294.
   2. In 5% of the random samples of size 20, the maximum absolute deviation 
      between the empirical distribution function and the theoretical
      distribution function will be 0.294 or greater. 
   3. If the null hypothesis is
         H0: random sample of size 20 is from a stated distribution
      If you observe a D statistic of 0.294 or greater for the sample,
      then reject the null hypothesis (at the (1-alpha)*100% = 95% confidence level) 
      that the sample came from the stated distribution. 
*/

/* call multiple times to create graph of Kolmogorov D distribution */
D = do(0.01, 0.99, 0.01);             /* sequence of D values */
CDF = j(1, ncol(D), .);
do i = 1 to ncol(D);
   CDF[i] = KolmogorovCDF(n, D[i]);   /* CDF for each D */ 
end;
title "Kolmogorov D CDF, n=20";
call series(D, CDF) grid={x y};       /* create a graph of the CDF */


/* create CDF and PDF for various values of n */
G = {. . .};
create KolmogorovCDF from G[colname={"n" "D" "CDF"}]; 
nn = T( {10 20 30 50 80 100} );
D =  T( do(0.005, 0.3, 0.005) || 
        do(0.31, 0.6, 0.01)   ||
        do(0.65, 0.95, 0.05) );
cdf = 0#D;
do j = 1 to nrow(nn);
   n = j(nrow(D), 1, nn[j]);
   do i = 1 to nrow(D);
      cdf[i] = KolmogorovCDF(nn[j], D[i]);
   end;
   G = n || D || cdf;
   append from G;
end;
close;

/************************************************/
/* A critical value is the root of the function
   f(D) = KolmogorovCDF(n, D) - (1-alpha)
   for any specified n and alpha value.
*/
start KolmogorovQuantile(D) global (n, alpha);
   return  KolmogorovCDF(n, D) - (1-alpha);
finish;

/* Example: Find critical value for n=20 and alpha=0.05 */
n = 20;
alpha = 0.05;
DCrit = froot("KolmogorovQuantile", {0.01, 0.99});
print n alpha DCrit[format=7.5];

/************************************************/
/* Create a table of critical values for the following vector of n and alpha values.
   This duplicates part of 
   "Table 1: Exact Critical Values" on p. 352 of Facchinetti (2009).
*/
nVals = do(20, 35, 5);
alphaVals = {0.001 0.01 0.05 0.10};
CritVal = j(ncol(nVals), ncol(alphaVals), .);
do i = 1 to ncol(nVals);
   n = nVals[i];
   do j = 1 to ncol(alphaVals);
      alpha = alphaVals[j];
      CritVal[i,j] = froot("KolmogorovQuantile", {0.01, 0.99});
   end;
end;
print CritVal[r=nVals c=alphaVals format=7.5 L="Critical Values (n,alpha)"];

QUIT;

/* Overlay the family of CDF plots. Add labels for easy viz */
data labels;
length labl $5;
input n D Y labl $;
datalines;
10  0.2  0.21 n=10
20  0.19 0.52 n=20
30  0.20 0.82 n=30
50  0.17 0.9  n=50
80  0.15 0.95 n=80
100 0.07 0.97 n=100
;

data KCDF;
set KolmogorovCDF labels;
run;

title "Distribution of Kolmogorov Dn";
proc sgplot data=KCDF noautolegend;
   series x=D y=CDF / group=n;
   text x=D y=Y text=labl / group=n position=right strip;
   xaxis grid max=1 offsetmax=0;
   yaxis grid label="Probability Pr(Dn < D)";
run;

/* use forward difference approximation to plot the PDF */
data KolmogorovPDF;
set KolmogorovCDF;
by n;
dcdf = dif(CDF);
dD = dif(D);
if first.n then do;
   dcdf = .; dD = .; PDF=0;
end;
else PDF = dcdf / dD;   /* finite difference derivative */
run;

title "Density of Kolmogorov Dn";
proc sgplot data=KolmogorovPDF;
   series x=D y=PDF / group=n;
   keylegend / sortorder=reverseauto;
   xaxis grid max=1 offsetmax=0;
   yaxis grid label="Density";
run;

/* for featured image */
title "Density of Kolmogorov Dn";
proc sgplot data=KolmogorovPDF;
   series x=D y=PDF / group=n;
   keylegend / sortorder=reverseauto position=E;
   xaxis grid max=0.4;
   yaxis grid label="Density";
run;
title;

