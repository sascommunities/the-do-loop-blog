/* SAS program to accompany the article 
   "The Poisson-binomial distribution for hundreds of parameters"
   by Rick Wicklin, published 07OCT2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/10/07/poisson-binomial-hundreds-of-parameters.html

   This program shows how to use the refined normal approximation (RNA)
   to approximate the PDF, CDF, and quantiles of 
   the Poisson-binomial distribution. 

   Reference: Hong, Y., (2013) 
   "On computing the distribution function for the Poisson binomial distribution,"
   Computational Statistics & Data Analysis, 59:41–51. 
   https://www.researchgate.net/publication/257017356_On_computing_the_distribution_function_for_the_Poisson_binomial_distribution
   
   The RNA is EQN 14 from Hong (2013, p. 10).
   You can use the RNA when N >= 500 and when 
   mu - 5*sigma is far from 0 and 
   mu + 5*sigma is far from N.
*/

proc iml;
/* For comparison: compute the Poisson-binomial PDF by using
   the RF1 recurrence relation. Hong (2013) says that RF1 
   can be used up to N=1000. See
   https://blogs.sas.com/content/iml/2020/09/30/pdf-cdf-quantile-poisson-binomial.html
*/
start PDFPoisBinom(_p);
   p = colvec(_p); N = nrow(p);
   R = j(1, N+1, 1);        /* R = previous row of probability matrix */
   S = j(1, N+1, 0);        /* S = current row of probability matrix */
   R[2:N+1] = cuprod(1-p);  /* first row (k=0) is cumulative product */
   pdf0 = R[N+1];           /* pdf(0) is last column of row */ 
   pdf = j(1, N);
   do k = 1 to N;           /* for each row k=1..N */
      S[,] = 0;             /* zero out the current row */
      do j = k to N;        /* apply recurrence relation from left to right */
         S[j+1] = (1-p[j])*S[j] + p[j]*R[j];  
      end;
      pdf[k] = S[N+1];      /* pdf(k) is last column of row */
      R = S;                /* store this row for the next iteration */
   end;
   return (pdf0 || pdf);    /* return PDF as a row vector */
finish;

/***********************************************************/
/* Compute the CDF by using the Refined Normal Approximation 
  (RNA), which is EQN 14 from Hong (2013, p. 10)           */
/***********************************************************/

start CDFPoisBinomRNA(_p, k=);
   p = colvec(_p); N = nrow(p);
   mu = sum(p);                              /* mean of P-B distribution */
   sigma = sqrt(sum(p#(1-p)));               /* std dev of P-B */
   skew = sum(p#(1-p)#(1-2*p)) / sigma##3;   /* skewness of P-B */
   if IsSkipped(k) then 
      k = 0:N;                               /* Default: Compute entire CDF */
   x = (k + 0.5 - mu) / sigma;
   /* Compute refined normal approximation by using Eqn 14 of Hong (2013) */
   CDF = cdf("Normal",x) + skew*(1-x##2)#pdf("Normal",x)/6;
   /* trap and map into [0,1]: 
      https://blogs.sas.com/content/iml/2020/10/05/trap-and-map-invalid-values.html */
   CDF = (CDF <> 0) >< 1;                    /* force CDF into [0,1] */
   return(CDF);
finish;

start PDFPoisBinomRNA(_p, k=);
   p = colvec(_p); N = nrow(p);
   if IsSkipped(k) then 
      k = 0:N;                               /* Default: Compute entire PDF */
   if k[1]=0 then do;
      CDF = CDFPoisBinomRNA(p, k);
      PDF = dif(CDF`);
      PDF[1] = 0;
   end;
   else do;
      CDF = CDFPoisBinomRNA(p, (k[1]-1) || k);
      PDF = dif(CDF`, 1, 1);   /* delete the first row (missing value) NEW FEATURE */
   end;
   return( PDF` );
finish;

/* The quantile is the smallest value for which the CDF is 
   greater than or equal to the given probability. 
   Consequently, the quantile function is a step function. 
*/
start QuantilePoisBinomRNA(p, _probs);     /* p = PB parameters; _prob = vector of probabilities for quantiles */
   cdf = CDFPoisBinomRNA(p);               /* index is one-based */
   probs = colvec(_probs);                 /* find quantiles for these probabilities */
   q = j(nrow(probs), 1, .);               /* allocate result vector */
   do i = 1 to nrow(probs);
      idx = loc( probs[i] <= CDF );        /* find all x such that p <= CDF(x) */ 
      q[i] = idx[1] - 1;                   /* smallest index. Subtract 1 b/c PDF is for k=0,1,...,N  */
   end;
   return ( q );
finish;
store module=(PDFPoisBinom CDFPoisBinomRNA PDFPoisBinomRNA QuantilePoisBinomRNA);

QUIT;

/********************************************************************************/

proc iml;
load module=(PDFPoisBinom CDFPoisBinomRNA PDFPoisBinomRNA QuantilePoisBinomRNA);

/* show a distribution for large N (large number of p) */
call randseed(12345);
N = 500;
p = randfun(N,"Uniform");        /* uniform distribution of p */
PDF_RF = PDFPoisBinom(p);

/* Note: PB is discrete distrib, but draw as series plot.
   Most of the PDF is almost zero, so "zoom in" */
kk = 150:350;
PDF = PDF_RF[,kk+1];

title "Poisson-Binomial Distribution";
title2 "N = 500, p_j ~ U(0,1)";
footnote J=L "mu=255.2, sigma=9.0"; 
call series(kk, PDF) grid={x y} 
            label={'Number of Successes' 'Density'};
footnote;

*--------------------------------------;

/* Compute the refined normal approximation and compare
   to the exact PDF */
PDF_RNA = PDFPoisBinomRNA(p);   /* approximation */
diff = PDF_RF - PDF_RNA;
print (max(abs(diff)))[L="Max Diff between PDFs"];

title "Difference Between Exact and Refined-Normal Approximation";
title2 "N = 500, p_j ~ U(0,1)";
call series(kk, Diff[kk]) grid={x y}
           label={'Number of Successes' 'Difference in Density'};

/* Note: You can also compute just a portion of the PDF */
PDF_Part = PDFPoisBinomRNA(p, kk);
title "Compute Only Part of the PDF";
call series(kk, PDF_Part) grid={x y}
           label={'Number of Successes' 'Difference in Density'};


*--------------------------------------;
/* For N = 500, 600, ..., 1000, compare the time to 
   compute the exact PDF to the time required to compute
   the RNA. 

   Repeat each computation numRep times and average results 
*/
numRep = 7;
NN = T( do(500, 1000, 100) );
results = j(nrow(NN),4,.);
results[,1] = NN;
do i = 1 to nrow(NN);
   N = NN[i];
   p = randfun(N,"Uniform");
   t0 = time();
   do j = 1 to numRep;
      PDF_RF = PDFPoisBinom(p);
   end;
   t1 = time() - t0;
   results[i,2] = t1 / numRep;
   t0 = time();
   do j = 1 to numRep;
      PDF_RNA = PDFPoisBinomRNA(p);
   end;
   t1 = time() - t0;
   results[i,3] = t1 / numRep;
   results[i,4] = max(abs(PDF_RF-PDF_RNA));
end;
print results[c={N RF RNA 'MaxDiff'}];

title "Maximum Difference between Densities";
title2 "Recursive Formula (RF) Versus Refined Normal Approximation (RNA)";
call series(NN, results[,4]) grid={x y};

SpeedUp = results[,2] / results[,3];
title "Approximate Speedup";
title2 "time(RNA) / time(RF)";
call series(NN, SpeedUp) grid={x y};

/* Conclusion: use RF (exact method) when N < 500 and use 
   RNA (approximate method) when N >= 500 
*/

*--------------------------------------------;

/* You can compute quantiles: */
probs = {0.05, 0.25, 0.75, 0.95};
qntl = QuantilePoisBinomRNA(p, probs);
print probs qntl;

*--------------------------------------------;
/* The accuracy of the RNA depends on whether the mean of the 
   PB distribution is close to 0 or N. If the skewness of the 
   distribution of parameters is extreme, the approximation
   might not be very good.  To demonstrate, randomly generate
   p[j] ~ Beta(0.5, 6), which has skewness of 7. Most of the 
   p[j] are very small (close to 0), so 
*/
N = 500;
sp = randfun(N,"Beta", 0.05, 6);  /* skewed distribution (skew(p)=7.1). Most near 0. */
/* Compute skewness of Beta distribution:
alpha=0.05; beta = 6; sk = 2*(beta-alpha)#sqrt(alpha+beta+1)/(alpha+beta+2)/sqrt(alpha*beta);
print sk; 
*/
k = 0:15;            /* b/c p[j] so small, E(num successes) <= 15 */
PDF_RF = PDFPoisBinom(sp)[,k+1];   /* exact */
PDF_RNA = PDFPoisBinomRNA(sp, k);  /* approximate */
Diff = PDF_RF - PDF_RNA;
maxDiff = max(abs(Diff));
print maxDiff;
/* graph it */
xx = k || k;
PDF = PDF_RF || PDF_RNA;
group = j(1,ncol(PDF_RF),"RF ") || j(1,ncol(PDF_RNA), "RNA");
title "Compare Exact and Approximate Methods";
title2 "Skewness(p) =7.1";
call series(xx,PDF)group=group;


QUIT;
