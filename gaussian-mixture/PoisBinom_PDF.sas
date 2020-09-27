

/* SAS program to accompany the article 
   "Density, CDF, and quantiles for the Poisson-binomial distribution"
   by Rick Wicklin, published 30Sep2020 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2020/09/30/pdf-cdf-quantile-poisson-binomial.html

   This program shows how to compute the PDF, CDF, and quantiles of 
   the Poisson-binomial distribution. The PDF is computed by using
   a recurrence relation.

   Reference: Hong, Y., (2013) 
   "On computing the distribution function for the Poisson binomial distribution,"
   Computational Statistics & Data Analysis, 59:41–51. 
   https://www.researchgate.net/publication/257017356_On_computing_the_distribution_function_for_the_Poisson_binomial_distribution

   Hong (2013) compares various methods and says that RF1 is 
   better than RF2 and can be used up to N=1000.
   Let's implement the RF1 method. 
*/
title;footnote;

/* This program generates the ENTIRE matrix of probabilities. 
   Used to illustrate how the algorithm works. 
   Do NOT compute the entire matrix in practice!
*/
proc iml;
p = {0.2, 0.2, 0.3, 0.3, 0.4, 0.6, 0.7, 0.8, 0.8, 0.9};  /* 10 trials */
N = nrow(p);

/* Initialize first row, which is cumulative product of (1-p[j]) */
R = j(1, N+1, 1);        /* previous row of recursion matrix, M */
S = j(1, N+1, 0);        /* current row of M */
R[2:N+1] = cuprod(1-p);  /* k=0: row is cumulative product */
xi0 = R[N+1];            /* xi[0] is Nth column of row */ 
xi = j(1, N);
M = j(N+1, N+1,0); M[1,] = R;
do k = 1 to N;           /* for each row k=1..N */
   S[,] = 0;             /* zero out the current row */
   do j = k to N;        /* apply recurrence relation */
      S[j+1] = (1-p[j])*S[j] + p[j]*R[j];  
   end;
   xi[k] = S[N+1];       /* xi[k] is Nth column of row */
   R = S;                /* store row and go to next iteration */
   M[k+1,] = R;
end;

print M[F=Best5. r=('k=0':'k=10') c=('j=0':'j=10')];

/* To illustrate first step: */
do i = 2 to N;
   M[i, i:N+1] = .;
end;
M[N+1,N+1]=.;
print M[F=Best5. r=('k=0':'k=10') c=('j=0':'j=10')];
QUIT;


/* Start the real program: Compute PDF, CDF, and quantiles */
proc iml;
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

p = {0.2, 0.2, 0.3, 0.3, 0.4, 0.6, 0.7, 0.8, 0.8, 0.9};  /* 10 trials */
pdf = PDFPoisBinom(p);
print pdf[c=('0':'10') F=Best5.];

N = nrow(p);
x = 0:N;
*call scatter(x, pdf);

/* write to SAS data set and plot distribution as a bar chart */
create PBPDF var {x pdf}; append; close;
submit;
title "Poisson-Binomial PDF";
title2 "N=10";
proc sgplot data=PBPDF;
   vbarbasic x / response=pdf;
   yaxis grid label="Probability";
   xaxis label="Number of Successes" integer type=linear;             /* force TYPE=LINEAR */
run;
endsubmit;

/*
Compute the CDF and quantiles of discrete distributions
https://blogs.sas.com/content/iml/2017/11/22/cdf-quantiles-discrete-distribution.html
*/
/* The CDF is the cumulative sum. It is a step function */
start CDFPoisBinom(p);
   pdf = PDFPoisBinom(p);
   return ( cusum(pdf) );
finish;

cdf = CDFPoisBinom(p);
print cdf[c=('0':'10') F=Best5.];
/*
call scatter(x, cdf) other="refline 0 1 / axis=y;";
*/

/* write to SAS data set and plot distribution as a step function */
create PBCDF var {x cdf}; append; close;
submit;
title "Poisson-Binomial CDF";
title2 "N=10";
proc sgplot data=PBCDF;
   step x=x y=cdf / markers markerattrs=(symbol=CircleFilled);
   yaxis grid label="Cumulative Probability";
   xaxis label="Number of Successes" integer type=linear;
run;
endsubmit;

/* 
The quantile is the smallest value for which the CDF is 
greater than or equal to the given probability. 
Consequently, the quantile function is a step function. 
*/
start QuantilePoisBinom(p, _probs);        /* p = PB parameters; _prob = vector of probabilities for quantiles */
   cdf = CDFPoisBinom(p);                  /* index is one-based */
   probs = colvec(_probs);                 /* find quantiles for these probabilities */
   q = j(nrow(probs), 1, .);               /* allocate result vector */
   do i = 1 to nrow(probs);
      idx = loc( probs[i] <= CDF );        /* find all x such that p <= CDF(x) */ 
      q[i] = idx[1] - 1;                   /* smallest index. Subtract 1 b/c PDF is for k=0,1,...,N  */
   end;
   return ( q );
finish;
 
probs = {0.05, 0.25, 0.75, 0.95};
qntl = QuantilePoisBinom(p, probs);
print probs qntl;

/* compare with empirical density of sample */
/* Random sample from the Poisson-Binomial distribution. The 
   probability if success for the i_th Bernoulli trial is p_i.
   Input: m = number of observations in sample
          p = column vector of N probablities 
   Output: m realizations of X ~ PoisBinom(p) 
*/
start RandPoisBinom(m, p);
   N = nrow(p)*ncol(p); /* number of indep trials */
   b = j(N, m);         /* each column is a binary indicator var */
   call randgen(b, "Bernoulli", p); 
   PB = b[+, ];         /* numSuccesses = sum each column */
   return PB;           /* m realizations of X ~ PoisBinom(p) */
finish;

PB = RandPoisBinom(5000, p);

/* Tabulate counts when there are unobserved categories
   https://blogs.sas.com/content/iml/2015/10/07/tabulate-counts-ref.html */
/* output levels and frequencies for categories in x, including all 
   levels in the reference set */
start TabulateLevels(OutLevels, OutFreq, x, refSet);
   call tabulate(levels, freq, x);        /* compute the observed frequencies */
   OutLevels = union(refSet, levels);     /* union of data and reference set */ 
   idx = loc(element(OutLevels, levels)); /* find the observed categories */
   OutFreq = j(1, ncol(OutLevels), 0);    /* set all frequencies to 0 */
   OutFreq[idx] = freq;                   /* overwrite observed frequencies */
finish;

call TabulateLevels(Levels, Freq, PB`, 0:N);
empirDensity = Freq/sum(Freq);

comp = x` || pdf` || empirDensity` || (pdf` - empirDensity`);
print comp[c={'x' 'pdf' 'empirDensity' 'Diff'}];

QUIT;
title; footnote;


/******************************************/
/* Alternate computation: The RF2 method of
   Chen, Dempster, and Liu (1994) 
   for computing Xi_k (Eqn 10 of Chen, 2003)
   This is not numerically stable for large N.
*/
proc iml;
p = {0.2, 0.2, 0.3, 0.3, 0.4, 0.6, 0.7, 0.8, 0.8, 0.9};  /* 10 trials */
N = nrow(p);

xi0 = prod(1-p);
*print xi0;
pRatio = p/(1-p);
t = j(1, N);
do i = 1 to N; 
   pRatioPow = pRatio##i;
   t[i] = sum(pRatioPow);
end;
*print t;

xi = j(1, N);
k = 1;
prevXi = xi0;
xi[1] = (1/k)*sum( (-1)##(0) # t[1] # prevXi );
do k = 2 to N;
   ll = 1:k;
   prevXi = xi[,k-1:1] || xi0;
   xi[k] = (1/k)*sum( (-1)##(ll-1) # t[,ll] # prevXi );
end;
pdf = xi0 || xi;
x = 0:N;
print (x`||pdf`)[L="RF2"];
cdf = cusum(pdf);
print (x`||cdf`)[L="RF2"];
QUIT;
