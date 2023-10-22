/* SAS program to accompany the articles 
   "The generalized birthday problem"
   by Rick Wicklin, published 23OCT2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/10/23/generalized-birthday-problem.html

   and
   "Quantiles of the generalized birthday problem"
   by Rick Wicklin, published 25OCT2023 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2023/10/25/quantiles-birthday-problem.html

   This program shows how to compute exact or approximate probabilities
   for the generalized birthday problem: https://en.wikipedia.org/wiki/Birthday_problem

   The approximation is Eqn 7.5 on p. 857 of
   Diaconis, P. and Mosteller F. (1989). Methods for studying coincidences. 
   JASA, 84, 853--861.
   https://media.socastsrm.com/wordpress/wp-content/blogs.dir/2323/files/2019/06/DC.pdf

   Typically, we only care about the CDF and quantile function for the 
   generalized birthday problem. For example, R supports the
   pbirthday and qbirthday functions:
   https://rdocumentation.org/packages/stats/versions/3.6.2/topics/birthday
   I used the R functions as models for the CDF and Quantile functions.

   For completeness, this article also supplies functions for the 
   PDF and for random-number generation.
*/

/* DEFINE THE FUNCTIONS FOR THE GENERALIZED BIRTHDAY PROBLEM */

/* The generalized birthday problem in IML:
   What is the probability that at least k=nDup people 
   in a group of N people share a characteristic
   (such as a birthday) if there are C=nCat different
   equally likely categories.
   In the classical problem, nDup=2 and nCat=365 equally likely birthdays.
   We want the probability of a shared birthday when there is 
   a group of N people. The "birthday paradox" is that in a group of 
   N=23, the probability exceeds 0.5.
*/

proc iml;
/* evaluate y=exp(x)-1 accurately when x is near 0.
   See https://blogs.sas.com/content/iml/2023/10/16/expxm1.html */
start expm1(x);
   if abs(x) > 0.03 then return( exp(x)-1.0 );
   /* approximate by using a Taylor polynomial */
   Taylor = x*( 1 + x/2*(1 + x/3*(1 + x/4*(1 + x/5*(1 + x/6)))) );
   return( Taylor );
finish;

/* The generalized birthday problem:
   In a room that contains N people who have equal probability of 
   belonging to nCat different categories, 
   return the probability that nDup people share a category.
   The classic birthday problem is nCat=365 and nDup=2.
*/
start CDFBirthday(N, nCat=365, nDup=2);
   if nDup < 2 then return(1.0);          /* everyone shares a BDay with themselves */
   if nDup > N then return(0.0);          /* can't have nDup matches if there aren't that many people! */
   if N > nCat*(nDup-1) then return(1.0); /* Ex: If N>365, there must be a shared BDay */
   /* For NDup=2, use an explicit formula. See
      https://blogs.sas.com/content/iml/2012/04/09/vectorized-computations-and-the-birthday-matching-problem.html
   */
   if nDup = 2 then do;
      i = 1:N;             /* enumerate the N individuals */
      iQ = 1 - (i-1)/nCat; /* individual probabilities */
      Q = prod(iQ);        /* compute probability of "no match" for N people */
      return(1 - Q);       /* prob of a match */
   end;
   /* Otherwise, use formula (Eqn 7.5 on p. 857) of
      Diaconis, P. and Mosteller F. (1989)
      See also the pbirthday function in R, which implements the formula. */
   k = nDup; c=NCat;                 /* for brevity, change notation */
   LHS = N * exp(-N/(c*k)) / (1 - N/(c*(k+1)))##(1/k);
   /* take log of both sides and solve for log(-log(1-p)) */
   LmL1mp = k*log(LHS) - (k-1)*log(c) - lgamma(k+1);
   /* now solve for p in y=log(-log(1-p)).
      Note: use of EXPM1 function in case of small probabilities */
   return( -expm1(-exp(LmL1mp)) );   /* = 1.0 - exp(-exp(LmL1mp)) */
finish;

/* for a discrete distribution, PDF(k) = CDF(k) - CDF(k-1) */
start PDFBirthday(N, nCat=365, nDup=2);
   CDF = CDFBirthday(N, nCat, nDup);
   if N <= 1 then    
      return CDF;
   return( CDF - CDFBirthday(N-1, nCat, nDup) );
finish;

/* Definition of quantile for a discrete distribution:
   For a probability p, return smallest integer, N, such that 
   CDF(N, nCat, nDup) >= p
*/
start QuantileBirthday(p, nCat=365, nDup=2);
   if p <= 0 then  
      return(1);
   if (p >= 1) then
      return(1 + nCat*(nDup-1));
   /* Otherwise, use formula (Eqn 7.5 on p. 857) of
      Diaconis, P. and Mosteller F. (1989)
      See also the qbirthday function in R, which implements the formula.
   */
   k = nDup; c=NCat;  /* for brevity, change notation */
   /* approx value of N */
   N = exp(((k-1)*log(c) + lgamma(k+1) + log(-log1px(-p)))/k);
   N = ceil(N);
   pN = CDFBirthday(N, c, k);
   if pN < p then do;
      /* keep increasing N until we satisfy inequality */
      N = N+1;
      do while (CDFbirthday(N, c, k) < p);
         N = N+1;
      end;
   end;
   else if (CDFBirthday(N-1, c, k) >= p) then do;
      /* decrease N while maintain the inequality */
      N = N - 1;
      do while (CDFBirthday(N-1, c, k) >= p);
         N = N-1;
      end;
   end;
   return N;
finish;

/* Return a random variate, N, from the birthday distribution.
   N represents the size of a group in which there are nDup
   shared birthdays */
start RandBirthday(sampleSize, nCat=365, nDup=2);
   randN = j(sampleSize, 1);
   u = j(sampleSize, 1);
   call randgen(u, "Uniform");
   do i = 1 to nrow(u);
      randN[i] = QuantileBirthday(u[i], nCat, nDup);
   end;
   return randN;
finish;
store module=(expm1 CDFBirthday PDFBirthday QuantileBirthday RandBirthday);
QUIT;
