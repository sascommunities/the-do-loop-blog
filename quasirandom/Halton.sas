/* A set of SAS IML functions that generate quasirandom points 
   in d dimensions by using Halton sequences. See
   https://blogs.sas.com/content/iml/2025/11/17/quasi-monte-carlo.html
*/
proc iml;
/* Compute LOG_b(x) for a vector of x values and for an integer base, b > 0.
   If an input value is not positive, this function returns a missing value.
   Computers represent numbers internally in base b=2, so use the change-of-base formula
   log_b(x) = log2(x) / log2(b)
   to compute the logarithm in any base, b.
*/
start logbase(x, base=10);
   if all(x)>0 then 
      return log2(x) / log2(base);
   /* otherwise, return missing values for x <= 0 */
   y = j(nrow(x), ncol(x), .);
   idx = loc(x>0);
   if ncol(idx)>0 then 
      y[idx] = log2(x[idx]) / log2(base);
   return y;
finish;

/* This function returns the number of digits in the base-b representation of an integer, x.
   If n >= 0 is a base-10 integer, it has
      k = ceil(log_b(n+1))
   digits when represented in base b. See
   https://blogs.sas.com/content/iml/2015/08/31/digits-in-integer.html
*/
start numDigBase(x, base=10);
   n = round(x);              /* ensure argument is an integer */
   return ceil( logbase(n+1,base) );
finish;

/* Convert integer x > 0 to a row vector in base b.
   If x is a column vector, return a matrix where each row is the base-b representation of x[i].
   The most significant bit is to the left; the least significant bit is to the right.
   For example, n=15 and base=3 gives (120)_3 because 15 = 1*3##2 + 2*3##1 + 0*3##0.
*/
start convertToBase(x, base);
   n = round(x);     /* ensure inputs are integers */
   numDig = numDigBase(max(n), base);
   /* For an explanation, see https://blogs.sas.com/content/iml/2015/08/31/digits-in-integer.html */
   c = j(nrow(n), numDig, 0);
   do i = 1 to numDig;
      a = mod(n, base);
      n = floor( n/base );
      c[ , numDig - i + 1] = a;
   end;
   return c;
finish;

/* Convert row vectors from any base to a fraction in base 10.
   The row vector represents coefficients for a fraction in base b.
   For example:
   Base 2:  c = {1 1 0 1}   represents 1/2 + 1/4 + 0/8 + 1/16
   Base 3:  c = {0 1 2 1}   represents 0/3 + 1/9 + 2*27 + 1/81.
   The input, c, can also be a matrix for which each row c[i,] is a vector of coefficients.
   For example, in base 3, define:
   c = {0 1 2 1,   
        1 2 1 0 };
   See https://blogs.sas.com/content/iml/2011/11/16/converting-from-base-2-to-base-10.html 
*/
start ConvertFromBase(c, base);
   pow = (ncol(c)-1):0;
   factor = base##pow;
   base10 = (c # factor)[ ,+];  /* c[k-1]*b##(k-1) + ... + c[1]*b##1 + c[0]*b##0 */ 
   return base10;
finish;

/* Convert row vectors from any base to a fraction in base 10.
   The row vector represents coefficients for a radix fraction in base b.
   For example:
   Base 2:  c = {1 1 0 1}   represents 1/2 + 1/4 + 0/8 + 1/16
   Base 3:  c = {0 1 2 1}   represents 0/3 + 1/9 + 2/27 + 1/81.
   The input, c, can also be a matrix for which each row c[i,] is a vector of coefficients.
*/
start ConvertFracFromBase(c, base);
   pow = -( 1:ncol(c) );         /* decreasing powers */
   factor = base##pow;           /* b^{-1}, b^{-2}, ... */
   fract10 = (c # factor)[ ,+];  /* fractions */
   return fract10;
finish;

/* reverse the order of the columns of a matrix */
start reverseCols(M);
   return M[, ncol(M):1];
finish;

/* generate a Halton sequence for a specified base
   INPUT: _n : a column vector of integers
          base : a positive scalar integer, often chosen to be a prime number 
   OUTPUT: a vector of length nrow(_n) that contains a Halton sequence
*/
start HaltonSeq1( _n, base );
   n = colvec(_n);
   C1 = ConvertToBase(n, base);          /* convert integers to base-b digits */
   C = reverseCols(C1);                  /* reverse digits to get fraction coefficients */
   return ConvertFracFromBase(C, base);  /* convert base-b fraction to base-10 decimal value */
finish;

/* Create a function to generate D Halton sequences as columns of a matrix.
   INPUT: N - number of points in each sequence
          bases - vector of D prime numbers (for example; {2, 3, 5})
   OUTPUT: an N x D matrix of quasirandom points in [0, 1]^D
   The i_th column is the Halton sequence for the digits 1:numPts in base b=bases[i].
*/
start HaltonSeqD( N, bases );
   D = nrow(colvec(bases));
   X = j(N, D, .);
   v = T(1:N);
   do i = 1 to D;
      X[,i] = HaltonSeq1(v, bases[i]);
   end;
   return X;
finish;

/* Create a function that generates numPts quasirandom points in dimension d 
   by using Halton sequences.
   INPUT:  numPts : a positive integer. The number of rows to output.
           dim    : 1 <= dim is the number of columns to output.
   OUTPUT: X is a matrix. Each column is a Halton sequence for the 
           digits 1:numPts in a prime-number base.
*/
start QRand_Halton( numPts, dim );
   primes = {   /* vector of 168 primes less than 1000 */
      2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
      59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
      127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
      191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
      257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
      331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
      401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
      467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557,
      563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619,
      631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
      709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787,
      797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
      877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953,
      967, 971, 977, 983, 991, 997 };
   primeVec = primes[1:dim];
   return HaltonSeqD( numPts, primeVec );
finish;

store module=(
  logbase numDigBase
  ConvertToBase ConvertFromBase
  reverseCols ConvertFracFromBase HaltonSeq1
  HaltonSeqD QRand_Halton
);
QUIT;
