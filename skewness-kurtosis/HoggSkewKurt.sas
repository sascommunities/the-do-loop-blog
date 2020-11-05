
proc iml;
/* Assume 0<= p <= 1.
   For a discussion of integer and fractional parts, see
   https://blogs.sas.com/content/iml/2020/02/10/fractional-part-of-a-number-sas.html
*/
/* ASSUME x is sorted in ascending order.
   Compute the mean of the lower p_th percentile of the data.
   Use linear interpolation if p*N is not an integer. Assume 0<= p <= 1. 
*/
start meanLowerPctl(p, x);
   N = nrow(x);  
   k = floor(p*N);                 /* index less than p*N */
   r = p*N - k;                    /* r >= 0 */
   if k<1 then m = x[1];           /* use x[1] as the mean */
   else if k>=N then m = mean(x);  /* use mean of all data */
   else do;                        /* use interpolation */
      sum = sum(x[1:k]) + r*x[k+1];
      m = sum / (k+r);
   end;
   return( m );
finish;

/* ASSUME x is sorted in ascending order.
   Compute the mean of the upper p_th percentile of the data.
*/
start meanUpperPctl(p, x);
   q = 1-p; 
   N = nrow(x); 
   k = ceil(q*N);                  /* index greater than p*N */
   r = k - q*N;                    /* r >= 0 */
   if k<=1 then m = mean(x);       /* use mean of all data */
   else if k>=N then m = x[N];     /* use x[1] as the mean */
   else do;                        /* use interpolation */
      sum = sum(x[k+1:N]) + r*x[k];
      m = sum / (N-k+r);
   end;
   return( m );
finish;

/* run the functions on a small test example */
x = {2, 4, 5, 7, 8, 8, 9, 9, 12, 16};  /* NOTE: sorted in increasing order */
obsNum = T(1:nrow(x));
Pctl = obsNum / nrow(x);
meanLow = j(nrow(Pctl), 1, .);
meanUp = j(nrow(Pctl), 1, .);
do i = 1 to nrow(Pctl);
   meanLow[i] = meanLowerPctl(Pctl[i], x);
   meanUp[i] = meanUpperPctl(Pctl[i], x);
end;
print Pctl meanLow meanUp;


/* plot the mean of the lower (or upper) tail as a function of p */
p = do(0, 1, 0.01);
mLow = j(1, ncol(p), .);
mUp  = j(1, ncol(p), .);
do i = 1 to ncol(p);
   mLow[i] = meanLowerPctl(p[i], x);
   mUp[i]  = meanUpperPctl(p[i], x);
end;

title "Mean of Lower p Percentile of Data";
call series(p, mLow) grid={x y} 
                   xvalues=do(0,1,0.1) 
                   label={"Percentile" "Mean"};
title "Mean of Upper p Percentile of Data";
call series(p, mUp) grid={x y} 
                   xvalues=do(0,1,0.1) 
                   label={"Percentile" "Mean"};

/* notice that this mean is NOT a linear function between data points! 
   For example, the mean between the first two data values is
   (x[1] + r*x[2]) / (1+r)
   which is not linear in r.
*/

/* estimate of Hogg skewness */
start HoggSkew(_x);
   x = _x;
   call sort(x);
   L05 = meanLowerPctl(0.05, x);
   U05 = meanUpperPctl(0.05, x);
   m25 = mean(x, 'trim', 0.25);
   Skew_Hogg = (U05 - m25) / (m25 - L05); /* Right Skewness: Hogg (1974, p. 918) */
   return Skew_Hogg;
finish;

/* estimate of Hogg kurtosis */
start HoggKurt(_x);
   x = _x;
   call sort(x);
   L20 = meanLowerPctl(0.20, x);
   U20 = meanUpperPctl(0.20, x);
   L50 = meanLowerPctl(0.50, x);
   U50 = meanUpperPctl(0.50, x);
   Kurt_Hogg = (U20 - L20) / (U50 - L50); /* Kurtosis: Hogg (1974, p. 913) */
   return Kurt_Hogg;
finish;


/* generate random sample from the exponential distribution */
call randseed(12345, 1);
N = 1E5;
x = randfun(N, "Expon");

Skew_Hogg = HoggSkew(x); /* find the Hogg skewness for exponential data */
print Skew_Hogg;

Kurt_Hogg = HoggKurt(x); /* find the Hogg kurtosis for exponential data */
print Kurt_Hogg;

*-----------------------------------------------------;

/* compare with exact value of the Hogg measures for the exponential DISTRIBUTION */
/* find the expected value of the truncated exponential distribution on the interval [a,b] 
   For details, see
   https://blogs.sas.com/content/iml/2020/11/04/expected-value-of-tail.html
*/
start Integrand(x);
   return x*pdf("Expon",x);
finish;
 
/* if f(x) is the PDF and F(x) is the CDF of a distribution, the expected value on [a,b] is
   (\int_a^b  x*f(x) dx) / (CDF(B) - CDF(a))
*/
start Expo1Moment(a,b);
   call quad(numer, "Integrand", a||b );
   /* define CDF(.M)=0 and CDF(.P)=1 */
   cdf_a = choose(a=., 0, cdf("Expon", a));  /* CDF(a) */
   cdf_b = choose(b=., 1, cdf("Expon", b));  /* CDF(b) */
   ExpectedValue = numer / (cdf_b - cdf_a);
   return ExpectedValue;
finish;
 

/* Hoggs skewness for exponential distribution 
   (U05 - m25) / (m25 - L05)
*/
p = 0.05;
/* expected value of lower p_th percentile */
qLow = quantile("Expon", p);
L05 = Expo1Moment(0, qLow);
/* expected value of upper p_th percentile */
qUp = quantile("Expon", 1-p);
U05 = Expo1Moment(qUp, .I); /* .I = infinity */
/* expected value of middle 50% */
q25 = quantile("Expon", 0.25);
q75 = quantile("Expon", 0.75);
m25 = Expo1Moment(q25, q75);
/* exact value for exponential distribution */
HoggSkew = (U05 - m25) / (m25 - L05);
print HoggSkew;


/* Hoggs kurtosis for exponential distribution */
q20 = quantile("Expon", 0.2);
L20 = Expo1Moment(0, q20);
q80 = quantile("Expon", 1-0.2);
U20 = Expo1Moment(q80, .I);     /* .I = infinity */
q50 = quantile("Expon", 0.5);
L50 = Expo1Moment(0, q50);
U50 = Expo1Moment(q50, .I);     /* .I = infinity */
HoggKurt = (U20 - L20) / (U50 - L50);
print HoggKurt;

QUIT;
title;
