/* SAS program to accompany articles 
   by Rick Wicklin, published on The DO Loop blog:
   https://blogs.sas.com/content/iml/2026/05/26/create-ecdf.html
   https://blogs.sas.com/content/iml/2026/06/08/confidence-bands-ecdf.html
   https://blogs.sas.com/content/iml/2026/06/15/confidence-band-qqplot.html

   This program defines three functions for creating ECDFs and QQ plots.
   NOTE: In SAS Viya, the ECDF function is built into SAS IML and does 
   not need to be loaded. For the Viya documentation, see
   https://go.documentation.sas.com/doc/en/pgmsascdc/v_073/casimllang/casimllang_common_sect105.htm
*/

proc iml;
/* Compute the empirical distribution function (ECDF) if a column vector of data.
   Optionally, evaluate the ECDF at values in a column vector, t. 
   If t is not specified, use t=x. So, a common syntax is 
   y = ecdf(x);
   Then use the STEP statement in PROC SGPLOT to plot y versus x.
   See https://blogs.sas.com/content/iml/2026/05/26/create-ecdf.html
   */
start ECDF(_x, _t=);
    x = colvec(_x);      /* ensure column vector */
    /* Tabulate the values of x.
       levels = unique sorted elements of x
       freq = number of duplicates for each level */
    call tabulate(levels, freq, x);
    levels = colvec(levels);
    freq = colvec(freq);
    y = cusum(freq) / sum(freq);

    if isSkipped(_t) then
       t = x;
    else
       t = colvec(_t);      /* ensure column vector */
    /* build a step function by assigning t to different bins, using levels as cutpoints */
    bins = bin(t, levels);

    /* first, assign ECDF[i] for nonmissing bins */
    ecdf = j(nrow(t), 1, .);    /* initialize ECDF to missing */
    idx = loc(bins ^= .);
    if ncol(idx) > 0 then
        ecdf[idx] = y[ bins[idx] ];

    /* Three reasons bins[i] could be a missing value:
       1. if t[i] is missing. No need to change ECDF[i] b/c it was initialized to missing.
       2. if t[i] < min(x). Set ECDF(t[i]) = 0.
       3. if t[i] > max(x). Set ECDF(t[i]) = 1.
    */
    LT_idx = loc(t^=. & t<min(x));
    if ncol(LT_idx) > 0 then 
       ecdf[LT_idx] = 0;
    GT_idx = loc(t^=. & t>=max(x));
    if ncol(GT_idx) > 0 then 
       ecdf[GT_idx] = 1;
    return ecdf;
finish;

/* Compute the Kolmogorov-Smirnov confidence bands for a given ECDF.
   The CL parameter specifies the confidence level, which defaults to 0.95.
   The CL parameter should be in the range [0.75, 0.99]. 
   If x is sorted, a common syntax is 
   y = ECDF(x);
   CL = ECDF_KSCL(y);      *95% CL;
   CL_Lower = CL[,1];
   CL_Upper = CL[,2];
   Then use STEP statement in PROC SGPLOT to plot y, CL_Lower, and CL_Upper agains x.
   See https://blogs.sas.com/content/iml/2026/06/08/confidence-bands-ecdf.html
    */
start ECDF_KSCL(_ECDF, CL=0.95);
   ECDF = colvec(_ECDF);      /* ensure column vector */
   n = countn(ECDF);          /* count the nonmissing values */
   KS_band = j(nrow(ECDF), 2, .);
   /* Enter a pretabulated set of critical values for confidence levels */
   dAlpha = 0.01;
   alpha = do(0.01, 0.25, dAlpha);
   conf_level = round(1 - alpha, dAlpha);
   ks_crit = { 1.6276236 1.517427 1.4490862 1.3985734 1.3580986 
               1.3241093 1.2946745 1.2686236 1.2451911 1.2238479 
               1.2042125 1.1860005 1.1689938 1.1530214 1.1379465 
               1.1236583 1.110065 1.0970904 1.0846701 1.0727492 
               1.0612805 1.0502232 1.0395418 1.0292048 1.0191847 };
   j = loc( conf_level = round(CL, dAlpha) );
   if ncol(j)=0 then do;
      print "ERROR in ECDF_KSCL: The confidence level must be between 0.75 and 0.99";
      return( KS_band );
   end;
   /* Compute the standard error for the K-S bands */
   SE_KS = ks_crit[j] / sqrt(n);

   /* construct upper and lower 95% CI bands centered on the ECDF.
      The >< operator returns the minimum, ensuring max value is 1.
      The <> operator returns the maximum, ensuring min value is 0.
      See https://blogs.sas.com/content/iml/2026/02/04/clip-values.html */
   CL_Upper = (ECDF + SE_KS) >< 1;
   CL_Lower = (ECDF - SE_KS) <> 0;
   KS_band = CL_Lower || CL_Upper;
   return KS_Band;
finish;

/* Comnpute the quantiles for a normal Q-Q plot, along with 
   95% confidence bands. The bands are the transformation of the K-S confidence 
   bands for the ECDF under the inverse CDF (quantile) function. 
   The CL parameter specifies the confidence level, which defaults to 0.95.
   The CL parameter should be in the range [0.75, 0.99].
   See https://blogs.sas.com/content/iml/2026/06/15/confidence-band-qqplot.html
*/
start QQ_KSCL(_x, CL=0.95);
   x = colvec(_x);
   call sort(x);
   n = countn(x);

   /* Compute theoretical plotting positions (Blom, 1958) 
      See https://blogs.sas.com/content/iml/2011/10/28/modeling-the-distribution-of-data-create-a-qq-plot.html */
   v = (T(1:n) - 0.375) / (n + 0.25);  
   q = quantile("Normal", v);

   /* Compute the ECDF and the 95% K-S probability bands */
   y = ECDF(x);
   KS_band = ECDF_KSCL(y);   /* 95% CL */

   /* Transform probability bands to quantile bounds.
      The domain of the QUANTILE function is the open interval (0,1),
      so clip the K-S values.
   */
   F_Lower = choose(KS_band[,1] > 0, KS_band[,1], .);
   F_Upper = choose(KS_band[,2] < 1, KS_band[,2], .);
   /* For a normal Q-Q plot, use the quantile of the normal distribution.
      Use the quantile function for other distributions (eg, "Expoential") to 
      construct confidence intervals for other Q-Q plots. */
   Q_Lower = quantile("Normal", F_Lower);
   Q_Upper = quantile("Normal", F_Upper);

   /* return the 4-column matrix {"x" "q" "Q_Lower" "Q_Upper"} */
   return( x || q || Q_Lower || Q_Upper );
finish;
store module=(ECDF ECDF_KSCL QQ_KSCL);
QUIT;

