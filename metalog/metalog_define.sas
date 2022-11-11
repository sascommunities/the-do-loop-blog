*proc iml;

/* SAS/IML implementation of the metalog ("metalogistic") 
   distribution (Keelin, 2016). 
   http://metalogdistributions.com/images/The_Metalog_Distributions_-_Keelin_2016.pdf
   and
   http://metalogdistributions.com/home.html
*/

/**********************************/
/* HELPER FUNCTIONS AND UTILITIES */
/**********************************/

/* SYSVER macro is supported only in PROC IML, I think.
   In SAS 9.4, this defines a function that emulates the PrintToLog call.
   In SAS Viya, the macro is empty and does nothing.
   After running this macro, call it as follows:
   call PrintToLog("This is a note", 0);
   call PrintToLog("This is a warning", 1);
   call PrintToLog("This is an error", 2);
*/
%macro Define_ML_ErrorFunc;
%if %sysevalf(&sysver = 9.4) %then %do;
start PrintToLog(msg,errCode);
   if errCode=0 then prefix = "NOTE: ";
   else if errCode=1 then prefix = "WARNING: ";
   else prefix = "ERROR: ";
   stmt = '%put ' + prefix + msg + ';';
   call execute(stmt);
finish;
store module=(PrintToLog);
%end;
%mend;
%Define_ML_ErrorFunc;

/* call with STRING: %StopOnError("My error msg"); */
%macro StopOnError(msg);
  do;
    msg = &msg;
    run PrintToLog(msg,2);
    STOP;
  end;
%mend;

/* Test whether a specified order is valid */
start Metalog_Order(_order);
   if type(_order) ^= 'N' then 
      %StopOnError("The order of the metalog model must be a number greater than or equal to 2.");
   order = int(_order);
   if order < 2 then 
      %StopOnError("The order of the metalog model must be a number greater than or equal to 2.");
   if order > 32 then 
      %StopOnError("The order of the metalog model must be less than or equal to 32.");
   return( order );
finish;

/******************************/
/* CONSTRUCTORS */
/******************************/


/*****************************************************/
/* FUNCTIONS THAT APPLY TO ALL METALOG DISTRIBUTIONS 
   The functions are 
   Metalog_Design
   Metalog_IsFeasible
*/
/*****************************************************/

/* When creating a ML distrib from (x,cdf) data, verify data are okay.
   This function OVERWRITES the X and CDF vectors if there are missing values.
   Replaces the vectors with the NONMISSING cases. */
start Metalog_ValidateData(x,y) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_ValidateData";
   if nrow(x)^=nrow(y) then 
      %StopOnError("Invalid data. The X and CDF vectors must be the same size.");
   idx = loc(x^=. & y^=.);
   if ncol(idx)<3 then
      %StopOnError("Invalid data. At least three nonmissing values are required.");
   if ncol(idx) < nrow(x) then do;
      x = x[idx]; y = y[idx];      /* all pairs are now nonmissing */
   end;
   if any(y <= 0 | y>=1) then
      %StopOnError("The input values must be probabilities in (0, 1).");
   return 1;
finish;

/* Create the metalog design matrix (Eqn 8 on p. 254 of Keelin 2016).
   y = cdf in (0,1)
   k = order of model = # of parameters = # of columns in design matrix

   Example syntax:
   y = T(do(0.1, 0.9, 0.1));
   M = Metalog_Design(y, 5);
   print y M;
*/
start Metalog_Design(_y, k) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_Design";
   y = colvec(_y);
   if any(y^=. & (y <= 0 | y>=1)) then 
      %StopOnError("(Metalog_Design) The input values must be probabilities in (0, 1).");
   if k<2 then
      %StopOnError("(Metalog_Design) The order of the metalog model must be greater than or equal to 2");
   M = j(nrow(y), k, 1);
   Ly = quantile('Logistic', y);  /* = LOGIT(y) = log(y/(1-y)) */
   M[,2] = Ly;
   if k=2 then return( M );
   shift = y - 0.5;
   M[,3] = shift # Ly;
   if k=3 then return( M );
   M[,4] = shift;
   if k=4 then return( M );
   do n = 5 to k;
      if mod(n,2) = 1 then do;  /* n is odd */
         pow = (n-1)/2;
         M[,n] = shift##pow;
      end;
      else do;                  /* n is even */
         pow = n/2 - 1;
         M[,n] = shift##pow # Ly;
      end;
   end;
   return( M );
finish;


/* Check the parameters for feasibility.
   Formulas for k=2,3,4 are from http://metalogdistributions.com/equations/feasibility.html

   The parameters must result in a monotone nondecreasing CDF or, 
   equivalently, a PDF for which PDF >= 0.
   The PDF is evaluated at quantiles between 0.001 and 0.999.
   The feasibility condition does not change for semi-, bounded, or unbounded. 

   Example syntax:
   a = {30, 10.5, 9.6, -20.5, -21.7};  * 5th order model;
   isValid = Metalog_IsFeasible(a); 
   print isValid;
*/
start Metalog_IsFeasible(_a) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_IsFeasible";
   a = colvec(_a);
   k = nrow(a);
   if k=2 then return( a[2]>0 );
   if k=3 then do;
      if a[2]<=0 then return( 0 );
      return( abs(a[3]) / a[2] <= 1.66711 );
   end;
   if all(a[2:k]=0) then return( 0 );
   if k=4 then do;
      if      a[2]<0 then return( 0 );
      else if a[2]=0 then return( a[3]=0 & a[4]>0 );
      else do; /* a[2] > 0 */
         a3a2 = abs(a[3]) / a[2];
         a4a2 = a[4] / a[2];
         b=4.5; c=8.5; d=1.93;
         if a4a2 < -4 then return( 0 );
         if -4.0 <= a4a2 & a4a2 <= 4.5 then 
            return( a3a2 <= d/c*sqrt(c##2-(a4a2-b)##2) );
         if  4.5 <  a4a2 & a4a2 <= 7.0 then 
            return( a3a2 <= 0.0216*(a4a2-4.5)+1.930 );
         if  7.0 <  a4a2 & a4a2 <=10.0 then 
            return( a3a2 <= 0.0040*(a4a2-7.0)+1.984 );
         if 10.0 <  a4a2 & a4a2 <=30.0 then 
            return( a3a2 <= 0.0002*(a4a2-10.)+1.996 );
         if 30.0 <  a4a2 then 
            return( a3a2 <= 2.0 );
      end;
   end;
   /* general case. No formula. Evaluate PDF and test if all PDF > 0 */
   dt = 0.001;
   p = 0.0001 // T( do(dt, 1-dt, dt) ) // 0.9999;
   pdf = Metalog_PDF(p, a, 0);   /* do not check feasibility in routine */
   return ( all(pdf>=0) );
finish;


/*****************************************************/
/* FUNCTIONS THAT ARE DIFFERENT FOR 
   UNBOUNDED, SEMIBOUNDED, and BOUNDED distributions. 
   The unbounded function do not have any identifying subset of letters. 
   For example:
   Metalog_ParamEst is for UNBOUNDED distributions whereas
   Metalog_SL_ParamEst is for semibounded (lower),
   Metalog_SU_ParamEst is for semibounded (upper), and
   Metalog_B_ParamEst  is for bounded metalog.

 */
/*****************************************************/

/************************/
/* ECDF                 */
/************************/
/* Returns the empirical cumulative probabilities for a set of data values.
   The i_th data value is assigned a rank, R[i]. 
   Tne rank is transformed into a value in (0,1) by using one of the 
   following methods:

   https://blogs.sas.com/content/iml/2011/10/28/modeling-the-distribution-of-data-create-a-qq-plot.html
   DEFAULT => VW      R[i] / (n + 1)          [van der Waerden]
              L. Makkonen (2008) https://doi.org/10.1080/03610920701653094
              This is the expected value of the order statistic for a uniform distrib.
   OTHER OPTIONS:
              Blom   (R[i] - 3/8) / (n + 1/4)
              Tukey  (R[i] - 1/3) / (n + 1/3)
              Hazen  (R[i] - 1/2) / n

   Example syntax:
   x = {1,2,3,8,3.3,9,2,4,3,2};
   ECDF = Metalog_ECDF(x);
   print x ECDF;
*/
start Metalog_ECDF(_x, Method="VW") global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_ECDF";
   x = colvec(_x);
   R = rank(x);
   N = countn(x);
   uMethod = upcase(Method);
   if      uMethod="VW" then 
      f = R / (N + 1);
   else if uMethod="BLOM" then 
      f = (R - 3/8) / (N + 1/4);
   else if uMethod="TUKEY" then 
      f = (R - 1/3) / (N + 1/3);
   else if uMethod="HAZEN" then 
      f = (R - 1/2) / N;
   else
      f = R / (N + 1);    /* use VW for invalid */
   return f;
finish;

/* Metalog_SL_ECDF is Metalog_ECDF( log(x-bL) ) */
start Metalog_SL_ECDF(_x, bL) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SL_ECDF";
   if any(_x <= bL) then
      %StopOnError("(Metalog_SL_ECDF) All X values must be less than the bL parameter.");
   z = log(_x - bL);
   cdf = Metalog_ECDF(z);
   return( cdf );
finish;

/* Metalog_SU_ECDF is Metalog_ECDF( -log(bU-x) ) */
start Metalog_SU_ECDF(_x, bU) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SU_ECDF";
   if any(_x >= bU) then
      %StopOnError("(Metalog_SU_ECDF) All X values must be greater than the bU parameter.");
   z = -log(bU -_x);
   cdf = Metalog_ECDF(z);
   return( cdf );
finish;

/* Metalog_B_ECDF is Metalog_ECDF(z) for z = log((x-bL)/(bU-x)) */
start Metalog_B_ECDF(_x, B) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_B_ECDF";
   bL = B[1]; bU = B[2];
   if any(_x <= bL | _x>= bU) then
      %StopOnError("(Metalog_B_ECDF) All X values must be in the interval (bL, bU).");
   z = log( (_x - bL)/(bU - _x) );
   cdf = Metalog_ECDF(z);
   return( cdf );
finish;
/************************/



/************************/
/* PARAMETER ESTIMATION */
/************************/

/* Fit parameters of metalog distribution to data. The input vectors are 
   pairs (x, y) where x contains the data values and y contains 
   empirical probabilities or expert estimates. By default, will not
   return estimates if they result in an invalid distribution.

   Example syntax:
   x = {8.0,  20.0, 30.0, 50.0, 85.0};
   y = {0.01, 0.1,  0.5,  0.90, 0.99};
   a = Metalog_ParamEst(x, y, 5); 
   print a;
*/
start Metalog_ParamEst(_x, _y, k) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_ParamEst";
   x = colvec(_x); y = colvec(_y);
   if nrow(x)^=nrow(y) then 
      %StopOnError("(Metalog_ParamEst) Invalid data. The X and Y vectors must be the same size.");
   idx = loc(x^=. & y^=.);
   if ncol(idx)<3 then 
      %StopOnError("(Metalog_ParamEst) Invalid data. At least three nonmissing values are required.");
   if ncol(idx) < nrow(x) then do;
      x = x[idx]; y = y[idx];      /* all pairs are now nonmissing */
   end;
   if any(y <= 0 | y>=1) then 
      %StopOnError("(Metalog_ParamEst) The input values must be probabilities in (0, 1).");

   M = Metalog_Design(y, k);
   /* Solve LS system a = GINV(M`*M)*(M`*x). Test whether M`*M is full rank */
   call appcort(a, numLinDep, M`*M, M`*x); /* find LS soln (uses QR) */
   if numLinDep > 0 then 
      call PrintToLog("(Metalog_ParamEst) The least squares system is singular", 1);
   return( a );
finish;

/* Metalog_SL_ParamEst = Metalog_ParamEst( log(x-bL) ) */
start Metalog_SL_ParamEst(_x, _y, order, _bL) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SL_ParamEst";
   bL = _bL[1];
   if any(_x <= bL) then
      %StopOnError("(Metalog_SL_ParamEst) All X values must be less than the bL parameter.");
   z = log(_x - bL);
   a = Metalog_ParamEst( z, _y, order );
   return( a );
finish;

/* Metalog_SU_ParamEst = Metalog_ParamEst( -log(bU-x) ) */
start Metalog_SU_ParamEst(_x, _y, order, _bU) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SU_ParamEst";
   if nrow(_bU)*ncol(_bU) > 1 then bU = _bU[2];
   else bU = _bU;
   if any(_x >= bU) then
      %StopOnError("(Metalog_SU_ECDF) All X values must be greater than the bU parameter.");
   z = -log(bU-_x);
   a = Metalog_ParamEst( z, _y, order );
   return( a );
finish;

/* Metalog_B_ParamEst = Metalog_ParamEst(z) for z = log((x-bL)/(bU-x)) */
start Metalog_B_ParamEst(_x, _y, order, B) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_B_ParamEst";
   bL = B[1]; bU = B[2];
   if any(_x <= bL | _x>= bU) then
      %StopOnError("(Metalog_B_ParamEst) All X values must be in the interval (bL, bU).");
   z = log( (_x - bL)/(bU - _x) );
   a = Metalog_ParamEst( z, _y, order );
   return( a );
finish;

/************************/
/* QUANTILES            */
/************************/

/* Given probabilities, p, find the quantile values that are associated
   with the metalog distribution with parameters a.

   Example syntax:
   p = {0.01, 0.1, 0.25, 0.5, 0.75, 0.90, 0.99};
   a = {30, 10.5, 9.6, -20.5, -21.7};  * 5th order model;
   x = Metalog_Quantile(p, a); 
   print x y;
*/
start Metalog_Quantile(_p, _a, CheckFeas=0) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_Quantile";
   a = colvec(_a);   p = colvec(_p);
   k = nrow(a);
   if k<2 then 
      %StopOnError("(Metalog_Quantile) Invalid coefficients. At least 2 coefficients are required.");
   if CheckFeas then do;
      isFeas = Metalog_isFeasible(a);
      if ^isFeas then 
         %StopOnError("(Metalog_Quantile) Invalid coefficients. The model is not a valid distribution.");
   end;
   M = Metalog_Design(p, k);
   Q = M*a;
   return( Q );
finish;

/* quantile for SL is bL + exp(M) */
start Metalog_SL_Quantile(_p, _a, _bL) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SL_Quantile";
   bL = _bL[1];
   Mq = Metalog_Quantile(_p, _a);
   Mlog = bL + exp(Mq); 
   return( Mlog );
finish;

/* quantile for SL is bU - exp(-M) */
start Metalog_SU_Quantile(_p, _a, _bU) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SU_Quantile";
   if nrow(_bU)*ncol(_bU) > 1 then bU = _bU[2];
   else bU = _bU;
   Mq = Metalog_Quantile(_p, _a);
   Mnlog = bU - exp(-Mq); 
   return( Mnlog );
finish;

/* quantile for SB is (bL + bU*exp(M)) / (1 + exp(M)) */
start Metalog_B_Quantile(_p, _a, B) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_B_Quantile";
   bL = B[1]; bU = B[2];
   Mq = Metalog_Quantile(_p, _a);
   Mlogit = (bL + bU*exp(Mq)) / (1 + exp(Mq)); 
   return( Mlogit );
finish;


/************************/
/* PDF                  */
/************************/

/* Given probabilities, p, find the quantile values (x) that are associated
   with the metalog distribution with parameters a. Return the 
   density of the metalog distribution at each value of x, as given by
   Eqn 9 on p. 254 of Keelin 2016.

   Example syntax:
   p = {0.01, 0.1, 0.25, 0.5, 0.75, 0.90, 0.99};
   a = {30, 10.5, 9.6, -20.5, -21.7};  * 5th order model;
   PDF = Metalog_PDF(p, a); 
   x = Metalog_Quantile(p, a);   * for reference, get quantile;
   print p x PDF;
*/
start Metalog_PDF(_y, _a, CheckFeas=0) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_PDF";
   a = colvec(_a);
   k = nrow(a);
   if k<2 then 
      %StopOnError("(Metalog_PDF) Invalid coefficients. At least 2 coefficients are required.");
   y = colvec(_y);
   if any(y^=. & (y <= 0 | y>=1)) then 
      %StopOnError("(Metalog_PDF) The input values must be probabilities in (0, 1).");
   y1my = y#(1-y);         /* y1my is short for y(1 MINUS y) */
   m = a[2] / y1my;
   if k>2 then do;
      Ly = quantile('Logistic', y);  /* = LOGIT(y) = log(y/(1-y)) */
      shift = y - 0.5;
      m = m + a[3]*(shift/y1my + Ly);
   end;
   if k>3 then do;
      m = m + a[4];
   end;

   do n = 5 to k;
      if mod(n,2) = 1 then do;  /* n is odd */
         pow = (n-3)/2;
         m = m + a[n]*(n-1)/2 * shift##pow;
      end;
      else do;                  /* n is even */
         pow = n/2 - 1;
         m = m + a[n]*(shift##pow / y1my + 
                 pow * shift##(pow-1) # Ly);
      end;
   end;
   bFeas = any(m<=0);
   if CheckFeas & bFeas then 
      %StopOnError("(Metalog_PDF) Invalid coefficients. The model is not a valid distribution.");
   /* if we don't check feasibility, protect against m=0 */
   if bFeas then do;
      mm = j(nrow(m), 1, .);
      idx = loc(m ^= 0);
      if ncol(idx)>0 then 
         mm[idx] = 1/m[idx];
      return ( mm );
   end;
   return( 1/m );
finish;


/* Metalog_SL_PDF = m # exp( -M ) where m is PDF and M is quantile */
start Metalog_SL_PDF(_p, _a, _bL) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SL_PDF";
   bL = _bL[1];
   Mq = Metalog_Quantile(_p, _a);
   m = Metalog_PDF(_p, _a);
   mPDFlog = m # exp( -Mq );
   return( mPDFlog );
finish;

/* Metalog_SU_PDF = m # exp( M ) where m is PDF and M is quantile */
start Metalog_SU_PDF(_p, _a, _bU) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SU_PDF";
   if nrow(_bU)*ncol(_bU) > 1 then bU = _bU[2];
   else bU = _bU;
   Mq = Metalog_Quantile(_p, _a);
   m = Metalog_PDF(_p, _a);
   mPDFlog = m # exp(Mq);
   return( mPDFlog );
finish;

/* Metalog_B_PDF = m # L( M ) where m is PDF and M is quantile 
   and L(M) = (1 + eM)##2 / ((bU-bL)*eM) where eM = exp(M) */
start Metalog_B_PDF(_p, _a, B) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_B_PDF";
   bL = B[1]; bU = B[2];
   Mq = Metalog_Quantile(_p, _a);
   m = Metalog_PDF(_p, _a);
   eM = exp(Mq);
   mlogit = m # (1 + eM)##2 / ((bU-bL)*eM);
   return( mlogit );
finish;


/************************/
/* RAND                 */
/************************/

/* Given parameters, a, for the metalog distribution, return a
   random sample of size n from the distribution.

   Example syntax:
   a = {30, 10.5, 9.6, -20.5, -21.7};  * 5th order model;
   x = Metalog_Rand(10, a); 
   print x;
*/
start Metalog_Rand(n, _a) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_Rand";
   if n<1 then 
      %StopOnError("(Metalog_Rand) Invalid size. The sample size must be a positive integer.");
   u = randfun(n, "uniform");
   a = colvec(_a);
   x = Metalog_Quantile(u, a);
   return( x );
finish;

start Metalog_SL_Rand(n, _a, _bL) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SL_Rand";
   bL = _bL[1];
   u = randfun(n, "uniform");
   x = Metalog_SL_Quantile(u, _a, bL);
   return( x );
finish;

start Metalog_SU_Rand(n, _a, _bU) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SU_Rand";
   if nrow(_bU)*ncol(_bU) > 1 then bU = _bU[2];
   else bU = _bU;
   u = randfun(n, "uniform");
   x = Metalog_SU_Quantile(u, _a, bU);
   return( x );
finish;

start Metalog_B_Rand(n, _a, B) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_B_Rand";
   u = randfun(n, "uniform");
   x = Metalog_B_Quantile(u, _a, B);
   return( x );
finish;


/**************************************************/
/* PDFX and CDFX: Given x, find PDF(x) and CDF(x) */
/**************************************************/

/* Given parameters, a, for the metalog distribution, return the
   density as a function of data values, x.

   Example syntax:
   x = {10, 20, 25, 30, 40, 50, 85};
   a = {30, 10.5, 9.6, -20.5, -21.7};  * 5th order model;
   PDF = Metalog_PDFX(x, a); 
   print x PDF;
*/
start Metalog_PDFX(_x, _a) global(g_MLCoeff, g_MLTargetQntl);
   x = colvec(_x);   a = colvec(_a);
   y = Metalog_CDFX(x, a);
   PDF = Metalog_PDF(y, a);
   if any(PDF<=0) then 
      %StopOnError("(Metalog_PDFX) Invalid coefficients. The model is not a valid distribution.");
   return( PDF );
finish;


/* Metalog_QuantileRoot is a PRIVATE helper function */
start Metalog_QuantileRoot(y) global(g_MLCoeff, g_MLTargetQntl);
   return( Metalog_Quantile(y, g_MLCoeff) - g_MLTargetQntl );
finish;

/* There is no formula to get CDF(x) as a function of x and a.
   However, we can solve the implicit equation by solving for 
   y such that y = CDF(x). */

/* Given parameters, a, for the metalog distribution, return the
   CDF as a function of data values, x.

   Example syntax:
   x = {10, 20, 25, 30, 40, 50, 85};
   a = {30, 10.5, 9.6, -20.5, -21.7};  * 5th order model;
   CDF = Metalog_CDFX(x, a); 
   print x CDF;
*/
start Metalog_CDFX(_x, _a) global(g_MLCoeff, g_MLTargetQntl);
   a = colvec(_a);
   if ^Metalog_IsFeasible(a) then
      %StopOnError("(Metalog_CDFX) Invalid coefficients. The model is not a valid distribution.");
   x = colvec(_x);
   g_MLCoeff = a;
   /* for efficiency, presearch and use smaller intervals */
   yGrid = {0.00001,0.001} // T(do(0.1, 0.99, 0.01)) // {0.999,0.99999};
   eps = 1E-12 / 3;
   yGrid = yGrid + eps; /* to help with robustness, perturb the yGrid values */
   xGrid = Metalog_Quantile(yGrid, a);
   /* evaluate CDF(x[i]) for each x[i] */
   b = bin(x, xGrid);
   y = j(nrow(x), 1, .);
   do i = 1 to nrow(x);
      if b[i]^=. then do;
         g_MLTargetQntl = x[i];
         /* decrease the LHS of the interval
            so that the root is in the interior of the interval */
         if b[i]=1 then 
           interval = yGrid[1]/2 || yGrid[2];
         else   /* b[i] > 1 so use next lowest and next highest bin boundary */
           interval = yGrid[b[i]-1] || yGrid[b[i]+1];
         *PRINT i (x[i])[L='x_i'] (b[i])[L='b_i'] interval;
         y[i] = froot("Metalog_QuantileRoot", interval );
      end;
   end;
   return( y );
finish;

%let Metalog_modules =
Metalog_ValidateData
Metalog_Order
Metalog_Design
Metalog_IsFeasible
Metalog_ECDF
Metalog_SL_ECDF
Metalog_SU_ECDF
Metalog_B_ECDF
Metalog_ParamEst
Metalog_SL_ParamEst
Metalog_SU_ParamEst
Metalog_B_ParamEst
Metalog_Quantile
Metalog_SL_Quantile
Metalog_SU_Quantile
Metalog_B_Quantile
Metalog_PDF
Metalog_SL_PDF
Metalog_SU_PDF
Metalog_B_PDF
Metalog_Rand
Metalog_SL_Rand
Metalog_SU_Rand
Metalog_B_Rand
Metalog_PDFX
Metalog_QuantileRoot
Metalog_CDFX
;

store module = (&Metalog_modules);

