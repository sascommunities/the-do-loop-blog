/*
Copyright © 2022, SAS Institute Inc., Cary, NC, USA.  All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
*/

/*********************************************/
/* High-level public API for Metalog package */
/*********************************************/

/* SAS/IML implementation of the metalog ("metalogistic") distribution 
   http://metalogdistributions.com/images/The_Metalog_Distributions_-_Keelin_2016.pdf
   and
   http://metalogdistributions.com/home.html

   This implementation uses an object-oriented approach. The metalog
   distribution is created as an object (technically, a named list)
   by using on of the following calls:
      ML_CreateFromData(x <,order> <,bounds>)
      ML_CreateFromCDF(x,cdf <,order> <,bounds>)
      ML_CreateFromCoef(coef <,bounds>)
   These constructors support the following optional parameters
      order  : (Default=5) Number of terms in metalog
      bounds : (Default={. .}) bounds for bounded or semibounded distrib
   They return a named list, L, which includes the following information
      L$'x'          : original data (size N) or missing
      L$'cdf'        : ECDF or user-specified (expert)
      L$'order'      : number of terms in metalog
      L$'bounds'     : Lower and upper bounds for bounded or semibounded distrib
      L$'boundType'  : type of ML distribution: 'U', 'SL', 'SU', or 'B'
      L$'a'          : vector (size=order) of metalog coefficients
      L$'isFeasible' : 0/1 flag to indicate feasibility of coefficients

   Naming convention: 
      - The ML_Create* methods return a named list.
      - Other methods with ML_ prefix take a named list as the input argument.
      - Low-level computational routines begin with Metalog_ and 
        take a series of matrix/vector/scalar arguments. The low-level functions
        are private and not part of the API.
*/

/******************************/
/* CONSTRUCTORS               */
/******************************/

/* Create a metalog object from (x,cdf) pairs, where x is a vector of 
   data values and cdf is a vector of the associate cumulative probabilities.
   Example syntax:
   x = {14,18,22.8,24.6,26.1,31,38,41};
   cdf = {0.1,0.25,0.4,0.5,0.6,0.75,0.9,0.95};
   L = ML_CreateFromCDF(x, cdf);   *default order and unbounded;
   order = 4;
   bounds = {0 .};
   L2 = ML_CreateFromCDF(x, cdf, order, bounds);   
*/
start ML_CreateFromCDF(_x, _cdf, _order=5, _bounds={. .});
   x = colvec(_x);
   cdf = colvec(_cdf);
   isValid = Metalog_ValidateData(x,cdf);  /* deletes missing values in x,cdf */
   order = Metalog_Order(_order);
   /* Ex: if symmetric-percentile triplet, must use order<= 3 */
   if nrow(x) < order then do;
      if ^IsSkipped(_order) then 
         run PrintToLog("(ML_CreateFromCDF): Invalid order. Using number of valid observations.",1);
      order = nrow(x);
   end;
   L = [#'x' = x,
        #'cdf' = cdf,
        #'order' = order,
        #'bounds'= ML_Bounds(_bounds),
        #'boundType' = ML_BoundType(_bounds),
        #'a'     = .,
        #'isFeasible' = .
       ];
   L$'a' = ML_ParamEst(L);
   L$'isFeasible' = Metalog_IsFeasible(L$'a');
   return( L );
finish;

/* Create a metalog object from x, which is a vector of 
   data values. The empirical CDF is used to fit the metalog distribution.
   Tied values are ranked arbitrarily. You can specify one of the 
   following scores to convert ranks to probabilities:
     VW:  The van der Waerden score is R[i] / (n+1). This is the default.
     Blom: Blom's score is (R[i] - 3/8) / (n+1/4).
     Tukey: The Tukey score is (R[i] - 1/3) / (n + 1/3)
     Hazen: The Haven score is (R[i] - 1/2) / n

   Example syntax:
   x = {14,18,22.8,24.6,26.1,31,38,41};
   L = ML_CreateFromData(x);   *default order and unbounded;
   order = 4;
   bounds = {0 .};
   method="Blom"
   L2 = ML_CreateFromData(x, order, bounds, method);
*/
start ML_CreateFromData(_x, order=5, bounds={. .}, method="VW");
   ECDF = Metalog_ECDF(_x, method);
   return ML_CreateFromCDF(_x, ECDF, order, bounds);
finish;

/* Create a metalog object from coefficients. The vector of 
   coefficients determines the order. Not all coefficients lead to 
   valid distributions.

   Example syntax:
   a   = {20, 1, 1, 2};
   L = ML_CreateFromCoef(a);

   bounds = {0 .};
   L2 = ML_CreateFromCoef(x, bounds);
*/
start ML_CreateFromCoef(_a, _bounds={. .});
   a = colvec(_a);
   b = Metalog_IsFeasible(a);
   L = [#'x' = .,
        #'cdf' = .,
        #'order' = nrow(a),
        #'bounds'= ML_Bounds(_bounds),
        #'boundType' = ML_BoundType(_bounds),
        #'a'     = a,
        #'isFeasible' = b
       ];
   return( L );
finish;

/**********************************/
/* HELPER FUNCTIONS AND UTILITIES */
/**********************************/

/* Print a summary of the metalog object */
start ML_Summary(L);
   if ^ML_ValidateObject(L) then return;
   Summary = strip(char(L$'order')) //
             L$'boundType'          //
             ML_BoundString(L)      //
             char(L$'isFeasible',1);
   Parameter = T('a1':('a'+strip(char(nrow(L$'a')))));
   Estimate = L$'a';
   print Summary[r={'Order' 'Type' 'Bounds' 'Is Feasible'} L="Model Summary"];
   print Estimate[r=Parameter];
finish;

/* Validate order in the range 2 < order <= 32; return order */
start ML_Order(L);
   if type(L)='L' then do;
      if ^ML_ValidateObject(L) then return(.);
      order=L$'order';
   end;
   else if type(L)='N' then order=L;
   else return( . );
   order = Metalog_Order(order);
   return( order );
finish;

/* Validate bounds; return bounds */
start ML_Bounds(L);
   if type(L)='L' then do;
      if ^ML_ValidateObject(L) then return({. .});
      bounds=L$'bounds';
   end;
   else if type(L)='N' then bounds=L;
   else return( {. .} );
   bounds = rowvec(bounds);
   if type(bounds) ^= 'N' | ncol(bounds) ^= 2 then do;
     msg = "The bounds for the metalog model must be a two-element numerical vector.";
     run PrintToLog(msg,2);
     STOP;
   end;
   if bounds[1]^=. & bounds[2]^=. then 
      if bounds[1] >= bounds[2] then do;
        msg = "The lower bound must be less than the upper bound.";
        run PrintToLog(msg,2);
        STOP;
      end;
   return( bounds );
finish;

/* return type of model based on the bounds:
   'U' = unbounded; 'B' = bounded
   'SU' = semi-bounded with upper bound; 'SL' = semi-bounded with lower bound
*/ 
start ML_BoundType(L);
   if type(L)='L' then do;
      if ^ML_ValidateObject(L) then return(' ');
      bounds=L$'bounds';
   end;
   else if type(L)='N' then bounds=L;
   else return( ' ' );
   return( Metalog_BoundType(bounds) );
finish;

/* Return the parameter estimates for the ML model.
   This is a column vectors with L$'order' elements */
start ML_Coef(L);
   if type(L)^='L' then return(.);
   if ^ML_ValidateObject(L) then return(.);
   return( L$'a' );
finish;

/* the feasibility conditions are the same regardless of type: 
   the PDF of the model must be strictly positive (except at the 
   endpoints of a semi-bounded or bounded domain, where the PDF is 0) */
start ML_IsFeasible(L);
   if type(L)^='L' then return(.);
   if ^ML_ValidateObject(L) then return(.);
   return( L$'isFeasible' );
finish;

/* INTERNAL HELPERS */
/* Test whether a named list object is a valid ML object */
start ML_ValidateObject(L);
   if type(L)^='L' then return( 0 );
   listNames = ListGetAllNames(L);
   obsNames = unique(upcase(listNames));
   targetNames = unique({x cdf order bounds boundType a isFeasible});
   if ncol(obsNames) ^= ncol(targetNames) then return(0);
   if any(obsNames ^= targetNames) then return( 0 );
   return( 1 );
finish;

start ML_BoundString(L);
   if type(L)='L' then do;
      if ^ML_ValidateObject(L) then return(' ');
      bounds=L$'bounds';
   end;
   else if type(L)='N' then bounds=L;
   else return(' ');
   LB = putn(bounds[1], "BEST9.");
   UB = putn(bounds[2], "BEST9.");
   return( cats('[',LB,',',UB,']') );
finish;

/******************************/
/* COMPUTATIONAL ROUTINES     */
/******************************/

/* Given a vector of cumulative probabilities, p, return the corresponding
   quantiles (x) for the metalog model

   Example syntax:
   x = {14,18,22.8,24.6,26.1,31,38,41};
   L = ML_CreateFromData(x);   
   p = {0.01, 0.1, 0.25, 0.5, 0.75, 0.90, 0.99};
   Q = ML_Quantile(L, p);
*/
start ML_Quantile(L, p);
   type = ML_BoundType(L);
   if type='B' then
      M = Metalog_B_Quantile(p, L$'a', L$'bounds');
   else if type='SL' then
      M = Metalog_SL_Quantile(p, L$'a', L$'bounds');
   else if type='SU' then
      M = Metalog_SU_Quantile(p, L$'a', L$'bounds');
   else  /* default is unbounded */
      M = Metalog_Quantile(p, L$'a');
   return(M);
finish;

/* Given a vector of cumulative probabilities, p, return the corresponding
   density (PDF) for the metalog model

   Example syntax:
   x = {14,18,22.8,24.6,26.1,31,38,41};
   L = ML_CreateFromData(x);   
   p = {0.01, 0.1, 0.25, 0.5, 0.75, 0.90, 0.99};
   PDF = ML_PDF(L, p);
   Q = ML_Quantile(L, p);   * get quantiles for reference;
   print p Q PDF;
*/
start ML_PDF(L, p);
   type = ML_BoundType(L);
   if type='B' then
      f = Metalog_B_PDF(p, L$'a', L$'bounds');
   else if type='SL' then
      f = Metalog_SL_PDF(p, L$'a', L$'bounds');
   else if type='SU' then
      f = Metalog_SU_PDF(p, L$'a', L$'bounds');
   else  /* default is unbounded */
      f = Metalog_PDF(p, L$'a');
   return(f);
finish;

/* Generate n random variates from the metalog distribution */
start ML_Rand(L, n);
   type = ML_BoundType(L);
   if type='B' then
      X = Metalog_B_Rand(n, L$'a', L$'bounds');
   else if type='SL' then
      X = Metalog_SL_Rand(n, L$'a', L$'bounds');
   else if type='SU' then
      X = Metalog_SU_Rand(n, L$'a', L$'bounds');
   else  /* default is unbounded */
      X = Metalog_Rand(n, L$'a');
   return(X);
finish;

/* Fit the metalog model to the data */
start ML_ParamEst(L);
   if all(L$'x' = .) | all(L$'cdf' = .) then 
      return( j(L$'order', 1, .) );
   type = ML_BoundType(L);
   if type='B' then
      a = Metalog_B_ParamEst(L$'x', L$'cdf', L$'order', L$'bounds');
   else if type='SL' then
      a = Metalog_SL_ParamEst(L$'x', L$'cdf', L$'order', L$'bounds');
   else if type='SU' then
      a = Metalog_SU_ParamEst(L$'x', L$'cdf', L$'order', L$'bounds');
   else  /* default is unbounded */
      a = Metalog_ParamEst(L$'x', L$'cdf', L$'order');
   return(a);
finish;


store module=(
ML_ValidateObject
ML_Order
ML_BoundString
ML_Bounds
ML_BoundType
ML_Coef
ML_IsFeasible
ML_CreateFromCDF
ML_CreateFromData
ML_CreateFromCoef
ML_Summary
ML_ParamEst
ML_Quantile
ML_PDF
ML_Rand
);

/***********************************************/
/***********************************************/
/***********************************************/


/***************************************************************/
/* Low-level private computational methods for Metalog package */
/***************************************************************/

/* SAS/IML computations for the metalog ("metalogistic") 
   distribution (Keelin, 2016). 
   http://metalogdistributions.com/images/The_Metalog_Distributions_-_Keelin_2016.pdf
   and
   http://metalogdistributions.com/home.html
*/

/**********************************/
/* HELPER FUNCTIONS AND UTILITIES */
/**********************************/

/* accurate computation of logit(y) for column vector y
   LOGIT(y) = log(y/(1-y)) = QUANTILE('Logistic',y) */
start Metalog_Logit(y);
   if any(y>0.9999) then do; /* use SQUANTILE for extreme right probabilities */
      logit = j(nrow(y),1,.);
      idx = loc(y>0.9);
      if ncol(idx)>0 then
         logit[idx] = squantile('Logistic', 1-y[idx]);
      idx = loc(y<=0.9);
      if ncol(idx)>0 then
         logit[idx] = quantile('Logistic', y[idx]);
   end;
   else
      logit = quantile('Logistic', y);
   return( logit );
finish;

/* Test whether a specified order is valid */
start Metalog_Order(_order);
  if type(_order) ^= 'N' then do;
     msg = "The order of the metalog model must be a number greater than or equal to 2.";
     run PrintToLog(msg,2);
     STOP;
   end;
   order = int(_order);
   if order < 2 then do;
     msg = "The order of the metalog model must be a number greater than or equal to 2.";
     run PrintToLog(msg,2);
     STOP;
   end;
   if order > 32 then do;
     msg = "The order of the metalog model must be less than or equal to 32.";
     run PrintToLog(msg,2);
     STOP;
   end;
   return( order );
finish;


/* return type of model based on the bounds:
   'U' = unbounded; 'B' = bounded
   'SU' = semi-bounded with upper bound; 'SL' = semi-bounded with lower bound
*/ 
start Metalog_BoundType(bounds);
   Lbound = (bounds[1]^=.);
   Ubound = (bounds[2]^=.);
   if  Lbound &  Ubound then return( 'B' );
   if ^Lbound &  Ubound then return( 'SU' );
   if  Lbound & ^Ubound then return( 'SL' );
   return( 'U' );
finish;


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
   if nrow(x)^=nrow(y) then do;
     msg = "Invalid data. The X and CDF vectors must be the same size.";
     run PrintToLog(msg,2);
     STOP;
   end;
   idx = loc(x^=. & y^=.);
   if ncol(idx)<3 then do;
     msg = "Invalid data. At least three nonmissing values are required.";
     run PrintToLog(msg,2);
     STOP;
   end;
   if ncol(idx) < nrow(x) then do;
      x = x[idx]; y = y[idx];      /* all pairs are now nonmissing */
   end;
   if any(y <= 0 | y>=1) then do;
     msg = "The input values must be probabilities in (0, 1).";
     run PrintToLog(msg,2);
     STOP;
   end;
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
   free missIdx;
   if k<2 then do;
     msg = "(Metalog_Design) The order of the metalog model must be greater than or equal to 2.";
     run PrintToLog(msg,2);
     STOP;
   end;
   if any(y <= 0 | y>=1) then do;
      missIdx = loc(y <= 0 | y>=1);
      if ncol(missIdx)>0 then 
         y[missIdx] = .;
   end;
   M = j(nrow(y), k, 1);
   if ncol(missIdx)>0 then M[missIdx,1] = .;
   Ly = Metalog_Logit(y);  /* = LOGIT(y) = log(y/(1-y)) */
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
start Metalog_SL_ECDF(_x, _bL) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SL_ECDF";
   bL = _bL[1];
   x = colvec(_x);
   /* CDF is 0 for all quantiles less than bL */
   if any(x <= bL) then do;
      cdf = j(nrow(x),1,0);
      idx = loc(x > bL);
      if ncol(idx)=0 then 
         return(j(nrow(x),1,0));
      z = log(x[idx] - bL);
      cdf[idx] = Metalog_ECDF(z);
   end;
   else do;
      z = log(x - bL);
      cdf = Metalog_ECDF(z);
   end;
   return( cdf );
finish;

/* Metalog_SU_ECDF is Metalog_ECDF( -log(bU-x) ) */
start Metalog_SU_ECDF(_x, _bU) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SU_ECDF";
   bU = _bU[1];
   x = colvec(_x);
   /* CDF is 1 for all quantiles greater than bU */
   if any(x >= bU) then do;
      cdf = j(nrow(x),1,1);
      idx = loc(x < bU);
      if ncol(idx)=0 then 
         return(j(nrow(x),1,1));
      z = -log(bU - x[idx]);
      cdf[idx] = Metalog_ECDF(z);
   end;
   else do;
      z = -log(bU - x);
      cdf = Metalog_ECDF(z);
   end;
   return( cdf );
finish;

/* Metalog_B_ECDF is Metalog_ECDF(z) for z = log((x-bL)/(bU-x)) */
start Metalog_B_ECDF(_x, B) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_B_ECDF";
   bL = B[1]; bU = B[2];
   x = colvec(_x);
   /* adjust CDF for quantiles outside of [bL,bU] */
   if any(x <= bL | x>= bU) then do;
      cdf = j(nrow(x),1,0);
      idx = loc(x >= bU);
      if ncol(idx)>0 then 
         cdf[idx] = 1;
      idx = loc(x > bL & x < bU);
      if ncol(idx)>0 then do;
         z = log( (x[idx] - bL)/(bU - x[idx]) );
         cdf[idx] = Metalog_ECDF(z);
      end;
   end;
   else do;
      z = log( (x - bL)/(bU - x) );
      cdf = Metalog_ECDF(z);
   end;
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
   if nrow(x)^=nrow(y) then do;
     msg = "(Metalog_ParamEst) Invalid data. The X and Y vectors must be the same size.";
     run PrintToLog(msg,2);
     STOP;
   end;
   idx = loc(x^=. & y>0 & y<1);
   if ncol(idx)<3 then do;
     msg = "(Metalog_ParamEst) You must specify at least three valid probabilities.";
     run PrintToLog(msg,1);
     return(j(k,1,.));
   end;
   if ncol(idx) < nrow(x) then do;
      x = x[idx]; y = y[idx];      /* all pairs are now valid */
   end;
 
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
   x = colvec(_x);
   y = colvec(_y);
   if any(x <= bL) then do;
      run PrintToLog("(Metalog_SL_ParamEst) Excluding values less than or equal to the lower bound.",0);
      idx = loc(x > bL);
      if ncol(idx)=0 then 
         return(j(order,1,.));
      x = x[idx];
      y = y[idx];
   end;
   z = log(x - bL);
   a = Metalog_ParamEst( z, y, order );
   return( a );
finish;

/* Metalog_SU_ParamEst = Metalog_ParamEst( -log(bU-x) ) */
start Metalog_SU_ParamEst(_x, _y, order, _bU) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SU_ParamEst";
   if nrow(_bU)*ncol(_bU) > 1 then bU = _bU[2];
   else bU = _bU;
   x = colvec(_x);
   y = colvec(_y);
   if any(x >= bU) then do;
      run PrintToLog("(Metalog_SU_ParamEst) Excluding values greater than or equal to the upper bound.",0);
      idx = loc(x < bU);
      if ncol(idx)=0 then 
         return(j(order,1,.));
      x = x[idx];
      y = y[idx];
   end;
   z = -log(bU-x);
   a = Metalog_ParamEst( z, y, order );
   return( a );
finish;

/* Metalog_B_ParamEst = Metalog_ParamEst(z) for z = log((x-bL)/(bU-x)) */
start Metalog_B_ParamEst(_x, _y, order, B) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_B_ParamEst";
   bL = B[1]; bU = B[2];
   x = colvec(_x);
   y = colvec(_y);
   if any(x <= bL | x>= bU) then do;
      run PrintToLog("(Metalog_B_ParamEst) Excluding values that are not strictly inside the bounding interval.",0);
      idx = loc(x > bL & x < bU);
      if ncol(idx)=0 then 
         return(j(order,1,.));
      x = x[idx];
      y = y[idx];
   end;
   z = log( (x - bL)/(bU - x) );
   a = Metalog_ParamEst( z, y, order );
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
   if k<2 then do;
     msg = "(Metalog_Quantile) Invalid coefficients. At least 2 coefficients are required.";
     run PrintToLog(msg,2);
     STOP;
   end;
   if CheckFeas then do;
      isFeas = Metalog_isFeasible(a);
      if ^isFeas then do;
        msg = "(Metalog_Quantile) Invalid coefficients. The model is not a valid distribution.";
        run PrintToLog(msg,2);
        STOP;
      end;
   end;
   M = Metalog_Design(p, k);
   /* M could contain rows of missing values if p contains invalid elements. */
   if all(M ^= .) then
      return(M*a);
   /* handle rows of missing values */
   Q = j(nrow(p), 1, .);
   nonMissIdx = loc(M[,1] ^= .);
   if ncol(NonMissIdx)=0 then return(Q);
   Q[nonMissIdx] = M[nonMissIdx,]*a;
   return( Q );
finish;

/* quantile for SL is bL + exp(M) */
start Metalog_SL_Quantile(_p, a, _bL) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SL_Quantile";
   bL = _bL[1];
   p = colvec(_p);
   cutoff = 1E-14;
   if any(p <= cutoff) then do;   /* for p=0 */
      Mlog = j(nrow(p),1,bL);
      idx = loc(p > cutoff);
      if ncol(idx)>0 then do;
         Mq = Metalog_Quantile(p[idx], a);
         Mlog[idx] = bL + exp(Mq); 
      end;
   end;
   else do;
      Mq = Metalog_Quantile(p, a);
      Mlog = bL + exp(Mq); 
   end;
   return( Mlog );
finish;

/* quantile for SL is bU - exp(-M) */
start Metalog_SU_Quantile(_p, a, _bU) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SU_Quantile";
   if nrow(_bU)*ncol(_bU) > 1 then bU = _bU[2];
   else bU = _bU;
   p = colvec(_p);
   cutoff = 1 - 1E-14;
   if any(p >= cutoff) then do;   /* for p=1 */
      Mnlog = j(nrow(p),1,bU);
      idx = loc(p < cutoff);
      if ncol(idx)>0 then do;
         Mq = Metalog_Quantile(p[idx], a);
         Mnlog[idx] = bU - exp(-Mq); 
      end;
   end;
   else do;
      Mq = Metalog_Quantile(_p, a);
      Mnlog = bU - exp(-Mq); 
   end;
   return( Mnlog );
finish;

/* quantile for SB is (bL + bU*exp(M)) / (1 + exp(M)) */
start Metalog_B_Quantile(_p, a, B) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_B_Quantile";
   bL = B[1]; bU = B[2];
   p = colvec(_p);
   eps = 1E-14;
   if any(p <= eps | p >= 1-eps) then do;   /* for p=0 or p=1 */
      Mlogit = j(nrow(p),1,bL);
      idx = loc(p >= 1-eps);
      if ncol(idx)>0 then 
         Mlogit[idx] = bU;
      idx = loc(p > eps & p < 1-eps);
      if ncol(idx)>0 then do;
         Mq = Metalog_Quantile(p[idx], a);
         Mlogit[idx] = (bL + bU*exp(Mq)) / (1 + exp(Mq)); 
      end;
   end;
   else do;
      Mq = Metalog_Quantile(_p, a);
      Mlogit = (bL + bU*exp(Mq)) / (1 + exp(Mq)); 
   end;
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
   if k<2 then do;
     msg = "(Metalog_PDF) Invalid coefficients. At least 2 coefficients are required.";
     run PrintToLog(msg,2);
     STOP;
   end;
   y = colvec(_y);
   if any(y^=. & (y <= 0 | y>=1)) then do;
     msg = "(Metalog_PDF) The input values must be probabilities in (0, 1).";
     run PrintToLog(msg,2);
     STOP;
   end;
   y1my = y#(1-y);         /* y1my is short for y(1 MINUS y) */
   m = a[2] / y1my;
   if k>2 then do;
      Ly = Metalog_Logit(y);  /* = LOGIT(y) = log(y/(1-y)) */
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
   if CheckFeas & bFeas then do;
     msg = "(Metalog_PDF) Invalid coefficients. The model is not a valid distribution.";
     run PrintToLog(msg,2);
     STOP;
   end;
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
start Metalog_SL_PDF(_p, a, _bL) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SL_PDF";
   bL = _bL[1];
   cutoff = 1E-14;
   if any(_p <= cutoff | _p >= 1-cutoff) then do;   /* for p=0 */
      MPDFlog = j(nrow(_p)*ncol(_p),1,.);
      idx = loc(_p <= cutoff);
      if ncol(idx)>0 then 
      	 MPDFlog[idx] = 0;
      idx = loc(_p > cutoff & _p < 1-cutoff);
      if ncol(idx)>0 then do;
          p = _p[idx];
          Mq = Metalog_Quantile(p, a);
          m = Metalog_PDF(p, a);
          mPDFlog[idx] = m # exp( -Mq );
      end;
   end;
   else do;
      Mq = Metalog_Quantile(_p, a);
      m = Metalog_PDF(_p, a);
      mPDFlog = m # exp( -Mq );
   end;
   return( mPDFlog );
finish;

/* Metalog_SU_PDF = m # exp( M ) where m is PDF and M is quantile */
start Metalog_SU_PDF(_p, a, _bU) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SU_PDF";
   if nrow(_bU)*ncol(_bU) > 1 then bU = _bU[2];
   else bU = _bU;
   cutoff = 1E-14;
   if any(_p <= cutoff | _p >= 1-cutoff) then do;   /* for p=1 */
      MPDFlog = j(nrow(_p)*ncol(_p),1,.);
      idx = loc(_p >= 1-cutoff);
      if ncol(idx)>0 then 
      	 MPDFlog[idx] = 0;
      idx = loc(_p > cutoff & _p < 1-cutoff);
      if ncol(idx)>0 then do;
         p = _p[idx];
         Mq = Metalog_Quantile(p, a);
         m = Metalog_PDF(p, a);
         mPDFlog[idx] = m # exp(Mq);
      end;
   end;
   else do;
      Mq = Metalog_Quantile(_p, a);
      m = Metalog_PDF(_p, a);
      mPDFlog = m # exp(Mq);
   end;
   return( mPDFlog );
finish;

/* Metalog_B_PDF = m # L( M ) where m is PDF and M is quantile 
   and L(M) = (1 + eM)##2 / ((bU-bL)*eM) where eM = exp(M) */
start Metalog_B_PDF(_p, a, B) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_B_PDF";
   bL = B[1]; bU = B[2];
   eps = 1E-14;
   if any(_p <= eps | _p >= 1-eps) then do;   /* for p=0 or p=1 */
      mlogit = j(nrow(_p)*ncol(_p), 1, 0);
      idx = loc(_p > eps & _p < 1-eps);
      if ncol(idx)>0 then do;
         p = _p[idx];
         Mq = Metalog_Quantile(p, a);
         m = Metalog_PDF(p, a);
         eM = exp(Mq);
         mlogit[idx] = m # (1 + eM)##2 / ((bU-bL)*eM);
      end;
   end;
   else do;
      Mq = Metalog_Quantile(_p, a);
      m = Metalog_PDF(_p, a);
      eM = exp(Mq);
      mlogit = m # (1 + eM)##2 / ((bU-bL)*eM);
   end;
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
   if n<1 then do;
     msg = "(Metalog_Rand) Invalid size. The sample size must be a positive integer.";
     run PrintToLog(msg,2);
     STOP;
   end;
   u = randfun(n, "uniform");
   a = colvec(_a);
   x = Metalog_Quantile(u, a);
   return( x );
finish;

start Metalog_SL_Rand(n, a, _bL) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SL_Rand";
   bL = _bL[1];
   u = randfun(n, "uniform");
   x = Metalog_SL_Quantile(u, a, bL);
   return( x );
finish;

start Metalog_SU_Rand(n, a, _bU) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_SU_Rand";
   if nrow(_bU)*ncol(_bU) > 1 then bU = _bU[2];
   else bU = _bU;
   u = randfun(n, "uniform");
   x = Metalog_SU_Quantile(u, a, bU);
   return( x );
finish;

start Metalog_B_Rand(n, a, B) global(_debug);
   if type(_debug)='N' then if _debug>0 then print "In Metalog_B_Rand";
   u = randfun(n, "uniform");
   x = Metalog_B_Quantile(u, a, B);
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
   if any(PDF<=0) then do;
     msg = "(Metalog_PDFX) Invalid coefficients. The model is not a valid distribution.";
     run PrintToLog(msg,2);
     STOP;
   end;
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
   if ^Metalog_IsFeasible(a) then do;
     msg = "(Metalog_CDFX) Invalid coefficients. The model is not a valid distribution.";
     run PrintToLog(msg,2);
     STOP;
   end;
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

store module = (
Metalog_ValidateData
Metalog_Logit
Metalog_Order
Metalog_BoundType
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
);

