*proc iml;

/* SAS/IML implementation of the metalog ("metalogistic") distribution 
   http://metalogdistributions.com/images/The_Metalog_Distributions_-_Keelin_2016.pdf
   and
   http://metalogdistributions.com/home.html

   This implementation uses an object-oriented approach. The metalog
   distribution is created as an object (actually, a named list)
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
        take a series of matrix/vector/scalar arguments
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
   if type(bounds) ^= 'N' | ncol(bounds) ^= 2 then 
      %StopOnError("The bounds for the metalog model must be a two-element numerical vector.");
   if bounds[1]^=. & bounds[2]^=. then 
      if bounds[1] >= bounds[2] then 
         %StopOnError("The lower bound must be less than the upper bound.");
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
   Lbound = (bounds[1]^=.);
   Ubound = (bounds[2]^=.);
   if  Lbound &  Ubound then return( 'B' );
   if ^Lbound &  Ubound then return( 'SU' );
   if  Lbound & ^Ubound then return( 'SL' );
   return( 'U' );
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



%let ML_modules =
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
ML_Rand;

store module=(&ML_modules);


