/* SAS program to accompany the article 
   "Standardized regression coefficients in PROC GLIMMIX"
   by Rick Wicklin, published 17MAY2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/05/17/standardized-coefficients-glimmix.html

   This program shows how to understand the standardized regression coeffcients for
   PROC GLMMIX, which are obtained by using the STDCOEF option on the MODEL stmt.
*/

/* output the standardized coefficients */
proc glimmix data=Sashelp.Class noitprint;
   model Weight = Height Age / ddfm=kr s STDCOEF; /* request standardization of estimates */
   ods select StandardizedCoefficients;
run;


/* How are the standardized coeffs in PROC GLIMMIX created? 
   1. The response variable is not standardized
   2. The standardization for the explanatory variables is
      (X - mean(X))/sqrt(CSS(X))
*/

/* Macro to standardize a set of variables by centering by the sample mean and 
   dividing by the sqrt(CSS) 
   For an explanation, see 
   https://blogs.sas.com/content/iml/2021/05/12/nonstandard-standardize-vars.html

   Example of calling syntax:
   %StdizeMeanSqrtCSS(Height Age,  <== list of variable in input data set
                      2,           <== number of variable in first arg
                      Sashelp.Class,  <== input data set
                      Want);       <== output data: same as input but the specified vars
                                       are standardized by X -> (X-mean(X))/sqrt(CSS(X))
*/
%macro StdizeMeanSqrtCSS(varNames, numVars, inData, outData);
/* output the MEAN and CSS statistics for each var */
proc means data=&inData noprint;
   var &varNames;
   output out=Stats0(drop=_TYPE_) mean=loc1-loc&numVars
                                  css=scale1-scale&numVars;
run;
/* transform CSS -> sqrt(CSS) */
data Stats;
set Stats0;
array scale[*] scale1-scale&numVars;
do i=1 to dim(scale);
   scale[i] = sqrt(scale[i]);
end;
run;

/* the statistics are in one long row. Make one row for each statistic */
data Parms; 
length _TYPE_ $14;
array x[*] &varNames;       
array loc[*] loc1-loc&numVars;
array scale[*] scale1-scale&numVars;
set Stats;
_TYPE_ = "LOCATION";
do i=1 to dim(loc);   x[i] = loc[i];   end;  output;
_TYPE_ = "SCALE";
do i=1 to dim(scale); x[i] = scale[i]; end;  output;
keep _TYPE_ &varNames;
run;

/* use parameters in data set to standardize X1-X3 */
proc stdize data=&inData out=&outData method=in(Parms);
   var &varNames;
run;
%mend;

/* Standardize 2 variables (Height and Age) in the Sashelp.Class data.
   Write the standardized variables to the Want data set. */
%StdizeMeanSqrtCSS(Height Age, 2, Sashelp.Class, Want);

/* run PROC GLIMMIX on the standardized data */
proc glimmix data=Want noitprint;          
   model Weight = Height Age / ddfm=kr s;
   ods select ParameterEstimates;
run;

