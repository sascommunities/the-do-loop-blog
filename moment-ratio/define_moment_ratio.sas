/***********************************************************************
 Programs from 
 Appendix E: Constructing a Moment-Ratio Diagram in SAS
 Wicklin, Rick, 2013, Simulating Data with SAS, SAS Institute Inc., Cary NC.

 If you run the code in this file, it 
 1. Creates an annotation data set named "Anno"
 2. Defines a macro %PlotMRDiagram(DS, annoDS, Transparency=0)
    that you can use to overlay a scatter plot on a simple moment-ratio diagram.
    DS = data set that has variables Skew and Kurt, which are the skewness and kurtosis
         for a distribution, for a data set, or for a set of random samples from a distribution.
    annoDS = pass in the "Anno" data set from (1)
    Transparency = a value between 0 and 1 (inclusive). If Transparency=0, the markers in the
         scatter plot are opague. If Transparency=1, the markers are invisible.

 USAGE:
 %PlotMRDiagram(MySkewKurt, Anno, Transparency=0.8)

 ***********************************************************************/

/* Create the annotation data set for the moment-ratio diagram */
/* All computations are in terms of EXCESS kurtosis.
   Use Y2AXIS to add FULL kurtosis axis. */
%let xL = -2.4; %let xR =  2.4;    /* range of skewness is [xL, xR] */
%let yB = -2.0; %let yT = 10.0;    /* range of kurtosis is [yB, yT] */
%let yB2 = %sysevalf(&yB+3);       /* range for full kurtosis       */
%let yT2 = %sysevalf(&yT+3);

/* Macro to convert between full kurtosis and excess kurtosis */
%macro Ex2Full(k); ((&k)+3)  %mend;
%macro Full2Ex(k); ((&k)-3)  %mend;


/* annotation data set for BOUNDARY of feasible (skew, kurt) region */
data Boundary(drop=x);
length function $12 Curve $12 LineColor $20;
retain DrawSpace "DataValue"
       LineColor "Black"
       Curve "Boundary";
function = "POLYLINE"; 
x1=&xL; y1=1+x1**2; y1=%Full2Ex(y1); output;
function = "POLYCONT";
do x = &xL to &xR by 0.1;
   x1=x; y1=1+x1**2;  y1=%Full2Ex(y1);
   output;
end;
run;


/* annotation data set for BETA region of moment-ratio diagram */
data Beta(drop=x);
length function $12 Curve $12 Label $24 LineColor $20;
retain DrawSpace "DataValue"
       Display "All"
       LineColor "LightGray" 
       Curve "Beta";
function = "TEXT"; Anchor="Left ";
label = "Beta";
x1=&xL; y1=2+x1**2; y1=%Full2Ex(y1); output;
Transparency = 0.5; FillTransparency = 0.5;
FillColor = "LightGray";
function = "POLYGON"; 
x1=&xL; y1=1+x1**2; y1=%Full2Ex(y1); output;
function = "POLYCONT";
n=1;
do x = &xL to &xR+0.05 by 0.1;
   x1=x; y1=1+x1**2; y1=%Full2Ex(y1);
   output;
end;
do x = &xR to &xL-0.05 by -0.1;
   x1=x; y1=3+1.5*x1**2; y1=%Full2Ex(y1);
   output;
end;
run;


/* annotation data set for LOGNORMAL curve of moment-ratio diagram */
data LogNormal;
length function $12 Label $24 Curve $12 LineColor $20;
retain DrawSpace "DataValue"
       LineColor "DarkGreen"
       Curve "LogNormal";
function = "TEXT"; Anchor="Right";
label = "LogN";
drop var var0;
var0 = 0.355; var=var0;
x1= (exp(var)+2)*sqrt(exp(var)-1);
y1 = exp(4*var) + 2*exp(3*var) + 3*exp(2*var) - 6;
output;
function = "POLYLINE"; Anchor=" ";
output;
function = "POLYCONT";
do var = Var0 to 0.005 by -0.005;
   x1= (exp(var)+2)*sqrt(exp(var)-1);
   y1 = exp(4*var) + 2*exp(3*var) + 3*exp(2*var) - 6;
   output;
end;
x1=0; y1=0;output;
run;


/* annotation data set for GAMMA curve of moment-ratio diagram */
data GammaCurve;
retain DrawSpace "DataValue";
length function $12 Label $24 Curve $12 LineColor $20;
retain LineColor "Magenta" Curve "Gamma";
drop a a0 dx N i;
function = "TEXT"; Anchor="Right";
label = "Gam";
a0 = (2/&xR)**2;
a = a0;
x1 = 2/sqrt(a); y1= 6/a; output;
Transparency = 0.5;
function = "POLYLINE"; Anchor=" ";  output;
dx=0.1; N=floor(&xR/dx);
function ="POLYCONT";
do i = 2 to N;
   a = (2/(2/sqrt(a) - 0.1))**2;
   x1 = 2/sqrt(a); y1= 6/a;
   if (&xL<= x1 <=&xR) & (&yB <= y1 <= &yT) then output;
end;
x1=0; y1=0; output;
run;


/* Annotation data set for special points in the moment-ratio diagram:
   Points for the exponential, normal, Gumbel, and t distribution with 
   DOF 5,6,7,8,10,11,12. Also text that marks the invalid region. */
data MRPoints(drop=nu);
retain DrawSpace "DataValue";
length function $12 Label $24 Curve $12;
function = "TEXT"; 
Label="E"; x1 = 2;    y1= 6;   Curve="Exponential"; output;
Label="N"; x1 = 0;    y1= 0;   Curve="Normal";      output;
Label="G"; x1 = 1.14; y1= 2.4; Curve="Gumbel";      output;
Label="SU"; x1 = 0.75; y1= 5; Curve="SU";      output;
Label="SB"; x1 = 2; y1= 4; Curve="SB";      output;
Curve="T";
do nu=5 to 8;
   Label=cats("T",nu); x1=0; y1 = 6/(nu-4); output;
end;
do nu=10 to 12;                              /* plot ellipses marks */
   Label="."; x1=0; y1 = 6/(nu-4); output;
end;
Curve="Region";
do x1=-2,2;
  function="TEXT";     Label="Invalid"; y1=-1; output;
  function="TEXTCONT"; Label="Region";         output;
end;
run;

/* Concatenate the pieces of the moment-ratio diagram */
data Anno;
  set Beta Boundary LogNormal GammaCurve MRPoints;
run;

/* Main macro: Plot the moment-ratio diagram for a given data set and a
   given annotation data set. By default, the (skew, kurt) scatter plot
   is fully opaque, but you can set the transparency in the macro call. */
%macro PlotMRDiagram(DS, annoDS, Transparency=0);
/* given a data set that contains variables KURT and SKEW, this macro
   adds the FULLKURT variable and labels for the three variables */
data &ds;
   set &ds;
   FullKurt = Kurt+3;
   label Kurt="Excess Kurtosis"
      FullKurt="Full Kurtosis"
      Skew="Skewness";
run;

proc sgplot data=&DS sganno=&annoDS noautolegend;
   scatter x=Skew y=Kurt /     transparency=&Transparency;
   scatter x=Skew y=FullKurt / y2axis transparency=1;  /* invisible */
   refline 0 / axis=x transparency=0.2;
   xaxis grid values=(-2.5 to 2.5 by 0.5);
   yaxis  reverse grid values=(-2 to 10);
   y2axis reverse grid values=( 1 to 13);
run;
%mend PlotMRDiagram;
