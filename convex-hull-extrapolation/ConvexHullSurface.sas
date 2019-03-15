/* SAS program to accompany the articles
   "Interpolation vs extrapolation: the convex hull of multivariate data"
   by Rick Wicklin, published 18MAR2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/03/18/interpolation-extrapolation-convex-hull.html

and
   "Truncate response surfaces"
   by Rick Wicklin, published 20MAR2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/03/18/truncate-response-surfaces.html ?

   This program computes and visualizes the convex hull of bivarate data.
   It then fits a quadratic response surface to the data.
   Although the observed responses are all positive, the predicted surface
   has negative predicted values, which occur when you extrapolate the 
   model outside of the convex hull of the data. If you interpolate the 
   model inside the convex hull, the predicted responses are positive,
   as expected.

   You can create graphs that visualize only the positive predicted values 
   by replacing the negative predicted values with missing values.
*/

data Sample;
input X Y Z @@;
datalines;
10 90 22  22 76 13  22 75  7 
24 78 14  24 76 10  25 63  5
26 62 10  26 94 20  26 63 15
27 94 16  27 95 14  29 66  7
30 69  8  30 74  8
;

title "Response Values for Bivariate Data";
proc sgplot data=Sample;
scatter x=x y=y / markerattrs=(size=12 symbol=CircleFilled)
        colorresponse=Z colormodel=AltThreeColorRamp;
xaxis grid; yaxis grid;
run;

proc means data=Sample;
run;

/* Compute convex hull of the data */
options ls=200;
proc iml;
use Sample;
read all var {x y} into points;
close;

/* get indices of points in the convex hull in counter-clockwise order */
indices = cvexhull( points );
/* positive indices are on boundary; negative indices are inside */
print (indices`)[L="indices"];

b = (indices > 0);     /* binary indicator variable for sorted vertices */
call sortndx(ndx, abs(indices)); /* get sort index  */
onBoundary = b[ndx];   /* binary indicator data in original oder */

title "Convex Hull of Bivariate Data";
call scatter(points[,1], points[,2]) group=onBoundary 
             option="markerattrs=(size=12 symbol=CircleFilled)";

/* you can extract the convex hull, in CC order */
hullNdx = indices[loc(b)];
convexHull = points[hullNdx, ];

/* In SAS/IML 14.1, you can use the polygon package to visualize the 
   convex hull:
   https://blogs.sas.com/content/iml/2016/04/27/packages-share-sas-iml-programs-html
*/
package help polygon;
package load polygon;
run PolyDraw(convexHull, points||onBoundary) GRID={x y}
    MARKERATTRS="size=12 symbol=CircleFilled";

/* or you can write the data to a SAS data set and use the 
   POLYGON statement in PROC SGPLOT */
p = points || onBoundary;
create TheData from p[colname={x y "onBoundary"}]; append from p; close;

poly = j(nrow(convexHull), 1, 1) || convexHull;
create Hull from poly[colname={ID cX cY}]; append from poly; close;
quit;

/* combine convex hull and data. Overlay the graphs */
data All; set TheData Hull; run;

proc sgplot data=All noautolegend;
   polygon x=cX y=cY ID=id / fill;
   scatter x=x y=y / group=onBoundary markerattrs=(size=12 symbol=CircleFilled); 
   xaxis grid; yaxis grid;
run;

/**************************************************/
/* END post on convex hull.                       */
/* BEGIN post on truncating the response surface. */
/**************************************************/


data Sample;
input X Y Z @@;
datalines;
10 90 22  22 76 13  22 75  7 
24 78 14  24 76 10  25 63  5
26 62 10  26 94 20  26 63 15
27 94 16  27 95 14  29 66  7
30 69  8  30 74  8
;

/* show the predicted response surfaces. Notice negative values */
ods graphics / width=400px height=400px ANTIALIASMAX=10000;
proc rsreg data=Sample plots=surface(fill=pred overlaypairs);
   model Z = Y X;
run;

proc rsreg data=Sample plots=surface(3d fill=Pred gridsize=80);
   model Z = Y X;
   ods select Surface;
   ods output Surface=Surface; /* use ODS OUTPUT to save surface data to a data set */
run;

/* rename vars and set negative responses to missing */
data Surf2;
set Surface(rename=(
       Predicted0_1_0_0 = Pred  /* rename the long ODS names */
       Factor1_0_1_0_0  = GY    /* 'G' for 'gridded' */
       Factor2_0_1_0_0  = GX))  
    Sample(in=theData);         /* combine with original data */
if theData then Type = "Data   ";
else            Type = "Gridded";
if Pred < 0 then Pred = .;      /* replace negative predictions with missing values */
label GX = 'X'  GY = 'Y';
run;

proc means data=Surf2 N NMiss Min Max;
var Pred GX GY;
run;

/* View the underlying GTL. Create new GTL templates to plot modified data */
proc template;
source Stat.Rsreg.Graphics.Contour;
source Stat.Rsreg.Graphics.Surface;
run;

proc template;                        /* surface plot with continuous color ramp */
define statgraph SurfaceTmplt;
dynamic _X _Y _Z _Title;              /* dynamic (gridded) variables */
 begingraph;
 entrytitle _Title;                   /* specify title at run time (optional) */
  layout overlay3d;
    surfaceplotparm x=_X y=_Y z=_Z /  /* specify variables at run time */
       name="surface" 
       surfacetype=fillgrid
       colormodel=threecolorramp      /* or =twocolorramp */
       colorresponse=_Z;              /* prior to 9.4m2, use SURFACECOLORGRADIENT= */
    continuouslegend "surface";
  endlayout;
endgraph;
end;
run;

proc template;                        /* surface plot with continuous color ramp */
define statgraph ContourTmplt;
dynamic _GX _GY _Z                    /* gridded variables */
        _X _Y _Title;                 /* dynamic variables */
 begingraph;
 entrytitle _Title;                   /* specify title at run time (optional) */
  layout overlay / aspectratio=1;
    contourplotparm x=_GX y=_GY z=_Z /  /* specify variables at run time */
       name="contour" 
       contourtype=labeledlinegradient
       gridded=TRUE
       NHINT=10
       colormodel=threecolorramp;     /* or =twocolorramp */
    scatterplot x=_X y=_Y / markerattrs=
            GRAPHDATADEFAULT(color=GraphOutlier:ContrastColor);
    continuouslegend "contour";
  endlayout;
endgraph;
end;
run;

proc sgrender data=Surf2 template = ContourTmplt;
 dynamic _GX='GX' _GY='GY' _Z='Pred' _X='X' _Y='Y'
         _Title="Truncated Response Contour for Z";
run;

proc sgrender data=Surf2 template = SurfaceTmplt;
where type = "Gridded";
 dynamic _X='GX' _Y='GY' _Z='Pred' 
         _Title="Truncated Response Surface for Z";
run;
