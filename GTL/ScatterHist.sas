/* SAS program to accompany the article 
   "How to create a scatter plot with marginal histograms in SAS"
   by Rick Wicklin, published 20May2011 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2011/05/20/how-to-create-a-scatter-plot-with-marginal-histograms-in-sas.html

   This program shows how to use the Graph Template Language (GTL)
   to create a scatter plot that has marginal histograms on the 
   X2 and Y2 axis.
*/

proc template;
  define statgraph scatterhist;
  dynamic XVAR YVAR TITLE TRANSVAL;
  begingraph / designwidth=600px designheight=400px;
    entrytitle TITLE;
    layout lattice / rows=2 columns=2 rowweights=(.2 .8) columnweights=(.8 .2)
                     rowdatarange=union columndatarange=union rowgutter=0 columngutter=0;
    /* histogram at X2 axis position */
    layout overlay / walldisplay=(fill)
                     xaxisopts=(display=none) yaxisopts=(display=none offsetmin=0);
        histogram XVAR / binaxis=false;
    endlayout;
    /* NOBS count cell */  
    layout overlay; 
      entry " "; *'NOBS = ' eval(n(XVAR));    
    endlayout;
    /* scatter plot with conditional axis types */
    layout overlay;         
      scatterplot y=YVAR x=XVAR / datatransparency=TRANSVAL markerattrs=(symbol=circlefilled);
    endlayout;
    /* histogram at Y2 axis position */
    layout overlay / walldisplay=(fill)
                     xaxisopts=(display=none offsetmin=0) yaxisopts=(display=none);
      histogram YVAR / orient=horizontal binaxis=false;
    endlayout;
    endlayout;    
  endgraph;
  end;
run;

/* default: use solid markers */
proc sgrender data=SasHelp.Class template=scatterhist;  
  dynamic YVAR="Weight" XVAR="Height"
          TITLE="Height-Weight Relationship";
run;

/* for larger data sets, use the TRANSVAL parameter to set transparency */
data SimNormal;
call streaminit(1);
do i = 1 to 1000;
   x = rand("Normal"); 
   y = x + rand("Normal"); 
   output;
end;
run;

proc sgrender data=SimNormal template=scatterhist;  
  dynamic YVAR="X" XVAR="Y" TRANSVAL=0.6
          TITLE="Bivariate Normal Data";
run;

