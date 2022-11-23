*proc iml;
/******************************/
/* GRAPHICAL ROUTINES         */
/******************************/

/* plot the CDF for the ML object */
start ML_PlotCDF(L, _p=);
   if IsSkipped(_p) then 
      p = {0.0001 0.0005 0.001} || do(0.01, 0.99, 0.01) || {0.999 0.9995 0.9999};
   else do;
      p = colvec(_p);
      call sort(p);
   end;
   x = ML_Quantile(L, p);
   call series(x, p) grid={x y} label={'x' 'Cumulative Probability'};
finish;

/* plot the PDF for the ML object */
start ML_PlotPDF(L, _p=);
   if IsSkipped(_p) then 
      p = {0.0001 0.0005 0.001} || do(0.01, 0.99, 0.01) || {0.999 0.9995 0.9999};
   else do;
      p = colvec(_p);
      call sort(p);
   end;
   x = ML_Quantile(L, p);
   f = ML_PDF(L, p);
   call series(x, f) grid={x y} label={'x' 'Density'};
finish;

/* plot the ECDF for data and overlay the CDF for the ML model */
start ML_PlotECDF(L, _p=);
   if IsSkipped(_p) then 
      p = {0.0001 0.0005 0.001} || do(0.01, 0.99, 0.01) || {0.999 0.9995 0.9999};
   else do;
      p = colvec(_p);
      call sort(p);
   end;
   q = ML_Quantile(L, p);
   x = L$'x';
   ECDF = L$'cdf';
   create Temp var {"x" "ECDF" "q" "p"}; append; close;
   order=L$'order';
   submit order;
   proc sgplot data=Temp;
      scatter x=x y=ECDF / markerattrs=(color=gray symbol=Circle) legendlabel="Data";
      series x=q y=p / legendlabel="Metalog_&order";
      yaxis grid label="Cumulative Probability";
      xaxis grid label="x";
      keylegend / location=inside across=1 opaque;
   run;
   endsubmit;
finish;


%let ML_modules_plot =
ML_PlotCDF
ML_PlotPDF
ML_PlotECDF
;

store module=(&ML_modules_plot);
 
