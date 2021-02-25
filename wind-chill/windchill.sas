/* SAS program to accompany the article 
   "The wind chill chart"
   by Rick Wicklin, published 24Feb2021 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2021/02/24/wind-chill-chart.html

   This program shows how to construct a wind chill chart from a 
   statistical formula that predicts the wind chill from current
   values of temperature and wind speed. Both US units and SI units
   are shown. The formulas are available from 
   https://www.weather.gov/safety/cold-wind-chill-chart    
   and
   https://en.wikipedia.org/wiki/Wind_chill
*/


/* Wind chill formula for US customary units (British Imperial units).
   The formula gives the apparent temperature (F) from 
   ambient temperature (F) and wind speed (mph).
   Source NOAA: https://www.weather.gov/safety/cold-wind-chill-chart    
*/
data WindChillUS;
label T = "Ambient Temperature (F)"
      V = "Wind Speed (mph)"
      WCT = "Wind Chill Temperature (F)";
do T = 40 to -40 by -5;
   do V = 5 to 40 by 5;
      WCT = 35.74 + 0.6215*T - 35.75*V**0.16 + 0.4275*T*V**0.16;
      output;
   end;
end;
run;


ods graphics / width=640px height=480px;
title "Wind Chill Chart (US Units)";
title2 "Dependence on Temperature";
footnote J=L "Source NOAA: https://www.weather.gov/safety/cold-wind-chill-chart";
proc sgplot data=WindChillUS;
   where V in (5 10 20 30);
   series x=T y=WCT / group=V;
   keylegend / position=E;
   xaxis grid values=(-40 to 40 by 5) valueshint; 
   yaxis grid values=(-80 to 40 by 5) valueshint;
run;

title2 "Dependence on Wind Speed";
proc sgplot data=WindChillUS;
   label T = "Temperature (F)";
   where T in (30 20 10 0);
   series x=V y=WCT / group=T lineattrs=(thickness=2);
   keylegend / position=E ;
   xaxis grid values=(-40 to 40 by 5) valueshint; 
   yaxis grid values=(-80 to 40 by 5) valueshint;
run;

title "Wind Chill Chart (US Units)";
footnote J=L "Source NOAA: https://www.weather.gov/safety/cold-wind-chill-chart";

/* use threshold to decide between white text and dark text */
data WindChillUS2;
set WindChillUS;
textColor = (WCT > -19.5);  /* use to color text */
run;

ods graphics / width=640px height=360px;
title2;
proc sgplot data=WindChillUS2;
   format WCT 3.0;
   heatmapparm x=T y=V colorresponse=WCT / name="w"
      colormodel=(DarkPurple Purple DarkBlue DarkCyan LightGray LightYellow LightOrange);
   /* trick: white text on a dark background and dark text on a light background */
   text x=T y=V text=WCT / colorresponse=textColor colormodel=(white black) textattrs=(weight=bold);
   gradlegend "w";
   xaxis values=(-40 to 40 by 5) valueshint; 
   yaxis values=(-100 to 40 by 5) valueshint;
run;

/*********************************/

/* Wind chill formula for SI units.
   The formula gives the apparent temperature (C) from 
   ambient temperature (C) and wind speed (kph).
   Source Environment Canada: https://en.wikipedia.org/wiki/Wind_chill
   and the formulas C = 5/9*F-32 and km = miles/1.609
*/
data WindChillSI;
label TC = "Ambient Temperature (C)"
      VK = "Wind Speed (kph)"
      WCT = "Wind Chill Temperature (C)";
do TC = 6 to -39 by -3;
   do VK = 5 to 65 by 5;
      WCT = 13.12 + 0.6215*TC - 11.37*VK**0.16 + 0.3965*TC*VK**0.16;
      output;
   end;
end;

/* use threshold to decide between white text and dark text */
data WindChillSI2;
set WindChillSI;
textColor = (WCT > -29.5);  /* use to color text */
run;

ods graphics / width=640px height=480px;
title "Wind Chill Chart (SI Units)";
proc sgplot data=WindChillSI2;
   format WCT 3.0;
   heatmapparm x=TC y=VK colorresponse=WCT / name="w"
      colormodel=(DarkPurple Purple DarkBlue DarkCyan LightGray LightYellow LightOrange);
   /* trick: white text on a dark background and dark text on a light background */
   text x=TC y=VK text=WCT / colorresponse=textColor colormodel=(white black) textattrs=(weight=bold);
   gradlegend "w";
   xaxis values=(-39 to 6 by 3) valueshint; 
   yaxis values=(5 to 100 by 5) valueshint;
run;

title;footnote;
ods graphics/reset;
