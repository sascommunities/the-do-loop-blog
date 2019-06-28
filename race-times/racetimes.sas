/* SAS program to accompany the article 
   "Visualize race times in SAS"
   by Rick Wicklin, published 01JUL2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/07/01/visualize-race-times.html

   This program shows how to create various graphs of race times for
   a 5K race.
   Data from an NC girls high school race in 2019. 
   https://nc.milesplit.com/meets/331151/results/616383/raw#.XRIMt4hKhhE
*/

data RaceTimes;
input time hhmmss.; /* read race time */
format time TIME10.;
Rank = _N_;
datalines;
18:33.24 
19:02.89 
19:22.91 
19:26.34 
19:35.62 
19:41.41 
19:55.91 
19:56.84 
19:57.09 
19:57.86 
20:03.08 
20:27.56 
20:27.84 
20:28.75 
20:31.58 
20:39.31 
20:41.52 
20:42.81 
20:46.88 
20:47.93 
20:48.86 
20:58.78 
20:59.33 
21:00.41 
21:07.35 
21:08.61 
21:11.68 
21:21.08 
21:22.00 
21:36.38 
21:36.74 
21:37.71 
21:42.86 
21:43.60 
21:55.75 
22:00.82 
22:01.33 
22:04.06 
22:13.05 
22:25.08 
22:42.61 
22:44.95 
22:45.55 
22:55.90 
23:00.73 
23:16.65 
23:17.02 
23:18.92 
23:20.28 
23:23.99 
23:29.83 
23:36.73 
23:43.08 
24:23.81 
24:42.94 
24:47.06 
25:42.27 
25:43.67 
25:45.64 
25:52.34 
25:58.02 
26:03.45 
26:05.64 
29:23.10 
31:01.10 
31:26.50 
;

title "Girls High School 5K";
proc sgplot data=RaceTimes noautolegend;
   scatter x=rank y=time;
   xaxis grid values=(1 10 to 70 by 10) valueshint;
   yaxis grid label="Race Time (minutes)" type=time
   values=(%eval(60*18) to %eval(60*32) by 60);
   format time TIME4.;
run;

proc sgplot data=RaceTimes noautolegend;
   histogram time / binwidth=60 binstart=%eval(60*17+30) showbins;
   fringe time / lineattrs=(thickness=2 color=black) transparency=0.6;
   yaxis grid offsetmin=0.05 ;
   xaxis grid label="Race Time (minutes)" /* type=time */
         values=(%eval(60*18) to %eval(60*32) by 60);
   format time TIME4.;
run;

ods graphics on;
proc univariate data=RaceTimes;
   var time;
   cdfplot time / odstitle="Cumulative Distribution of Race Times"
                  odstitle2="Created by Using PROC UNIVARIATE"
                  ;
   ods select cdfplot;
   format time TIME4.;
run;

title "Rank versus Time";
proc sgplot data=RaceTimes noautolegend;
   scatter x=time y=rank;
   yaxis grid;
   xaxis grid label="Race Time (minutes)" type=time
   values=(%eval(60*18) to %eval(60*32) by 60);
   format time TIME4.;
run;



