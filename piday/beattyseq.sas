/* SAS program to accompany the article 
   "The Beatty sequence for pi"
   by Rick Wicklin, published 9MAR2022 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2022/03/09/beatty-sequence-pi.html

   This program shows how to compute, visualize, and validate
   the Beatty sequence for pi.

   For an overview, see
   https://www.johndcook.com/blog/2015/03/14/beatty-sequence/
   Beatty's Thm: 
   Every positive integer is either part of the sequence 
   [ n pi ] or the sequence 
   [ n pi/(pi – 1) ] where n ranges over positive integers, 
   and no positive integer is in both sequences.

   In general, if r > 1 is irrational, then define s = r/(r-1). 
   The sequences 
    B[n; r] = floor(n*r)  and
    C[n; r] = floor(n*r/(r-1)) 
   partition the positive integers into disjoint sets.
*/
%let maxN = 355;            /* because pi ~ 355/113 */
data Beatty;
length Slope $9;
pi = constant('pi');
Sequence = 'B'; Slope = 'pi';
do n = 1 to &maxN;
   s = floor(pi*n);         /* the Beatty sequence */
   output;
end;
Sequence = 'C'; Slope = 'pi/(pi-1)';
do n = 1 to &maxN;
   s = floor(pi/(pi-1)*n);  /* the complementary sequence */
   output;
end;
/* find the lessor of the maximum values of the sequences */
MB = floor(pi* &maxN);
MC = floor(pi/(pi-1)* &maxN);         
call symputx('maxS', min(MB, MC));    /* for maxN=355, maxS=520 */
drop MB MC pi;
run;

%put &=maxS;   /* display value in SAS log */
proc sort data=Beatty;
   by s;
run;

ods graphics / reset;
ods graphics / width=400px height=560px;
title "Beatty and Complementary Sequences for pi";
proc sgplot data=Beatty;
   scatter x=n y=s / group=Sequence markerattrs=(symbol=CircleFilled);
   xaxis grid;
   yaxis max=&maxS grid LABEL="Sequence Value";
run;

proc print data=Beatty(obs=12) noobs;
   var s Sequence;
run;


data SeqViz;
set Beatty;
z = 1;
run;
ods graphics / width=640px height=200px;
title "Beatty and Complementary Sequences for pi";
title2 "First 60 Integers";
proc sgplot data=SeqViz noautolegend;
   where s<=60;
   yaxis display=none;
   xaxis grid gridattrs=(color=CXF0F0F0) values=( 0 to 21 by 3
                      25 to 43 by 3
                      47 to 60 by 3 ) valueshint;
   scatter x=s y=z / group=Sequence datalabel=Sequence
           datalabelpos=top datalabelattrs=(size=12)
           markerattrs=(symbol=SquareFilled size=12);
run;

title2 "Another 60 Integers";
proc sgplot data=SeqViz noautolegend;
   where 302 <= s <= 361;
   yaxis display=none;
   xaxis grid gridattrs=(color=CXF0F0F0) values=(304,307,
                      311 to 329 by 3,
                      333 to 354 by 3, 358,361) valueshint;
   scatter x=s y=z / group=Sequence datalabel=Sequence
           datalabelpos=top datalabelattrs=(size=12)
           markerattrs=(symbol=SquareFilled size=12);
run;
/*  The continued fraction expansion of pi is:
3/1
22/7
333/106
355/113
103993/33102
104348/33215
208341/66317
312689/99532
*/


data Converge;
set Beatty(where=(s<= &maxS));
NumInBeatty + (Sequence='B');
Ratio = s / NumInBeatty;
run;

ods graphics / width=640px height=480px;
title "Ratio of Counts: s / #(BettyNumbers <= s)";
proc sgplot data=converge noautolegend;
   scatter x=s y=Ratio / group=Sequence;
   refline 3.14159 / axis=y label='pi';
   refline 3 22 333 355 / axis=x;
   yaxis min=3.1 max=3.18;
   xaxis max = 400;
run;

/* Do the two sequences intersect? Is the union all integers in [1, &maxS]? */
proc iml;
use Beatty;
   read all var 's' into B where(Sequence='B');
   read all var 's' into C where(Sequence='C');
close;

/* is there any intersection? */
intersect = xsect(B, C);
if IsEmpty(intersect) then 
   print "B and C do not intersect";
else print "The intersection is:", intersect;

max = min( max(B), max(C) );

/* use the UNION function */
union = union(B, C);
union = union[ , 1:max];  /* only test for values in [1, &maxS] */
if all(union = (1:max)) then 
   print "All integers are found";
else 
   print "Some integer is missing";
print max;
QUIT;


