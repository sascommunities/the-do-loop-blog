
proc iml;
load module=(expm1 CDFBirthday PDFBirthday QuantileBirthday RandBirthday);

/*********************/
/* CDF and PDF       */
/*********************/

/* some tests for evaluating the birthday CDF */
Prob2BDay  = CDFbirthday(23);           * N=23, nCat=365, nDup=2;
Prob3BDay = CDFBirthday(23,  365, 3);   * N=23, nCat=365, nDup=3;
EquiProb3 = CDFBirthday(88,  365, 3);   * N=88, nCat=365, nDup=3;
EquiProb4 = CDFBirthday(187, 365, 4);   * N=187, nCat=365, nDup=4;
print Prob2BDay Prob3BDay EquiProb3 EquiProb4;

/* if N=23, what is the probability that three or more share a birth month?
   What about a birth day-of-month (eg, the 12th or some month) */
Prob3Month = CDFBirthday(23, 12, 3);        * N=23, nCat=12, nDup=3;
Prob3DayOfMonth = CDFBirthday(23, 30, 3);   * N=23, nCat=30, nDup=3;
print Prob3Month Prob3DayOfMonth;

/* graph the distribution and density */
N = T(1:187);
CDF3 = j(nrow(N), 1, .);
do i = 1 to nrow(N);
   CDF3[i] = CDFBirthday(N[i], 365, 3); * nCat=365 nDup=3;
end;
title "Probability of Three or More Matching Birthdays";
title2 "Assume 365 Equally Likely Birthdays";
call series(N, CDF3) grid={x y} label={"Number of People (N)" "Prob(Match)"};

/* compute and graph probability for multiple values of nDup */
N = T(1:187);
CDF = j(nrow(N), 3, .);
do nDup = 2 to 4;
   if nDup=2 then max=100; else max=nrow(N);
   do i = 1 to max;
      CDF[i,nDup-1] = CDFBirthday(N[i], 365, nDup); * nCat=365;
   end;
end;

create CDF from N CDF[c={'N' 'CDF2' 'CDF3' 'CDF4'}];
append from N CDF;
close;

N = T(1:60);
PDF = j(nrow(N), 1);
do i = 1 to nrow(N);
   PDF[i] = PDFBirthday(N[i], 365, 2);    * C=365, D=2;
end;
 
title "Probability Density for Two Shared Birthdays";
title2 "Assume 365 Equally Likely Days";
call series(N, PDF) grid={x y} label={"Number of People (N)" "Density"};

/******************************************************/
/* run the tests from the R documentation */
/******************************************************/

*classical probability of 2 or more matching birthdays;
p = CDFbirthday(23);
print p (0.5072972)[L="R Answer"];

*probability of 3 or more matching birthdays;
p = CDFbirthday(23, 365, 3);
print p (0.01441541)[L="R Answer"];

*prob of 4 or more shared birthdays in 150 people;
p = CDFbirthday(150, 365, 4);
print p (0.2690146)[L="R Answer"];

*100 or more coincident birthdays in 1000 people: very rare;
p = CDFbirthday(1000, 365, 100);
print p[F=E13.] (1.531434e-113)[L="R Answer" F=E13.];

QUIT;

title "Probability of Multiple Matching Birthdays";
title2 "Assume 365 Equally Likely Birthdays";
proc sgplot data=CDF;
   series x=N y=CDF2 / curvelabel="2 Matches";
   series x=N y=CDF3 / curvelabel="3 Matches";
   series x=N y=CDF4 / curvelabel="4 Matches";
   xaxis grid label="Number of People (N)" values=(0 to 180 by 30) valueshint;
   yaxis grid label="Prob(Match)" values=(0 to 1 by 0.1);
run;



proc iml;
load module=(expm1 CDFBirthday PDFBirthday QuantileBirthday RandBirthday);

/*********************/
/* QUANTILE and RAND */
/*********************/
 
*classical birthday problem;
N2 = QuantileBirthday(0.5); /* find N s.t. Pr(2 BDays) >= 0.5 */
N3 = QuantileBirthday(0.5, 365, 3); /* find N s.t. Pr(3 BDays) >= 0.5 */
N4 = QuantileBirthday(0.5, 365, 4); /* find N s.t. Pr(4 BDays) >= 0.5 */
print N2 N3 N4;

*Example from Diaconis & Mosteller p. 858;
*Three BDays are on the 16th, but different months;
*What size group gives equal prob of this occurence?;
N3_30 = QuantileBirthday(0.5, 30, 3); *assumes 30 days in month;
print N3_30;

prob = do(0.05, 0.95, 0.0025);
Q = j(1, ncol(prob));
do i = 1 to ncol(prob);
   Q[i] = QuantileBirthday(prob[i], 365, 2); *default nCat=365, nDup=2;
end;
title "Quantile Function for Classic Birthday Problem";
title2 "Assume 365 Equally Likely Birthdays";
call series(prob, Q) grid={x y} 
                     xvalues=do(0.1,0.9,0.1)
                     yvalues=do(5,50,5)
                     label={"Prob(Matching)" "Number of People (N)"};


/******************************************/
/* run the tests from the R documentation */
/******************************************/                  

*classical birthday problem;
N = QuantileBirthday(0.5);
print N (23)[L="R Answer"];

*0.9 probability of three or more shared birthdays;
N = QuantileBirthday(0.9, 365, 3);
print N (135)[L="R Answer"];

*Example from Diaconis & Mosteller p. 858;
*Three BDays are on the 16th, but different months;
*What size group gives equal prob of this occurence?;
N = QuantileBirthday(0.5, 30, 3); *assumes 30 days in month;
print N (18)[L="R Answer"];

*other examples from R doc;
N = QuantileBirthday(0.5, 365, 4);  *exact value=187;
print N (187)[L="R Answer"];

N = QuantileBirthday(0.5, 365, 10); *exact value=1181;
print N (1179)[L="R Answer"];

*group size for 0.5 prob of sharing a 4-digit PIN number;
N = QuantileBirthday(0.5, 1E4);
print N (119)[L="R Answer"];

* test the RandBirthday function;
call randseed(54321, 1);
Sizes = RandBirthday(5);
print Sizes;


NN = RandBirthday(500);
title "Random Sizes of Rooms That Have Two Shared Birthdays";
call histogram(NN);

QUIT;
