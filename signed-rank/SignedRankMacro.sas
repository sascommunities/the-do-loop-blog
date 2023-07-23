/*
Definition of the %SIGNEDRANK macro.
Cite as:
Leuchs AK, Neuhäuser M. (2010) "A SAS/IML algorithm for exact nonparametric paired tests."
GMS Med Inform Biom Epidemiol. 6(1)
Freely available from: http://www.egms.de/en/journals/mibe/2010-6/mibe000104.shtml

SYNTAX:
%signedrank(data, label_diff, test, alternative, round=4);

*data: data set 
*var_diff: name of variable (difference) to evaluate 
*test (optional): Specification of the test 
    - 'Signed' for Wilcoxon's signed rank test 
    - 'Pratt' for modified test according to Pratt 
    - 'Original' for test based on original data 
    DEFAULT: TEST='Signed' 
*alternative (optional): Specification of the alternative 
    - 'less' one-sided test with H1: mu<0 
    - 'greater' one-sided test with H1: mu>0 
    - 'two' two-sided test with H1: mu ^= 0
    - 'all' all three tests 
    DEFAULT: ALTERNATIVE='two' 
*round (optional): specifies the number of decimal places to round on 
 (rounding only for the test based on original data) 
 DEFAULT: ROUND=4 

This version has been slightly modified and reformatted by 
Rick Wicklin 19JUL2023.
See 
https://blogs.sas.com/content/iml/2023/07/24/exact-signed-rank-pratt.html
1. Make test and alternative parameters optional
2. Use %UPCASE to make parameters insensitive to case
3. Explain why the test statistic in this macro is different from PROC UNIVARIATE
*/
 
%MACRO signedrank(data, var_diff, test, alternative, round=4); 
* Wicklin Modification: define default values for skipped parameters;
%if not %length(&test)        %then %let test='Signed';
%if not %length(&alternative) %then %let alternative='Two';

data data;
   set &data(rename=(&var_diff=obs_diff));
   diff_abs=abs(obs_diff);
   keep diff_abs obs_diff;
run;

proc sql noprint;
   select count(*) into :diff0 from data where diff_abs=0;
   select count(*) into :n_all from data;
quit;

*abort if all differences are zero;
%IF &n_all=&diff0 %THEN %DO;
%put NOTE: All differences are zero;
%RETURN;
%END;

*assigning ranks for nonzero values;
proc rank data=data out=data;
   where obs_diff ^= 0;
   var diff_abs;
   ranks rank_diff;
run;

*signed ranks;
data data;
   set data;
   rank_sign=rank_diff*sign(obs_diff);
run;

proc iml;
use data;
*Wilcoxons signed rank test;
if upcase(&test)='SIGNED' then do;
   read all var {rank_sign} into d;     *signed ranks;
end;
*Pratts modification of Wilcoxons signed rank test;
if upcase(&test)='PRATT' then do;
   read all var {rank_sign} into d;
   *signed ranks;
   d=(abs(d)+&diff0)#sign(d);
   *assigning ranks according to Pratts modification;
end;
*test based on original data;
if upcase(&test)='ORIGINAL' then do;
   read all var {obs_diff} into d;
   *rounding;
   d=round((10**&round)*d);
end;
close data;

/* Note by Rick Wicklin 20JUL2023: 
   For &test='signed', the following test statistic 
   is different from the statistic by PROC UNIVARIATE.
   Let S be the statistic used by the UNIVARIATE procedure.
   Then S = Tplus - nt(nt+1)/4, where Tplus is the sum of the ranks where x>0
   and nt is the number of values that are not 0.
   The %signedrank macro reports Tplus, not S.
   As discussed in Wicklin (2023) "On the computation of the Wilcoxon 
   signed rank statistic," S and Tplus are equivalent test statistics. See
   https://blogs.sas.com/content/iml/2023/07/19/wilcoxon-signed-rank.html
*/
*computation of the statistic;
tstat=sum(d#(d>=0));

/* shift-algorithm */
start shift(d);
   n=NROW(d);
   *for the shift-algorithm the values in d need to be integer;
   potenz=0;

   do while (sum(d^=int(d))^=0);
      d=10*d;
      potenz=potenz+ 1;
   end;
   ad=abs(d);
   *absolute values;
   *determine largest common factor;
   ggt=0;
   k=min(ad + (d=0)#max(ad));   *smallest |d| ^= 0;
   do while (ggt=0 & k>=1);
      if d/k=int(d/k) then
         ggt=k;
      else
         k=k-1;
   end;
   *required constants;
   adshort=ad/ggt;
   lng=sum(adshort);
   values=((0:lng)*ggt)`;
   *possible values for statistic;
   *algorithm;
   if lng ^= 0 then do;
      prob=1 // j(lng, 1, 0);
      do k=1 to n;
         if adshort[k] ^= 0 then do;
            shift=j(adshort[k], 1, 0) // prob[1:lng+1 -adshort[k]];
            prob=prob + shift;
         end;
         else
            prob=2*prob;
         *difference can be 0;
      end;
      prob=prob/(2**n);
      *reverse the process which changes the observations/ranks 
                      to integers (this was needed for the shift algorithm);
      values=values/(10**potenz);
      *resulting distribution (1th column: statistic -- 2nd column: probability);
      dist=values || prob;
   end;
   d=d* 10**(-1*potenz);
   return (dist);
finish;

dist=shift(d);

*computation of p-values;
*one-sided;
pvalue_gr=sum(dist[, 2]#(dist[, 1]>=tstat));
pvalue_less=sum(dist[, 2]#(dist[, 1]<=tstat));
*two-sided;
tstatOUT=tstat;
*output-statistic;

if upcase(&test)='SIGNED' then do;
   ew=(&n_all-&diff0)*((&n_all-&diff0)+1)/4;
end;
else if upcase(&test)='PRATT' then do;
   ew=((&n_all*(&n_all+ 1)-(&diff0*(&diff0+1)))/4);
end;
else if upcase(&test)='ORIGINAL' then do;
   tstatOUT=tstat/(10**&round);
   ew=sum(abs(d))/2;
end;

if (tstat<ew) then do;
   d=-1*d;
   tstat=sum(d#(d>=0));
end;
lower=ew-(tstat-ew);
pvalue_two=sum(dist[, 2]#(dist[, 1]<=lower | dist[, 1]>=tstat));

*OUTPUT;
if upcase(&alternative)='ALL' then do;
   label={'n' 'n_nonzero' 'statistic' 'pvalue_gr' 'pvalue_less' 'pvalue_two'};
   out=&n_all || &n_all-&diff0 || tstatOUT || pvalue_gr || pvalue_less || pvalue_two;
end;

if UPCASE(&alternative)='GREATER' then do;
   label={'n' 'n_nonzero' 'statistic' 'pvalue_gr'};
   out=&n_all || &n_all-&diff0 || tstatOUT || pvalue_gr;
end;

if upcase(&alternative)='LESS' then do;
   label={'n' 'n_nonzero' 'statistic' 'pvalue_less'};
   out=&n_all|| &n_all-&diff0 || tstatOUT || pvalue_less;
end;

if upcase(&alternative)='TWO' then do;
   label={'n' 'n_nonzero' 'statistic' 'pvalue_two'};
   out=&n_all|| &n_all-&diff0 || tstatOUT || pvalue_two;
end;

create output from out [colname=label];
append from out;
QUIT;


data output;
   set output;
   test=&test;
run;

%IF %UPCASE(&alternative)='ALL' %THEN %DO;
proc report data=output nowd headline;
   column test n n_nonzero statistic pvalue_less pvalue_gr pvalue_two;
   define n / display 'n';
   define n_nonzero / display 'n (nonzero)';
   define statistic / display 'statistic';
   define pvalue_less / display 'p-value (less)';
   define pvalue_gr / display 'p-value (greater)';
   define pvalue_two / display 'p-value (two-sided)';
run;
%END;
%IF %UPCASE(&alternative)='LESS' %THEN %DO;
proc report data=output nowd headline;
   column test n n_nonzero statistic pvalue_less;
   define n / display 'n';
   define n_nonzero / display 'n (nonzero)';
   define statistic / display 'statistic';
   define pvalue_less / display 'p-value (less)';
run;
%END;
%IF %UPCASE(&alternative)='GREATER' %THEN %DO;
proc report data=output nowd headline;
   column test n n_nonzero statistic pvalue_gr;
   define n / display 'n';
   define n_nonzero / display 'n (nonzero)';
   define statistic / display 'statistic';
   define pvalue_gr / display 'p-value (greater)';
run;
%END;
%IF %UPCASE(&alternative)='TWO' %THEN %DO;
proc report data=output nowd headline;
   column test n n_nonzero statistic pvalue_two;
   define n / display 'n';
   define n_nonzero / display 'n (nonzero)';
   define statistic / display 'statistic';
   define pvalue_two / display 'p-value (two-sided)';
run;
%END;
%MEND;
