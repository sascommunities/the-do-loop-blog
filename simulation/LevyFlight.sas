/* SAS program to accompany the article 
   "Levy flight and vectorizing a simulation in SAS"
   by Rick Wicklin, published 04NOV2024 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2024/11/04/levy-flight-vectorize.html

   This program shows how to write a simulation of two models of animal locomotion:
   1. A Gaussian random walk, in which the animal stays near its initial location.
   2. A Levy flight, in which the animal makes a few moves near it's initial 
      loation, then moves a long distance in a random location. This process
      repeats.
*/

proc iml;
/* generate n random variates from folded normal (Gaussian) */
start RandAbsNormal(n, sigma); 
   return( abs( randfun(N, "Normal", 0, sigma) ) ); 
finish;

/* helper function: return n random variates from IGamma(alpha,beta) 
   See https://blogs.sas.com/content/iml/2021/01/27/inverse-gamma-distribution-sas.html */
start RandIGamma(n, alpha, beta);
   return (1 / randfun(n, 'Gamma', alpha, 1/(beta)));
finish;
/* generate n random variates from truncated Levy(0,c) with variates truncated at MAX.
   See https://blogs.sas.com/content/iml/2014/12/01/max-and-min-rows-and-cols.html
*/
start RandTruncLevy(n, c, max); 
   x = RandIGamma(n, 0.5, c/2);  /* x ~ Levy(0, c) */
   return (x >< max);
finish;
store module = _all_;

call randseed(12345);
n = 1000;                 /* number of steps */
sigma = 0.1;              /* scale parameter for N(0,sigma) */
c = 0.3*sigma;            /* Levy parameter that makes medians similar */
max = 3;                  /* no flights can be longer than MAX units */

GaussSteps = RandAbsNormal(n, sigma);   /* n values from folded normal */
LevySteps  = RandTruncLevy(n, c, max);  /* n values from truncated Levy */

Distribution = repeat( {'Gaussian', 'Levy'}, 1, nrow(GaussSteps) );
Steps = GaussSteps // LevySteps;
create Steps var {'Distribution' 'Steps'};
append;
close;
QUIT;

title "Gaussian and Truncated Levy Distributions";
proc sgpanel data=Steps;
panelby Distribution / columns=1;
histogram Steps / binwidth=0.1 binstart=0.05;
rowaxis grid;
colaxis grid values=(0 to 3 by 0.5) valueshint;
run;

proc iml;
load module = _all_;

/* generate directions uniformly at random. See
   https://blogs.sas.com/content/iml/2021/02/03/random-points-sphere.html
*/
start GetRandomUnitVectors(N, dim);
   x = randfun(N//dim, "Normal");   /* x ~ MVN(0, I(p)) */
   u = x / sqrt(x[ ,##]);           /* unit vectors in random directions */ 
   return ( u );
finish;

/* Let x0 be an initial position in R^d.
   Let the rows of the Nxd matrix, V, be vectors in R^d.
   Each row describes how to move relative to the previous position.
*/
start MakePath(x0, V);
   N = nrow(V);
   d = ncol(V);
   path = j(N,d,.); 
   do j = 1 to d;
      path[,j] = cusum(V[,j]);  /* component-wise sum of steps */
   end;
   return( x0 // path );
finish;

/* simulate taking N steps from an initial position in dim dimensions
   where step size is abs(N(0,sigma))
*/   
start SimGaussianWalk(N, x0, sigma);
   step = RandAbsNormal(N, sigma); /* Gaussian step sizes */
   u =  GetRandomUnitVectors(N, ncol(x0));
   s = step # u;                    /* steps in random directions */
   path = MakePath(x0, s);
   return( path );
finish;

/* simulate taking N steps from an initial position in dim dimensions
   where step size is Levy(0,c) where c is scale parameter
*/   
start SimLevyFlight(N, x0, c, max);
   step = RandTruncLevy(n, c, max); /* <= only this statement is different! */
   u =  GetRandomUnitVectors(N, ncol(x0));
   s = step # u;                    /* steps in random directions */
   path = MakePath(x0, s);
   return( path );
finish;


/* call each simulation with the same random number stream */
call randseed(12345, 1);
N = 100;                   /* number of steps in random walk */
x0 = {0 0};                /* initial position */
sigma = 0.1;               /* scale parameter for N(0,sigma) */
walk = SimGaussianWalk(N, x0, sigma);

/*
title "Vectorized Gaussian Random Walk";
call series(walk[,1], walk[,2]) grid={x y} option="markers";
*/

call randseed(12345, 1);   /* reset random number stream */
c = 0.3*sigma;             /* Levy parameter: makes medians similar */
max = 3;                   /* no flights longer than MAX units */
flight = SimLevyFlight(N, x0, c, max);

/*
title "Vectorized Levy Flight";
call series(flight[,1], flight[,2]) grid={x y} option="markers";
*/
/* write data set that contains both paths */
group = repeat( {'Gaussian', 'Levy'}, 1, N+1 );
numObs = T(0:N) // T(0:N);
x = walk[,1] // flight[,1];
y = walk[,2] // flight[,2];
create Compare var {'numObs' 'Group' 'x' 'y'};
append;
close;

QUIT;

/* overlay the paths */
ods graphics / width=480px height=480px;
title "Overlay Gaussian Walk and Levy Flight";
proc sgplot data=Compare aspect=1;
   series x=x y=y / group=group markers;
   xaxis grid values=(-2 to 12 by 2) ;
   yaxis grid values=(-2 to 12 by 2) ;
run;
