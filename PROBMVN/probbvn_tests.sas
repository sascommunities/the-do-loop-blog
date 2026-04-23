/* 2-D CDF test cases */
proc iml;
load module= _all_;
tol = 1e-7;

/* Validate 1-D rectangular Test Cases */
sigma = 1;
mu = 0;
p1 = probuvn_mod( -1, 1, sigma );
p2 = probuvn_mod( .M, 1, sigma );
p3 = probuvn_mod( -1, .I, sigma );
p4 = probuvn_mod( .M, .I, sigma );
correct = {0.6826895 0.8413447 0.8413447 1};
p = p1 || p2 || p3 || p4;
maxDiff = max(abs(p-correct));
if maxDiff > tol then
   print "--- 1-D Test Std ERROR ---" maxDiff;
else print "--- 1-D Test Std passes ---";

sigma = 2;
mu = 0.5;
p1 = probuvn_mod( -1.5, 2.5, sigma, mu );
p2 = probuvn_mod( .M, 2.5, sigma, mu );
p3 = probuvn_mod( -1.5, .I, sigma, mu );
p4 = probuvn_mod( .M, .I, sigma, mu );
correct = {0.6826895 0.8413447 0.8413447 1};
p = p1 || p2 || p3 || p4;
maxDiff = max(abs(p-correct));
if maxDiff > tol then
   print "--- 1-D Test Scale ERROR ---" maxDiff;
else print "--- 1-D Test Scale passes ---";


/* The test cases validate the correctness of the 
   bivariate normal CDF implementation. Compare the computed 
   probabilities against known values for specific regions in 
   the 2-D plane, such as half planes and quadrants. 
   The test cases cover various scenarios, including finite limits 
   and infinite limits, where missing values are used to indicate infinity. 
   Verify that the probabilities sum to 1 across a partition of the space.
*/
rho = 0.4;
x0 = 0.5; y0 = 0.8;
 
/* probability for half planes */
probLeft  = cdfbvn_std(x0, .I, rho); /* left half plane {(x,y) | x < x0} */
probRight = 1 - probLeft;            /* right half plane */
probLower = cdfbvn_std(.M, y0, rho); /* lower half plane {(x,y) | y < y0} */
probUpper = 1 - probLower;           /* upper half plane */
prob = probLeft || probRight || probLower || probUpper;
correct = {0.6914625 0.3085375 0.7881446 0.2118554};
maxDiff = max(abs(prob-correct));
if maxDiff > tol then
   print "--- 2-D Half-Plane Test ERROR ---" maxDiff;
else print "--- 2-D Half-Plane Test passes ---";

/* use half planes to compute quadrants */
P1 = cdfbvn_std(x0, y0, rho);        /* probability SW from (x0,y0) */
P2 = cdfbvn_std(.M, y0, rho) - P1;   /* probability SE from (x0,y0) */
P3 = cdfbvn_std(x0,.I, rho) - P1;    /* probability NW from (x0,y0) */
P4 = 1 - (P1 + P2 + P3);             /* probability NE from (x0,y0) */
prob = P1 || P2 || P3 || P4;
correct = {0.5896302 0.1985144 0.1018323 0.1100231};
maxDiff = max(abs(prob-correct));
if maxDiff > tol then
   print "--- 2-D Quadrant Test ERROR ---" maxDiff;
else print "--- 2-D Quadrant Test passes ---";

/* 2-D Test Cases: 9 tic-tac-toe regions */
rho = 0.4;
a = -0.5; b =  0.5;
c = -1.2; d =  0.8;
tol = 1e-12;
region = {"SW", "S", "SE", "W", "C", "E", "NW", "N", "NE"};

L = ( .M||.M ) //
    (  a||.M ) //
    (  b||.M ) //
    ( .M|| c ) //
    (  a|| c ) //
    (  b|| c ) //
    ( .M|| d ) //
    (  a|| d ) //
    (  b|| d );

U = (  a|| c ) //
    (  b|| c ) //
    ( .I|| c ) //
    (  a|| d ) //
    (  b|| d ) //
    ( .I|| d ) //
    (  a||.I ) //
    (  b||.I ) //
    ( .I||.I );

Prob = j(nrow(region), 1, .);
Correct = j(nrow(region), 1, .);
do i = 1 to nrow(region);
   Prob[i] = probbvn_std( L[i,], U[i,], rho );
end;

Correct[1] = cdfbvn_std(a, c, rho);                                            
Correct[2] = cdfbvn_std(b, c, rho) - cdfbvn_std(a, c, rho);                      
Correct[3] = CDF("Normal", c) - cdfbvn_std(b, c, rho);                         
Correct[4] = cdfbvn_std(a, d, rho) - cdfbvn_std(a, c, rho);                      
Correct[5] = cdfbvn_std(b, d, rho) - cdfbvn_std(a, d, rho)                       
           - cdfbvn_std(b, c, rho) + cdfbvn_std(a, c, rho);                      
Correct[6] = CDF("Normal", d) - CDF("Normal", c)                           
           - cdfbvn_std(b, d, rho) + cdfbvn_std(b, c, rho);                      
Correct[7] = CDF("Normal", a) - cdfbvn_std(a, d, rho);                         
Correct[8] = CDF("Normal", b) - CDF("Normal", a)                           
           - cdfbvn_std(b, d, rho) + cdfbvn_std(a, d, rho);                      
Correct[9] = 1 - CDF("Normal", b) - CDF("Normal", d)                        
           + cdfbvn_std(b, d, rho);                                            

Diff = abs(Prob - Correct);
maxDiff = max(Diff);
if maxDiff > tol then
   print "--- 2-D Region Test ERROR ---" maxDiff,
         region Prob Correct Diff;
else print "--- 2-D Region Test passes ---";

partitionErr = abs(Prob[+] - 1);
if partitionErr > tol then
   print "--- 2-D Partition Test ERROR ---" partitionErr;
else print "--- 2-D Partition Test passes ---";

QUIT;
