/* Successive Difference Replication Method 

Stephen Ash
"Using successive difference replication for estimating variances"
Survey Methodology, June 2014 47
Vol. 40, No. 1, pp. 47-59
Statistics Canada, Catalogue No. 12-001-X

RA (Row Assignment schemes) according to Ash (2014, p. 54-55)
Ash suggests two simple RAs:

Let k be the order (=number of rows) of a Hadamerd matrix 
   RA1: Modular cycles, given by the 5 steps on p. 54.
   RA2: Generate the simple RA 
        (1,2)(2,3)(3,4)...(k,1)
        Repeat this RA until you have obtained assignments for 
        all n rows.
*/
proc iml;
start RA2(n, k);
   w = ceil(n / k);  /* make w copies of 1:k */
   s1 = T(1:k);
   s2 = s1[2:k] // s1[1];
   R = repeat( s1||s2, w);
   return R[1:n, ];
finish;

/* Helper Functions for the RA1 assignment scheme 
   Generate one or more cycles. Given a vector 1:k,
   takes steps (mod k) starting at beginVal in 
   increments of Step.
   Stop when we return to the initial value and return
   the vector of values in the first column. The second
   column is the LAG of the first column.

   k = number of rows in Hadamard matrix
   beginVal = starting value, beginVal <= k
   Step = step size, Step < k

   Ex: RACycleLoop(8, 1, 1)  => first col 1:8
       RACycleLoop(8, 1, 2)  => {1,3,5,7}
       RACycleLoop(8, 1, 4)  => {1,5}
*/
start RACycleLoop(k, beginVal, Step);
   r1 = beginVal;
   currVal = beginVal;
   nextVal = .;
   
   do cnt = 1 to k until (nextVal = beginVal);
      nextVal = mod(currVal + Step, k);
      if nextVal = 0 then nextVal = k;
      currVal = nextVal;
      if nextVal ^= beginVal then 
         r1 = r1 // nextVal;
   end;
   r2 = r1[2:nrow(r1)] // r1[1];
   return (r1 || r2);
finish;

/* Helper: RACycle
   Input:
   k = number of rows in Hadamard matrix
   Step = step size
   Return: kx2 matrix. First column is a set of cycles with 
   increasing start values.  The second
   column is the LAG of each cycle.
   Ex: RACycle(8, 4) => {1,5,2,6,3,7,4,8}
*/
start RACycle(k, _Step);
   Step = _Step;
   /* if Step >= k, MOD by k. If Step=k, reset to Step = 1 */
   if Step >= k then 
      Step = mod(Step, k);
   if Step=0 then Step = 1;
   if mod(k, Step)=0 then do; /* Step divides k evenly */
      numLoops = Step;   
      sizeLoops = k / Step;   
      R = j(k, 2, .);
      do i = 1 to numLoops;
         i1 = 1 + sizeLoops*(i-1);
         i2 = sizeLoops*i;
         R[i1:i2, ] = RACycleLoop(k, i, Step);
      end;
   end;
   else 
      R = RACycleLoop(k, 1, Step);
   return R;
finish;

/* RA1 : Assignment scheme Ash (2014, p. 54)
   n = total number of pairs to generate
   k = order of Hadamard matrix 

   The algorithm returns an nx2 matrix where the first column 
   is a set of results from RACycle with increasing step sizes.
*/
start RA1(n, k);
   if n <= k then do;
      R = RACycle(k, 1);
      return R[1:n,];
   end;
   /* if n > k, create multiple pairs for d=1,2,..k-1 */
   G = j(n, 2, .);
   numCycles = ceil(n / k);
   sizeLoops = k;
   do i = 1 to numCycles;
      i1 = 1 + sizeLoops*(i-1);
      u = sizeLoops*i;
      i2 = min(n, u);
      gap = u - i2;
      R = RACycle(k, i);
      if gap=0 then 
         G[i1:i2, ] = R;
      else 
         G[i1:i2, ] = R[1:nrow(R)-gap, ];
   end;
   return G;
finish;

store module=(RACycleLoop RACycle RA1 RA2);
QUIT;


proc iml;
load module=(RACycleLoop RACycle RA1 RA2);

print "=== Test the RA2 method ===";
R = RA2(10, 4);
print R;
R = RA2(3, 4);
print R;

print "=== Test RACycleLoop ===";
k = 8;  /* order of Hadamard k <= n */
d = 1;  /* step */
beginVal = 1;
R = RACycleLoop(k, beginVal, d);  print k beginVal d, R;
beginVal = 2;
R = RACycleLoop(k, beginVal, d);  print k beginVal d, R;

d = 2;
beginVal = 1;
R = RACycleLoop(k, beginVal, d);  print k beginVal d, R;
beginVal = 2;
R = RACycleLoop(k, beginVal, d);  print k beginVal d, R;
d = 4;
beginVal = 1;
R = RACycleLoop(k, beginVal, d);  print k beginVal d, R;

print "=== Test RACycle ===";
k = 8;  /* order of Hadamard k <= n */
d = 1;
R = RACycle(k, d);  print k d, R;
d = 2;
R = RACycle(k, d);  print k d, R;
d = 3;
R = RACycle(k, d);  print k d, R;
d = 4;
R = RACycle(k, d);  print k d, R;
/* if d = k, use d=1.
   if d > k, subtract multiple of k. Ex: 11 becomes mod(11,8)=3 */
d = 11;   /* equivalent to d=3 */
R = RACycle(k, d);  print k d, R;

print "=== Test RA1 ===";
n = 7;
k = 8;
pairs = RA1(n, k);
print n k, pairs[r=(1:n)];

n = 14;
k = 16;
pairs = RA1(n, k);
print n k, pairs[r=(1:n)];
n = 23;
pairs = RA1(n, k);
print n k, pairs[r=(1:n)];

n = 20;
k = 16;
pairs = RA1(n, k);
print n k, pairs[r=(1:n)];

n = 23;
pairs = RA1(n, k);
print n k, pairs[r=(1:n)];
QUIT;

/* Ash, p. 52, shows an example that does not use
   the RA1 or RA2 scheme, but instead uses a custom assignment.
   Run the example for the RA1 and RA2 scheme.
*/
proc iml;
load module=(RACycleLoop RACycle RA1 RA2);

/* Ash uses n=14 and k=16 rows. The nonnormal Hadamard
   matrix is the Kronecker product of two 4x4 matrices 
   For this example, the RA1 and RA2 schemes give the same
   replicates.
*/
H4a = hadamard(4);
H4b = H4a;
H4b[2, ] = -H4b[2, ]; 
H4b[ ,2] = -H4b[ ,2];
H = H4a @ H4b;

print "=== Ash Example 1 using RA1 scheme ===";
n = 14;
k = nrow(H);
pairs = RA1(n, k);
print n k, pairs[r=(1:n) L='RA1'];

print "=== Show how to get replicates from the RA scheme ===";
Ha = H[pairs[,1], ];
Hb = H[pairs[,2], ];
Replicates = 1 + 2**(-3/2) * (Ha - Hb);
print replicates[r=(1:n) c=(1:k) F=3.1];


print "=== Ash Example 1 using RA2 scheme ===";
pairs = RA2(n, k);
print n k, pairs[r=(1:n) L='RA2'];

print "=== Show how to get replicates from the RA scheme ===";
Ha = H[pairs[,1], ];
Hb = H[pairs[,2], ];
Replicates = 1 + 2**(-3/2) * (Ha - Hb);
print replicates[r=(1:n) c=(1:k) F=3.1];

/***********************************/
print "=== Test a larger example that has n=3*k, which is divisible by k ===";
n = 3*k;
pairs1 = RA1(n, k);
pairs2 = RA2(n, k);
print n k, (pairs1 || pairs2)[r=(1:n) 
           c={'RA1_1' 'RA1_2' 'RA2_1' 'RA2_2'} L='RA1 vs RA2'];

QUIT;
