/* SAS program to accompany the article 
   "Implement a SMOTE simulation algorithm in SAS"
   by Rick Wicklin, published 05MAY2025 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2025/05/05/smote-simulation-sas.html

   This program shows how to implement a simple SMOTE algorithm in SAS 
   from "first principals." The algorithm assumes all variables are continuous.

   First, define and store HELPER FUNCTIONS: 
          NearestNbr RandRowsAndNbrs SMOTESimContin
   For convenience, put the definitions into a separate file.
   USAGE:
 
   %INCLUDE "SMOTE_define.sas";
   proc iml;
   load module=(NearestNbr RandRowsAndNbrs SMOTESimContin);
   ...use the functions...
*/

proc iml;
/* Compute indices (row numbers) of k nearest neighbors.
   INPUT:  X    an (N x d) data matrix
           k    specifies the number of nearest neighbors (k>=1) 
   OUTPUT: NN   an (N x k) matrix of row numbers. NN[,j] contains
                the row numbers (in X) of the j_th closest neighbors
           dist an (N x k) matrix. dist[,j] contains the distances
                between X and the j_th closest neighbors
   See https://blogs.sas.com/content/iml/2016/09/14/nearest-neighbors-sas.html
*/
start NearestNbr(NN, dist, X, k=1);
   N = nrow(X);
   NN = j(N, k, .);             /* j_th col contains obs numbers for j_th closest obs */
   dist = j(N, k, .);           /* j_th col contains distance for j_th closest obs */
   D = distance(X);
   D[do(1, N*N, N+1)] = .;      /* set diagonal elements to missing */
   do j = 1 to k;
      dist[,j] = D[ ,><];       /* smallest distance in each row */
      NN[,j] = D[ ,>:<];        /* column of smallest distance in each row */
      if j < k then do;         /* prepare for next closest neighbors */
         ndx = sub2ndx(dimension(D), T(1:N)||NN[,j]);
         D[ndx] = .;            /* set elements to missing */
      end;      
   end;
finish;
/* Randomly select points and a nearest neighbor to the selected point.
   In the data matrix, X, each point is a row. We will name each 
   point by using the row number: X1, X2, ..., XN
   
   Input: An Nxk integer matrix of nearest neighbor information
          NN[i,j] is the name (1-N) of j_th nearest neighbor to the i_th point
          
   Output: Return a Tx2 integer matrix, R. The first column of R 
          contains random rows of X. The second column contains random
          rows among the k nearest neighbors.  
          That is, for i=1..T, 
          1. Randomly choose a row of X
          2. Randomly choose one of the k nearest neighbors
   The columns of R can be used to select pairs of neighbors:  
   X[R[,1], ] is the set of selected points
   X[R[,2], ] is the set of selected nearest neighbors
*/
start RandRowsAndNbrs(T, NN);
   N = nrow(NN);
   k = ncol(NN);
   rand_row = randfun(T, "Integer", 1, N);  /* random row of nn matrix */
   rand_col = randfun(T, "Integer", 1, k);  /* random col of nn matrix */
   rand_subs = rand_row || rand_col;        /* => random [i,j] position of nn matrix */
   rand_idx = sub2ndx(N||k, rand_subs);     /* => random index into nn matrix */
   nn_row = nn[rand_idx];
   return rand_row || nn_row;
finish;

/* The simulation step of the SMOTE algorithm, assuming that all data are 
   continuous. 
   Input: X is an Nxd data matrix. Each row is a data point.
          T is the integer number of new synthetic points to generate
          k is the number of nearest neighbors to use
   Output: A Txd 
*/
start SMOTESimContin(X, T, k=5);
   run NearestNbr(NN, dist, X, k);
   R = RandRowsAndNbrs(T, NN);
   p1 = X[R[,1], ];
   p2 = X[R[,2], ];
   u = randfun(T, "Uniform");                  * random variates in U[(0,1) ;
   p_new = p1 + u#(p2 - p1);
   return p_new;
finish;

store module=(NearestNbr RandRowsAndNbrs SMOTESimContin);
quit;
