/* SAS/IML program to accompany the article
   "How to create biplots in SAS"
   by Rick Wicklin, published 13NOV2019 on The DO Loop blog:
   https://blogs.sas.com/content/iml/2019/11/30/create-biplots-in-sas.html

   This program definee three SAS/IML modules that you an use to create biplots:
   1. The CalcBiplot module returns the matrices that contain the projection of
      the observations (scores) and the projection of the variables (loadings)
      onto the first few principal components.
*/



/***************************************************/

proc iml;
/*
start CalcBiplot(outObs, OutVars, OutSingVals, mX, dim, FacType, StdMethod, scale=0 );
*/
start CalcBiplot(outObs,      /* (output) an (N x dim) matrix of observations for the biplot. 
                                 This is the first dim columns of G (see below) */
                 OutVars,     /* (output) a (p x dim) matrix of ordered pairs representing the endpts 
                                 of variable vectors. This is the first dim columns of H (see below) */
                 OutSingVals, /* (output) the vector of singular values of the X matrix. This tells you
                                 how much variance is accounted for in the biplot. */
                 mX,          /* (input) the Nxp data matrix */
                 dim,         /* (input) the dimensions for the biplot (usually dim=2) */
                 FacType,     /* (input) Determines the scalings of the observations and variable vectors.
                                  The SVD is X = ULV` = [ U*(L##alpha) ] * [ (L##(1-alpha))*V` ]
                                               =        G         *       H`
                              If FacType is a string, there are four valid values:
                              'GH'   alpha=0.  1) The length of the variable vectors (cols of H)
                                            are proportional to the variances of the variables.
                                            2) Euclidean distance between i_th and j_th rows of G
                                            is proportional to the Mahalanobis distance between i_th and j_th
                                            observations in the data set.
                              'COV'   alpha=0 but also multiply G by sqrt(N-1) and divide H by the same.
                                            1) The length of the variable vectors (cols of H)
                                            are EQUAL to the variances of the variables.
                                            2) Euclidean distance between i_th and j_th rows of G
                                            is EQUAL to the Mahalanobis distance between i_th and j_th
                                            observations in the data set.
                              'JK'   alpha=1.   1) Positions of points in biplot identical to points in score
                                            plot of first two principal components.
                                            2) Euclidean distance between i_th and j_th rows of G
                                            is EQUAL to the Euclidean distance between i_th and j_th
                                            observations in the data set.
                              'SYM'  alpha=0.5  Observations and variables are treated symmetrically.

                              If FacType is numerical, alpha = FacType and 0 <= FacType <= 1 */

                 StdMethod,   /* (input)    StdMethod a string with three valid values
                                 'GrandMean'   adjust X by the grand mean X[:]
                                 'Mean'        adjust each column of X by the mean of that column X[:,i]
                                 'Std'         standardize by centering X and dividing each column by its variance */
                 scale=1 );   /* (input)    If 0 or empty, use automatic scaling between vectors and obs
                                            If 1, do not scale
                                            Otherwise, scale vectors by the specified scale factor */
   if type(FacType)='N' then
        alpha = FacType;
   else if type(FacType)='C' then do;
     factorization = upcase(FacType);
     if      factorization='GH'  then alpha=0;
     else if factorization='SYM' then alpha=0.5;
     else if factorization='JK'  then alpha=1;
     else if factorization='COV' then alpha=0;
     else stop "ERROR: Invalid biplot scaling";
   end;
   else do;
      factorization='COV';
      alpha=0;
   end;

   /* https://blogs.sas.com/content/iml/2015/02/23/complete-cases.html */
   /* return rows that have no missing values */
   NonMissingRows = loc(countmiss(mX, "row")=0);
   if ncol(NonMissingRows)=0 then 
      stop "ERROR: Every row of the data matrix contains a missing value. Cannot continue.";
   bMissingValues = (ncol(NonMissingRows) ^= nrow(X));
   X = mX[NonMissingRows,];         /* use complete cases */
   N = nrow(X);

   std=upcase(StdMethod);
   if std = 'NONE' then
      ;
   else if std = 'GRANDMEAN' then
      X = X - X[:];             /* remove grand mean */
   else
      X = X - mean(X);          /* remove column means */
   if std = 'STD' then do;
      S = std(X);
      /* if any column is constant, can't standardize (this is PROC STDIZE behavior) */
      idx = loc( S>0 );
      if ncol(idx)=0 then return;         /* all constant */
      X[,idx] = X[,idx] # (1 / S[,idx]);  /* X * diag(1/S) */
   end;

   /* X = U diag(L) V` */
   call svd(U,L,V,X);
   OutSingVals = L; /* return full spectrum */

   U = U[,1:dim];
   V = V[,1:dim];
   L = L[1:dim];

   if alpha=0 then do;
      A = U;
      B = V # L`;                   /* V * diag(L) */
   end;
   else if alpha=1 then do;
      A = U # L`;                   /* U * diag(L); */
      B = V;
   end;
   else do;
      A = U # (L ## alpha)`;       /* U*diag(L ## alpha ) */
      B = V # (L ## (1-alpha))`;   /* V*diag(L ## (1-alpha) ) */
   end;
   if factorization = 'COV' then do;
      A = sqrt(N-1) # A ;
      B = B / sqrt(N-1);
   end;

   /* Optionally scale the vectors to aid visualization.
      If scale=0, perform automatic scaling based on ange of data. */
   if scale<=0 then 
      sc = max(sqrt(A[,##])) /  max(sqrt(B[,##]));
   else sc = scale;         /* use manual scaling provided by user */
   print ("Scale Factor: " + strip(char(sc, 5,2)));

   B = B # sc;                                                          

   /* set return arguments */
   if bMissingValues then do;
      mA = A;
      A = j( nrow(mX), ncol(mA), .);
      A[ NonMissingRows, ] = mA;
   end;
   OutObs  = A;
   OutVars = B;
finish CalcBiplot;

/* The WriteBiplot moduile writes biplot data for plotting:
   1. Create _Scores, which contains the projected observations
   2. Create _Vectors, which contain the projected variables
   3. Create macro variables &MinAxis and &MaxAxis.
      MinAxis is the minimum of all X and Y coordinates.
      MaxAxis is the maximum of all X and Y coordinates.
*/
start WriteBiplot(X,          /* The data matrix */
                  ID,         /* The ID variable to label rows of X.
                                 If you pass in an empty matrix, 
                                 obs numbers are used */
                  _Variable_, /* Names of columns of X */
                  FacType,    /* 'GH', 'COV', 'JK', or 'SYM' */
                  StdMethod,  /* 'None', 'Mean', or 'Std' */
                  Scale=1);   /* Additional scaling applied to vectors */
   dim = 2;
   run CalcBiplot(A, B, L, X, dim, FacType, StdMethod, Scale);
   percent = L##2 / L[##];     /* percent variation explained */
   values = putn(percent, "PERCENT7.2" );

   /* Columns of A are PCs; write to data set */
   PCNames = {'Prin1' 'Prin2'};
   labs = {'Dimension 1' 'Dimension 2'};
   labels = labs + " (" + values[1:dim]` + ")";
   Tbl = TableCreate(PCNames, A);
   call TableSetVarLabel(Tbl, PCNames, labels);  /* use table to supply labels */
   if ncol(ID)>0 then
      call TableAddVar(Tbl, "_ID_", ID);
   else 
      call TableAddVar(Tbl, "_ID_", T(1:nrow(A)));  /* use row numbers if no ID variable */
   call TableWriteToDataSet(Tbl, "_Scores");

   /* columns of B are scaled vectors; write to data set */
   create _Vectors from B[rowname=_Variable_ colname={'vx' 'vy'}];
   append from B[rowname=_Variable_];
   close;

   /* best to use same scale on both axes. Create macro variables. */
   minAxis = min(A, B);
   call symputx("MinAxis", minAxis);
   maxAxis = max(A, B);
   call symputx("MaxAxis", maxAxis);
finish WriteBiplot;


/* The Biplot modules creates biplot by using ODS graphics.
   By default, vectors are NOT scaled and observations are NOT labeled. 
   These are different defaults than the %BIPLOT macro, which 
   scales vectors and labels observations by default. */
start Biplot(X,          /* The data matrix */
             ID,         /* The ID variable to label rows of X.
                            If you pass in an empty matrix, 
                            obs numbers are used */
             _Variable_, /* Names of columns of X */
             FacType,    /* 'GH', 'COV', 'JK', or 'SYM' */
             StdMethod,  /* 'None', 'Mean', or 'Std' */
             Scale=1,    /* Additional scaling applied to vectors */
             labelPoints=0); /* should points be labeled by ID variable? */

   call WriteBiplot(X, ID, _Variable_, FacType, StdMethod, Scale);
   paramStr = "FacType=" + strip(FacType) + "; Std=" + strip(StdMethod) + ";";
   labelStr = choose(labelPoints, "datalabel=_ID_", " ");
   submit paramStr labelStr;
      data Biplot;
         set _Scores _Vectors;
      run;
      title2 "&paramStr";
      proc sgplot data=Biplot aspect=1 noautolegend;
         refline 0 / axis=x; refline 0 / axis=y;
         scatter x=Prin1 y=Prin2 / &labelStr; 
         vector  x=vx    y=vy    / datalabel=_Variable_
                                   lineattrs=GraphData2 datalabelattrs=GraphData2;
         xaxis grid offsetmin=0.1 offsetmax=0.1 min=&minAxis max=&maxAxis;
         yaxis grid min=&minAxis max=&maxAxis;
      run;
   endsubmit;
finish Biplot;

store module=(CalcBiplot WriteBiplot Biplot);
QUIT;
