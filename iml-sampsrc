/****************************************************************/
/*          S A S   S A M P L E   L I B R A R Y                 */
/*                                                              */
/*    NAME: LstBST.sas                                          */
/*   TITLE: Implement a binary search tree by using lists       */
/* PRODUCT: IML                                                 */
/*  SYSTEM: ALL                                                 */
/*                                                              */
/* SUPPORT: Rick Wicklin                UPDATE: July 2016       */
/*     REF:                                                     */
/*    MISC:                                                     */
/* Modules: BSTNewNode, BSTCreate, BSTInsert, BSTLookup,        */
/*          BSTGetKeyDepth, BSTDepth, BSTGetEdges, BSTPlot,     */
/*          BSTNewNode, BSTLookup, BSTCreate, BSTInsert         */
/****************************************************************/
%include sampsrc(LstQueue.sas);
proc iml;
load module=(QueueCreate QueuePush QueuePop);

/* A node is a three-item list:
   node[1] contains the KEY   value
   node[2] contains the LEFT  value (or empty if null)
   node[3] contains the RIGHT value (or empty if null) */
start BSTNewNode(value);
   node = ListCreate(3);     /* create list with 3 null items */
   call ListSetItem(node, 1, value);   /* set KEY value */
   return node;
finish;

/* pass a vector of key values to this routine to create a BST
   that has those values as keys */
start BSTCreate(x);
   bst =  BSTNewNode( x[1] );
   do i = 2 to nrow(colvec(x));
      run BSTInsert(bst, x[i]);
   end;
   return bst;
finish;

/* Insert a new branch for a key value in a BST. If the value
   already exists, do nothing (so there are never duplicates) */
start BSTInsert(root, value);
   if ListLen(root)=0 then do;   /* List empty. Set root node */
      root = BSTNewNode(value);
      return;
   end;
   /* otherwise, search tree to find value */
   found = BSTLookup(path, root, value);  /* if found, return */
   if ^found then                         /* else add to sub-path */
      call ListSetSubItem(root, path, BSTNewNode(value));
finish;


/* Search for a target value in a binary search tree.
   Input: root is the root node of a BST
          value is the target value
   Output: path contains the path to the node that contains the target
          value or the node where the target value can be inserted.
   Return: 1 if the target value is in the tree; 0 otherwise */
start BSTLookup(path, root, value);
   KEY = 1; LEFT = 2; RIGHT = 3;
   path = {};
   T = root;
   do while (1);
      if value = T$KEY then
         return 1;                  /* found it: return path to subitem */
      else if value < T$KEY then do;
         path = path || LEFT;       /* add to path */
         T = T$LEFT;                /* new root is left child */
      end;
      else do;
         path = path || RIGHT;      /* add to path */
         T = T$RIGHT;               /* new root is right child */
      end;
      if type(T)='U' then
         return 0;                  /* not found: return path to subitem */
   end;
finish;

/******************************************************/
/* iterative algorithm to find depth of a binary tree */
/******************************************************/
start BSTGetKeyDepth(root);
   KEY = 1; LEFT = 2; RIGHT = 3;
   v = {};        /* vector of key values (assume conformal matrices) */
   if ListLen(root)=0 then return v;
   Q = QueueCreate(root);
   depth = 0;

   do while (1);
      /* nodeCount (queue size) = number of nodes at current level */
      nodeCount = ListLen(Q);
      if nodeCount = 0 then return v;
      depth = depth + 1;
      /* pop all nodes of current level and push all nodes of next level */
      do while (nodeCount > 0);
         node = QueuePop(Q);
         v = v // (node$KEY || depth);
         child = node$LEFT;
         if type(child)^='U' then
            call QueuePush(Q, child);
         child = node$RIGHT;
         if type(child)^='U' then
            call QueuePush(Q, child);
         nodeCount = nodeCount - 1;
      end;
   end;
finish;

/* find depth of a binary tree */
start BSTDepth(root);
   v = BSTGetKeyDepth(root); /* (key, depth) */
   return max(v[,2]);
finish;

/* iterative algorithm to find edges that connect the nodes of a binary tree.
   Return edges as a two-column matrix (from, to) of nodes.
   The nodes are enumerated from top to bottom and from left to right. */
start BSTGetEdges(root);
   KEY = 1; LEFT = 2; RIGHT = 3;
   v = {};        /* vector of key values (assume conformal matrices) */
   if ListLen(root)=0 then return v;
   Q = QueueCreate(root);
   do while (1);
      /* nodeCount (queue size) = number of nodes at current level */
      nodeCount = ListLen(Q);
      if nodeCount = 0 then return v;
      /* pop all nodes of current level and push all nodes of next level */
      do while (nodeCount > 0);
         node = QueuePop(Q);
         parentID = node$KEY;
         child = node$LEFT;
         if type(child)^='U' then do;
            call QueuePush(Q, child);
            childID = child$KEY;
            v = v // (parentID || childID);
         end;
         child = node$RIGHT;
         if type(child)^='U' then do;
            call QueuePush(Q, child);
            childID = child$KEY;
            v = v // (parentID || childID);
         end;
         nodeCount = nodeCount - 1;
      end;
   end;
finish;

/* plot a binary search tree */
start BSTPlot(bst);
   v = BSTGetKeyDepth(bst);
   NodeID = v[,1];
   x = rank(NodeID);
   y = v[,2];
   create nodes var {"NodeID" "x" "y"}; append; close;

   edges = BSTGetEdges(bst);
   FromNode = edges[,1];
   ToNode = edges[,2];
   LinkID = T(1:nrow(FromNode));
   N = 2*nrow(FromNode);
   LX = j(N,1,.);
   LY = j(N,1,.);
   do i = 1 to nrow(FromNode);
      j = loc(NodeID = FromNode[i]);  /* index into NodeID */
      k = 2*i - 1;
      LX[k,] = X[j];
      LY[k,] = Y[j];
      j = loc(NodeId = ToNode[i]);
      LX[k+1,] = X[j];
      LY[k+1,] = Y[j];
      *print i (LX || LY);
   end;
   LinkId = LinkId || LinkID;
   create LinkCoords var {"LinkID" "LX" "LY"}; append; close;

   submit;
   data All;
   set nodes LinkCoords;
   run;

   proc sgplot data=All noautolegend nocycleattrs;
   format nodeID 3.;
   series x=LX y=LY / group=LinkID lineattrs=(color=black pattern=solid);
   scatter x=x y=y /  markerattrs=(symbol=circlefilled size=32)
                      filledoutlinedmarkers markerfillattrs=(color=white)
                      markeroutlineattrs=(thickness=2);
   scatter x=x y=y / markerchar=nodeID;
   yaxis reverse display=none;
   xaxis display=none;
   run;
   endsubmit;

   call delete({Nodes LinkCoords All});
finish;

store module=_all_;
quit;
