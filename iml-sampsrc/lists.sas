/****************************************************************/
/*          S A S   S A M P L E   L I B R A R Y                 */
/*                                                              */
/*    NAME: Lists.sas                                           */
/*   TITLE: Lists and Data Structures                           */
/* PRODUCT: IML                                                 */
/*  SYSTEM: ALL                                                 */
/*                                                              */
/* SUPPORT: Rick Wicklin                UPDATE: July 2016       */
/*     REF:                                                     */
/*    MISC:                                                     */
/* Modules:                                                     */
/*                                                              */
/****************************************************************/


/**************************************************
 List Utilities: The ListUtil Package
 **************************************************/

proc iml;
package load ListUtil;
package help ListUtil;   /* display overview in SAS log */
quit;


/**************************************************
 Create a Growable List of Matrices
 **************************************************/

proc iml;
/* create a list of matrices; use ListSetItem to fill */
L = ListCreate(3);            /* allocate list of 3 elements */
do n = 1 to ListLen(L);       /* for each element in list */
   A = j(n, n, n-1);          /* define n x n matrix */
   call ListSetItem(L, n, A); /* assign n_th element of L */
end;

sum = j(ListLen(L), 1);
do n = 1 to ListLen(L);       /* for each element in list */
   B = ListGetItem(L, n);     /* get n_th matrix of L */
   sum[n] = sum(B);           /* compute sum of elements */
end;

C = 1:3;
call ListAddItem(L, C);    /* add C as 4th element to L */
D = {4 3, 2 1};
call ListAddItem(L, D);    /* add 5th element to L */

package load Listutil;
run struct(L);

call ListDeleteItem(L, {1 3 5}); /* remove three elements */

quit;

/**************************************************
 Create a List of Items of Different Types
 **************************************************/

proc iml;
M = {1 2, 3 4};
C = "A":"G";
tbl = TableCreateFromDataSet("Sashelp", "Class", "obs=5");
sublist = ListCreate(3);
do i = 1 to ListLen(sublist);
   call ListSetItem(sublist, i, j(i, i, i##2));
end;

list = ListCreate();
call ListAddItem(list, {1.2 3.45 6.789});   /* add numeric vector */
call ListAddItem(list, {"Male" "Female"});  /* add character vector */
call ListAddItem(list, sublist);            /* add sublist */

package load ListUtil;                      /* load ListPrint module */
run ListPrint(list);

quit;

/**************************************************
 Create an Associative Array
 **************************************************/

proc reg data=sashelp.class plots=none;
   where sex="M";
   model weight = height;
   output out=Out p=Pred r=Res;
   ods output ParameterEstimates=PE;
quit;

proc iml;
use PE; read all var {"Variable" "Estimate"}; close;
use Out; read all var {"Weight" "Pred" "Res"}; close;

StructNames = {"Variable" "Estimate" "DepVar" "Predicted" "Residual"};
RegModel = ListCreate( StructNames );
call ListSetItem(RegModel, "Variable", Variable);
call ListSetItem(RegModel, "Estimate", Estimate);
call ListSetItem(RegModel, "DepVar", Weight);
call ListSetItem(RegModel, "Predicted", Pred);
call ListSetItem(RegModel, "Residual", Res);

package load ListUtil;                      /* load ListPrint module */
run ListPrint(RegModel);

/* Module that creates a plot of observed vs predicted response.
   Pass in a list that contains elements named "DepVar"and "Predicted" */
start PredPlot(L);
   Observed = ListGetItem(L, "DepVar");
   Predicted = ListGetItem(L, "Predicted");
   call Scatter(Observed, Predicted) procopt="noautolegend"
                    other="lineparm x=0 y=0 slope=1 / clip";
finish;

run PredPlot(RegModel);

quit;

/**************************************************
 Create a List of Lists
 **************************************************/

proc iml;
use Sashelp.Class;           /* read data */
read all var {Age Sex Name};
close Sashelp.Class;

Age = char(Age, 2);          /* convert to character vector */
ages = unique(Age);
L = ListCreate(ages);        /* outer list: elements named "11":"16" */

gender = unique(Sex);
K = ListCreate(gender);      /* inner list: elements named {"F" "M"} */

do i = 1 to ncol(ages);      /* For each age level... */
   idx = loc(Age=ages[i]);   /* Find observations for this age */
   do j = 1 to ncol(gender);              /* for each gender... */
      jdx = loc( Sex[idx]=gender[j] );    /* find this age and gender */
      if ncol(jdx)=0 then students = {};  /* no students found */
      else students = Name[idx[jdx]];     /* get student names */
      call ListSetItem(K, gender[j], students);  /* value of inner list */
   end;
   call ListSetItem(L, ages[i], K);  /* set sublist as value */
end;

quit;

/**************************************************
 Construct a Stack
 **************************************************/

proc iml;
/* implement a stack, which is a 1-D FILO structure */
start StackCreate( item= );
   S = ListCreate();              /* create empty list   */
   if ^IsSkipped(item) then       /* if item specified,  */
      call ListAddItem(S, item);  /*    add item to list */
   return S;
finish;

/* push an item onto the stack */
start StackPush(S, item);
   call ListAddItem(S, item);     /* add item to the end */
finish;

/* pop an item from the stack */
start StackPop(S);
   A = ListGetItem(S, ListLen(S), 'd'); /* get & remove last item */
   return A;
finish;

/* peek at the item at the top of the stack without removing it */
start StackPeek(S);
   A = ListGetItem(S, ListLen(S), 'c'); /* get last item */
   return A;
finish;

/* return 1 if stack is empty; 0 otherwise */
start StackIsEmpty(S);
   return (ListLen(S) = 0);
finish;

/* return number of elements in stack */
start StackLen(S);
   return ListLen(S);
finish;

store module=_all_;
quit;


/**************************************************
 Reverse the Words in a Sentence
 **************************************************/

proc iml;
load module=(StackCreate StackPush StackPop
             StackPeek StackIsEmpty StackLen);

/* Create sentence. Break into vector of words. */
str = "Now is the time for all good men to come to the aid of their party.";
n = countw(str, " .");         /* use blanks and period as delimiters */
words = scan(str, 1:n, " .");  /* character vector */

S = StackCreate();             /* create an empty stack */
do i = 1 to ncol(words);
   run StackPush(S, words[i]); /* push each element onto the stack */
end;

print (StackPeek(S))[L="Top of Stack"]; /* the last word is on top */

/* retrieve the data in reverse order */
w = j(1, StackLen(S), "     ");
do i = 1 to ncol(w);           /* pop each element; insert into stack */
   w[i] = StackPop(S);
end;
print w[L="Reversed Words"];

if StackIsEmpty(S) then
   print "Stack is empty";
else print (StackPeek(S))[L="Top of Stack"];

quit;

/**************************************************
 Implement a Postfix Calculator
 **************************************************/

proc iml;
load module=(StackCreate StackPush StackPop);

/* Given a binary operator, return  the expression
  (L op R) where op is in the set {+, -, *, /} */
start BinaryCalc(operator, L, R);
   if      operator="+" then return L + R;
   else if operator="-" then return L - R;
   else if operator="*" then return L * R;
   else if operator="/" then return L / R;
   else return .;
finish;

/* Input a space-separated string that represents a valid
   arithmetic operation in postfix notation. The string
   can contain numbers and the binary operators {+, -, *, /}.
   The string must represent a valid operation; no error
   checking is performed. */
start PostfixCalc(str);
   n = countw(str, " ");
   tokens = scan(str, 1:n, " ");  /* character vector */
   binOps = {"+","-","*","/"};    /* four binary operators */
   S = StackCreate();             /* create an empty stack */
   do i = 1 to ncol(tokens);
      token = tokens[i];          /* get the token */
      if element(token, binOps) then do; /* it's binary op */
         R = StackPop(S);         /* retrieve the previous */
         L = StackPop(S);         /*    two numbers        */
         result = BinaryCalc(token, num(L), num(R));
         run StackPush(S, char(result)); /* push result on stack */
      end;
      else
         run StackPush(S, token);  /* push number on stack */
   end;
   return num(StackPop(S));        /* return result as number */
finish;

/* examples of parsing postfix expressions */
str = {"2 2.8 7.2 + * 5 /",    /* 2*(2.8+7.2) / 5  =   4 */
       "4 5 7 2 + - *",        /* 4*(5 - (7+2))    = -16 */
       "4 -5 + 6 -2 -  *",     /* (4 + -5)*(6 - -2)=  -8 */
       "2 2 2 2 * * *"    };   /*  2**4            =  16 */
result = j(nrow(str), 1);
do i = 1 to nrow(str);
   result[i] = PostfixCalc(str[i]);
end;
print str result;

quit;

/**************************************************
 Construct a Binary Search Tree
 **************************************************/

/* L[i] is key value; L[2] is left child; L[3] is right child */
proc format;
    value BSTFmt  1='Key'  2='Left'  3='Right';
run;

proc iml;
/* A node is a three-element list:
   node[1] contains the KEY   value
   node[2] contains the LEFT  value (or empty if null)
   node[3] contains the RIGHT value (or empty if null) */
start BSTNewNode(value);
   node = ListCreate(3);     /* create list with 3 null elements */
   call ListSetItem(node, 1, value);   /* set KEY value */
   return node;
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
      if value = ListGetItem(T, KEY) then
         return 1;                  /* found it: return path to subitem */
      else if value < ListGetItem(T, KEY) then do;
         path = path || LEFT;       /* add to path */
         T = ListGetItem(T, LEFT);  /* new root is left child */
      end;
      else do;
         path = path || RIGHT;      /* add to path */
         T = ListGetItem(T, RIGHT); /* new root is right child */
      end;
      if type(T)='U' then
         return 0;                  /* not found: return path to subitem */
   end;
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

x = {5 3 1 9 1 6 4}`;
bst = BSTCreate(x);

found = BSTLookup(path, bst, 6);
print found[L="Was 6 found?"], path[L="Path from root" F=BSTFmt.];
found = BSTLookup(path, bst, 10);
print found[L="Was 10 found?"], path[L="Path from root" F=BSTFmt.];
quit;

/**************************************************
 Plot a Binary Search Tree
 **************************************************/

%include sampsrc(LstBST.sas);      /* define modules */
proc iml;
load module = _all_;                /* load modules */
x = {5 3 1 9 1 6 4}`;
bst = BSTCreate(x);
title "Diagram of Binary Search Tree";
call BSTPlot(bst);
quit;
