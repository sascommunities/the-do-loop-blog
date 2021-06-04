/****************************************************************/
/*          S A S   S A M P L E   L I B R A R Y                 */
/*                                                              */
/*    NAME: LstStack.sas                                        */
/*   TITLE: Implement a stack structure by using lists          */
/* PRODUCT: IML                                                 */
/*  SYSTEM: ALL                                                 */
/*                                                              */
/* SUPPORT: Rick Wicklin                UPDATE: July 2016       */
/*     REF:                                                     */
/*    MISC:                                                     */
/* Modules: StackCreate, StackPush, StackPop,                   */
/*          StackPeek, StackIsEmpty, StackLen                   */
/****************************************************************/
proc iml;
/* implement a stack, which is a 1-D FILO structure */
start StackCreate( item= );
   S = [];                        /* create empty list   */
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
