/****************************************************************/
/*          S A S   S A M P L E   L I B R A R Y                 */
/*                                                              */
/*    NAME: LstQueue.sas                                        */
/*   TITLE: Implement a queue structure by using lists          */
/* PRODUCT: IML                                                 */
/*  SYSTEM: ALL                                                 */
/*                                                              */
/* SUPPORT: Rick Wicklin                UPDATE: July 2016       */
/*     REF:                                                     */
/*    MISC:                                                     */
/* Modules: QueueCreate, QueuePush, QueuePop                    */
/*          QueuePeek, QueueIsEmpty                             */
/****************************************************************/
proc iml;
/* implement a queue, which is a 1-D FIFO structure */
start QueueCreate( item= );
   Q = [];
   if ^IsSkipped(item) then
      call ListAddItem(Q, item);
   return Q;
finish;

/* push item onto the back of the queue */
start QueuePush(Q, item);
   call ListAddItem(Q, item);
finish;

/* pop item from the front of the queue */
start QueuePop(Q);
   A = ListGetItem(Q, 1, 'd');
   return A;
finish;

/* peek at the item at the front of the queue without removing it */
start QueuePeek(Q);
   A = ListGetItem(Q, 1, 'c');
   return A;
finish;

/* return 1 if queue is empty; 0 otherwise */
start QueueIsEmpty(Q);
   return (ListLen(Q) > 0);
finish;

store module=_all_;
quit;
