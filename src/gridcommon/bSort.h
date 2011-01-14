#ifndef BSORT_H
#define BSORT_H

#include <stdio.h>
#include <string.h>

#include "bPoints.h"

namespace bStd { class bSort; };

class bStd::bSort {
   public:
      bSort();
      ~bSort();
      
      template <class T> void static qsort( const T*&, int );
      template <class T> void static qsort( T*&, int );
      template <class T> void static _qsort( T*&, int, int );
      template <class T> void static _qsort( T**&, int, int );
   
      template <class T, class U> void static qsort( T* &, U*&, int, bool =(false) );
      template <class T, class U> void static _qsort( T* &, U*&, int, int );
      template <class T, class U> void static _qsortRev( T* &, U*&, int, int );

      template <class T, class U> void static qsortDbl( T*&, U*&, int );
      template <class T, class U> void static _qsortDbl( T*&, U*&,int, int );

      static void qsortPnts( bStd::bPoints & );
      template <class T>  void static _qsortPnts(T* &,int,int);

      template <class T> static void qsortPnts( const bStd::bPoints &, T* );
      template <class T, class U> void static _qsortPnts( const T*, U*, int, int );

      template <class T> static void swap( T*&, int, int );
      template <class T> static void swap( T**&, int, int );
      template <class T, class U> static void swap( T*&, U*&, int, int );
      template <class T> static void swapPnts( T*&, int, int );

   private:
}; // end class


/*
   Initial Quick Sort calls
*/

/* QUICK SORT :: list pointer */
//~ template <class T> void bStd::bSort::qsort( const T* &list, int num ) { 
   //~ _qsort( list, 0, num - 1 );
//~ }
template <class T> void bStd::bSort::qsort( T* &list, int num ) {
   _qsort( list, 0, num - 1 );
   return;
}
template <class T> void bStd::bSort::_qsort( T* &list, int left, int right ) {
   if(left == right) { return; } // return if we're done
   else if(left+1 == right) {
      if( list[left] > list[right] ) { swap( list, left, right ); }
      return;
   }
   else {}
      
   int pivot = (right + left) >> 1;
   T pivotVal = list[pivot];
   //~ if( list[left] == pivotVal || list[right] == pivotVal ) { return; }
      
   swap( list, pivot, right );
   pivot = right;
      
   int store = left;
   for( int i=left; i < right; ++i ) {
      if( list[i] == pivotVal ) {
         swap( list, i--, --right );
      }
      else if( list[i] < pivotVal ) {
         swap( list, i, store );
         ++store;
      }
      else {}
   }
   
   int newStore = store;
   if( store < right ) {
      while( right < pivot ) {
         swap( list, store, right );
         ++store; ++right;
      }
   }
   
   swap( list, pivot, store );
   _qsort( list, left, newStore );
   _qsort( list, store, right );
   return;
}

template <class T> void bStd::bSort::_qsort( T** &list, int left, int right ) {
   if(left == right) { return; } // return if we're done
   else if(left+1 == right) { return; }
   else {}
      
   int pivot = (right + left) >> 1;
   T* pivotVal = list[pivot];
   //~ if( strcmp(list[left], pivotVal) == 0 || strcmp(list[right], pivotVal) == 0 ) { return; }
      
   swap( list, pivot, right );
   pivot = right;
      
   int store = left;
   for( int i=left; i < right; ++i ) {
      int cmp = strcmp( list[i], pivotVal ) < 0;
      if( cmp == 0 ) {
         swap( list, i--, --right );
      }
      else if( cmp < 0 ) {
         swap( list, i, store );
         ++store;
      }
      else {}
   }

   int newStore = store;
   if( store < right ) {
      while( right < pivot ) {
         swap( list, store, right );
         ++store; ++right;
      }
   }

   swap( list, pivot, store );
   _qsort( list, left, newStore );
   _qsort( list, store, right );
   return;
}

/****** QSORT :: indirect */
template <class T, class U> void bStd::bSort::qsort( T* &list, U* &o, int size, bool rev ) {
   rev ? _qsortRev( list, o, 0, size - 1 ) : _qsort( list, o, 0, size - 1);
   return;
}
template <class T, class U> void bStd::bSort::_qsort( T* &list, U* &o, int left, int right ) {
   if(left == right) { return; } // return if we're done
   else if( (left + 1) == right) {
      if( list[ o[left] ] > list[ o[right] ] ) { swap( o, left, right ); }
      return;
   }
   else {}
   
   int pivot = (right + left) >> 1; // get the pivot (i.e. the middle point)
   T pivotVal = list[ o[pivot] ]; // save the value
   //~ if( list[ o[left] ] == pivotVal || list[ o[right] ] == pivotVal ) { return; }

   swap( o, pivot, right ); // move the pivot element to the end  (temporarily)
   pivot = right;

   int store = left; // start at the left 
   for(int i = left; i < right; ++i) {
      if( list[ o[i] ] == pivotVal ) {
         swap( o, i--, --right );
      }
      else if( list[ o[i] ] < pivotVal) {    // if less than pivot
         swap(o, i, store);  // swap points
         ++store;  // move the store pointer
      }
      else {}
   }

   int newStore = store;
   if( store < right ) {
      while( right < pivot ) {
         swap( o, store, right );
         ++store;
         ++right;
      }
   }

   swap( o, pivot, store ); // put the pivot back in place
   _qsort( list, o, left, newStore ); // recursively iterate through both sides
   _qsort( list, o, store, right );

   return;
}
template <class T, class U> void bStd::bSort::_qsortRev( T* &list, U* &o, int left, int right ) {
   if( left == right ) { return; } // return if we're done
   else if( (left + 1) == right ) {
      if( list[ o[left] ] < list[ o[right] ] ) { swap( o, left, right ); }
      return;
   }
   else {}
   
   int pivot = (right + left) >> 1; // get the pivot (i.e. the middle point)
   T pivotVal = list[ o[pivot] ]; // save the value
   //~ if( list[ o[left] ] == pivotVal || list[ o[right] ] == pivotVal ) { return; }

   swap( o, pivot, right ); // move the pivot element to the end  (temporarily)
   pivot = right;

   int store = left; // start at the left 
   for(int i = left; i < right; ++i) {
      if( list[ o[i] ] == pivotVal ) {
         swap( o, i--, --right );
      }
      else if( list[ o[i] ] > pivotVal) {    // if less than pivot
         swap(o, i, store);  // swap points
         ++store;  // move the store pointer
      }
      else {}
   }

   int newStore = store;
   if( store < right ) {
      while( right < pivot ) {
         swap( o, store, right );
         ++store;
         ++right;
      }
   }

   swap( o, pivot, store ); // put the pivot back in place
   _qsortRev( list, o, left, newStore ); // recursively iterate through both sides
   _qsortRev( list, o, store, right );

   return;
}


/* QUICK SORT :: double sort -- list A first, then list B */
template <class T, class U> void bStd::bSort::qsortDbl( T* &lista, U* &listb, int num ) {
   _qsortDbl( lista, listb, 0, num - 1 );
   int left = 0;
   for( int i = 0; i < num - 1; ++i ) {
      left = i;
      while( i < num - 1 && lista[i] == lista[i + 1] ) { ++i; }
      //~ printf("second sort: %d - %d\n",left,i);
      _qsortDbl( listb, lista, left, i );
   }
   return;
}
template <class T, class U> void bStd::bSort::_qsortDbl( T* &lista, U* &listb, int left, int right ) {
   if(left == right) { return; } // return if we're done
   else if( (left + 1) == right) {
      if( lista[left] > lista[right] ) { swap( lista, listb, left, right ); }
      return;
   }
   else {}
      
   int pivot = (right + left) >> 1;
   T pivotVal = lista[pivot];
   //~ if( lista[left] == pivotVal || lista[right] == pivotVal ) { return; }

   swap( lista, listb, pivot, right );
   pivot = right;

   int store = left;
   for( int i=left; i < right; ++i ) {
      if( lista[i] == pivotVal ) {
         swap( lista, listb, i--, --right );
      }
      else if( lista[i] < pivotVal ) {
         swap( lista, listb, i, store );
         ++store;
      }
      else {}
      //~ else {}
      //~ if( lista[i] <= pivotVal ) {
         //~ swap( lista, listb, i, store );
         //~ ++store;
      //~ }
   }

   int newStore = store;
   if( store < right ) {
      while( right < pivot ) {
         swap( lista, listb, store, right );
         ++store; ++right;
      }
   }

   swap( lista, listb, pivot, store );
   _qsortDbl( lista, listb, left, newStore );
   _qsortDbl( lista, listb, store, right );
   return;
}

/* QUICK SORT :: Save Sort in Index */
template <class T> void bStd::bSort::qsortPnts(const bStd::bPoints &p, T *a) {
   _qsortPnts(p.pnts_,a, 0, p.numPnts_ - 1);
   return;
}
template <class T, class U> void bStd::bSort::_qsortPnts(const T* pnts, U *a, int left, int right) {
   if(left == right) { return; } // return if we're done
   else if( (left + 1) == right) {
      if( pnts[a[left]*3] > pnts[a[right]*3] ) { swap( a, left, right ); }
      return;
   }
   else {}
   
   int pivot = ((right - left) >> 1) + left; // get the pivot (i.e. the middle point)
   T pivotVal = pnts[a[pivot]*3]; // save the value

   swap(a, pivot, right); // move the pivot element to the end  (temporarily)
   pivot = right;

   int store = left; // start at the left 
   for(int i = left; i < right; ++i) {
      if( pnts[a[i]*3] == pivotVal ) {    // if less than pivot
         swap( a, i--, --right );
      }
      else if( pnts[a[i]*3] < pivotVal ) {    // if less than pivot
         swap( a, i, store );  // swap points
         ++store;  // move the store pointer
      }
      else {}
   }

   int newStore = store;
   if( store < right ) {
      while( right < pivot ) {
         swap( a, store, right );
         ++store; ++right;
      }
   }

   swap( a, pivot, store ); // put the pivot back in place
   _qsortPnts( pnts, a, left, newStore ); // recursively iterate through both sides
   _qsortPnts( pnts, a, store, right );

   return;
}

/* SORT ROUTINES */

/* QUICK SORT :: sort points according to their x coordinate */
template <class T> void bStd::bSort::_qsortPnts(T* &pnts, int left, int right) {
   if(left == right) { return; } // return if we're done
   else if( (left + 1) == right) {
      if( pnts[left*3] > pnts[right*3] ) { swap( pnts, left, right ); }
      return;
   }
   else {}
   
   int pivot = ((right - left) >> 1) + left;
   T pivotVal = pnts[pivot*3];
   
   swapPnts(pnts, pivot, right); // move the pivot element to the end  (temporarily)
   pivot = right;
   
   int store = left; // start at the left 
   for(int i = left; i < right; ++i) {
      if( pnts[i*3] == pivotVal ) {    // if less than pivot
         swap( pnts, i--, --right );
      }
      else if( pnts[i*3] < pivotVal) {    // if less than pivot
         swapPnts( pnts, i, store );  // swap points
         ++store;                   // move the store pointer
      }
      else {}
   }
   swapPnts(pnts, pivot, store); // put the pivot back in place

   //~ std::printf("l|r|p: %.2f | %.2f | %.2f\n", pnts[left*3], pnts[right*3], pnts[pivot*3]);

   int newStore = store;
   if( store < right ) {
      while( right < pivot ) {
         swap( pnts, store, right );
         ++store; ++right;
      }
   }

   swap( pnts, pivot, store ); // put the pivot back in place
   _qsortPnts( pnts, left, newStore ); // recursively iterate through both sides
   _qsortPnts( pnts, store, right );


   return;
}



/* SWAP POINTS */
template <class T> void bStd::bSort::swapPnts(T* &pnts, int p1, int p2) {
   if(p1 == p2) { return; }
   
   p1 *= 3;
   p2 *= 3;
   T temp;
   
      
   for(int i=0; i<3; ++i) {
      temp = pnts[p1];
      pnts[p1] = pnts[p2];
      pnts[p2] = temp;
      ++p1;
      ++p2;
   }
   
   return;
}

/* SWAP POINTS :: save by index */
template <class T> void bStd::bSort::swap( T* &a, int p1, int p2 ) {
   if(p1 == p2) { return; }
   T temp = a[p1];
   a[p1] = a[p2];
   a[p2] = temp;
   return;
}

/* SWAP POINTS :: two sets */
template <class T, class U> void bStd::bSort::swap( T* &a, U* &b, int p1, int p2 ) {
   if(p1 == p2) { return; }
   T tempt = a[p1];
   a[p1] = a[p2];
   a[p2] = tempt;

   U tempu = b[p1];
   b[p1] = b[p2];
   b[p2] = tempu;
   return;
}

template <class T> void bStd::bSort::swap( T** &l, int p1, int p2 ) {
   if(p1 == p2) { return; }
   T* tempt = l[p1];
   l[p1] = l[p2];
   l[p2] = tempt;
   return;
}


#endif
