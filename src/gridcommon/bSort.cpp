// Note: this file is simply a placeholder
// for the makefile. Due to the templated
// nature of the bSort class, all methods
// are in the header file (even if untemplated)

#include "bSort.h"
using namespace bStd;

bSort::bSort() {}
bSort::~bSort() {}

/* Sort Points */
void bStd::bSort::qsortPnts(bStd::bPoints &p) {
   _qsortPnts(p.pnts_, 0, p.numPnts_ - 1);
   return;
}

namespace bStd {
/*
   
template <> void bSort::_qsort<char>( char* &list, int left, int right ) {
   if(left == right) { return; } // return if we're done
   else if(left+1 == right) { return; }
   else {}
      
   int pivot = ((right - left) >> 1) + left;
   char pivotVal = list[pivot];
      
   //~ list[pivot] = list[right]; // swap
   //~ list[right] = pivotVal;
   bSort::swap( list, pivot, right );
   pivot = right;
      
   int store = left;
   //~ T temp = 0;
   printf("here.\n");
   for( int i=left; i < right; ++i ) {
      if( list[i] < pivotVal ) {
         //~ temp = list[store];
         //~ list[store] = list[i];
         //~ list[i] = temp;
         bSort::swap( list, i, store );
         ++store;
      }
   }
   //~ list[right] = list[store];
   //~ list[store] = pivotVal;
   bSort::swap( list, pivot, store );
   bSort::_qsort( list, left, store );
   bSort::_qsort( list, store, right );
   return;
}
*/

};
