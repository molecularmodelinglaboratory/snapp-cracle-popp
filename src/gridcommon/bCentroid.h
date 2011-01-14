#ifndef BCENTROID_H
#define BCENTROID_H

#include "bDelTess.h"
#include "bPoints.h"

namespace bStd { class bCentroid; };

class bStd::bCentroid {
   typedef unsigned int uint;
   typedef unsigned int ushort;
   
   private:
      char fileBase_[64];
      char filePath_[64];

      static std::deque<float>  centroids;
      static std::deque<std::string> residues;
      static std::deque<uint>   residuesPos;
      static std::deque<uint>   chainPos;
      static std::deque<uint>   chainSize;

   public:
      bCentroid();
      bCentroid( char[], char[] );
      bCentroid( const bCentroid & );
      ~bCentroid();
   
      void setSrc( char[], char[] );
   
      bool pdb2Centroids();
      static bool pdb2Centroids( char[], char[], char[] =(NULL) );
      static bool pdb2Centroids( bPoints**&, int&, char[], char[] );
      static bool _pdb2Centroids( char[], char[] );
   
      static uint dt2Centroids( bDelTess& );
      static uint dt2Centroids( bDelTess&, bPoints*& );
      static uint dt2Centroids( bDelTess&, char[], char[] );

};

#endif

