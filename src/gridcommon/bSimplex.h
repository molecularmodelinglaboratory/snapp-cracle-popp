#ifndef BSIMPLEX_H
#define BSIMPLEX_H

#include <stdio.h>
#include "bVBase.h"
#include "bHex.h"

namespace bStd { class bSimplex; };

class bStd::bSimplex : virtual bStd::bVBase {
   friend class bDelTess;
   friend class bCentroid;

   private:
      
      static short count_;
      unsigned short id_;
      //~ double vrtx_[20]; // vertex coordinates [p_h, p_x, p_y, p_z, p_q]
      double det_[5]; // determinant
      bSimplex *nghbr_[4]; // pointer to neighbors

   protected:

   public:
      // Constructors
      bSimplex();
      bSimplex( const bSimplex& );
      ~bSimplex();

      // Assignment Operators
      bSimplex& operator=( const bSimplex& );
      bSimplex& operator=( const double* );

      // Conditional Operators
      bool operator==( const bSimplex& );
      bool operator!=( const bSimplex& );
      short  inSphere( const double* );

      // Setup
      bool setup( const double* );
      bool orient( double* );

      // Calculations
      void detByMinors( double*, const double*, int );
      void swapRow( double[], int, int, int );
      void swapRow( double[], int, int );

      //Output
      void print( FILE* =(stdout) ) const;
      //~ void pymolSimplex( FILE*, char[], char[], int =(0) ) const;
};

#endif
