
#include <stdio.h>
#include <stdlib.h>

#include "bVBase.h"
#include "bSimplex.h"
#include "bHex.h"

using namespace std;
using namespace bStd;

/* Static Declarations */
short bSimplex::count_ = 0;

/**/
/****** Constructors */
/* Default */
bSimplex::bSimplex() {
   for( int i=0; i < 5; ++i ) {
      this->det_[i] = 0;
   }
   for( int i=0; i < 4; ++i ) { this->nghbr_[i] = NULL; }
}

/* Copy */
bSimplex::bSimplex( const bSimplex &rhs ) {
   this->id_ = rhs.id_;
   this->det_[0] = rhs.det_[0]; this->det_[1] = rhs.det_[1];
   this->det_[2] = rhs.det_[2]; this->det_[3] = rhs.det_[3]; this->det_[4] = rhs.det_[4];
   this->nghbr_[0] = (rhs.nghbr_[0] == NULL) ? NULL : rhs.nghbr_[0];
   this->nghbr_[1] = (rhs.nghbr_[1] == NULL) ? NULL : rhs.nghbr_[1];
   this->nghbr_[2] = (rhs.nghbr_[2] == NULL) ? NULL : rhs.nghbr_[2];
   this->nghbr_[3] = (rhs.nghbr_[3] == NULL) ? NULL : rhs.nghbr_[3];
}

/* Destructor */
bSimplex::~bSimplex() {
   this->nghbr_[0] = NULL; this->nghbr_[1] = NULL; this->nghbr_[2] = NULL; this->nghbr_[3] = NULL;
}

/**/
/****** Assignment Operators */
bSimplex& bSimplex::operator=( const bSimplex &rhs ) {
   if( this == &rhs ) { return *this; }

   this->id_ = rhs.id_;
   this->det_[0] = rhs.det_[0]; this->det_[1] = rhs.det_[1];
   this->det_[2] = rhs.det_[2]; this->det_[3] = rhs.det_[3]; this->det_[4] = rhs.det_[4];
   this->nghbr_[0] = (rhs.nghbr_[0] == NULL) ? NULL : rhs.nghbr_[0];
   this->nghbr_[1] = (rhs.nghbr_[1] == NULL) ? NULL : rhs.nghbr_[1];
   this->nghbr_[2] = (rhs.nghbr_[2] == NULL) ? NULL : rhs.nghbr_[2];
   this->nghbr_[3] = (rhs.nghbr_[3] == NULL) ? NULL : rhs.nghbr_[3];

   return *this;
}

bSimplex& bSimplex::operator=( const double* rhs ) {
   this->setup( rhs );
   return *this;
}

/**/
/****** Conditional Operators */
/* Equal */
bool bSimplex::operator==( const bSimplex &rhs ) {   
   bool isEqual = true;
   for( int i=0; i < 5 && isEqual; ++i ) {
      isEqual &= (this->det_[i] == rhs.det_[i]);
   }
   return isEqual;
}

/* Not Equal */
bool bSimplex::operator!=( const bSimplex &rhs ) {
   return !( *this == rhs );
}

short bSimplex::inSphere( const double* testPt ) {
   double test[5] = { 1.0, testPt[0], testPt[1], testPt[2], 0.0 };
   double r = testPt[0]; r *= r; test[4] += r;
          r = testPt[1]; r *= r; test[4] += r;
          r = testPt[2]; r *= r; test[4] += r;

   r = 0.0;
   test[0] *= this->det_[0]; r += test[0];
   test[1] *= this->det_[1]; r += test[1];
   test[2] *= this->det_[2]; r += test[2];
   test[3] *= this->det_[3]; r += test[3];
   test[4] *= this->det_[4]; r += test[4];

   if( r < 0 ) { return 1; }
   else if( r > 0 ) { return -1; }
   else {}
   return 0;
}

/**/
/****** Setup */
/* Construct Simplex, Given Coordinates */
bool bSimplex::setup( const double* pnts ) {
   // Temporary array for determinant calculation
   double mat[25];
   
   // Initialize homology coordinate
   mat[0] = 1.0; mat[5] = 1.0; mat[10] = 1.0; mat[15] = 1.0;
   
   // Copy over array
   mat[1] = pnts[0]; mat[2] = pnts[1]; mat[3] = pnts[2];
   mat[6] = pnts[3]; mat[7] = pnts[4]; mat[8] = pnts[5];
   mat[11] = pnts[6]; mat[12] = pnts[7]; mat[13] = pnts[8];
   mat[16] = pnts[9]; mat[17] = pnts[10]; mat[18] = pnts[11];
   
   // Calculate distance coord
   double a = 0.0;
   a = pnts[0]; a *= a; mat[4]  = a; a = pnts[1];  a *= a; mat[4]  += a; a = pnts[2];  a *= a; mat[4]  += a;
   a = pnts[3]; a *= a; mat[9]  = a; a = pnts[4];  a *= a; mat[9]  += a; a = pnts[5];  a *= a; mat[9]  += a;
   a = pnts[6]; a *= a; mat[14] = a; a = pnts[7];  a *= a; mat[14] += a; a = pnts[8];  a *= a; mat[14] += a;
   a = pnts[9]; a *= a; mat[19] = a; a = pnts[10]; a *= a; mat[19] += a; a = pnts[11]; a *= a; mat[19] += a;

   // Set last coord
   mat[20] = 1.0; mat[21] = 1.0; mat[22] = 1.0; mat[23] = 1.0; mat[24] = 1.0;

   bool valid = this->orient(mat);
   return valid;
}

/****** ORIENT ******
throw: 0, same point; 1, co-planar */
bool bSimplex::orient( double* pts ) {
   double minx = pts[1];
   double miny = pts[2];
   double minz = pts[3];

   bool swap = false;
   unsigned int row = 0;
   for( unsigned int i=1; i< 4; ++i ) {
      if( i == row ) { continue; }
      swap = false;

      int xX = i; xX *= 5; ++xX;
      int yX = xX; ++yX;
      int zX = yX; ++zX;

      if( miny > pts[yX] ) { swap = true; }
      else if( miny == pts[yX] ) {
         if( minx > pts[xX] ) { swap = true; }
         else if( minx == pts[xX] ) {
            if( minz > pts[zX] ) { swap = true; }
            else if( minz == pts[zX] ) {
               printf("...same points: [min:0] %.2f, %.2f, %.2f   [vtx:%d] %.2f, %.2f, %.2f\n", minx, miny, minz,i, pts[xX], pts[yX], pts[zX]);
               exit(1);
            }
            else {}
         }
         else {}
      }
      else {}

      if( swap ) {
         minx = pts[xX];
         miny = pts[yX];
         minz = pts[zX];
         row = i;
      }
   }

   if( row != 0 ) { this->swapRow( pts, 0, row, 5 ); }
   ++row;

   this->detByMinors( this->det_, pts, 5 );

   // check for swap
   if( this->det_[4] > 0 ) {
      this->swapRow( pts, 1, 2 );
      this->detByMinors( this->det_, pts, 5 );
   }

   // Check for co-planar
   bool notCoPlanar = true;
   if( this->det_[4] == 0 ) {
      notCoPlanar = false;
   }

   return notCoPlanar;
}

void bSimplex::detByMinors( double* det, const double* mat, int dim ) {
   if( dim < 3 ) {
      det[0] = (mat[0] * mat[3]); // ad - bc
      det[1] = -(mat[1] * mat[2]);
      return;
   }

   // SETUP MINORS
   int minDim = dim - 1;
   double minMat[ minDim * minDim ];
   double minDet[ minDim ];

   // CALCULATE MINORS -- RECURSIVE!!
   int minX = 0, matX = 0;
   for( int i=0; i < dim; ++i ) {
      matX = 0; minX = 0;
      for( int k=0; k < minDim; ++k ) {
         for( int m=0; m < dim; ++m ) {
            if( m == i ) { ++matX; continue; }
            minMat[ minX ] = mat[ matX ];
            ++minX; ++matX;
         }
      }
      detByMinors( minDet, minMat, minDim );
      det[i] = 0.0;
      for( int k=0; k < minDim; ++k ) {
         det[i] += minDet[k];
      }
   }

   // ADJUST MINOR DETERMINANTS (i.e., mult by -1 if needed)
   matX = dim; matX *= minDim;
   for( int i=0; i < dim; ++i ) {
      if( !(i & 1) ) { det[i] *= -1.0; }
      det[i] *= mat[ matX ];
      ++matX;
   }

   return;
}

/* Swap Row -- array */
void bSimplex::swapRow( double mat[], int a, int b, int dim ) {
   if( a == b ) { return; }
   double temp = 0.0;
   int aI = a * dim, bI = b * dim;
   for( int i=0; i < dim; ++i ) {
      temp = mat[ aI ];
      mat[ aI ] = mat[ bI ];
      mat[ bI ] = temp;
      ++aI;
      ++bI;
   }
   return;
}

void bSimplex::swapRow( double mat[], int a, int b ) {
   if( a == b ) { return; }
   register int aX = a, bX = b, cX = 0;
   aX *= 5; bX *= 5; ++aX; ++bX; cX = aX;
   double temp[4] = { mat[cX], mat[++cX], mat[++cX], mat[++cX] };
   mat[aX] = mat[bX];
   mat[++aX] = mat[++bX];
   mat[++aX] = mat[++bX];
   mat[++aX] = mat[++bX];
   mat[bX] = temp[3];
   mat[--bX] = temp[2];
   mat[--bX] = temp[1];
   mat[--bX] = temp[0];
   return;
}

/**/
/****** Output */
/* Print */
void bSimplex::print( FILE* op ) const {
   fprintf(op, "%xu |", this->id_);
   for( int i=0; i < 5; ++i ) {
      fprintf(op, "%#2.2f ", this->det_[i]);
   } fprintf(op, "|\n");

   return;
}





