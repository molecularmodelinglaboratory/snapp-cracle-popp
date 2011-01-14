
#include <stdio.h>
#include <math.h>
#include "bStamp.h"
using namespace std;
using namespace bStd;

/* */
/* Constructors and Destructors */
bStamp::bStamp() {
   this->stamp_ = NULL;
   this->slice_ = NULL;
   this->radius_ = 0;
   this->size_ = 0;
}

bStamp::bStamp( const bStamp& s ) {
   this->size_ = s.size_;
   this->radius_ = s.radius_;
   this->stamp_ = new ulong[ this->size_ ];
   this->slice_ = new ulong*[ this->size_ ];
   for( int i=0; i < this->size_; ++i ) {
      this->stamp_[i] = s.stamp_[i];
      this->slice_[i] = new ulong[8];
      for( int k=0; k < 8; ++k ) { this->slice_[i][k] = s.slice_[i][k]; }
   }
}

bStamp::~bStamp() {
   delete [] this->stamp_;
   this->stamp_ = NULL;
   for( int i=0; i < this->size_; ++i ) {
      delete [] this->slice_[i];
      this->slice_[i] = NULL;
   }
   delete [] this->slice_;
   this->slice_ = NULL;
}

bStamp& bStamp::operator=( const bStamp& rhs ) {
   if( this->stamp_ != NULL ) {
      delete [] this->stamp_;
      this->stamp_ = NULL;
      for( int i=0; i < this->size_; ++i ) {
         delete [] this->slice_[i];
         this->slice_[i] = NULL;
      }
      delete [] this->slice_;
      this->slice_ = NULL;
   }

   this->size_ = rhs.size_;
   this->radius_ = rhs.radius_;
   this->stamp_ = new ulong[ this->size_ ];
   this->slice_ = new ulong*[ this->size_ ];
   for( int i=0; i < this->size_; ++i ) {
      this->stamp_[i] = rhs.stamp_[i];
      this->slice_[i] = new ulong[8];
      for( int k=0; k < 8; ++k ) {
         this->slice_[i][k] = rhs.slice_[i][k];
      }
   }
   return *this;
}

/* Initialize stamp_ */
bool bStamp::initializeStamp( const int r, const int m) {
   if( this->stamp_ != NULL ) { delete [] this->stamp_; this->stamp_ = NULL; }
   if( this->slice_ != NULL ) {
      for( int i=0; i < this->size_; ++i ) {
         delete [] this->slice_[i];
         this->slice_[i] = NULL;
      }
      delete [] this->slice_;
      this->slice_ = NULL;
   }

   // radius_ from calling function: int radius_ = (fit_ / res_);
   this->radius_ = r;
   int diamet = 2 * this->radius_;

   // initialize bins
   this->size_ = 0;
   for(int i = 0; i <= this->radius_; ++i) { this->size_ += i; }

   // resize! (and fill in with ones)
   this->stamp_ = new ulong[ this->size_ ];
   this->slice_ = new ulong*[ this->size_ ];

   // loop through each outer layer of the grid
   int index = 0;
   int n = m;
   int p = 2*radius_ - 1;

   // setup the middle point
   float temp = ( (float)diamet - 1 ) / 2;
   float mid[3] = { temp, temp, temp };
   float scaledLength = (float)this->radius_;// - 0.5;

   for(int k = 0; k < this->radius_; ++k) {

      // loop through each unique position on the outer layer
      for(int i = k; i < this->radius_; ++i) {

         // initialize the valarray
         this->slice_[index] = new ulong[8];

         // calculate each of the eight points
         // note: these are mapped to the grid, not to the small box of the
         // exclusion. MEANING: we only need to add the x dimension (usually i
         // in our implementation) to the indirect arrays to get the
         // appropriate index on the grid.
         slice_[index][0] = i + (k*n);
         slice_[index][1] = (p-i) + (k*n);
         slice_[index][2] = (i*n) + k;
         slice_[index][3] = (i*n) + (p-k);
         slice_[index][4] = (p-i)*n + k;
         slice_[index][5] = (p-i)*n + (p-k);
         slice_[index][6] = (p-k)*n + i;
         slice_[index][7] = (p-k)*n + (p-i);
         
         // Setup the stamp
         stamp_[index] = 0x0;
         float gpt[3] = { (float)k, (float)i, 0.0 };
         for( int z = 0; z < radius_; ++z ) {
            gpt[2] = (float)z;

            // check the distance between mid pt and tmp pt
            // -- if w/in range, add symmetrically
            if( pointDistance( mid, gpt ) < scaledLength ) {
               stamp_[index] |= (ulong)pow( 2.0, z );
               stamp_[index] |= (ulong)pow( 2.0, (diamet - 1 - z) );
            }
         }
         ++index;
      }
   }

   return 1;
}

/* Print the Stamp */
void bStamp::print() {
   int index = 0;
   for( int i=0; i < this->radius_; ++i ) {
      for( int j=0; j < i; ++j ) { printf("  "); }
      for( int j=i; j < this->radius_; ++j ) {
         if( this->stamp_[index] > 0) { printf("1 "); }
         else { printf("0 "); }
         ++index;
      }
      printf("\n");
   }
   printf( "%d : %d\n", this->size_, index );
   return;
}

/* Calculate Point Distances */
float bStamp::pointDistance( const float* a, const float* b ) {
   float c[3] = { a[0], a[1], a[2] };
   c[0] -= b[0]; c[1] -= b[1]; c[2] -= b[2];
   c[0] *= c[0]; c[1] *= c[1]; c[2] *= c[2];
   float sum = c[0]; sum += c[1]; sum += c[2];
   return sqrt( sum );
}





