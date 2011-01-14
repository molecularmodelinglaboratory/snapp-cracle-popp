//~ #include <iostream>
//~ #include <fstream>
#include <cmath>
//~ #include <valarray>
//~ #include <string>
//~ #include <sstream>
//~ #include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bGrid.h"
#include "bPoints.h"
#include "bStamp.h"
#include "bHex.h"
using namespace std;
using namespace bStd;

const int DEBUG_COUNT=-1;

/* */
/* Constructors and Destructors */
bGrid::bGrid() : res_(0.5), fit_(4), thk_(2)  {
   if(!bHex::isHamRdy_) { bHex::setupHammingTable(); } 

   // Grid
   this->grd_ = NULL;

   // Dimensions
   this->min_[0] = 0.0; this->max_[0] = 0.0;
   this->min_[1] = 0.0; this->max_[1] = 0.0;
   this->min_[2] = 0.0; this->max_[2] = 0.0;

   // Flags
   this->haveFStmp_ = false;
   this->haveTStmp_ = false;
   this->isStamped_ = false;

   // Measurements
   this->size_ = 0;
   this->length_ = 0;
   this->height_ = 0;
   this->depth_ = 0;
}

bGrid::bGrid( const bGrid &rhs ) {
   this->copyBasicParam( rhs );

   // Copy data (if available)
   this->grd_ = new bHex[ this->size_ ];
   if( rhs.isStamped_ ) {
      for( int i=0; i < rhs.size_; ++i ) {
         this->grd_[i] = rhs.grd_[i];
      }
   }
   else {
      this->isStamped_ = false;
      for( int i=0; i < this->size_; ++i ) {
         this->grd_[i].resize( this->depth_ );
      }
   }

   // Copy fit stamp (if available), or create
   if( rhs.haveFStmp_ ) {
      (this->fStmp_) = (rhs.fStmp_);
      this->haveFStmp_ = rhs.haveFStmp_;
   }
   else {
      int radius = this->fit_; // exclusion
      this->haveFStmp_ = this->fStmp_.initializeStamp( radius, this->max_[0] );
   }

   // Copy thk stamp (if available), or create
   if( rhs.haveTStmp_ ) {
      this->tStmp_ = rhs.tStmp_;
      this->haveTStmp_ = rhs.haveTStmp_;
   }
   else {
      int radius = this->thk_ + this->fit_; // inclusion
      this->haveTStmp_ = this->tStmp_.initializeStamp( radius, this->max_[0] );
   }

   // Copy protein source
   this->prt_ = rhs.prt_; // point to same data
}

bGrid::~bGrid() {
   delete [] this->grd_;
   this->grd_ = NULL;
   this->prt_ = NULL;
}

/* deep assignment */   
bGrid& bGrid::operator=(const bGrid &rhs) {
   if(this == &rhs ) { return *this; }
   this->copyBasicParam( rhs );

   // Delete old data
   if( this->grd_ != NULL ) {
      delete [] this->grd_;
      this->grd_ = NULL;
   }

   // Copy data (if available)
   this->grd_ = new bHex[ this->size_ ];
   if( rhs.isStamped_ ) {
      for( int i=0; i < rhs.size_; ++i ) {
         this->grd_[i] = rhs.grd_[i];
      }
   }
   else {
      this->isStamped_ = false;
   }

   // Copy fit stamp (if available), or create
   if( rhs.haveFStmp_ ) {
      this->fStmp_ = rhs.fStmp_;
      this->haveFStmp_ = rhs.haveFStmp_;
   }
   else {
      int radius = this->fit_; // exclusion
      this->haveFStmp_ = this->fStmp_.initializeStamp( radius, this->max_[0] );
   }

   // Copy thk stamp (if available), or create
   if( rhs.haveTStmp_ ) {
      this->tStmp_ = rhs.tStmp_;
      this->haveTStmp_ = rhs.haveTStmp_;
   }
   else {
      int radius = this->thk_ + this->fit_; // inclusion
      this->haveTStmp_ = this->tStmp_.initializeStamp( radius, this->max_[0] );
   }

   // Copy protein source
   this->prt_ = NULL;
   this->prt_ = rhs.prt_; // point to same data

   return *this;
}

void bGrid::copyBasicParam( const bGrid &rhs ) {
   // Copy measurements
   this->size_ = rhs.size_;
   this->length_ = rhs.length_;
   this->height_ = rhs.height_;
   this->depth_ = rhs.depth_;

   // Copy parameters
   this->res_ = rhs.res_;
   this->fit_ = rhs.fit_;
   this->thk_ = rhs.thk_;
   this->offset_ = rhs.offset_;
   this->buffer_ = rhs.buffer_;

   // Copy min and max
   this->min_[0] = rhs.min_[0]; this->max_[0] = rhs.max_[0];
   this->min_[1] = rhs.min_[1]; this->max_[1] = rhs.max_[1];
   this->min_[2] = rhs.min_[2]; this->max_[2] = rhs.max_[2];

   return;
}

/* Validate Same Size Grids */
bool bGrid::sameDim( const bGrid &g ) {
   bool same = true;
   if( this->length_ != g.length_ ) { same = false; }
   else if( this->height_ != g.height_ ) { same = false; }
   else if( this->depth_ != g.depth_ ) { same = false; }
   else {}
   return same;
}

/* Remove Intersection */
bGrid& bGrid::operator-=( const bGrid &rhs ) {
   if( this == &rhs ) { return *this; }
   else if( !(this->sameDim( rhs )) ) { return *this; }
   else if( !rhs.isStamped_ ) { return *this; }
   else {}

   this->isStamped_ = true;
   for( int i=0; i < rhs.size_; ++i ) {
      this->grd_[i] ^= rhs.grd_[i];
   }

   return *this;
}

/* Add */
bGrid& bGrid::operator+=( const bGrid &rhs ) {
   if( this == &rhs ) { return *this; }
   else if( !(this->sameDim( rhs )) ) { return *this; }
   else if( !rhs.isStamped_ ) { return *this; }
   else {}

   this->isStamped_ = true;
   for( int i=0; i < rhs.size_; ++i ) {
      this->grd_[i] |= rhs.grd_[i];
   }

   return *this;
}

/* Intersection Only */
bGrid& bGrid::operator/=( const bGrid &rhs ) {
   if( this == &rhs ) { return *this; }
   else if( !(this->sameDim( rhs )) ) { return *this; }
   else if( !rhs.isStamped_ ) { return *this; }
   else {}

   this->isStamped_ = true;
   for( int i=0; i < rhs.size_; ++i ) {
      this->grd_[i] &= rhs.grd_[i];
   }

   return *this;
}

/* Overlay */
void bGrid::overlay( const bGrid &rhs ) {
   // This function essentially duplicates rhs, but without the actual grid.
   if( this == &rhs ) { return; }
   this->copyBasicParam( rhs );

   // Copy fit stamp (if available), or create
   if( rhs.haveFStmp_ ) {
      this->fStmp_ = rhs.fStmp_;
      this->haveFStmp_ = rhs.haveFStmp_;
   }
   else {
      int radius = this->fit_; // exclusion
      this->haveFStmp_ = this->fStmp_.initializeStamp( radius, this->max_[0] );
   }

   // Copy thk stamp (if available), or create
   if( rhs.haveTStmp_ ) {
      this->tStmp_ = rhs.tStmp_;
      this->haveTStmp_ = rhs.haveTStmp_;
   }
   else {
      int radius = this->thk_ + this->fit_; // inclusion
      this->haveTStmp_ = this->tStmp_.initializeStamp( radius, this->max_[0] );
   }

   // Copy protein source
   this->prt_ = rhs.prt_; // point to same data

   // Resize data
   if( this->grd_ != NULL ) {
      delete [] this->grd_;
      this->grd_ = NULL;
   }
   this->grd_ = new bHex[ this->size_ ];
   for( int i=0; i < this->size_; ++i ) { this->grd_[i].resize( this->depth_ ); }

   return;
}

/* */
/* Handle Command Line Arguements */
bool bGrid::cla(int numArg, char** argArray) {
   (void) printf( "[bGrid] Reading Parameters...\n" );
   // loop through all given parameters (i=1 to skip command)
   for(int i=1; i<numArg; ++i) {

      // check for a flag
      if(argArray[i][0] == '-') {
         switch(argArray[i][1]) {
            case 'd':
               (void) printf("hmm\n");
               break;
            case 'r':
               res_ = atof(argArray[++i]);
               (void) printf("\tres: %.2f\n",res_);
               break;
            case 'f':
               fit_ = atof(argArray[++i]);
               (void) printf("\tfit: %.2f\n",fit_);
               break;
            case 't':
               thk_ = atof(argArray[++i]);
               (void) printf("\tthk: %.2f\n",thk_);
               break;
            default:
               (void) printf("\tUnknown parameter: '%s'\n\n",argArray[i]);
               exit(1);
               break;
         }

      }
      // bad command line input...
      else {
         (void) printf("Usage: ./main.e [-i] <prot file> [<tet file>] [options]\n");
         (void) printf("\tr: resolution [1.0,0.5]\n");
         (void) printf("\tf: fit; defines exclusion radius\n");
         (void) printf("\tt: thickness; defines thickness of grid\n\n");
         exit(1);
      }
   }
   this->offset_ = this->fit_ + this->thk_;
   this->buffer_ = this->fit_;
   return 1;
}

void bGrid::setParam( float param[] ) {
   this->res_ = param[0];
   this->fit_ = param[1];
   this->thk_ = param[2];
   this->offset_ = this->fit_ + this->thk_;
   this->buffer_ = this->fit_;
   (void) printf("[bGrid] Parameters:\n");
   (void) printf("\tres: %.2f\n", this->res_ );
   (void) printf("\tfit: %.2f\n", this->fit_ );
   (void) printf("\tthk: %.2f\n", this->thk_ );
   (void) printf("\toffset: %.2f\n", this->offset_ );
   (void) printf("\tbuffer: %.2f\n", this->buffer_ );
   return;
}

/* */
/* Get point object */
void bGrid::setPointObject(bPoints& p) {
   this->prt_ = &p;
   this->prt_->setGridParam( this->res_, this->fit_, this->thk_, this->offset_, this->buffer_ );
   //~ this->prt_->fit_ = this->fit_;
   //~ this->prt_->thk_ = this->thk_;
   //~ this->prt_->res_ = this->res_;
   //~ this->prt_->offset_ = this->offset_;
   //~ this->prt_->buffer_ = this->buffer_;
   return;
}

/* Initialize Grid */
void bGrid::initializeGrid() {
   if( !(this->prt_->isInPos_) || !(this->prt_->isInRes_) ) { throw "[bGrid] Points not in docking space!"; }

   // So far, we have saved the fit, thk, and res
   // -- bPoints has done everything else
   // >> This should be the first bGrid called (after setPointObject)
   // >> ALL bGrid params must ALWAYS be in resolution
   this->fit_ /= this->res_;
   this->thk_ /= this->res_;
   this->offset_ /= this->res_;
   this->buffer_ /= this->res_;

   // save the min and max
   // -- remember: if it is in the docking space, it will already have offsets
   //              for fit and thk.
   if( !(this->prt_->haveMM_) ) { this->prt_->_findMinMax(); }
   this->min_[0] = prt_->min_[0]; this->max_[0] = prt_->max_[0];
   this->min_[1] = prt_->min_[1]; this->max_[1] = prt_->max_[1];
   this->min_[2] = prt_->min_[2]; this->max_[2] = prt_->max_[2];
   if( !(this->prt_->isInRes_) ) {
      this->min_[0] /= this->res_;   this->max_[0] /= this->res_;
      this->min_[1] /= this->res_;   this->max_[1] /= this->res_;
      this->min_[2] /= this->res_;   this->max_[2] /= this->res_;
   } 

   // adjust for grid usage
   this->min_[0] -= this->offset_; this->max_[0] += this->offset_;
   this->min_[1] -= this->offset_; this->max_[1] += this->offset_;
   this->min_[2] -= this->offset_; this->max_[2] += this->offset_;

   // Calculate dimensions
   this->length_ = (this->max_[0] + 1) + (this->buffer_);
   this->height_ = (this->max_[1] + 1) + (this->buffer_);
   this->depth_ = (this->max_[2] + 1);
   this->size_ = this->length_ * this->height_;
   (void) printf("[bGrid] dim: {%d x %d x %d, %d}\n", this->length_, this->height_, this->depth_, this->size_ );
   (void) printf("[bGrid] min: "); bPoints::printPoint( min_ );
   (void) printf("[bGrid] max: "); bPoints::printPoint( max_ );

   // Delete previous, resize, and reallocate
   if( this->grd_ != NULL ) {
      delete [] this->grd_;
      this->grd_ = NULL;
   }
   this->grd_ = new bHex[ this->size_ ];
   for( int i=0; i < this->size_; ++i ) { this->grd_[i].resize( this->depth_ ); }

   this->initializeStamps();
   return;
}

/* Initialize Exclusion and Inclusion Matrices (or grids...yes, I know) */
bool bGrid::initializeStamps() {
   int radius = this->fit_; // exclusion (fit)
   this->haveFStmp_ = this->fStmp_.initializeStamp( radius, this->length_ );

   radius = this->thk_ + this->fit_; // inclusion (fit + thk)
   this->haveTStmp_ = this->tStmp_.initializeStamp( radius, this->length_ );

   return this->haveFStmp_ & this->haveTStmp_;
}

/* */
/* Stamp Protein Points with [In|Ex]clusion Grids */
bool bGrid::stampPoints() {
   bGrid tmp;
   tmp.overlay( *this );
   tmp.fit_ = 5;// / 2;
   tmp.initializeStamps();

   // check for exclusion data
   if( this->haveFStmp_ ) {
      (void) printf("[bGrid] Stamping Protein...\n");
      this->isStamped_ = this->_stampPoints( this->fStmp_, this->prt_->pnts_, this->prt_->numPnts_ );
      if( this->prt_->haveTets_ ) {
         (void) printf("[bGrid] Stamping Centroids...\n");
         this->isStamped_ = this->_stampPoints( tmp.fStmp_, this->prt_->tets_->pnts_, this->prt_->tets_->numPnts_ );
      }
   }

   // check for inclusion data
   if( this->haveTStmp_ ) {
      (void) printf("[bGrid] Stamping Inclusions...\n");
      tmp.isStamped_ = tmp._stampPoints( this->tStmp_, tmp.prt_->pnts_, tmp.prt_->numPnts_ );
      *this -= tmp;
   }

   return isStamped_;
}

/* Stamp Grid Around a Point */
bool bGrid::_stampPoints( bStamp &stamp, float* pnts, int pntCnt ) {

   // Limit the number of points stamped
   if(DEBUG_COUNT != -1) { pntCnt = DEBUG_COUNT; }

   // adjust the point coordinates from the center to the corner
   int moveBy = stamp.radius_ - 1;
   for( int i=0; i < pntCnt * 3; ++i ) { pnts[i] -= moveBy; }

   // declare and initialize temporary values (needed here for scope)
   float oldPnts[3] = { 0, 0, 0 };
   int nDepth = 0;
   size_t nIndex = 0;
   size_t oIndex = 0;
   size_t chngIndex = 0;

   // loop through each point
   int index = 0;
   for( int i=0; i < pntCnt; ++i ) {

      float newPnts[3] = { pnts[index], pnts[++index], pnts[++index] };
      ++index;

      // calculate the new depth and index
      nIndex = (int)newPnts[1] * this->length_ + (int)newPnts[0];
      nDepth = (int)newPnts[2];

      // calculate the amount we need to change by
      chngIndex = nIndex - oIndex;

      // save the old values
      oIndex = nIndex;

      // adjust the stamp and stamp it
      for( int k=0; k < stamp.size_; ++k ) {
         stamp.slice_[k][0] += chngIndex; stamp.slice_[k][1] += chngIndex;
         stamp.slice_[k][2] += chngIndex; stamp.slice_[k][3] += chngIndex;
         stamp.slice_[k][4] += chngIndex; stamp.slice_[k][5] += chngIndex;
         stamp.slice_[k][6] += chngIndex; stamp.slice_[k][7] += chngIndex;
         this->grd_[ stamp.slice_[k][0] ].assign( stamp.stamp_[k], nDepth );
         this->grd_[ stamp.slice_[k][1] ].assign( stamp.stamp_[k], nDepth );
         this->grd_[ stamp.slice_[k][2] ].assign( stamp.stamp_[k], nDepth );
         this->grd_[ stamp.slice_[k][3] ].assign( stamp.stamp_[k], nDepth );
         this->grd_[ stamp.slice_[k][4] ].assign( stamp.stamp_[k], nDepth );
         this->grd_[ stamp.slice_[k][5] ].assign( stamp.stamp_[k], nDepth );
         this->grd_[ stamp.slice_[k][6] ].assign( stamp.stamp_[k], nDepth );
         this->grd_[ stamp.slice_[k][7] ].assign( stamp.stamp_[k], nDepth );
      }

      // prepare for the next set of coordinates
      oldPnts[0] = newPnts[0]; oldPnts[1] = newPnts[1]; oldPnts[2] = newPnts[2];
   } // LOOP: points

   // reset the coordinates to the center
   for( int i=0; i < pntCnt * 3; ++i ) { pnts[i] += moveBy; }

   // reset stamp to original position
   for( int i=0; i < stamp.size_; ++i ) {
      stamp.slice_[i][0] -= nIndex; stamp.slice_[i][1] -= nIndex;
      stamp.slice_[i][2] -= nIndex; stamp.slice_[i][3] -= nIndex;
      stamp.slice_[i][4] -= nIndex; stamp.slice_[i][5] -= nIndex;
      stamp.slice_[i][6] -= nIndex; stamp.slice_[i][7] -= nIndex;
   }

   return true;
}


/**/
void bGrid::_stampPoint( bStamp &stamp, float p[] ) {
   // adjust point
   int moveBy = stamp.radius_ - 1;
   p[0] -= moveBy; p[1] -= moveBy; p[2] -= moveBy;
   
   ulong adjust = ( (int)p[1] * this->length_ ) + (int)p[0];
   int depth = p[2];
   for( int k=0; k < stamp.size_; ++k ) {
      this->grd_[ stamp.slice_[k][0] + adjust ].assign( stamp.stamp_[k], depth );
      this->grd_[ stamp.slice_[k][1] + adjust ].assign( stamp.stamp_[k], depth );
      this->grd_[ stamp.slice_[k][2] + adjust ].assign( stamp.stamp_[k], depth );
      this->grd_[ stamp.slice_[k][3] + adjust ].assign( stamp.stamp_[k], depth );
      this->grd_[ stamp.slice_[k][4] + adjust ].assign( stamp.stamp_[k], depth );
      this->grd_[ stamp.slice_[k][5] + adjust ].assign( stamp.stamp_[k], depth );
      this->grd_[ stamp.slice_[k][6] + adjust ].assign( stamp.stamp_[k], depth );
      this->grd_[ stamp.slice_[k][7] + adjust ].assign( stamp.stamp_[k], depth );
   }

   p[0] += moveBy; p[1] += moveBy; p[2] += moveBy;

   this->isStamped_ = true;
   return;
}

/* */
/* Find Nearby Points */
uint bGrid::countNearby( uint a[], FILE* op ) {
   //~ int nearby = 0;
   bHex range( this->depth_ );
   range.assign( 0x1F, a[2] - 2 );
   bHex check( this->depth_ );
   uint sumTot = 0;
   //~ bPoints test;
   //~ test.sizeAndSpace( *prt_ );
   //~ test.isInRes_ = true;
   //~ test.isInPos_ = true;
   for( uint k=(a[1] + 2); k > (a[1] - 3); --k ) {
      int index = k; index *= this->length_; index += a[0] - 2;
      
      for( uint i=(a[0] - 2); i < (a[0] + 3); ++i ) {
         check.erase();
         check ^= range & this->grd_[ index ];
         ++index;
         sumTot += check.getNumFlip();
         //~ for( uint z=(a[2] - 2); z < (a[2] + 3); ++z ) {
            //~ float pnt[3] = { (float)i, (float)k, float(z) };
            //~ test.addPoint( pnt );
         //~ }
      }
   }
   
   //~ if( sumTot < 10 && sumTot > 0 ) {
      //~ char name[16];
      //~ sprintf( name, "test_%u%u%u", a[0], a[1], a[2]);
      //~ test.setToNormalSpace();
      //~ bPoints::pymolPoints( op, test.pnts_, test.numPnts_, name, "green", 0.05 );
   //~ }
   return sumTot;
}

/* remove internal
void bGrid::removeInternal( FILE* op ) {
   
   bGrid tmp( *this );
   uint x[2] = { min_[0] + this->thk_ + this->fit_ + 10, max_[0] - this->thk_ - this->fit_ - 10 };
   uint y[2] = { min_[1] + this->thk_ + this->fit_ + 10, max_[1] - this->thk_ - this->fit_ - 10 };
   uint z[2] = { min_[2] + this->thk_ + this->fit_ + 10, max_[2] - this->thk_ - this->fit_ - 10 };
   
   // Loop over grid
   tmp.pymol3dGrid(op, "woTest", "orange");
   this->pymol3dGrid(op, "thTest", "orange");
   for( uint k=y[0]; k < y[1]; k += 5 ) {
      for( uint i=x[0]; i < x[1]; i += 5 ) {
         for( uint m=z[0]; m < z[1]; m +=5 ) {
            uint pnt[3] = { (float)i, (float)k, (float)m };
            uint num = tmp.countNearby( pnt, op );
            if( num < 8 && num  > 0 ) {
               printf( "before: %u\t", num );
               // Erase points
               bHex range( this->depth_);
               range.assign( 0x1F, pnt[2] - 2 );
               range.flipBits();
               for( uint sk=(k + 2); sk > (k - 3); --sk ) {
                  int index = sk; index *= this->length_; index += pnt[0] - 2;
                  
                  for( uint si=(i - 2); si < (i + 3); ++si ) {
                     tmp.grd_[index] &= range;
                     //~ tmp.grd_[index] ^= range;
                     ++index;
                  }
               }
               printf( "after:%u\n", tmp.countNearby( pnt, op ) );
            }
         }
      }
   }
   //~ tmp.pymol3dGrid(op, "woInternal", "orange");
   bGrid diff( *this );
   //~ diff.pymol3dGrid(op, "diffTest", "purple");
   diff -= tmp;
   //~ diff.pymol3dGrid(op, "diff", "purple",true);
   return;
}


/* Validate Grid Point */
bool bGrid::isaGridPoint( float point[] ) const {
   int nIndex = (int)point[1];
       nIndex *= this->length_;
       nIndex += (int)point[0];
   int nDepth = (int)point[2];
   bool isa = false;

   //~ printf("stamp size: %d\t", stamp.size_);
   //~ printf("attempted depth: %d\t", nDepth);
   //~ printf("grid depth: %d\t", this->depth_);
   //~ printf("hex depth: %d [%d]\n", this->grd_[nIndex].numBits_, this->grd_[nIndex].size_);
   if(nIndex < 0 || nDepth < 0) { isa = false; }
   else if( nIndex > this->size_ ) { isa = false; }
   else if( nDepth > this->depth_ ) { isa = false; }
   else if( grd_[nIndex] & nDepth ) { isa = true; }
   else {};

   return isa;
}

/****** PRINT */
/* Write PyMol Grid */
void bGrid::pymol3dGrid( FILE* op, const char name[], const char color[], bool full ) {
   float size = 0.1;
   
   bHex** toPymol = &(this->grd_);
   bGrid tmp;
   if( !full ) {
      bool inDS = this->prt_->isInDockingSpace();
      if( !inDS ) { this->prt_->setToDockingSpace(); }
      
      tmp.overlay( *this );
      --tmp.fit_;
      tmp.thk_ = 2;
      tmp.initializeStamps();
      tmp.isStamped_ = tmp._stampPoints( tmp.tStmp_, tmp.prt_->pnts_, tmp.prt_->numPnts_ );
      tmp /= *this;
      toPymol = &(tmp.grd_);
      
      if( !inDS ) { this->prt_->setToNormalSpace(); }
   }
   
   this->pymol3dGrid( op, *toPymol, this->length_, this->height_, name, color, size );
   return;
}
void bGrid::pymol3dGrid( FILE* op, bHex* grd, int len, int hei, const char name[], const char color[], float size ) {

   (void) fprintf(op,"%s = [\n",name);

   int index = 0;
   for(int k=0; k< hei; ++k) {
      for(int i=0; i< len; ++i) {
         int   num = grd[index].getNumFlip();
         if( num == 0 ) { ++index; continue; }
         int* list = grd[index].getActive();
         ++index;

         for( int n=0; n < num; ++n ) {
            float a = (i * this->res_) + this->prt_->planeDisplacement_[0];
            float b = (k * this->res_) + this->prt_->planeDisplacement_[1];
            float c = (list[n] * this->res_) + this->prt_->planeDisplacement_[2];
            (void) fprintf(op,"\tSPHERE, %.2f, %.2f, %.2f, %.2f,\n",a,b,c,size);
         }
         list = NULL;
      }
   }
   (void) fprintf(op,"\t]\n");
   (void) fprintf(op,"cmd.load_cgo(%s,'%s')\n\n",name,name);
   (void) fprintf(op,"cmd.color(\"%s\", \"%s\")\n", color, name );

   //~ (void) fprintf(op,"cmd.center('%s')\n",name);
   //~ (void) fprintf(op,"cmd.hide(%s)\n",name);

   return;
}

void bGrid::printProtein() {
   // Easy loop through grid -- SAVE
   int index = 0;
   for(int k=0; k < this->height_; ++k) {
      for(int i=0; i < this->length_; ++i) {
         int   num = this->grd_[index].getNumFlip();
         if( num == 0 ) { continue; }
         int* list = this->grd_[index].getActive();
         ++index;

         int cnt = 0;
         for( int n=0; n < num; ++n ) {
            while( cnt < list[n] ) { ++cnt; (void) printf("0 "); }
            (void) printf("1 ");
         }
         while( cnt < this->depth_ ) { ++cnt; (void) printf("0 "); }
         (void) printf("\n");
         list = NULL;
      }
      (void) printf("\n");
   }
   return;
}

void bGrid::printFStmp() {
   this->fStmp_.print();
   return;
}
void bGrid::printTStmp() {
   this->tStmp_.print();
   return;
}


void bGrid::print2dProtein() {
   this->print2dGrid( this, 1 );
   return;
}
void bGrid::print2dGrid( bGrid *g, int showPrt) {

   int prtSize = prt_->capPnts_;
   if(DEBUG_COUNT != -1) { prtSize = DEBUG_COUNT; }
   else if(!showPrt) { prtSize = 0; }
   else {}

   int index = 0;
   for(int k=0; k < g->height_; ++k) {
      (void) printf("%3d:",k);
      for(int i=0; i < g->length_; ++i) {
         bool pointhere = false;
         int num = g->grd_[index].getNumFlip();
         ++index;

         if( showPrt ) {
            for( int m=0; m < prtSize; ++m ) {
               int b = m * 3;
               int xy[2] = {g->prt_->pnts_[b], g->prt_->pnts_[++b] };
               if( i == xy[0] && k == xy[1] ) {
                  pointhere = true;
                  m = prtSize;
               }
            }
         }

         if( pointhere ) { (void) printf("?"); }
         else if( num == 0 ) { (void) printf(" "); continue; }
         else { (void) printf("1"); }
      }
      (void) printf(":%d\n", g->length_);
   }
   (void) printf("::%d\t", g->height_);

   (void) printf("Correction");
   for( int i=0; i < 3; ++i ) {
      (void) printf(": %f", g->prt_->planeDisplacement_[i]);
   }
   (void) printf("\n\n");
   return;
}

/* */

void bGrid::clear() {
   for( int i=0; i < this->size_; ++i ) { this->grd_[i].erase(); }
   return;

   delete [] this->grd_;
   this->grd_ = NULL;
   this->grd_ = new bHex[ this->size_ ];
   for( int i=0; i < this->size_; ++i ) { this->grd_[i].resize( this->depth_ ); }
   return;
}

