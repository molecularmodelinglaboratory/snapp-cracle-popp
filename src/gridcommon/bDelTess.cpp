#include <deque>

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include "bAA.h"
#include "bCalc.h"
#include "bHex.h"
#include "bPoints.h"
#include "bSimplex.h"
#include "bDelTess.h"
#include "bSort.h"
#include "MersenneTwister.h"
#include "bSnapp.h"

//~ #include "bDocker.h"

using namespace std;
using namespace bStd;

int bDelTess::numInst_ = 0;

/****** Constructor */

/* Default */
bDelTess::bDelTess() :
   printEdg_(true), printVal_(false), printInv_(false), printTrm_(false), printRem_(false), printOni_(false), printInf_(false)
{
   this->erase();
}

/* Copy */
bDelTess::bDelTess( const bDelTess &rhs, const bool saveSrc ) :
   printEdg_(true), printVal_(false), printInv_(false), printTrm_(false), printRem_(false), printOni_(false), printInf_(false)
{
   if( rhs.src_ == NULL ) { this->erase(); }
   else {
      this->printEdg_ = rhs.printEdg_;

      // Copy self
      this->vrtx_ = rhs.vrtx_;
      this->simplex_ = rhs.simplex_;
      this->simplID_ = rhs.simplID_;
      this->edges_ = rhs.edges_;
      this->max_ = rhs.max_;
      this->min_ = rhs.min_;
      this->size_ = rhs.size_;
      this->score_ = rhs.score_;

      // Copy sources
      if( saveSrc && rhs.src_ != NULL ) {
         this->numChains_ = rhs.numChains_;
         this->src_ = new bPoints* [this->numChains_];
         this->src_[0] = new bPoints( *(rhs.src_[0]) );
         for( int i=1; i < this->numChains_; ++i ) { this->src_[i] = rhs.src_[i]; }
         this->toAdd_ = rhs.toAdd_;
      }
      else {
         this->src_ = NULL;
         this->numChains_ = 0;
         this->toAdd_.clear();
      }
      
      // Copy information;
      this->newID_ = rhs.newID_;
      this->newList_ = rhs.newList_;
      this->hullID_ = rhs.hullID_;
      this->open_ = rhs.open_;
      this->currPt_ = new double[3];
      this->currPt_[0] = rhs.currPt_[0]; this->currPt_[1] = rhs.currPt_[1]; this->currPt_[2] = rhs.currPt_[2];

      // Copy flags
      this->isReady_ = rhs.isReady_;
      this->isVerified_ = rhs.isVerified_;
      this->isSlim_ = rhs.isSlim_;
      this->isTrim_ = rhs.isTrim_;
      this->haveEdge_ = rhs.haveEdge_;
      this->haveType_ = rhs.haveType_;
      this->haveScore_ = rhs.haveScore_;

      // debugging
      //~ if( doDebug_ ) debug_ = fopen("debug.txt", "w");

      ++numInst_;
   }
}

/* Destructor */
bDelTess::~bDelTess() {

   // self
   this->vrtx_.clear();
   this->simplex_.clear();
   this->simplID_.clear();
   this->simplTy_.clear();
   this->edges_.clear();

   // sources
   if( this->src_ != NULL ) {
      if( this->src_[0] != NULL ) { delete this->src_[0]; }
      for( int i=0; i < this->numChains_; ++i ) { this->src_[i] = NULL; }
      delete [] this->src_;
      this->src_ = NULL;
   }
   this->toAdd_.clear();

   // Information (constantly changing)
   this->newID_.clear();
   this->newList_.clear();
   this->hullID_.clear();
   this->open_.clear();
   delete [] this->currPt_;
   this->currPt_ = NULL;

   // debugging
   //~ this->invalidS_.clear();
   //~ this->invalidP_.clear();
   //~ this->removed_.clear();
   //~ this->removedWhen_.clear();
   //~ if( doDebug_ ) fclose( debug_ );

   --numInst_;
}

bDelTess& bDelTess::operator=( const bDelTess& rhs ) {
   // Copy self
   this->vrtx_ = rhs.vrtx_;
   this->simplex_ = rhs.simplex_;
   this->simplID_ = rhs.simplID_;
   this->simplTy_ = rhs.simplTy_;
   this->edges_ = rhs.edges_;
   this->max_ = rhs.max_;
   this->min_ = rhs.min_;
   this->size_ = rhs.size_;
   this->score_ = rhs.score_;

   // Remove sources
   if( this->src_ != NULL ) {
      if( this->src_[0] != NULL ) { delete this->src_[0]; }
      for( int i=0; i < this->numChains_; ++i ) { this->src_[i] = NULL; }
      delete [] this->src_;
      this->src_ = NULL;
   }

   // Copy sources
   this->numChains_ = rhs.numChains_;
   this->src_ = new bPoints* [this->numChains_];
   this->src_[0] = new bPoints( *(rhs.src_[0]) );
   for( int i=1; i < this->numChains_; ++i ) { this->src_[i] = rhs.src_[i]; }
   this->toAdd_ = rhs.toAdd_;

   // Copy information;
   this->newID_ = rhs.newID_;
   this->newList_ = rhs.newList_;
   this->hullID_ = rhs.hullID_;
   this->open_ = rhs.open_;
   if( this->currPt_ != NULL ) { delete [] this->currPt_; this->currPt_ = NULL; }
   this->currPt_ = new double[3];
   this->currPt_[0] = rhs.currPt_[0]; this->currPt_[1] = rhs.currPt_[1]; this->currPt_[2] = rhs.currPt_[2];

   // Copy flags
   this->isReady_ = rhs.isReady_;
   this->isVerified_ = rhs.isVerified_;
   this->isSlim_ = rhs.isSlim_;
   this->isTrim_ = rhs.isTrim_;
   this->haveEdge_ = rhs.haveEdge_;
   this->haveType_ = rhs.haveType_;
   this->haveScore_ = rhs.haveScore_;

   // debugging
   //~ this->invalidS_ = rhs.invalidS_;
   //~ this->invalidP_ = rhs.invalidP_;
   //~ this->removed_ = rhs.removed_;
   //~ this->removedWhen_ = rhs.removedWhen_;
   return *this;
}



void bDelTess::erase() {
   // self
   this->vrtx_.clear();
   this->simplex_.clear();
   this->simplID_.clear();
   this->simplTy_.clear();
   this->edges_.clear();
   this->max_ = 0;
   this->min_ = 0;
   this->size_ = 0;
   this->score_ = 0.0;


   // sources
   this->numChains_ = 0;
   this->src_ = NULL;
   this->toAdd_.clear();

   // Information (constantly changing)
   this->newID_.clear();
   this->newList_.clear();
   this->hullID_.clear();
   this->open_.clear();
   this->currPt_ = new double[3];
   this->currPt_[0] = 0.0; this->currPt_[1]  = 0.0; this->currPt_[2]  = 0.0;

   // flags
   this->isReady_ = false;
   this->isVerified_ = false;
   this->isSlim_ = false;
   this->isTrim_ = false;
   this->haveEdge_ = false;
   this->haveType_ = false;
   this->haveScore_ = false;

   ++numInst_;
}

void bDelTess::reset() {
   this->clear();
   this->vrtx_.clear();
   this->toAdd_.clear();
   this->isReady_ = false;
   this->prep();
   return;
}

void bDelTess::clear() {
   // self
   while( this->vrtx_.size() > 4 ) { vrtx_.pop_back(); }
   this->simplex_.clear();
   this->simplID_.clear();
   this->simplTy_.clear();
   this->simTrm_.clear();
   this->edges_.clear();
   this->score_ = 0.0;

   // flags
   this->haveType_ = false;
   this->haveEdge_ = false;
   this->isTrim_ = false;
   this->isSlim_ = false;
   this->isVerified_ = false;
   this->haveScore_ = false;

   // Information (constantly changing)
   if( this->currPt_ != NULL ) { delete [] this->currPt_; this->currPt_ = NULL; }
   this->currPt_ = new double[3];
   this->newID_.clear();
   this->newList_.clear();
   this->hullID_.clear();
   this->open_.clear();
}



/**/
/****** Internal */
void bDelTess::addSrc( bPoints &ptSrc ) {
   // Ensure ptSrc is ready
   if( ptSrc.numPnts_ == 0 ) { printf("[bDelTess] Empty point set.\n"); }
   else if( !(ptSrc.haveMM_) ) { ptSrc._findMinMax(); }
   else {}

   bPoints** temp = NULL;
   if( this->src_ != NULL ) {
      if( this->src_[0] != NULL ) {
         delete this->src_[0];
         this->src_[0] = NULL;
      }

      if( this->numChains_ > 1 ) {
         temp = new bPoints* [this->numChains_];
         for( int i=1; i < this->numChains_; ++i ) {
            temp[i] = this->src_[i];
            this->src_[i] = NULL;
         }
      }
      delete [] this->src_;
      this->src_ = NULL;
   }
   else { ++this->numChains_; }
   ++this->numChains_;

   // New source
   this->src_ = new bPoints* [this->numChains_];
   this->src_[0] = NULL;
   if( temp != NULL ) {
      for( int i=1; i < (this->numChains_ - 1); ++i ) { 
         this->src_[i] = temp[i];
         temp[i] = NULL;
      }
      delete [] temp;
      temp = NULL;
   }
   this->src_[ (this->numChains_ - 1) ] = &ptSrc;

   return;
}

void bDelTess::addSrc( bPoints** ptSrc, int num ) {
   // Clear Source
   if( this->src_ != NULL ) {
      if( this->src_[0] != NULL ) {
         delete this->src_[0];
         this->src_[0] = NULL;
      }
      for( int i=1; i < this->numChains_; ++i ) {
         if( this->src_[i] != NULL ) { this->src_[i] = NULL; }
      }
      delete [] this->src_;
      this->src_ = NULL;
   }

   // Reset srouces
   this->numChains_ = num;
   ++(this->numChains_);
   this->src_ = new bPoints*[this->numChains_];
   this->src_[0] = NULL;
   for( int i=1; i < this->numChains_; ++i ) { this->src_[i] = ptSrc[i - 1]; }

   return;
}

bool bDelTess::relink( bPoints* o, bPoints* n ) {
   bool found = false;
   for( int i=1; i < this->numChains_; ++i ) {
      if( (this->src_[i]) == o ) {
         this->src_[i] = NULL;
         this->src_[i] = n;
         i = this->numChains_;
         found = true;
      }
   }
   return found;
}

bool bDelTess::prep() {
   bool prep = false;
   if( this->numChains_ <= 1 ) { return prep; }
   
   // initial max and min
   if( !this->src_[1]->haveMM_ ) { this->src_[1]->_findMinMax(); }
   int max[3] = { this->src_[1]->max_[0], this->src_[1]->max_[1], this->src_[1]->max_[2] };
   int min[3] = { this->src_[1]->min_[0], this->src_[1]->min_[1], this->src_[1]->min_[2] };
   
   // find min and max extremes
   for( int i=2; i < this->numChains_; ++i ) {
      if( !this->src_[i]->haveMM_ ) { this->src_[1]->_findMinMax(); }
      if( max[0] < this->src_[i]->max_[0] ) { max[0] = this->src_[i]->max_[0]; }
      if( max[1] < this->src_[i]->max_[1] ) { max[1] = this->src_[i]->max_[1]; }
      if( max[2] < this->src_[i]->max_[2] ) { max[2] = this->src_[i]->max_[2]; }
      if( min[0] > this->src_[i]->min_[0] ) { min[0] = this->src_[i]->min_[0]; }
      if( min[1] > this->src_[i]->min_[1] ) { min[1] = this->src_[i]->min_[1]; }
      if( min[2] > this->src_[i]->min_[2] ) { min[2] = this->src_[i]->min_[2]; }
   }

   // find absolute max and min
   this->max_ = max[0]; this->min_ = min[0];
   if( this->max_ < max[1] ) { this->max_ = max[1]; }
   if( this->max_ < max[2] ) { this->max_ = max[2]; }
   if( this->min_ > min[1] ) { this->min_ = min[1]; }
   if( this->min_ > min[2] ) { this->min_ = min[2]; }

   // calculate infinity points
   int len = this->max_ - this->min_;
   int mid = this->min_ + len / 2;
   float infPt[12] = { mid - len * 4, mid - len * 1, mid - len * 3,
                       mid + len * 1, mid - len * 2, mid + len * 4,
                       mid - len * 1, mid + len * 4, mid - len * 1,
                       mid + len * 4, mid - len * 1, mid + len * 0
   };

   // Save infinity points
   if( this->src_[0] != NULL ) {
      delete this->src_[0];
      this->src_[0] = NULL;
   }
   this->src_[0] = new bPoints(4);
   this->src_[0]->addPoints( infPt, 4 );

   // Add first points
   ptData pt( 0, 0 );
   this->vrtx_.clear();
   this->vrtx_.push_front( pt ); pt.pos_ = 1;
   this->vrtx_.push_front( pt ); pt.pos_ = 2;
   this->vrtx_.push_front( pt ); pt.pos_ = 3;
   this->vrtx_.push_front( pt );

   // Randomize
   this->toAdd_.clear();
   this->randomizePts();

   return true;
}

/* Uses the Mersenne Twister algorithm to get a random number
   -- uses the nextInt method from Java to constrain w/in [0-range]
*/

/* Randomize Point Order */
void bDelTess::randomizePts() {
   rndNum_.seed((unsigned int)time(NULL)); // seeds the Mersenne Twister RNG
   this->toAdd_.clear();

   // create unaltered point list
   deque<ptData> ptList;
   this->size_ = 0;
   ptData p( 0, 0 );
   for( int i=1; i < this->numChains_; ++i ) {
      this->size_ += this->src_[i]->numPnts_;
      p.chn_ = i;
      for( int k=0; k < this->src_[i]->numPnts_; ++k ) { 
         p.pos_ = k;
         ptList.push_back( p );
      }
   }

   // generate random numbers and add index 
   int rndm = 0;
   int last = this->size_ - 1;
   for( int i=0; i < (this->size_ - 1); ++i ) {
      rndm = bCalc::getRandomNumber( last );
      while( rndm >= (int)ptList.size() ) { rndm = bCalc::getRandomNumber( last ); }
      toAdd_.push_back( ptList[rndm] );
      ptList.erase( ptList.begin() + rndm );
      --last;
   }

   toAdd_.push_back( ptList.front() );
   return;
}

/**/

/***** TESSELLATE METHODS */
bool bDelTess::tessellate_full( bPoints &ptSrc, const char* path, const char* base ) {
   this->addSrc( ptSrc );
   this->prep();

   bool isValid = false;
   uint cnt = 0;
   do{
      isValid = this->tessellate();
      isValid &= this->verify();
      if( !isValid ) { this->reset(); }
   } while( !isValid && ++cnt < 10 );

   if( isValid ) {
      this->removeExcess();
      this->trim( 12.0 );
      this->findEdges();

      // Create Filename
      char basename[96];
      memset( basename, '\0', sizeof(char) * 96 );
      if( path == NULL ) { strcpy( basename, this->src_[1]->pntPath_ ); }
      else { strcpy( basename, path ); }
      if( base == NULL ) { 
         strcat( basename, this->src_[1]->pntBase_ );
         for( int i=2; i < this->numChains_; ++i ) {
            if( this->src_[i]->pntBase_ == NULL ) { continue; }
            else if( this->src_[i]->pntBase_[0] == '\0' ) { continue; }
            else { sprintf( basename, "%s_%s", basename, this->src_[i]->pntBase_ ); }
         }
      }
      else { strcat( basename, base ); }

      // Output files
      char fullname[96];
      strcpy( fullname, basename );
      strcat( fullname, "_gdDT.del4" );
      FILE* op = fopen( fullname, "w" );
      this->printTetPos( op );
      fclose( op );

      strcpy( fullname, basename );
      strcat( fullname, "_gdDT.tet" );
      op = fopen( fullname, "w" );
      this->printTetRes( op );
      fclose( op );
   }
   return isValid;
}


bool bDelTess::tessellate_full( bPoints **ptSrc, int num, const char* base, const char* path, float t, bool doprint ) {
   this->addSrc( ptSrc, num );
   this->prep();
   bool isValid = false;
   uint cnt = 0;
   do{
      isValid = this->tessellate();
      isValid &= this->verify();
      if( !isValid ) { this->reset(); }
   } while( !isValid && ++cnt < 10 );

   if( isValid ) {
      this->trim( t );
      this->removeExcess();
      this->findEdges();
      this->score();
   }

   // Create Filename
   if( path == NULL ) {
      path = this->src_[1]->pntPath_;
      base = this->src_[1]->pntBase_;
   }
   int psize = strlen( path );
   int bsize = strlen( base );
   int bnsize = psize + bsize + 1;
   char basename[ bnsize ];
   memset( basename, '\0', bnsize );
   memmove( basename, path, psize );
   memmove( basename + psize, base, bsize );

   // Output files
   int fsize = bnsize + 6;
   char fullname[ fsize ];
   memset( fullname, '\0', fsize );
   memmove( fullname, basename, --bnsize );

   // Write Tet by Res (e.g., 'AAAA 0')
   memmove( fullname + bnsize, ".tet", 4 );
   FILE* op = fopen( fullname, "w" );
   this->printTetRes( op );
   fclose( op );
   if( doprint ) {
      op = fopen( fullname, "r" );
      if( op ) { printf("wrote: %s\n", fullname); }
      fclose(op);
   }

   // Write Tet by Pos (e.g., '1 2 3 4')
   memmove( fullname + bnsize, ".del4", 5 );
   op = fopen( fullname, "w" );
   this->printTetPos( op );
   fclose( op );
   if( doprint ) {
      op = fopen( fullname, "r" );
      if( op ) { printf("wrote: %s\n", fullname); }
      fclose(op);
   }

   // Write Point Coord (e.g., '1.0 2.0 3.0 1:0 G')
   memmove( fullname + bnsize, ".pnt", 5 );
   op = fopen( fullname, "w" );
   this->printTetPnt( op );
   fclose( op );
   if( doprint ) {
      op = fopen( fullname, "r" );
      if( op ) { printf("wrote: %s\n", fullname); }
      fclose(op);
   }

   // Write Tess Catalog (e.g., everything else)
   memmove( fullname + bnsize, ".cat", 5 );
   op = fopen( fullname, "w" );
   this->printTetCat( op );
   fclose( op );
   if( doprint ) {
      op = fopen( fullname, "r" );
      if( op ) { printf("wrote: %s\n", fullname); }
      fclose(op);
   }

   return isValid;
}


bool bDelTess::tessellate_std( bPoints** ptSrc, int num, int focus ) {
   this->addSrc( ptSrc, num );
   this->prep();

   bool isValid = false;
   uint cnt = 0;
   do{
      isValid = this->tessellate();
      isValid &= this->verify();
      if( ++cnt > 10 ) { break; }
      if( !isValid ) { this->reset(); }
   } while( !isValid );

   if( isValid ) {
      if( focus != 0 ) { this->focus( focus ); }
      this->trim( 12.0 );
      this->removeExcess();
      this->findEdges();
   }

   return isValid;
}


bool bDelTess::tessellate_l1o( bPoints** ptSrc, int num, int chn, int pos ) {
   // Save the l1o point
   ptData l1o( pos, chn );

   // Setup the tessellation
   this->addSrc( ptSrc, num );
   this->prep();

   // Remove the point from the randomization
   for( uint i=0; i < this->toAdd_.size(); ++i ) {
      if( this->toAdd_[i].chn_ == chn ) {
         if( this->toAdd_[i].pos_ == pos ) {
            this->toAdd_.erase( this->toAdd_.begin() + i );
         }
      }
   }

   // Perform the tessellation
   bool isValid = this->tessellate();
   if( !isValid ) {
      this->reset();
      isValid = this->tessellate();
   }

   // Add the point back in
   this->toAdd_.push_back( l1o );

   return isValid;
}

bool bDelTess::tessellate_l1o_fin( bPoints* ptSrc, int num, int focus ) {

   // Set the new source (diff addr, but same coord save one)
   this->src_[num] = ptSrc;

   // Move the point over
   this->vrtx_.push_back( this->toAdd_.front() );
   this->toAdd_.pop_front();

   // Finish the tessellation
   bool isValid;
   this->resetHandlers();
   isValid = this->checkPoint( this->vrtx_.size() - 1 );
   this->resetHandlers();
   this->newList_.clear();

   // Verify and polish
   isValid = this->verify();
   if( isValid ) {
      if( focus != 0 ) { this->focus( focus ); }
      this->trim( 12.0 );
      this->removeExcess();
   }

   return isValid;
}



/****** TESSELLATE ******/
bool bDelTess::tessellate() {
   //~ printf("[bDelTess] Tessellating protein...\n");

   // Reset simplex counter
   bSimplex::count_ = 0;

   // Create first simplex
   bSimplex tet;
   int pts[4] = { 0, 1, 2, 3 };
   this->addSimplex( pts );

   // Prep
   if( this->currPt_ != NULL ) { delete [] this->currPt_; this->currPt_ = NULL; }
   this->currPt_ = new double[3];
   this->resetHandlers();

   // Add each point
   bool isValid = true;
   while( !this->toAdd_.empty() ) {

      // Add to list, remove from queue
      this->vrtx_.push_back( this->toAdd_.front() );
      this->toAdd_.pop_front();

      // Reset Handlers
      isValid = this->checkPoint( this->vrtx_.size() - 1 );
      this->resetHandlers();
      this->newList_.clear();

   } // end while loop

   return isValid;
}

/****** Simplex Handling */
/* Remove */
/* Check All Simplexes Against Single Point */
bool bDelTess::checkPoint( int ptId ) {
   if( this->currPt_ == NULL ) { this->currPt_ = new double[3]; }
   
   bool valid = true;
   int chn = this->vrtx_[ ptId ].chn_;
   int pos = this->vrtx_[ ptId ].pos_;
   
   this->currPtId_ = ptId;
   this->currPt_[0] = this->src_[chn]->pnts_[ ( 0 + (pos * 3) ) ];
   this->currPt_[1] = this->src_[chn]->pnts_[ ( 1 + (pos * 3) ) ];
   this->currPt_[2] = this->src_[chn]->pnts_[ ( 2 + (pos * 3) ) ];
   
   for( uint i=0; i < this->simplex_.size(); ++i ) {
      if( this->simplID_[i].empty() ) { continue; } // deleted simplex (open)
      if( this->simplID_[i] & this->currPtId_ ) { continue; } // simplex uses point

      int isIn = this->simplex_[i].inSphere( this->currPt_ );
      if( isIn == 1 ) {} // outside
      else if( isIn == 0 || isIn == -1 ) {
         this->outbreak( NULL, &(this->simplex_[i]) );
         valid = this->treat();
         break; // should have found all sick simplexes
      }
      else {}
   }
   return valid;
}

/* Check All Points Against Single Simplex */
void bDelTess::checkSimplex( bSimplex &sx ) {
   if( this->currPt_ == NULL ) { this->currPt_ = new double[3]; }
   for( uint i=4; i < this->vrtx_.size(); ++i ) {
      if( this->simplID_[ sx.id_ ] & (int)i ) { continue; } // simplex uses point

      int chn = this->vrtx_[ i ].chn_;
      int pos = this->vrtx_[ i ].pos_;
      this->currPtId_ = i;
      this->currPt_[0] = this->src_[chn]->pnts_[ ( 0 + (pos * 3) ) ];
      this->currPt_[1] = this->src_[chn]->pnts_[ ( 1 + (pos * 3) ) ];
      this->currPt_[2] = this->src_[chn]->pnts_[ ( 2 + (pos * 3) ) ];
      int isIn = sx.inSphere( this->currPt_ );
      if( isIn == 1 ) {} // outside
      else if( isIn == -1 || isIn == 0 ) {
         this->outbreak( NULL, &(sx) );
         this->treat();
         break; // if infected, it'll already be deleted
      }
      else {}
   }
   return;
}

/* Outbreak */
void bDelTess::outbreak( bSimplex *prev, bSimplex *curr ) {
   if( this->sick_ & curr->id_ ) { return; }

   this->sick_ |= curr->id_;
   int numNull = 0;
   for( int i=0; i < 4; ++i ) { // check neighbors

      // Check simplex validity
      if( curr->nghbr_[i] == NULL ) {
         if( numNull > 0 ) {
            bHex test( this->simplID_[ curr->id_ ] );
            test.filter( (ulong)0xF, 0 );
            if( test.getNumFlip() < 4 ) {
               printf("[outbreak] Too many null neighbors.\n");
               exit(1);
            }
         }
         else {
            this->infectedHull( curr );
            ++numNull;
         }
         continue;
      } // no neighbor
      else if( curr->nghbr_[i] == prev ) {
         continue;
      } // previous simplex
      else {}

      // Check simplex point containment
      int isIn = curr->nghbr_[i]->inSphere( this->currPt_ );
      if( this->simplID_[ curr->nghbr_[i]->id_ ] & this->currPtId_ ) {
         isIn = 0;
      }

      // Check
      if( isIn == 1 ) { this->noInfection( curr, curr->nghbr_[i] ); } // outside
      else if( isIn == -1 || isIn == 0 ) { this->outbreak( curr, curr->nghbr_[i] ); } // inside | edge
      else {} // shouldn't get here
   } // end neighbor loop

   return;
}

/* Infected Hull */
void bDelTess::infectedHull( bSimplex *curr ) {

   bHex similarity( this->simplID_[ curr->id_ ] ); // copy the infected id
   similarity.filter( (ulong)0xF, 0 ); // identify hull points (infinity points)
   similarity |= this->currPtId_; // add current point
   if( similarity.getNumFlip() < 4 ) {
      printf("[infectedHull] Will create invalid simplex.\n");
      exit(1);
   } // throw if we don't have a simplex
   else if( similarity.getNumFlip() == 5 ) {
      for( int i=0; i < 4; ++i ) {
         similarity ^= i;
         this->newID_.push_back( similarity );
         similarity |= i;
      }
   } // handle the very first simplex
   else {
      this->newID_.push_back( similarity ); // add to list to create
   } // typical handling
   return;
}

/* No Infection */
void bDelTess::noInfection( bSimplex *prev, bSimplex *curr ) {
   this->well_ |= curr->id_;
   bHex similarity( this->simplID_[ prev->id_ ] ); // copy previous id
   similarity &= this->simplID_[ curr->id_ ]; // identify common points
   similarity |= this->currPtId_; // add the current point
   if( similarity.getNumFlip() != 4 ) {
      printf("[noInfection] Will create invalid simplex.\n");
      exit(1);
   } // throw if we don't have a simplex
   this->newID_.push_back( similarity ); // add to list to create
   return;
}

/* Treat */
bool bDelTess::treat() {
   bool valid = true;
   if( (valid = this->addNewSimplexes()) ) {
      this->removeSimplexes();
      this->identifyNeighbors();
   }
   else { valid = this->cleanIteration(); }
   return valid;
}

bool bDelTess::cleanIteration() {
   bool valid = true;
   for( uint k=0; k < this->newList_.size(); ++k ) {
      this->delSimplex( this->newList_[k] );
   }
   if( this->toAdd_.empty() ) { valid = false; }
   else { 
      this->toAdd_.push_back( this->vrtx_.back() ); // place point at end of list
      this->vrtx_.pop_back(); // remove from cu
   }
   return valid;
}

/* Add Simplexes in NewID */
bool bDelTess::addNewSimplexes() {
   for( uint i=0; i < this->newID_.size(); ++i ) {
      int pos;
      pos = this->addSimplex( this->newID_[i] );
      if( pos == -1 ) { return false; }
      this->newList_.push_back( pos );
   }
   return true;
}

int  bDelTess::addSimplex( bHex &id ) {
   int *plist = id.getActive();
   int pos;
   pos = this->addSimplex( plist );
   return pos;
}

int  bDelTess::addSimplex( const int pt[] ) {
   // pt holds the indexes as listed by vrtx_, not src_

   // SIMPLEX :: variables
   bHex id( this->size_ );
   double ptSet[ 12 ];

   // Mask points to be included
   // Save points at infinity at the end
   // Create the simplex ID
   int sX = 0;
   int pX = 0;
   for(int i=0; i < 4; ++i ) {
      id |= pt[i];
      pX = i * 3;
      sX = this->vrtx_[ pt[i] ].pos_ * 3;
      ptSet[pX]   = this->src_[ this->vrtx_[pt[i]].chn_ ]->pnts_[sX];
      ptSet[++pX] = this->src_[ this->vrtx_[pt[i]].chn_ ]->pnts_[++sX];
      ptSet[++pX] = this->src_[ this->vrtx_[pt[i]].chn_ ]->pnts_[++sX];
   }

   // Create the simplex
   bSimplex s;
   bool valid = s.setup( ptSet );

   int pos = 0;
   if( valid ) {
      // Add to data structures
      if( this->open_.empty() ) { // add a new simplex if no open spots
         pos = this->simplex_.size();
         this->simplex_.push_back( s );
         this->simplID_.push_back( id );
      }
      else { // overwrite an old simplex
         pos = this->open_.front();
         this->simplex_[ this->open_.front() ] = s;
         this->simplID_[ this->open_.front() ] = id;
         this->open_.pop_front();
      }
      this->simplex_[ pos ].id_ = pos;
   }
   else { pos = -1; }

   this->isVerified_ = false; // might cause problems -- move back to id ver
   return pos;
}



/* Remove Sick Simplexes */
void bDelTess::removeSimplexes() {
   int *sick = this->sick_.getActive();
   int sickFlip = this->sick_.getNumFlip();

   // Nullify sick simplexes
   for( int i=0; i < sickFlip; ++i ) { this->delSimplex( sick[i] ); }
   sick = NULL;
   return;
}

void bDelTess::delSimplex( const int pos ) {
   this->open_.push_back( pos );
   this->simplID_[pos] = 0; // Clear ID
   for( int i=0; i < 4; ++i ) { // Remove all references to this simplex
      if( this->simplex_[pos].nghbr_[i] == NULL ) { continue; }
      for( int k=0; k < 4; ++k ) {
         if( this->simplex_[pos].nghbr_[i]->nghbr_[k] == NULL ) { continue; }
         if( this->simplex_[pos].nghbr_[i]->nghbr_[k] == &(this->simplex_[pos]) ) {
            this->simplex_[pos].nghbr_[i]->nghbr_[k] = NULL;
            k = 4;
         }
         else {}
      }
      this->simplex_[pos].nghbr_[i] = NULL;
   }
   return;
}

/* Identify Neighbors for New Simplexes */
void bDelTess::identifyNeighbors() {
   
   // Get array
   int *well = this->well_.getActive();
   int wellFlip = this->well_.getNumFlip();

   // Point to neighbors
   for( uint i=0; i < this->newList_.size(); ++i ) {
      
      int nghbrCnt = 0;
      bHex similarity( this->simplID_[ this->newList_[i] ] );
      
      // Check for neighbors amongst other new simplexes
      for( uint k=(i + 1); k < this->newList_.size() && nghbrCnt < 4; ++k ) {
         similarity &= this->simplID_[ this->newList_[k] ]; // find common points
         if( similarity.getNumFlip() == 3 ) {
            while( nghbrCnt < 4 && this->simplex_[ this->newList_[i] ].nghbr_[ nghbrCnt ] != NULL ) { ++nghbrCnt; }
            if( nghbrCnt >= 4 ) { throw "[bDelTess] Too many neighbors."; }
            this->simplex_[ this->newList_[i] ].nghbr_[ nghbrCnt ] = &( this->simplex_[ this->newList_[k] ] );
            ++nghbrCnt;
            for( int m=0; m < 4; ++m ) {
               if( this->simplex_[ this->newList_[k] ].nghbr_[m] == NULL ) {
                  this->simplex_[ this->newList_[k] ].nghbr_[m] = &(this->simplex_[ this->newList_[i] ]);
                  m = 4;
               } // save neighbor in empty slot
            } // loop through neighbor's neighbors
         } // three common points
         similarity = this->simplID_[ this->newList_[i] ]; // reset similarity
      }
      
      // Check for neighbors amongst existing simplexes
      for( int k=0; k < wellFlip && nghbrCnt < 4; ++k ) {

         similarity &= this->simplID_[ well[k] ]; // find common points
         if( similarity.getNumFlip() == 3 ) { 
            while( this->simplex_[ this->newList_[i] ].nghbr_[ nghbrCnt ] != NULL ) { ++nghbrCnt; }
            if( nghbrCnt > 4 ) { throw "[bDelTess] Too many neighbors."; }
            this->simplex_[ this->newList_[i] ].nghbr_[ nghbrCnt ] = &( this->simplex_[ well[k] ] );
            for( int m=0; m < 4; ++m ) {
               if( this->simplex_[ well[k] ].nghbr_[m] == NULL ) {
                  this->simplex_[ well[k] ].nghbr_[m] = &(this->simplex_[ this->newList_[i] ]);
                  m = 4;
               } // save neighbor in empty slot
               else {
               }
            } // loop through neighbor's neighbors
            ++nghbrCnt;
         } // three common points
         similarity = this->simplID_[ this->newList_[i] ]; // reset similarity
      }
   } // end neighbors

   well = NULL;
   return;
}


/* Immunize */
void bDelTess::immunize() {
   // Checks the validity of all new simplexes (and repairs)
   while( ! this->newList_.empty() ) {
      this->resetHandlers();
      this->checkSimplex( this->simplex_[ this->newList_.front() ] );
      this->newList_.pop_front(); // remove checked simplex
      this->reconcile( this->newList_, this->sick_ ); // remove newly deleted
   } // check every new simplex -- even those just made here
   return;
}

/* Reconcile First List */
void bDelTess::reconcile( dei &list, bHex &toRemove ) {
   int *removal = toRemove.getActive();
   int numRemove = toRemove.getNumFlip();
   for( int i=0; i < numRemove; ++i ) {
      for( uint k=0; k < list.size(); ++k ) {
         if( list[k] == removal[i] ) {
            list.erase( list.begin() + k );
         }
      }
   }
}

/* Reset Handlers */
void bDelTess::resetHandlers() {
   int multiplier = this->vrtx_.size() > 80 ? 10 : 6;
   int newSize = this->vrtx_.size() * multiplier;
   this->sick_.resize( newSize );
   this->well_.resize( newSize );
   this->sick_ = 0;
   this->well_ = 0;
   this->newID_.clear();
   this->hullID_.clear();
   return;
}

bool bDelTess::isNew( bSimplex *curr ) {
   return !(this->sick_ & curr->id_);
   return true;
}

void bDelTess::dumbCheck() {
   printf(" == BEGIN DUMB CHECK == ");
   char dumb[128];
   bool doprint = false;
   for( uint i=0; i < this->simplex_.size(); ++i ) {
      memset( dumb, '\0', 128 );
      doprint = false;
      sprintf(dumb, "[%d] ", i); if( this->simplID_[i].empty() ) { continue; }
      for( uint k=0; k < 4; ++k ) {
         if( this->simplex_[i].nghbr_[k] == NULL ) { continue; }
         sprintf(dumb, "%u", this->simplex_[i].nghbr_[k]->id_);
         if( this->simplID_[ this->simplex_[i].nghbr_[k]->id_ ].empty() ) { doprint = true; sprintf( dumb, "(NULL), "); }
         else { sprintf(dumb, ", "); }
      }
      if( doprint ) { printf("%s\n",dumb); }
   }
   printf(" == END DUMB CHECK ==\n");
   return;
}

/* Verify Tessellation */
bool bDelTess::verify() {
   bool validTess = true;

   // Test simplex validity
   bool missingPoints = false;
   bool validSimplex = true;
   //~ vaf test( 3 );
   double test[3] = { 0.0, 0.0, 0.0 };
   bHex missingPt( this->simplID_[0] );
   for( uint i=0; i < this->simplex_.size(); ++i ) {
      if( skip( i ) ) { continue; }
      if( this->simplID_[i].getNumFlip() < 4 ) {
         printf("[bDelTess] Found simplex with out enough bits flipped.\n");
         validTess = false;
      }

      missingPt |= this->simplID_[i]; // accumulate points to test for missing

      for( uint k=4; k < this->vrtx_.size(); ++k ) {
         if( this->simplID_[i] & (int)k ) { continue; } // is point in simplex?

         
         int chn = this->vrtx_[ k ].chn_;
         int pos = this->vrtx_[ k ].pos_;
         //~ printf("checking [%d:%d]\n",chn,pos);
         test[0] = this->src_[chn]->pnts_[ ( 0 + (pos * 3) ) ];
         test[1] = this->src_[chn]->pnts_[ ( 1 + (pos * 3) ) ];
         test[2] = this->src_[chn]->pnts_[ ( 2 + (pos * 3) ) ];
         short inSph = this->simplex_[i].inSphere( test );

         if( inSph == 1 ) { validSimplex= true; }
         else if( inSph == -1 ) {
            validSimplex = false;
         }
         else if( inSph == 0 ) {
            validSimplex = false;
         }
         else { printf("\tscrewed"); }

         if( !validSimplex ) {
            validTess = false;
            break;
         }
         
      } // point loop
   } // simplex loop

   if( validTess ) {

      // Find missing points
      missingPt.assign( (ulong)0xF, 0 );
      missingPt.flipBits();

      if( !(missingPt.empty()) ) {
         // List the missing
         int *mP = missingPt.getActive();
         int nF = missingPt.getNumFlip();
         for( int i=0; i < nF; ++i ) {
            if( mP[i] > (int)(this->vrtx_.size() - 1) ) { i = nF; }
            else {
               missingPoints = true;
            }
         }
         mP = NULL;
      }

   }

   validTess &= !missingPoints;
   this->isVerified_ = validTess;
   return validTess;
}

void bDelTess::retry() {
   bool hadEdge = this->haveEdge_;
   bool hadType = this->haveType_;
   bool wasTrim = this->isTrim_;
   bool wasSlim = this->isSlim_;
   this->reset();
   this->tessellate();
   if( wasTrim ) this->trim();
   if( wasSlim ) this->removeExcess();
   if( hadEdge ) this->findEdges();
   if( hadType ) this->findTypes();

   return;
}

/****** To Work With */
void bDelTess::focus( bPoints* focus ) {

   // Identify focus chain
   int chn = 0;
   for( int i=1; i < this->numChains_; ++i ) {
      if( this->src_[i] == focus ) { chn = i; i = this->numChains_; }
   }
   if( chn == 0 ) { return; }
   else { this->focus( chn ); }
   return;
}

void bDelTess::focus( int chn ) {
   if( ! this->isVerified_ ) { return; printf("[bDelTess::focus] Please validate tessellation.\n"); }
   if( this->numChains_ < 3 || chn >= this->numChains_ ) { return; }
   // Remove non-focus simplexes
   bool hasFocus = false;
   int* alist;
   for( uint i=0; i < this->simplID_.size(); ++i ) {
      if( this->skip(i) ) { continue; }
      
      hasFocus = false;
      alist = this->simplID_[i].getActive();
      if( this->vrtx_[ alist[0] ].chn_ == chn ) { hasFocus = true; }
      else if( this->vrtx_[ alist[1] ].chn_ == chn ) { hasFocus = true; }
      else if( this->vrtx_[ alist[2] ].chn_ == chn ) { hasFocus = true; }
      else if( this->vrtx_[ alist[3] ].chn_ == chn ) { hasFocus = true; }
      else { }
         
      if( !hasFocus ) { this->delSimplex( i ); }
      alist = NULL;
   }

   return;
}

void bDelTess::intersect() {
   if( this->numChains_ < 3 ) { return; }

   bool isIntersect = false;
   int* alist;
   for( uint i=0; i < this->simplID_.size(); ++i ) {
      if( this->skip(i) ) { continue; }
      
      isIntersect = false;
      alist = this->simplID_[i].getActive();
      int chn = this->vrtx_[ alist[0] ].chn_;
      if( this->vrtx_[ alist[1] ].chn_ != chn ) { isIntersect = true; }
      else if( this->vrtx_[ alist[2] ].chn_ != chn ) { isIntersect = true; }
      else if( this->vrtx_[ alist[3] ].chn_ != chn ) { isIntersect = true; }
      else { }
         
      if( !isIntersect ) { this->delSimplex( i ); }
      alist = NULL;
   }
   return;
}

void bDelTess::removeExcess() {
   // This function will remove any infinity simplexes and empty simplexes
   
   if( ! this->isVerified_ ) { return; printf("[bDelTess::excess] Please validate tessellation.\n"); }
   if( this->isSlim_ ) { return; }
   
   //~ triangles_.clear();
   
   //~ printf("[bDelTess] Removing empty and invalid simplexes...\n");
   bool remove;
   uint swapPos = 0;
   
   bSimplex tmpS;
   bHex tmpH;
   //~ uint refSaver[ simplID_.size() ][4];
   for( uint i=0; i < this->simplID_.size(); ++i ) {
      remove = false;
      //~ refSaver[i] = i;
      if( !this->skip( i ) ) { 
         //~ refSaver[i][0] = this->simplex_[i].nghbr_[0]->id_;
         //~ refSaver[i][1] = this->simplex_[i].nghbr_[1]->id_;
         //~ refSaver[i][2] = this->simplex_[i].nghbr_[2]->id_;
         //~ refSaver[i][3] = this->simplex_[i].nghbr_[3]->id_;
         
         continue;
      }

      //~ delSimplex( i );
   //~ printf("...%d\n",i);
      // If the swap pos and current pos are different, change them
      // Also ensure that each pointer is updated
      //~ printf("OLD [(%d) -> %d]\t", swapPos, i); printNghbr( swapPos );
      //~ printf("   \t"); printNghbr( i );
      //~ printf("%4d: %p\t%p\t%p\t%p\n", swapPos,
      delSimplex( i );

      if( swapPos != i ) {
         //~ tmpS = this->simplex_[i];
         //~ tmpH = this->simplID_[i];
         this->simplex_[i] = this->simplex_[swapPos];
         this->simplID_[i] = this->simplID_[swapPos];
         //~ this->simplex_[swapPos] = tmpS;
         //~ this->simplID_[swapPos] = tmpH;
         //~ this->simplex_[swapPos].id_ = swapPos;
         this->simplex_[i].id_ = i;
         for( uint k=0; k < 4; ++k ) {
            if( this->simplex_[i].nghbr_[k] == NULL ) { continue; }
            for( uint m=0; m < 4; ++m ) {
               if( this->simplex_[i].nghbr_[k]->nghbr_[m] == NULL ) { continue; }
               if( this->simplex_[i].nghbr_[k]->nghbr_[m] == &( this->simplex_[swapPos] ) ) {
                  this->simplex_[i].nghbr_[k]->nghbr_[m] = &( this->simplex_[i] );
               }
            }
         }
      }
      //~ printf("New [%d -> (%d)]\t", swapPos, i); 
      //~ printNghbr( swapPos );
      //~ printf("   \t");
      //~ printNghbr( i );
      // Nullify the old swap ( i now holds the swap )
      //~ delSimplex( swapPos );
      //~ printf("Del\t"); printNghbr( swapPos );
      //~ printf("\n");
      //~ this->simplID_[swapPos] = 0x0;
      //~ this->simplex_[swapPos].nghbr_[0] = NULL;
      //~ this->simplex_[swapPos].nghbr_[1] = NULL;
      //~ this->simplex_[swapPos].nghbr_[2] = NULL;
      //~ this->simplex_[swapPos].nghbr_[3] = NULL;
      ++swapPos;
   }

   //~ printf("num to remove: %d\ttotal num: %u\n", pos, this->simplID_.size());
   //~ printf("before: [%3d] ", this->simplex_[swapPos].id_); printNghbr( swapPos );
   //~ bSimplex* del[4];
   for( uint i=0; i < swapPos; ++i ) {
   //~ while( !(this->simplID_.empty()) && this->simplID_.front().empty() ) {
      //~ printf("removing...%d\n", this->simplID_.size());
      //~ printf("removing %d...\n", this->simplex_[i].id_);
      //~ for( int k=0; k < 4; ++k ) {
         //~ del[k] = this->simplex_[i].nghbr_[k];
         //~ printf("{ ");
         //~ for( int m=0; m < 4; ++m ) {
            //~ if( del[k]->nghbr_[m] == NULL ) { printf("NULL, "); }
            //~ else{ printf("%4d (%p), ", del[k]->nghbr_[m]->id_, del[k]->nghbr_[m]); }
         //~ }
         //~ printf(" }, ");
      //~ } printf("\n");
      //~ this->simInf_.push_back( this->simplID_.front() );
      this->simplID_.pop_front();
      this->simplex_.pop_front();
      
      //~ for( int k=0; k < 4; ++k ) {
         //~ printf("{ ");
         //~ for( int m=0; m < 4; ++m ) {
            //~ if( del[k]->nghbr_[m] == NULL ) { printf("NULL, "); }
            //~ else{ printf("%4d (%p), ", del[k]->nghbr_[m]->id_, del[k]->nghbr_[m]); }
         //~ }
         //~ printf(" }, ");
      //~ } printf("\n\n");
      //~ --swapPos;
   }
   //~ printf(" after: [%3d] ", this->simplex_[0].id_); printNghbr( 0 );

   for( uint i=0; i < this->simplex_.size(); ++i ) {
      this->simplex_[i].id_ = (short)i;
   }
   //~ printf(" after: [%3d] ", this->simplex_[0].id_); printNghbr( 0 );
   //~ printf("swapPos: %d\n", swapPos);

   this->haveType_ = false;
   this->isSlim_ = true;
   //~ printf("Info:\n");
   //~ for( int i=0; i < this->numChains_; ++i ) { 
      //~ if( this->src_[i] == 0 ) { printf("[%d] NULL\n",i); }
      //~ else { printf("[%d] %d :: %s\n", i, this->src_[i]->numPnts_, this->src_[i]->aaSeq_); }
   //~ }
   //~ printf("\n");

   return;
}

void bDelTess::printNghbr( int pos ) {
   printf("\n%4d [%p]=> ", this->simplex_[pos].id_, &(this->simplex_[pos]) );
   for( int i=0; i < 4; ++i ) {
      if( this->simplex_[pos].nghbr_[i] == NULL ) { printf(" %87dNULL, ",0); if( i==1 ) printf("\n\t\t   "); continue; }
      else { printf("%4d [%p] { ", this->simplex_[pos].nghbr_[i]->id_, this->simplex_[pos].nghbr_[i]); }
      for( int k=0; k < 4; ++k ) {
         if( this->simplex_[pos].nghbr_[i]->nghbr_[k] == NULL ) { printf(" %11dNULL, ",0); }
         else { printf("%4d [%p], ", this->simplex_[pos].nghbr_[i]->nghbr_[k]->id_, this->simplex_[pos].nghbr_[i]->nghbr_[k]); }
      }
      printf("}, ");
      if( i == 1 ) { printf("\n\t\t   "); }
   }
   printf("\n");
}

void bDelTess::chkNghbr() {
   char out[2048];
   printf("here\n");
   bool doprint = false;
   for( uint pos=0; pos < this->simplex_.size(); ++pos ) {
      memset( out, '\0', 2048 );
      doprint = true;
      sprintf(out, "%4d => ",pos);
      for( int i=0; i < 4; ++i ) {
         if( this->simplex_[pos].nghbr_[i] == NULL ) { sprintf(out, "%sNULL, ",out); continue; }
         else { sprintf(out,"%s%4d { ", out,this->simplex_[pos].nghbr_[i]->id_); }
         bool haveref = false;
         for( int k=0; k < 4; ++k ) {
            if( this->simplex_[pos].nghbr_[i]->nghbr_[k] == NULL ) { sprintf(out, "%sNULL, ",out); continue; }
            else { sprintf(out, "%s%4d, ", out, this->simplex_[pos].nghbr_[i]->nghbr_[k]->id_); }
            if( this->simplex_[pos].id_ == this->simplex_[pos].nghbr_[i]->nghbr_[k]->id_ ) { haveref = true; }
         }
         doprint &= haveref;
         sprintf(out, "%s}, ",out);
      }
      printf("%s\n",out);
   }
   printf("gone\n");
   return;
}

void bDelTess::trim( float threshold ) {
   if( ! this->isVerified_ ) { return; printf("[bDelTess::trim] Please validate tessellation.\n"); }
   if( this->isTrim_ ) { return; }
   if( this->src_[1]->isInRes_ ) { threshold /= this->src_[1]->res_; }

   //~ printf("[bDelTess] Trimming long edges...\n");
   int nPt = this->size_;
   float** dist = new float*[nPt];
   for( int i=0; i < nPt; ++i ) {
      dist[i] = new float[nPt];
      for( int k=0; k < nPt; ++k ) { dist[i][k] = 0.0; }
   }
   float pRef[3];
   float pCmp[3];

   // Measure pt distance in all vertices
   for( uint i=0; i < this->simplID_.size(); ++i ) {
      if( skip( i ) ) { continue; }
      //~ if( ! nullNgh(i) ) { continue; }

      //~ this->trimTet( i, threshold );
      int* list = this->simplID_[i].getActive();
      for( int y=0; y < 4; ++y ) {
         int yP   = list[y]; yP -= 4;
         int yPos = this->vrtx_[ list[y] ].pos_; yPos *= 3;
         int yChn = this->vrtx_[ list[y] ].chn_;

         for( int z=(y+1); z < 4; ++z ) {
            int zP   = list[z]; zP -= 4;
            int zPos = this->vrtx_[ list[z] ].pos_; zPos *= 3;
            int zChn = this->vrtx_[ list[z] ].chn_;

            if( dist[ yP ][ zP ] == 0.0 ) {
               pRef[0] = this->src_[yChn]->pnts_[yPos + 0];
               pRef[1] = this->src_[yChn]->pnts_[yPos + 1];
               pRef[2] = this->src_[yChn]->pnts_[yPos + 2];
               pCmp[0] = this->src_[zChn]->pnts_[zPos + 0];
               pCmp[1] = this->src_[zChn]->pnts_[zPos + 1];
               pCmp[2] = this->src_[zChn]->pnts_[zPos + 2];
               dist[yP][zP] = bPoints::pointDistance( pRef, pCmp ); 
            } // TEST if distance calculated

            if( dist[yP][zP] > threshold ) {
               this->delSimplex( i );
               y = 4;
               z = 4;
            } // DELETE simplex if edges are too long
         } // LOOP upper diagonal
      } // LOOP through vertices
   } // LOOP through simplexes
   
   //~ simplID_[0] = 0;
   this->isTrim_ = true;
   this->isSlim_ = false;
   return;
}


void bDelTess::onion() {
   // Remove the outermost layer (i.e., anything with a NULL neighbor)
   if( ! this->isVerified_ ) { return; printf("[bDelTess::onion] Please validate tessellation.\n"); }

   this->simOni_.clear();
   deque<uint> toDelete;
   int num = 0;
   
   // Loop through all simplexes
   for( uint i=0; i < this->simplex_.size(); ++i ) {
      if( skip( (int)i ) ) { continue; }
      ++num;

      for( uint k=0; k < 4; ++k ) {
         if( this->simplex_[i].nghbr_[k] == NULL || this->simplID_[ this->simplex_[i].nghbr_[k]->id_ ] & (ulong) 0xF ) {
            toDelete.push_back(i);
            break;
         }
      }
   }

   for( uint i=0; i < toDelete.size(); ++i ) {
      this->delSimplex( toDelete[i] );
   }

   this->isTrim_ = false;
   this->haveEdge_ = false;
   this->haveType_ = false;
   this->haveScore_ = false;
   return;
}

/* Find Edges */
void bDelTess::findEdges() {
   if( ! this->isVerified_ ) { throw "[bDelTess::edges] Please validate tessellation"; }

   //~ printf("[bDelTess] Finding edges...\n");
   if( !(this->edges_.empty()) ) { this->edges_.clear(); }

   bHex edge( this->vrtx_.size() );
   for( uint i=0; i < this->vrtx_.size(); ++i ) {
      this->edges_.push_back( edge );
   } // initialize edges

   for( uint i=0; i < this->simplID_.size(); ++i ) {
      if( skip( i ) ) { continue; }

      bHex currID( this->simplID_[i] ); // save ID (temp)
      int *active = currID.getActive(); // get points
      if( currID.getNumFlip() < 4 ) { currID.print(); throw "[bDelTess] Found simplex without enough bits flipped."; continue; }
      for( uint k=0; k < 4; ++k ) {
         currID ^= active[k]; // remove current point
         this->edges_[ active[k] ] |= currID; // save
      }// loop through points
      active = NULL;
   } // loop through all ids

   this->haveEdge_ = true;
   return;
}

void bDelTess::findTypes() {
   if( ! this->isVerified_ ) { throw "[bDelTess::type] Please validate tessellation"; }
   if( ! this->isSlim_ ) { this->removeExcess(); }
   if( this->haveType_ ) { return; }
   if( !(this->simplTy_.empty()) ) { this->simplTy_.clear(); }
   for( uint i=0; i < this->simplID_.size(); ++i ) {
      if( skip( i ) ) { this->simplTy_.push_back( -1 ); continue; }
      this->simplTy_.push_back( this->findType( i ) );
   }
   this->haveType_ = true;
   return;
}

/* Type of Tetrahedral */
int bDelTess::findType( int pos ) {
   if( this->haveType_ && this->isVerified_ ) { return this->simplTy_[pos]; }

   int* list = this->simplID_[pos].getActive();
   int type = 0;

   // Type considerations -- flatten loop for speed
   register int* c = new int[4];
   c[0] = this->vrtx_[ list[0] ].chn_; c[1] = this->vrtx_[ list[1] ].chn_;
   c[2] = this->vrtx_[ list[2] ].chn_; c[3] = this->vrtx_[ list[3] ].chn_;
   register int* p = new int[4];
   p[0] = this->vrtx_[ list[0] ].pos_; p[1] = this->vrtx_[ list[1] ].pos_;
   p[2] = this->vrtx_[ list[2] ].pos_; p[3] = this->vrtx_[ list[3] ].pos_;
   bSort::qsortDbl( c, p, 4 );

   //~ register int p_two = this->vrtx_[ list[1] ].pos_;
   register bool consec = false;
   int p_curr = 0;
   int p_next = 0;

   // Check 0-1
   if( c[0] == c[1] ) {
      p_curr = this->src_[ c[0] ]->pos_ == NULL ? p[0] : this->src_[ c[0] ]->pos_[ p[0] ];
      p_next = this->src_[ c[1] ]->pos_ == NULL ? p[1] : this->src_[ c[1] ]->pos_[ p[1] ];
      ++p_curr;
      if( p_curr == p_next ) { ++type; consec = true; }
   }
   
   // Check 1-2
   if( c[1] == c[2] ) {
      p_curr = this->src_[ c[1] ]->pos_ == NULL ? p[1] : this->src_[ c[1] ]->pos_[ p[1] ];
      p_next = this->src_[ c[2] ]->pos_ == NULL ? p[2] : this->src_[ c[2] ]->pos_[ p[2] ];
      ++p_curr;
      if( p_next == p_curr ) {
         ++type;
         if( consec ) { ++type; }
         consec = true;
      }
      else { consec = false; }
   }
   else { consec = false; }
   
   // Check 2-3
   if( c[2] == c[3] ) {
      p_curr = this->src_[ c[2] ]->pos_ == NULL ? p[2] : this->src_[ c[2] ]->pos_[ p[2] ];
      p_next = this->src_[ c[3] ]->pos_ == NULL ? p[3] : this->src_[ c[3] ]->pos_[ p[3] ];
      ++p_curr;
      if( p_curr == p_next ) {
         ++type;
         if( consec ) { ++type; }
         consec = true;
      }
      else { consec = false; }
   }
   else { consec = false; }
   if( type == 5 ) --type;

   delete [] p;
   delete [] c;
   p = NULL;
   c = NULL;
   
   //~ printf("tet [%d]: %d\n", pos, type);

   return type;
}


float bDelTess::score() {
   if( ! this->isVerified_ ) { printf("[bDelTess::score] Please validate tessellation.\n"); }
   if( !(this->isSlim_) ) { this->removeExcess(); }
   if( !(this->haveType_) ) { this->findTypes(); }
   if( this->haveScore_ ) { return this->score_; }
   if( !(this->simplSc_.empty()) ) { this->simplSc_.clear(); }
   this->score_ = 0.0;
   char* res = new char[5];
   res[0] = '\0';
   int* list = NULL;
   int cnt = 0;
   for( uint i=0; i < this->simplID_.size(); ++i ) {
      if( this->skip(i) ) { continue; }
      
      list = this->simplID_[i].getActive();
      res[0] = this->src_[ this->vrtx_[list[0]].chn_ ]->aaSeq_[ this->vrtx_[list[0]].pos_ ];
      res[1] = this->src_[ this->vrtx_[list[1]].chn_ ]->aaSeq_[ this->vrtx_[list[1]].pos_ ];
      res[2] = this->src_[ this->vrtx_[list[2]].chn_ ]->aaSeq_[ this->vrtx_[list[2]].pos_ ];
      res[3] = this->src_[ this->vrtx_[list[3]].chn_ ]->aaSeq_[ this->vrtx_[list[3]].pos_ ];

      bSort::qsort( res, 4 );
      this->simplSc_.push_back( bSnapp::score( res, this->simplTy_[cnt] ) );
      this->score_ += this->simplSc_.back();
      //~ printf("tet [%u]:\n%-5d\n%-6s\n%-6.2f\n", i, this->simplTy_[cnt], res, this->simplSc_.back());
      ++cnt;
      list = NULL;

      delete [] res;
      res = NULL;
      res = new char[5];


   }
   delete [] res;
   res = NULL;
   
   this->haveScore_ = true;
   return this->score_;
}


/****** Checking */
bool bDelTess::skip( int pos ) const {
   //~ if( this->isSlim_ && this->isVerified_ ) { return false; }
   //~ if( this->isVerified_ ) { return false; }
   if( this->simplID_[pos].empty() ) { return true; }
   else if( (this->simplID_[pos] & (ulong)0xF) ) { return true; }
   else {}
   return false;
}


bool bDelTess::nullNgh( int pos ) {
   bool nullNgh = false;
   if( this->simplex_[pos].nghbr_[0] == NULL ) { nullNgh = true; }
   else if( this->simplex_[pos].nghbr_[1] == NULL ) { nullNgh = true; }
   else if( this->simplex_[pos].nghbr_[2] == NULL ) { nullNgh = true; }
   else if( this->simplex_[pos].nghbr_[3] == NULL ) { nullNgh = true; }
   else if( skip( this->simplex_[pos].nghbr_[0]->id_ ) ) { nullNgh = true; }
   else if( skip( this->simplex_[pos].nghbr_[1]->id_ ) ) { nullNgh = true; }
   else if( skip( this->simplex_[pos].nghbr_[2]->id_ ) ) { nullNgh = true; }
   else if( skip( this->simplex_[pos].nghbr_[3]->id_ ) ) { nullNgh = true; }
   else {}
   return nullNgh;
}

/****** Output */

/* Print Tetrahedral Composition */
/* Print Tetrahedral Composition */

void bDelTess::print( FILE* op ) const {
   printf("[bDelTess]\n");
   printf("\tnum chains: %d\n", this->numChains_);
   printf("\tnum points: %d\n", this->size_);
   printf("\t  is valid: %u\n", this->isVerified_);
   printf("\thave edges: %u\n", this->haveEdge_);
   printf("\thave types: %u\n", this->haveType_);
   printf("\thave score: %u\n", this->haveScore_);
   return;
}

void bDelTess::printTetRes( FILE* op ) {
   if( ! this->isVerified_ ) { throw "[bDelTess::printByRes] Please validate tessellation"; }
   //~ if( !this->haveType_ ) { printf("Unable to print tet by res -- find type first\n"); return; }
   
   for( int i=1; i < this->numChains_; ++i ) {
      if( this->src_[i]->aaSeq_ == NULL ) { return; }
   }
   if( !this->haveType_ ) { this->findTypes(); }

   char* aaList = new char[5];
   memset( aaList, '\0', 5 );
   //~ aaList[4] = '\0';
   int list[4] = {0};
   for( uint i=0; i < this->simplID_.size(); ++i ) {
      if( skip( i ) ) { continue; }

      this->simplID_[i].getActive(list);
      aaList[0] = this->src_[ this->vrtx_[list[0]].chn_ ]->aaSeq_[ this->vrtx_[list[0]].pos_ ];
      aaList[1] = this->src_[ this->vrtx_[list[1]].chn_ ]->aaSeq_[ this->vrtx_[list[1]].pos_ ];
      aaList[2] = this->src_[ this->vrtx_[list[2]].chn_ ]->aaSeq_[ this->vrtx_[list[2]].pos_ ];
      aaList[3] = this->src_[ this->vrtx_[list[3]].chn_ ]->aaSeq_[ this->vrtx_[list[3]].pos_ ];
      bSort::qsort( aaList, 4 );

      fprintf(op, "%s %d %.2f\n", aaList, this->simplTy_[i], this->simplSc_[i]);
   }
   delete [] aaList; aaList = NULL;
   return;
}

void bDelTess::printTetPos( FILE* op ) const {
   if( ! this->isVerified_ ) { throw "[bDelTess::printByTet] Please validate tessellation"; }
   register int* plist = new int[4];
   register int* clist = new int[4];
   int list[4] = {0};
   for( uint i=0; i < this->simplID_.size(); ++i ) {
      if( skip( i ) ) { continue; }

      this->simplID_[i].getActive(list);
      plist[0] = this->vrtx_[ list[0] ].pos_;
      plist[1] = this->vrtx_[ list[1] ].pos_;
      plist[2] = this->vrtx_[ list[2] ].pos_;
      plist[3] = this->vrtx_[ list[3] ].pos_;
      clist[0] = this->vrtx_[ list[0] ].chn_;
      clist[1] = this->vrtx_[ list[1] ].chn_;
      clist[2] = this->vrtx_[ list[2] ].chn_;
      clist[3] = this->vrtx_[ list[3] ].chn_;

      //~ // setup continuous numbering
      for( int i=0; i < 4; ++i ) {
         if( clist[i] > 1 ) {
            for( int k=1; k < clist[i]; ++k ) {
               plist[i] += this->src_[k]->numPnts_;
            }
         }
      }

      bSort::qsort( plist, 4 );

      fprintf(op, "%d %d %d %d\n", 
         plist[0],
         plist[1],
         plist[2],
         plist[3] );
   }
   delete [] plist; plist = NULL;
   delete [] clist; clist = NULL;
   return;
}

void bDelTess::printTetPnt( FILE* op ) const {
   for( int i=1; i < this->numChains_; ++i ) {
      uint pX = 0;
      bPoints* curr = this->src_[i];
      for( int k=0; k < curr->numPnts_; ++k ) {
         fprintf(op, "%.3f ", curr->pnts_[ pX ]); ++pX;
         fprintf(op, "%.3f ", curr->pnts_[ pX ]); ++pX;
         fprintf(op, "%.3f %d:%d %c\n", curr->pnts_[ pX ], i, k, curr->aaSeq_[k]);
      }
   }
   return;
}

void bDelTess::printTetCat( FILE* op ) const {
   if( ! this->isVerified_ ) { throw "[bDelTess::printByRes] Please validate tessellation"; }

   // SNAPP Score
   float avgSc = this->score_; avgSc /= (float) this->size_;
   fprintf(op, "SNAPP    %.3f %.3f\n", this->score_, avgSc);

   // Type Count
   int typeCnt[5] = {0};
   for( uint i=0; i < this->simplTy_.size(); ++i ) {
      ++typeCnt[ this->simplTy_[i] ];
   }
   fprintf(op, "TYPE_CNT %d %d %d %d %d\n", typeCnt[0],typeCnt[1],typeCnt[2],typeCnt[3],typeCnt[4]);
   
   // Residue Count
   for( int i=1; i < this->numChains_; ++i ) {
      int resCnt[26] = {0};
      char* seq = this->src_[i]->aaSeq_;
      for( int k=0; k < this->src_[i]->numPnts_; ++k ) {
         ++resCnt[ seq[k] - 'A' ];
      }
      int cnt = 5;
      for( int k=0; k < 26; ++k ) {
         if( ++cnt == 6 ) {
            fprintf(op, "\nCHN_%-5d",i);
            cnt = 0;
         }
         fprintf(op, "%c:%-5d", k + 'A', resCnt[k]);
      }
   }
   return;
}

/* PyMol */
void bDelTess::pymol( FILE* op, char givenName[], char givenColor[], char label[] ) {

   //~ char colorRed[5] = "ruby";
   char colorG90[7] = "gray90";
   char colorG40[7] = "gray40";
   char colorBlu[5] = "blue";
   char colorGre[6] = "green";
   char colorPnk[5] = "pink";
   char colorBlk[6] = "black";

   //~ char* name = new char[13];
   int len = strlen( givenName );
   char* name = new char[ len + 5 ];
   memset( name, '\0', len + 5 );
   memmove( name, givenName, len );

   //~ bDocker::_pymolHeader( op, 1 );
   this->pymolPseudoatoms( op, givenName, givenColor );

   // Print points causing invalid simplexes
   for( int m=0; m < (int)this->pntInv_.size(); ++m ) {
      fprintf(op, "cmd.show(\"dots\",\"%s//%d/%d/\")\n", givenName, this->pntInv_[m].chn_, this->pntInv_[m].pos_ );
   }

   /* Print Edges */ if( this->printEdg_ ) {
      this->pymolEdges( op, givenName, givenColor, label );
   }

   /* Print Invalid */ if( this->printInv_ ) {
      memmove( name + len, "_inv", 5 );
      this->pymolSimplexes( op, name, colorPnk, this->simInv_ );
      memset( name + len, '\0', 5 );
   }

   /* Print Trimmed */ if( this->printTrm_ && this->simTrm_.size() > 0 ) {
      printf("Printing invalid simplexes...\n");
      memmove( name + len, "_trm", 5 );
      this->pymolSimplexes( op, name, colorBlu, this->simTrm_ );
      memset( name + len, '\0', 5 );
   }
   
   /* Print Onions */ if( this->printOni_ && this->simOni_.size() > 0 ) {
      memmove( name + len, "_oni", 5 );
      this->pymolSimplexes( op, name, colorGre, this->simOni_ );
      memset( name + len, '\0', 5 );

   }

   /* Print Simplexes */ if( this->printVal_ ) {
      printf("Printing valid simplexes %u[%ld]...\n", printVal_, this->simplex_.size() );
      memmove( name + len, "_val", 5 );
      this->pymolSimplexes( op, name, colorG90, this->simplID_, true );
      memset( name + len, '\0', 5 );
   }

   /* Print Infinity */ if( this->printInf_ ) {
      printf("Printing infinity simplexes [%ld]...\n", this->simInf_.size() );
      memmove( name + len, "_inf", 5 );
      this->pymolSimplexes( op, name, colorG40, this->simInf_ );
      memset( name + len, '\0', 5 );
   }

   /* Print Removed */ if( this->printRem_ ) {
      memmove( name + len, "_rem", 5 );
      this->pymolSimplexes( op, name, colorBlk, this->simRem_ );
      memset( name + len, '\0', 5 );
   }

   delete [] name;
   name = NULL;
   //~ printf("done\n");
   return;
}

/* PyMol Points */
void bDelTess::pymolPseudoatoms( FILE* op, char name[], char color[] ) {
   for( int i=1; i < this->numChains_; ++i ) {
      int index = 0;
      for( int k=0; k < this->src_[i]->numPnts_; ++k ) {
         fprintf(op, "cmd.pseudoatom(\"%s\",pos=[", name );
         fprintf(op, "%.2f,", this->src_[i]->pnts_[index]); ++index;
         fprintf(op, "%.2f,", this->src_[i]->pnts_[index]); ++index;
         fprintf(op, "%.2f]", this->src_[i]->pnts_[index]); ++index;
         if( this->src_[i]->aaSeq_ != NULL ) { fprintf(op, ", resn=\"%c\"", this->src_[i]->aaSeq_[k]); }
         fprintf(op, ", resi=\"%d\"", k);
         fprintf(op, ", segi=\"1\"");
         fprintf(op, ", chain=\"%d\",name=\"SCC\")\n", i );
      }
      
   }
   fprintf(op, "cmd.color(\"%s\",\"%s\")\n", color, name);
}

/* PyMol Edges */
void bDelTess::pymolEdges( FILE* op, char name[], char color[], char label[] ) {
   if( ! this->isVerified_ ) { throw "[bDelTess::pymol:edges] Please validate tessellation"; }
   for( uint i=4; i < this->edges_.size(); ++i ) {
      int* active = this->edges_[i].getActive();
      for( int k=0; k < this->edges_[i].getNumFlip(); ++k ) {
         fprintf(op,"cmd.bond(\"%s//%d/%d/\",\"%s//%d/%d/\")\n",
            name, this->vrtx_[i].chn_, this->vrtx_[i].pos_,
            name, this->vrtx_[ active[k] ].chn_, this->vrtx_[ active[k] ].pos_ );
      }
      active = NULL;
   }

   if( label != NULL ) {
      int i = 0;
      while( this->edges_[i].getNumFlip() < 2 ) ++i;
      int* tolabel = this->edges_[i].getActive();
      bPoints* a = this->src_[ this->vrtx_[ tolabel[0] ].chn_ ];
      bPoints* b = this->src_[ this->vrtx_[ tolabel[1] ].chn_ ];
      int ap = this->vrtx_[ tolabel[0] ].pos_; ap *= 3;
      int bp = this->vrtx_[ tolabel[1] ].pos_; bp *= 3;
      float c[3] = { a->pnts_[ ap ], a->pnts_[ ++ap ], a->pnts_[ ++ap ] };
      c[0] += b->pnts_[ bp ]; c[1] += b->pnts_[ ++bp ]; c[2] += b->pnts_[ ++bp ];
      c[0] /= 2; c[1] /= 2; c[2] /= 2;
      fprintf(op,"cmd.pseudoatom(\"%s\",segi=\"LABEL\",pos=[%.2f,%.2f,%.2f],label=\"%s\")\n", name, c[0], c[1], c[2], label );
      a = NULL;
      b = NULL;
   }

   fprintf(op,"cmd.color(\"%s\",\"%s\")\n", color, name ); // color the simplex

   return;
}

void bDelTess::pymolSimplexes( FILE* op, char name[], char color[], deH &ids, bool onlyVal ) {
   //~ if( ! this->isVerified_ ) { throw "[bDelTess::pymol:edges] Please validate tessellation"; }
   char n[32];
   sprintf(n,"%s",name);
   for( uint i=0; i < ids.size(); ++i ) {
      if( onlyVal && this->skip(i) ) { continue; }
      int* active = ids[i].getActive();
         sprintf(name, "%s_%03d",n,i);
      
      ushort chn[4] = { 0, 0, 0, 0 };
      ushort pos[4] = { 0, 0, 0, 0 };
      char res[4]; memset( res, '\0', 4 );
      for( uint k=0; k < 4; ++k ) {
         chn[k] = this->vrtx_[ active[k] ].chn_;
         pos[k] = this->vrtx_[ active[k] ].pos_;
         uint nX = pos[k] * 3;
         if( this->src_[chn[k]]->haveSeq_ ) memmove( res, bAA::aa1to3(this->src_[ chn[k] ]->aaSeq_[ pos[k] ]), 3 );
         fprintf(op,"cmd.pseudoatom(\"%s\", pos=[%.2f, %.2f, %.2f], segi=\"%d\", chain=\"%d\", resi=\"%d\", resn=\"%s\", name=\"%d\")\n",
            name, this->src_[ chn[k] ]->pnts_[nX], this->src_[ chn[k] ]->pnts_[nX + 1], this->src_[ chn[k] ]->pnts_[nX + 2],
            i, chn[k], pos[k], res, k );
         for( uint m=0; m < k; ++m ) {
            fprintf(op,"cmd.bond(\"%s/%d///%d\",\"%s/%d///%d\")\n",
               name, i, k, name, i, m );
         }
      }
   }
   fprintf(op, "cmd.color(\"%s\",\"%s\")\n", color, name);
   return;
}
