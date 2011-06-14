//~ #include <iostream>
#include <fstream>
#include <cmath>
//~ #include <valarray>
#include <string>
#include <sstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <sys/stat.h>
//~ #include "bVBase.h"
#include "bPoints.h"
#include "bAA.h"
#include "bSys.h"
using namespace std;
using namespace bStd;

const int DEBUG_COUNT=-1;

/*--------------------------------*/
/***** CONSTRUCTOR & DESTRUCTORS */
bPoints::bPoints( int c )
   : capPnts_(c), pnts_( NULL ), tets_( NULL ), pos_( NULL ), aaSeq_( NULL )
{

   // set flags
   this->isInRes_ = 0;
   this->isInPos_ = 0;
   this->numPnts_ = 0;
   this->havePnts_ = false;
   this->haveTets_ = false;
   this->haveMM_ = false;
   this->haveSeq_ = false;

   // allocate memory
   this->pntPath_ = NULL; //new char[96];
   this->pntBase_ = NULL; //new char[32];

   this->min_[0] = 0.0; this->max_[0] = 0.0;
   this->min_[1] = 0.0; this->max_[1] = 0.0;
   this->min_[2] = 0.0; this->max_[2] = 0.0;

   this->planeDisplacement_[0] = 0.0;
   this->planeDisplacement_[1] = 0.0;
   this->planeDisplacement_[2] = 0.0;

   this->fit_ = 0.0;
   this->thk_ = 0.0;
   this->res_ = 0.0;
   this->offset_ = 0.0;
   
   // check for point capacity
   if( this->capPnts_ != 0 ) { this->resize( this->capPnts_ ); }
}

bPoints::bPoints( const bPoints &rhs )
   : capPnts_(0), pnts_( NULL ), tets_( NULL ), pos_( NULL ), aaSeq_( NULL )
{
   // copy flags
   this->havePnts_ = rhs.havePnts_;
   this->haveTets_ = rhs.haveTets_;
   this->isInRes_ = rhs.isInRes_;
   this->isInPos_ = rhs.isInPos_;
   this->haveMM_ = rhs.haveMM_;
   this->haveSeq_ = rhs.haveSeq_;

   // Copy meta data
   for(int i=0; i<3; ++i) {
      this->min_[i] = rhs.min_[i];
      this->max_[i] = rhs.max_[i];
      this->planeDisplacement_[i] = rhs.planeDisplacement_[i];
   }

   // Copy data
   if( rhs.capPnts_ > 0 ) {
      this->pnts_ = new float[ rhs.capPnts_ * 3 ];
      this->pos_ = new uint[ rhs.capPnts_ ];
      for( int i=0; i < rhs.capPnts_ * 3; ++i ) { this->pnts_[i] = rhs.pnts_[i]; }
      for( int i=0; i < rhs.capPnts_; ++i ) { this->pos_[i] = rhs.pos_[i]; }
   }
   //~ else { this->pnts_ = NULL; }

   // Counters
   this->numPnts_ = rhs.numPnts_;
   this->capPnts_ = rhs.capPnts_;

   if( rhs.tets_ != NULL ) {
      this->tets_ = new bPoints;
      this->tets_ = rhs.tets_;
   }
   else { this->tets_ = NULL; }

   if( rhs.aaSeq_ != NULL ) {
      int len = strlen( rhs.aaSeq_ );
      this->aaSeq_ = new char[len + 1];
      memset( this->aaSeq_, '\0', len + 1 );
      memmove( this->aaSeq_, rhs.aaSeq_, len );
   }
   
   // Copy file info
   if( rhs.pntBase_ != NULL && rhs.pntPath_ != NULL ) {
      this->pntBase_ = new char[32];
      this->pntPath_ = new char[96];
      memset( this->pntBase_, '\0', 32 );
      memset( this->pntPath_, '\0', 96 );
      memmove( this->pntBase_, rhs.pntBase_, strlen(rhs.pntBase_) );
      memmove( this->pntPath_, rhs.pntPath_, strlen(rhs.pntPath_) );
   }
   else {
      this->pntBase_ = NULL;
      this->pntPath_ = NULL;
   }

   // Copy grid parameters
   this->fit_ = rhs.fit_;
   this->thk_ = rhs.thk_;
   this->res_ = rhs.res_;
   this->offset_ = rhs.offset_;
}

bPoints::~bPoints() {
   if( this->pnts_ != NULL ) {
      delete [] this->pnts_;
      this->pnts_ = NULL;
   }

   if( this->tets_ != NULL ) {
      delete this->tets_;
      this->tets_ = NULL;
   }

   if( this->pos_ != NULL ) {
      delete [] this->pos_;
      this->pos_ = NULL;
   }

   if( this->aaSeq_ != NULL ) {
      delete [] this->aaSeq_;
      this->aaSeq_ = NULL;
   }

   if( this->pntPath_ != NULL ) {
      delete [] this->pntPath_;
      this->pntPath_ = NULL;
   }

   if( this->pntBase_ != NULL ) {
      delete [] this->pntBase_;
      this->pntBase_ = NULL;
   }
}

/*---------------------------*/
/***** OVERLOADED OPERATORS */
bPoints& bPoints::operator=( const bPoints &rhs ) {
   if(this == &rhs ) { return *this; }// check for self assignment
   
   // copy flags
   this->havePnts_ = rhs.havePnts_;
   this->haveTets_ = rhs.haveTets_;
   this->isInRes_ = rhs.isInRes_;
   this->isInPos_ = rhs.isInPos_;
   this->haveMM_ = rhs.haveMM_;
   this->haveSeq_ = rhs.haveSeq_;

   // Copy meta data
   for(int i=0; i<3; ++i) {
      this->min_[i] = rhs.min_[i];
      this->max_[i] = rhs.max_[i];
      this->planeDisplacement_[i] = rhs.planeDisplacement_[i];
   }

   // Copy data
   if( rhs.capPnts_ > 0 ) {
      if( this->capPnts_ == rhs.capPnts_ ) { this->clear(); }
      else { this->resize( rhs.capPnts_ ); }
      for( int i=0; i < rhs.capPnts_ * 3; ++i ) {
         this->pnts_[i] = rhs.pnts_[i];
      }
   }
   else { this->pnts_ = NULL; }

   if( rhs.aaSeq_ != NULL ) {
      int len = strlen( rhs.aaSeq_ );
      this->aaSeq_ = new char[len + 1];
      memset( this->aaSeq_, '\0', len + 1 );
      memmove( this->aaSeq_, rhs.aaSeq_, len );
   }

   // Counters
   this->numPnts_ = rhs.numPnts_;
   this->capPnts_ = rhs.capPnts_;

   if( this->tets_ != NULL ) {
      delete this->tets_;
      this->tets_ = NULL;
   }
   if( rhs.tets_ != NULL ) {
      this->tets_ = new bPoints;
      *(this->tets_) = *(rhs.tets_);
   }

   // Copy file info
   if( rhs.pntBase_ != NULL && rhs.pntPath_ != NULL ) {
      this->pntBase_ = new char[32];
      this->pntPath_ = new char[96];
      memset( this->pntBase_, '\0', 32 );
      memset( this->pntPath_, '\0', 96 );
      memmove( this->pntBase_, rhs.pntBase_, strlen(rhs.pntBase_) );
      memmove( this->pntPath_, rhs.pntPath_, strlen(rhs.pntPath_) );
   }
   else {
      this->pntBase_ = NULL;
      this->pntPath_ = NULL;
   }

   this->fit_ = rhs.fit_; // copy grid parameters
   this->thk_ = rhs.thk_;
   this->res_ = rhs.res_;
   this->offset_ = rhs.offset_;

   return *this;
}

/*---------------------------*/
/***** OBJECT MANIPULATION */
void bPoints::clear() {
   this->resize( this->capPnts_ );
   this->haveMM_ = false;
   return;
}

int  bPoints::size( int num ) {
   return this->numPnts_;
}
void bPoints::resize( int num ) {
   if( this->capPnts_ == num && this->numPnts_ != 0 ) { return; }

   if( this->pnts_ != NULL ) {
      delete [] this->pnts_;
      this->pnts_ = NULL;
   }
   if( this->pos_ != NULL ) {
      delete [] this->pos_;
      this->pos_ = NULL;
   }
   
   this->numPnts_ = 0;
   this->capPnts_ = num;
   this->pnts_ = new float[this->capPnts_ * 3];
   this->pos_ = new uint[ this->capPnts_ ];
   
   for( int i=0; i < this->capPnts_ * 3; ++i ) { this->pnts_[i] = 0.0; } // init to 0
   for( int i=0; i < this->capPnts_; ++i ) { this->pos_[i] = 0; } // init to 0

   this->haveMM_ = false;
   return;
}

void bPoints::resizeCopy( int num ) {
   if( this->capPnts_ == num ) { return; }
   else if( this->pnts_ == NULL || this->numPnts_ == 0 ) { this->resize( num ); return; }
   else {}

   int sizePnt = this->capPnts_ * 3;
   float tempPnt[ sizePnt ];

   // Save, Erase, Resize, Rewrite :: Pnts
   for( int i=0; i < sizePnt; ++i ) { tempPnt[i] = this->pnts_[i]; } // save
   delete [] this->pnts_;
   this->pnts_ = NULL;
   this->pnts_ = new float[ num * 3 ];
   for( int i=0; i < sizePnt; ++i ) { this->pnts_[i] = tempPnt[i]; } // reload
   for( int i=sizePnt; i < num * 3; ++i ) { this->pnts_[i] = 0; } // init to 0
   this->capPnts_ = num;

   // Save position data
   uint* tmp = new uint[ num ];
   for( int i=0; i < this->numPnts_; ++i ) { tmp[i] = this->pos_[i]; }
   for( int i=this->numPnts_; i < num; ++i ) { tmp[i] = 0; }
   delete [] this->pos_;
   this->pos_ = NULL;
   this->pos_ = tmp;

   return;
}

/*--------------------*/
/***** DATA SOURCING */
bool bPoints::haveFile() {
   char* file = NULL;
   bool haveFile = false;
   if( this->pntPath_ != NULL && this->pntBase_ != NULL ) {
      bSys::multicat( file, this->pntPath_, this->pntBase_, ".out" );
      haveFile = bSys::fileExists( file );
      printf("pntfile: %s [%u]\n", file, haveFile);
   }
   delete [] file;
   file = NULL;
   return haveFile;
}

bool bPoints::cla(int numArg, char** argArray) {
   // Outdated -- see backups
   // loop through all given parameters (i=1 to skip command)
   return false;
}

bool bPoints::setSrc( const char base[], const char path[], bool doWrite ) {
   bool validPath = false;

   // Initialize memory
   this->pntPath_ = new char[96];
   this->pntBase_ = new char[32];
   memset( this->pntPath_, '\0', 96 );
   memset( this->pntBase_, '\0', 32 );

   // Copy the path and base
   strcpy( this->pntPath_, path );
   strcpy( this->pntBase_, base );

   // Remove extension if needed
   uint i = strcspn( this->pntBase_, "." );
   if( i != strlen(this->pntBase_) ) { bSys::removeExt( this->pntBase_ ); }

   // Save the lengths, and add the last slash if needed
   int ppl = strlen( this->pntPath_ );
   int pbl = strlen( this->pntBase_ );
   if( this->pntPath_[ ppl - 1 ] != '/' ) { strcat( this->pntPath_, "/" ); ++ppl; }

   // Construct the needed file name
   char file[96];
   memset( file, '\0', 96 );
   memmove( file, this->pntPath_, ppl );
   memmove( file + ppl, this->pntBase_, pbl );
   if( memcmp( file + strlen(file) - 4, ".out", 4 ) != 0 ) strcat( file, ".out" );

   // Check if the file exists, otherwise try a different path
   if( bSys::fileExists(file) ) { validPath = true; }
   else {
      strcat( this->pntPath_, base );
      strcat( this->pntPath_, "/" );
      ppl = strlen( this->pntPath_ );
      memset( file, '\0', 96 );
      memmove( file, this->pntPath_, ppl );
      memmove( file + ppl, this->pntBase_, pbl );
      if( memcmp( file + strlen(file) - 4, ".out", 4 ) != 0 ) strcat( file, ".out" );
      if( bSys::fileExists(file) ) { validPath = true; }
   }

   // Save the filenames
   if( doWrite ) printf("[bPoints] Files:\n\tBase: %s\n\tPath: %s\n", this->pntBase_, this->pntPath_ );

   //~ if( validPath ) this->setTets();
   return validPath;
}

bool bPoints::setTets() {
   if( this->haveTets_ ) return true;

   bool haveTetFile = false;
   char* file = NULL;
   bSys::multicat( file, this->pntPath_, this->pntBase_, "_CENTROIDS.out" );

   if( bSys::fileExists( file ) ) {
      if( this->tets_ != NULL ) { delete [] this->tets_; this->tets_ = NULL; }
      haveTetFile = true;
      this->tets_ = new bPoints;
      this->tets_->pntPath_ = new char[96];
      this->tets_->pntBase_ = new char[32];
      memset( this->tets_->pntPath_, '\0', 96 );
      memset( this->tets_->pntBase_, '\0', 32 );
      memmove( this->tets_->pntPath_, this->pntPath_, strlen(this->pntPath_) );
      memmove( this->tets_->pntBase_, this->pntBase_, strlen(this->pntBase_) );
      bSys::multicat( this->tets_->pntBase_, "_CENTROIDS" );
      this->tets_->isInPos_ = this->isInPos_;
      this->tets_->isInRes_ = this->isInRes_;
   }
   else { haveTetFile = false; }

   delete [] file;
   file = NULL;
   return haveTetFile;
}

void bPoints::setGridParam( float r, float f, float t, float o, float b ) {
   float mod = isInRes_ ? r : 1;
   this->res_ = r;
   this->fit_ = f * mod;
   this->thk_ = t * mod;
   this->offset_ = o * mod;
   this->buffer_ = b * mod;
   return;
}

void bPoints::_saveFilename( char saveAs[], char file[] ) {
   strcpy(saveAs,file);
   return;
}
void bPoints::_saveFilenameNoExt( char saveAs[], char file[] ) {
   strncpy(saveAs,file,strrchr(file,'.')-file);
   return;
}

/*-----------------*/
/***** DATA INPUT */
/* Read in Protein from File */
bool bPoints::readPoints() {

   char* file = NULL;
   if( !(this->havePnts_) ) {
      bSys::multicat( file, this->pntPath_, this->pntBase_, ".out" );
      if( bSys::fileExists( file ) ) {
         this->numPnts_ = this->readPoints( file );
         this->havePnts_ = this->numPnts_ > 0 ? true : false;
      }
      else { printf("[bPoints] can't open: %s\n", file); throw "[bPoints] Unable to open protein file."; }
   }
   delete [] file;
   file = NULL;

   if( DEBUG_COUNT != -1) {
      this->numPnts_ = DEBUG_COUNT;
      if( this->haveTets_ ) this->tets_->numPnts_ = DEBUG_COUNT;
   }

   this->capPnts_ = this->numPnts_;
   this->getSeq();
   this->_findMinMax();
   return this->havePnts_;
}

/* Read Points from a File */
int bPoints::readPoints( char* file ) {
   // open the file
   ifstream ip;
   ip.open(file, ifstream::in);
   if(!ip) { return 0; }

   // read in the file... not the most efficient, but we need the size
   // -- we read in w/o doing anything to get a cound of the residues
   deque<string> data;
   string bffr;
   uint i = 0;
   for( i = 0; getline( ip, bffr ); ++i ) { data.push_back(bffr); }
   ip.close();

   // resize the array and prepare the stream
   // -- this is an unfortunate side effect of using the valarray
   // => it clears everything upon resize =/
   //~ pntSet.resize(data.size()*3);
   this->pnts_ = new float[ i * 3 ];
   this->pos_ = new uint[ i ];

   // loop through each line and save the coordinates
   int index = 0;
   for( i = 0; i < data.size(); ++i ) {
      istringstream ss(data[i]);
      float flt;

      // loop through each value on the line
      for(int j = 0; ss >> flt; ++j) {
         this->pnts_[index] = flt;
         ++index;
      }
   }

   for( uint k=0; k < i; ++k ) { this->pos_[k] = k; }

   this->numPnts_ = i;
   this->capPnts_ = i;
   this->_findMinMax();

   return i;
}

/* STATIC Read Points */
int bPoints::_readPoints(char* file, float* &pntSet) {
   // open the file
   ifstream ip;
   ip.open(file, ifstream::in);
   if(!ip) { return 0; }

   // read in the file... not the most efficient, but we need the size
   // -- we read in w/o doing anything to get a cound of the residues
   deque<string> data;
   string bffr;
   uint i = 0;
   for( i = 0; getline( ip, bffr ); ++i ) { data.push_back(bffr); }
   ip.close();

   // resize the array and prepare the stream
   // -- this is an unfortunate side effect of using the valarray
   // => it clears everything upon resize =/
   //~ pntSet.resize(data.size()*3);
   pntSet = new float[ i * 3 ];

   // loop through each line and save the coordinates
   int index = 0;
   for( i = 0; i < data.size(); ++i ) {
      istringstream ss(data[i]);
      float flt;

      // loop through each value on the line
      for(int j = 0; ss >> flt; ++j) {
         pntSet[index] = flt;
         ++index;
      }
   }

   return i;
}

void bPoints::addPoint( float newPt[], int pos ) {
   // Take another look at this later...seems redundant with addPointAtPos (as coded)
   // Also, double check for out-of-bounds
   if( this->capPnts_ == 0 ) { this->resize( pos ); }
   else if( this->numPnts_ == this->capPnts_ ) { this->resizeCopy( this->capPnts_ + 4 ); }
   else {}

   this->pos_[ this->numPnts_ ] = (pos == -1) ? this->numPnts_ : pos;

   pos = this->numPnts_; pos *= 3;
   this->pnts_[pos] = newPt[0];
   this->pnts_[++pos] = newPt[1];
   this->pnts_[++pos] = newPt[2];
   ++this->numPnts_;
   this->haveMM_ = false;
   return;
}

void bPoints::addPointAtPos( float newPt[], int pos ) {
   if( this->capPnts_ == 0 ) { this->resize(4); }
   else if( this->numPnts_ == this->capPnts_ ) { this->resizeCopy( this->capPnts_ + 4 ); }
   else {}

   if( pos > this->capPnts_ ) { this->resizeCopy( pos + 4 ); }
   pos -= 1; // position to index
   this->pos_[ pos ] = pos;
   pos *= 3; // 

   this->pnts_[pos] = newPt[0];
   this->pnts_[++pos] = newPt[1];
   this->pnts_[++pos] = newPt[2];
   ++this->numPnts_;
   this->haveMM_ = false;
   return;
}

void bPoints::addPoints( float newPts[], int num, bool clear ) {
   //~ printf("size: %d, %d => %d\n", this->numPnts_, this->capPnts_, num);
   if( clear || capPnts_ == 0 ) { this->resize( num ); }
   else if( this->numPnts_ + num > capPnts_ ) {
      //~ printf("copy resize\n");
      this->resizeCopy( this->numPnts_ + num );
   }
   else{}


   int pos = this->numPnts_ * 3;
   int last = num * 3;
   //~ for( int i=0; i < num; ++i ) { 
      //~ printf("[%d] %u ... %d\n", i, this->pos
      //~ this->pos_[ this->numPnts_ + i ] = this->numPnts_ + i; }
   //~ printf("size: %d, %d => %d {%d, %d}\n", this->numPnts_, this->capPnts_, num, pos, last);
   for( int i=0; i < last; ++i ) {
      //~ printf("[%d] %6.2f, %6.2f\n", i, this->pnts_[pos+i], newPts[i]);
      this->pnts_[ pos + i ] = newPts[i];
   }
   this->numPnts_ += num;
   return;
}

void bPoints::getSeq() {
   if( !this->pntPath_ ) { throw "[getSeq] Path not set.\n"; return; }

   // Check for file
   char file[96];
   memset( file, '\0', 96 );
   strcpy( file, this->pntPath_ );
   strcat( file, this->pntBase_ );
   strcat( file, ".seq" );
   if( !(bSys::fileExists( file )) ) { return; }

   // Setup String
   this->aaSeq_ = new char[ this->numPnts_ + 1 ];
   memset( this->aaSeq_, '\0', this->numPnts_ + 1 );

   // Open File
   ifstream ip;
   ip.open( file, ifstream::in );
   string bffr;

   for( int i=0; getline( ip, bffr ); ++i ) {
      istringstream ss( bffr );
      char code[4];
      ss >> code;
      this->aaSeq_[i] = bAA::aa3to1( code );
   }

   ip.close();
   this->haveSeq_ = true;
   return;
}
void bPoints::setSeq( char seq[] ) {
   if( this->aaSeq_ != NULL ) {
      delete [] this->aaSeq_;
      this->aaSeq_ = NULL;
   }

   int len = strlen( seq );
   this->aaSeq_ = new char[len + 1];
   strcpy( this->aaSeq_, seq );

   return;
}

/****** SIZE AND SPACE */
void bPoints::sizeAndSpace( const bPoints &rhs, int nP ) {
   if(this == &rhs ) { return; }// check for self assignment
   if( nP > 0 ) this->resize( nP );

   // Data
   this->planeDisplacement_[0] = rhs.planeDisplacement_[0];
   this->planeDisplacement_[1] = rhs.planeDisplacement_[1];
   this->planeDisplacement_[2] = rhs.planeDisplacement_[2];
   this->fit_ = rhs.fit_; // copy grid parameters
   this->thk_ = rhs.thk_;
   this->res_ = rhs.res_;
   
   // Move if necessary
   if( this->pnts_ != NULL ) {
      if( rhs.isInPos_ && rhs.isInRes_ ) {
         if( ! isInPos_ ) {
            _translatePointPlane();
         }
         if( ! isInRes_ ) {
            _changePointSpace();
         }
      }
      else {
         if( isInRes_ ) {
            _changePointSpace();
         }
         if( isInPos_ ) {
            _translatePointPlane();
         }
      }
      this->_findMinMax();
   }
   else {
      // just initializing -- no points yet
      this->havePnts_ = false;
      this->haveTets_ = false;

      // we assume that we're already in docking space
      this->isInRes_ = rhs.isInRes_;
      this->isInPos_ = rhs.isInPos_;
   }

   // Copy data
   return;
}

void bPoints::prep() {
   if( this->planeDisplacement_[0] == (this->planeDisplacement_[1] == (this->planeDisplacement_[2] == 0.0)) ) {
      // calculate the amount required to move into the positive plane for each axis
      // adjust min and max to encompass the grid, not just the protein
      // NOTE: we don't adjust fit or thk b/c bGrid should not have put them in resolution yet!!
      this->_findMinMax(); // get the min and max of the point set
      this->planeDisplacement_[0] = this->min_[0] - (this->offset_ + this->buffer_);
      this->planeDisplacement_[1] = this->min_[1] - (this->offset_ + this->buffer_);
      this->planeDisplacement_[2] = this->min_[2] - (this->offset_ + this->buffer_);
   }
   if( this->haveTets_ ) { this->tets_->sizeAndSpace( *this ); }

}


/* Move all points into the resolution space for docking

>> Description of problem

  -2  -1   0   1   2   3 (in angstroms)
           0 1 2 3 4 5 6 (in angstroms/[res_=0.5])
   |---|---|---|---|---|
    *[-1.75]     *[1.5]
1) <--(f+t)              : pointDisplacement = min - buffer
   pD = -1.75 - ([fit=1] + [thk=1]) = -3.75
2) (pD)-->   (pD)-->     : shift all points by pD
3)                       : put in resolution
  -2  -1   0   1   2   3   4   5     (in angstroms)
           0 1 2 3 4 5 6 7 8 9 A B   (in angstroms/[res_=0.5])
   |---|---|---|---|---|---|---|---|
                   *[2~4]       *[5.25~10.5]
4)                       : place point in a box of size 2*buffer
  -2  -1   0   1   2   3   4   5     (in angstroms)
           0 1 2 3 4 5 6 7 8 9 A B C D E F  (in angstroms/[res_=0.5])
   |---|---|---|---|---|---|---|---|---|---|
                   *_           *
           |---|---|_|---|---|  _
             t   f     |---|---|_|---|---|
5)                       : min = (int)pointMin - (t+f)
                         : max = (int)pointMax + (t+f) + 1
   min = 2 - (1+1) | 4 - (2+2) = 0
   max != [5 + (1+1) | 10 + (2+2) = 14] ~> 15!
   
   It is important to wait until in resolution to get the min/max.
   Also important to add one to the value of max
   Remember: max is the last needed index!
*/
bool bPoints::isInDockingSpace() { return this->isInRes_ & this->isInPos_; }
bool bPoints::setToDockingSpace() { return setToDockingSpace( *this ); }
bool bPoints::setToDockingSpace( std::deque<bPoints*> &p, int nC ) {
   bool worked = false;
   for(int i=0; i < nC; ++i) { // loop through all point sets
      worked &= setToDockingSpace( *p[i] );
   }
   return worked;
}
bool bPoints::setToDockingSpace( bPoints &p ) {
   // make sure we have the protein data
   if( ! p.isInPos_ ) {

     // move the data points for protein and tetrahedrals
      p._translatePointPlane();
   }
   
   if( ! p.isInRes_ ) {
      // change to resolution space and then find the min and max
      // -- we're in the positive plane
      // -- must change to resolution space BEFORE getting max and min
      p._changePointSpace();
   }
   p._findMinMax();

   // Also change tets if we have 'em
   if( p.haveTets_ ) { p.tets_->setToDockingSpace(); }

   return (p.havePnts_ & p.haveTets_);
}

/* Set to Normal Space */
bool bPoints::setToNormalSpace() { return setToNormalSpace( *this ); }
bool bPoints::setToNormalSpace( std::deque<bPoints*> &p, int nC ) {
   for(int i=0; i < nC; ++i) { // loop through all point sets
      setToNormalSpace( *p[i] );
   }
   return true;
}

bool bPoints::setToNormalSpace( bPoints &p ) {
   if(p.isInRes_) {
      p._changePointSpace();
      
   }
   if(p.isInPos_) {
      p._translatePointPlane();
   }
   p._findMinMax();

   // Also change tets if we have 'em
   if( p.haveTets_ ) { p.tets_->setToNormalSpace(); }

   return true;
}

/* Translate points back and forth while in normal space */
bool bPoints::_translatePointPlane() {
   // skip if we're in resolution
   if( isInRes_ ) return isInPos_;

   // if we're in the positive plane, go to the original & vice versa
   int direction = 1;
   if( isInPos_ ) { direction = -1; }
   for( int i=0; i < 3; ++i ) {
      this->planeDisplacement_[i] *= direction;
   } // redirect displacement

   // move each axis
   int index = 0;
   for( int i=0; i < this->numPnts_; ++i ) {
      for( int k=0; k < 3; ++k ) {
         this->pnts_[index] -= planeDisplacement_[k];
         ++index;
      } // axis loop
   } // point loop

   for( int i=0; i < 3; ++i ) {
      this->planeDisplacement_[i] *= direction;
   } // reset displacement
   
   // flip the bit and return
   this->haveMM_ = false;
   isInPos_ ^= 1;
   return isInPos_;
}

/* Diffract points from normal to resolution space and back */
bool bPoints::_changePointSpace() {

   // doint do this unless we're in the positive plane
   if(!this->isInPos_) { return this->isInRes_; }

   // go to resolution space, unless we're there already
   float diffraction = 1 / this->res_;
   if(this->isInRes_) { diffraction = this->res_; }

   // diffract / refract....whatever...
   for( int i=0; i < this->numPnts_ * 3; ++i ) {
      this->pnts_[i] *= diffraction;
   }

   // flip the bit and return
   this->haveMM_ = false;
   this->isInRes_ ^= 1;
   return isInRes_;
}

/*-------------------*/
/***** MEASUREMENTS */
/* Find Min and Max */
void bPoints::_findMinMax() {
   if( this->haveMM_ ) return;

   for( int i=0; i < 3; ++i ) {
      this->min_[i] = this->pnts_[i];
      this->max_[i] = this->pnts_[i];
   } // reset max and min
   
   int index = 0;
   for( int i=0; i < this->numPnts_; ++i ) {
      for( int k=0; k < 3; ++k ) {
         if( this->min_[k] > this->pnts_[index] )
            this->min_[k] = this->pnts_[index];
         if( this->max_[k] < this->pnts_[index] )
            this->max_[k] = this->pnts_[index];
         ++index;
      } // loop through axes
   } // loop through points

   this->haveMM_ = true;
   return;
}

/* FindNearbyPoints */
bool bPoints::findNearbyPoints( const bPoints &searchPnts, bPoints &nearby, bool extend ) {
   bool foundNearby = false;
   
   // Calculate the box dimensions
   if( !(this->haveMM_) ) { this->_findMinMax(); }
   float ma[3], mi[3];
   float bnd = extend ? 3.0 : 2.2;
   for( int i=0; i < 3; ++i ) {
      mi[i] = this->min_[i]; mi[i] -= (this->fit_ * bnd); //1.8);
      ma[i] = this->max_[i]; ma[i] += (this->fit_ * bnd); //1.8);
   }

   // Identify which points in the search space are nearby
   // -- Test bounds of each point (per axis) and save the index if valid
   deque<int> found;
   int index = 0;
   for( int i=0; i < searchPnts.numPnts_; ++i ) {
      bool isValid = true;
      index = i * 3;
      for(int k=0; k < 3 && isValid; ++k) {
         if( searchPnts.pnts_[index + k] > ma[k] ) { isValid = false; continue; }
         else if( searchPnts.pnts_[index + k] < mi[k] ) { isValid = false; continue; }
         else {}
      } // loop axis
      if( isValid ) { found.push_back( i ); }
   } // loop points
   
   if( found.size() > 0 ) {
      foundNearby = true;
      // Prepare nearby
      nearby.pntPath_ = new char[96];
      nearby.pntBase_ = new char[32];
      memset( nearby.pntPath_, '\0', 96 );
      memset( nearby.pntBase_, '\0', 32 );
      memmove( nearby.pntPath_, searchPnts.pntPath_, strlen(searchPnts.pntPath_) );
      memmove( nearby.pntBase_, searchPnts.pntBase_, strlen(searchPnts.pntBase_) );
      nearby.sizeAndSpace( *this, found.size() );
      nearby.aaSeq_ = new char[ found.size() + 1 ];
      memset( nearby.aaSeq_, '\0', found.size() + 1 );

      // Save search points
      for( uint k=0; k < found.size(); ++k ) {
         int fpos = found[k] * 3;
         
         nearby.aaSeq_[k] = searchPnts.aaSeq_[ found[k] ];
         
         float npt[3];
         npt[0] = searchPnts.pnts_[ fpos ];
         npt[1] = searchPnts.pnts_[ ++fpos ];
         npt[2] = searchPnts.pnts_[ ++fpos ];
         nearby.addPoint( npt, searchPnts.pos_[ found[k] ] );
         
      } // loop points
   }

   return foundNearby;
}

/* Calculate Point Distances */
float bPoints::pointDistance(float* a, float* b) {
   float c[3] = { 0, 0, 0 };
   float sum = 0;
   for( int i=0; i < 3; ++i ) {
      c[i] = a[i] - b[i];
      c[i] *= c[i];
      sum += c[i];
   }
   return sqrt( sum );
}

double bPoints::rmsd( const bPoints &rhs ) {
   return rmsd( *this, rhs );
}

double bPoints::rmsd( const bPoints &lhs, const bPoints &rhs ) {
   if( lhs.numPnts_ != rhs.numPnts_ ) { return -1.0; }
   register double rmsd = 0.0;
   register double x, y, z;
   int pX = 0;
   for( int i=0; i < lhs.numPnts_; ++i ) {
      x = lhs.pnts_[ pX ];   x -= rhs.pnts_[ pX ];
      y = lhs.pnts_[ ++pX ]; y -= rhs.pnts_[ pX ];
      z = lhs.pnts_[ ++pX ]; z -= rhs.pnts_[ pX ];
      x *= x; y *= y; z *= z;
      rmsd += x; rmsd += y; rmsd += z;
   }
   rmsd /= lhs.numPnts_;
   rmsd = sqrt( rmsd );
   return rmsd;
}

/*--------------*/
/****** ROTATE */
bool bPoints::rotateTheta( float point[], int theta ) {
   // Rotate Theta -- inclination
   register double t = theta * _PI_ / 180; // adjust to radians
   double rotMat[9] = {0.0};
   rotMat[0] = cos(t); // initialize rotation matrix
   rotMat[1] = -sin(t);
   rotMat[3] = -rotMat[1];
   rotMat[4] = rotMat[0];
   rotMat[8] = 1;

   rotate( point, rotMat );
   return true;
}
bool bPoints::rotatePhi( float point[], int phi) {
   // Rotate Phi -- Azimuth 
   register float t = phi * _PI_ / 180; // adjust to radians
   double rotMat[9] = {0.0};
   rotMat[0] = cos(t); // initialize rotation matrix
   rotMat[2] = sin(t);
   rotMat[4] = 1;
   rotMat[6] = -rotMat[2];
   rotMat[8] = rotMat[0];

   rotate( point, rotMat );
   return true;
}

void bPoints::rotate( float point[], const double rotMat[] ) {
   double ptCopy[3] = { point[0], point[1], point[2] };
   double temp[3] = { 0.0 };
   for( int i=0; i < 3; ++i ) {
      int index = i * 3;
      temp[0] = rotMat[ index ];
      temp[1] = rotMat[ ++index ];
      temp[2] = rotMat[ ++index ];
      temp[0] *= ptCopy[0];
      temp[1] *= ptCopy[1];
      temp[2] *= ptCopy[2];
      point[i] = temp[0];
      point[i] += temp[1];
      point[i] += temp[2];
   }
   return;
}

/****** Conversion Functions */
bool bPoints::c2s( float* pt ) {
   float t[3] = { pt[0], pt[1], pt[2] };
   t[0] *= t[0]; t[1] *= t[1]; t[2] *= t[2];
   t[0] = t[0] + t[1] + t[2]; t[0] = sqrt( t[0] );
   pt[2] = atan2( pt[2], pt[0] );
   pt[1] = asin( pt[1] / t[0] );
   pt[0] = t[0];
   pt[1] *= 180 / _PI_;
   pt[2] *= 180 / _PI_;
   return true;
}

bool bPoints::s2c( float* pt ) {
   pt[1] *= _PI_ / 180;
   pt[2] *= _PI_ / 180;
   float t[3] = { pt[0], pt[1], pt[2] };
   pt[0] = t[0] * cos(t[1]) * cos(t[2]);
   pt[1] = t[0] * sin(t[1]);
   pt[2] = t[0] * cos(t[1]) * sin(t[2]);
   return true;
}

/****** PyMol */
void bPoints::pymol( FILE *op, char name[], char color[], int chain, char chnLabel[] ) const {
   this->pymolConnectedPseudoatoms( op, name, color, chain, chnLabel );
   return;
}
void bPoints::pymol_dc( FILE *op, char name[], char color[], int chain, char chnLabel[] ) const {
   this->pymolConnectedPseudoatoms( op, name, color, chain, chnLabel );
   return;
}
void bPoints::pymolPoints( FILE *op, float* pnts, int numPnts, char name[] ) {
   char color[6] = "green";
   float size = 0.2;
   pymolPoints( op, pnts, numPnts, name, color, size );
   return;
}
void bPoints::pymolPoints( FILE *op, float* pnts, int numPnts, char name[], char color[], float size) {
   fprintf(op,"%s = [\n",name);
   for( int i=0; i < numPnts; ++i ) {
      int index = i*3;
      fprintf(op,"\tSPHERE, %.2f, %.2f, %.2f, %.2f,\n",
         pnts[index], pnts[index+1], pnts[index+2], size);
   }
   fprintf(op,"\t]\n");
   fprintf(op,"cmd.load_cgo(%s,'%s')\n",name,name);
   fprintf(op,"cmd.color( \"%s\", \"%s\" )\n", color, name );
   return;
}
void bPoints::pymolConnectedPoints( FILE *op, float* pnts, int numPnts, char name[], char color[], float size ) {
   fprintf(op,"%s = [\n",name);
   fprintf(op,"\tCOLOR, %s\n", color);
   for(int i=0; i<numPnts; ++i) {
      int index = i*3;
      fprintf(op,"\tSPHERE, %.2f, %.2f, %.2f, %.2f,\n",
         pnts[index], pnts[index+1], pnts[index+2], size);
   }
   fprintf(op,"\tLINEWIDTH, 2.0,\n");
   fprintf(op,"\tBEGIN, LINES,\n");
   fprintf(op,"\tCOLOR, %s\n", color);
   for(int i=1; i<numPnts; ++i) {
      int index = i*3;
      fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",pnts[index],pnts[index+1],pnts[index+2]);
      fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",pnts[index-3],pnts[index-2],pnts[index-1]);
   }
   fprintf(op,"\tEND\n");
   fprintf(op,"\t]\n");
   fprintf(op,"cmd.load_cgo(%s,'%s')\n",name,name);
   return;
}


/* PYMOL CONNECTED :: Multiple bPoints */
void bPoints::pymolConnectedPseudoatoms( FILE *op, bPoints** src, int numSrc, char name[], char color[], char** chnLabel ) {
   char* label;
   for( int i=0; i < numSrc; ++i ) {
      label = (chnLabel == NULL) ? NULL : chnLabel[i];
      src[i]->pymolConnectedPseudoatoms( op, name, color, i, label );
      label = NULL;
   }
   return;
}

/* PYMOL CONNECTED :: self */
void bPoints::pymolConnectedPseudoatoms( FILE *op, char name[], char color[], int chain, char chnLabel[] ) const {
   char* labels = this->aaSeq_;
   this->pymolConnectedPseudoatoms( op, this->pnts_, this->numPnts_, name, color, labels, chain, chnLabel );
   labels = NULL;
   return;
}

/* PYMOL CONNECTED */
void bPoints::pymolConnectedPseudoatoms( FILE *op, float* pnts, int numPnts, char name[], char color[], char label[], int chain, char chnLabel[] ) {
   pymolPseudoatoms( op, pnts, numPnts, name, color, label, chain, chnLabel );
   for( int i=1; i < numPnts; ++i ) {
      fprintf(op,"cmd.bond(\"%s//%d/%d/\",\"%s//%d/%d/\")\n", name, chain, i, name, chain, i-1 );
   }
   return;
}

/* PYMOL :: Multiple bPoints */
void bPoints::pymolPseudoatoms( FILE *op, bPoints** src, int numSrc, char name[], char color[], char** chnLabel ) {
   char* label;
   for( int i=0; i < numSrc; ++i ) {
      label = (chnLabel == NULL) ? NULL : chnLabel[i];
      src[i]->pymolPseudoatoms( op, name, color, i );
      label = NULL;
   }
   return;
}

/* PYMOL :: self */
void bPoints::pymolPseudoatoms( FILE* op, char name[], char color[], int chain, char chnLabel[] ) const {
   char* labels = (this->aaSeq_ == NULL) ? NULL : this->aaSeq_;
   this->pymolPseudoatoms( op, this->pnts_, this->numPnts_, name, color, labels, chain, chnLabel );
   labels = NULL;
   return;
}
void bPoints::pymolPseudoatoms( FILE *op, float* pnts, int numPnts, char name[], char color[], char label[], int chain, char chnLabel[] ) {

   // print points
   int index = 0;
   for( int i=0; i < numPnts; ++i ) {
      fprintf(op,"cmd.pseudoatom(\"%s\",pos=[%.2f,", name, pnts[index]); ++index;
      fprintf(op,"%.2f,",pnts[index]); ++index;
      fprintf(op,"%.2f],",pnts[index]); ++index;
      fprintf(op,"chain=\"%d\",", chain);
      fprintf(op,"resi=\"%d\",", i);
      if( label != NULL ) fprintf(op,"resn=\"%s\",", bAA::aa1to3(label[i]) );
      fprintf(op,"name=\"SCC\")\n");
   }
   
   if( chnLabel != NULL ) {
      float mid[3] = { pnts[0], pnts[1], pnts[2] };
      mid[0] += pnts[3]; mid[1] += pnts[4]; mid[2] += pnts[5];
      mid[0] /= 2; mid[1] /= 2; mid[2] /= 2;
      fprintf(op,"cmd.pseudoatom(\"%s\",segi=\"LABEL\",pos=[%.2f,%.2f,%.2f],label=\"%s\")\n", name, mid[0], mid[1], mid[2], chnLabel );
   }
   
   fprintf(op,"cmd.color(\"%s\",\"%s//%d//\")\n", color, name, chain); // color peptide
   if( strcmp( color, "purple" ) == 0 ) {
      fprintf(op,"cmd.color(\"white\",\"%s//%d/0/\")\n", name, chain); // color first residue
   }
   else {
      fprintf(op,"cmd.color(\"purple\",\"%s//%d/0/\")\n", name, chain); // color first residue
   }
   //~ fprintf(op,"cmd.clip(\"slab\",100)\n");
   return;
}

/****** Print functions */
void bPoints::printPoint( const int* pt ) {
   for(int i=0; i<3; ++i) { printf("%d | ",pt[i]); } printf("\n");
   return;
}
void bPoints::printPoint( const float* pt ) {
   for(int i=0; i<3; ++i) { printf("%.2f | ",pt[i]); } printf("\n");
   return;
}
void bPoints::print( FILE* op ) const {
   fprintf(op, "numPoints: %d\n", this->numPnts_);
   fprintf(op, "min: "); this->printPoint( this->min_ );
   fprintf(op, "max: "); this->printPoint( this->max_ );
   if( this->aaSeq_ != NULL ) { fprintf(op, "Seq: %s\n", this->aaSeq_); }
   this->printPoints( this->pnts_, this->numPnts_, op );
   return;
}
//~ void bPoints::printPoints( const float* pnts, int num ) { printPoints( stout, pnts, num ); }
void bPoints::printPoints( const float* pnts, int num, FILE* op ) {
   for( int i=0; i < num; ++i ) {
      fprintf(op, "\t%d: ",i);
      for(int j=0; j<3; ++j) {
         int index = (i*3)+j;
         fprintf(op, "%#4.2f | ",pnts[index]);
      }
      fprintf(op, "\n");
   }
   fprintf(op, "\n");
   return;
}





/* END bPoints */


/* BEGIN bPVar */
bPVar::bPVar () :
   pos_(0), var_(NULL), res_(NULL)
{
}
bPVar::bPVar( double pnt[], uint p, char v, char r ) :
   pos_(p), var_(v), res_(r)
{
   this->pnt_[0] = pnt[0]; this->pnt_[1] = pnt[1]; this->pnt_[2] = pnt[2];
}
bPVar::~bPVar()
{
}

void bPVar::set( double pnt[], uint p, char v, char r ) {
   this->pos_ = p;
   this->var_ = v;
   this->res_ = r;
   this->pnt_[0] = pnt[0]; this->pnt_[1] = pnt[1]; this->pnt_[2] = pnt[2];
   return;
}
void bPVar::setPos( uint p ) { this->pos_ = p; return; }
void bPVar::setVar( char v ) { this->var_ = v; return; }
void bPVar::setRes( char r ) { this->res_ = r; return; }
void bPVar::setPnt( double pnt[] ) {
   this->pnt_[0] = pnt[0];
   this->pnt_[1] = pnt[1];
   this->pnt_[2] = pnt[2];
   return;
}

uint bPVar::pos() { return this->pos_; }
char bPVar::var() { return this->var_; }
char bPVar::res() { return this->res_; }
void bPVar::pnt( double pnt[] ) {
   pnt[0] = this->pnt_[0]; pnt[1] = this->pnt_[1]; pnt[2] = this->pnt_[2];
   return;
}

/* END bPVar */


