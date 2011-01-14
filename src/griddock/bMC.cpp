
#include <deque>
#include <list>
#include <fstream>
#include <sstream>

#include <string.h>
#include <stdio.h>

#include "bMC.h"
#include "bCalc.h"
#include "bHex.h"
#include "bDelTess.h"
#include "bGrid.h"
#include "bPoints.h"
#include "bSort.h"
#include "bList.h"


using namespace std;
using namespace bStd;

bMC::bMC( bGrid* g, bPoints* s, bConf* p ) : 
   orgGrd_(g), orgSrc_(s), orgPep_(p), havePerturbSchema_(false), perturbSchema_(NULL), perturbDone_(NULL), numOpt_(0), threshold_(0.0), numIter_(1000)
{
   this->perturbRangeMin_[0] = -2;
   this->perturbRangeMin_[1] = -2;
   this->perturbRangeMin_[2] = -2;
   this->perturbRangeMax_[0] = 2;
   this->perturbRangeMax_[1] = 2;
   this->perturbRangeMax_[2] = 2;
   
   this->currPnt_ = 0;
   this->currPtb_ = 0;
   
   //~ if( c == NULL ) { myList_ = true; this->conf_ = new bList; }
   //~ this->confSc_ = NULL;
   //~ this->confOr_ = NULL;
}

bMC::bMC( const bMC& rhs ) : 
   orgGrd_(NULL), orgSrc_(NULL), orgPep_(NULL), havePerturbSchema_(false), perturbSchema_(NULL), perturbDone_(NULL), numOpt_(0)
{
   this->orgGrd_ = rhs.orgGrd_;
   this->orgSrc_ = rhs.orgSrc_;
   this->orgPep_ = rhs.orgPep_;
   //~ this->confSc_ = new float[ this->conf_->size() ];
   //~ this->confOr_ = new int[ this->conf_->size() ];
   //~ for( uint i=0; i < this->conf_->size(); ++i ) {
      //~ this->confSc_[i] = rhs.confSc_[i];
      //~ this->confOr_[i] = rhs.confOr_[i];
   //~ }
}

bMC::~bMC() {
   if( this->perturbSchema_ != NULL ) {
      for( int i=0; i < this->numOpt_; ++i ) {
         if( this->perturbSchema_[i] != NULL ) {
            delete [] this->perturbSchema_[i];
            this->perturbSchema_[i] = NULL;
         }
      }
      delete [] this->perturbSchema_;
      this->perturbSchema_ = NULL;
   }
   
   if( this->perturbDone_ != NULL ) {
      for( int i=0; i < this->orgPep_->pp_->numPnts_; ++i ) {
         if( this->perturbDone_[i] != NULL ) {
            delete this->perturbDone_[i];
            this->perturbDone_[i] = NULL;
         }
      }
      delete [] this->perturbDone_;
      this->perturbDone_ = NULL;
   }
   
   //~ if( this->confSc_ != NULL ) { delete [] this->confSc_; this->confSc_ = NULL; }
   //~ if( this->confOr_ != NULL ) { delete [] this->confOr_; this->confOr_ = NULL; }
   
   //~ this->gpnts.clear();
   this->orgGrd_ = NULL;
   this->orgSrc_ = NULL;
   this->orgPep_ = NULL;
}

void bMC::clear() {
   //~ if( this->perturbSchema_ != NULL ) {
      //~ for( int i=0; i < this->numOpt_; ++i ) {
         //~ if( this->perturbSchema_[i] != NULL ) {
            //~ delete [] this->perturbSchema_[i];
            //~ this->perturbSchema_[i] = NULL;
         //~ }
      //~ }
      //~ delete [] this->perturbSchema_;
      //~ this->perturbSchema_ = NULL;
   //~ }
   
   //~ if( this->confSc_ != NULL ) { delete [] this->confSc_; this->confSc_ = NULL; }
   //~ if( this->confOr_ != NULL ) { delete [] this->confOr_; this->confOr_ = NULL; }


   for( int i=0; i < this->orgPep_->pp_->numPnts_; ++i ) { this->perturbDone_[i]->erase(); }
   this->orgGrd_ = NULL;
   this->orgSrc_ = NULL;
   this->orgPep_ = NULL;
}

void bMC::setup( const bGrid &grd, const bPoints &src, bConf &pep ) {
   this->orgGrd_ = &grd;
   this->orgSrc_ = &src;
   this->orgPep_ = &pep;
   //~ pep.print();
   return;
}

bool bMC::metropolis() {
   //~ this->orgPep_->pp_->print();
   const bPoints* sor = this->orgSrc_;
   bPoints* pep = this->orgPep_->pp_;
   int      num = this->orgPep_->pp_->numPnts_;
   float    scr = this->orgPep_->sc_;
   bList*   conf = this->orgPep_->new_chld();
   this->perturbTested_ = 0;
   ushort allTested = (1 >> num) - 1;

   // Declare and initialize the original set
   //~ bConf clean;
   //~ *clean.pp_ = *pep;

   // Setup stock
   bPoints* nearby = new bPoints;
   pep->findNearbyPoints( *sor, *nearby, true );
   bPoints* dtsrc[2] = { nearby, pep };
   bDelTess** dtstock = new bDelTess*[ num ];
   for( int i=0; i < num; ++i ) {
      dtstock[i] = new bDelTess;
      dtstock[i]->tessellate_l1o( dtsrc, 2, 2, i );

      //~ dtstock[i]->verify();
      //~ dtstock[i]->findEdges();
      //~ char name[32]; sprintf( name, "p%d", i);
      //~ char color[32] = "red";
      //~ FILE *op = fopen("d2.py", "a");
      //~ dtstock[i]->pymolDelTess( op, name, color );
      //~ fclose(op);
      //~ printf("size: %d | %d | %d\n", dtstock[i]->src_[0]->numPnts_,dtstock[i]->src_[1]->numPnts_,dtstock[i]->src_[2]->numPnts_);

   }
//~ exit(1);
   // Setup Perturb
   int min[3] = { -3, -3, -3 }, max[3] = { 3, 3, 3 };
   this->setupPerturb( min, max );

   // Perform set number of iterations
   this->numIter_ = this->numOpt_;
   this->numIter_ *= num; // temporary!!
   this->numIter_ *= .2;
   //~ numIter_ = 26;

   //~ bList conf;
   int ppt = 0;
   bool valid = false;
   for( uint i=0; i < this->numIter_; ++i ) {
      bConf* curr = new bConf( *(pep) );

      // Perturb
      if( (ppt = this->perturb( *curr )) != -1) {
         // Copy and finish tessellation
         //~ printf("point: %d\n",ppt);
         //~ printf("[%d] size: %d | %d | %d\n", ppt, dtstock[ppt]->src_[0]->numPnts_,dtstock[ppt]->src_[1]->numPnts_,dtstock[ppt]->src_[2]->numPnts_);
         bDelTess dt( *(dtstock[ ppt ]) );
         
         valid = dt.tessellate_l1o_fin( curr->pp_, 2, 2 );
         //~ dt.score();
         //~ conf->print();
         //~ curr->print();

         
      //~ char name[32]; sprintf( name, "pf%d_%d",ppt,i);
      //~ char color[32] = "red";
      //~ FILE *op = fopen("d2.py", "a");

      //~ dt.verify();
      //~ dt.isVerified_ = true;
      //~ dt.findEdges();
      //~ dt.pymolDelTess( op, name, color );
      //~ dtstock[ppt]->verify();
      //~ dtstock[ppt]->findEdges();
      //~ dtstock[ppt]->pymolDelTess( op, name, color );
      //~ fclose(op);
         // Test and score
         if( !valid ) {
            --i;
            continue;
         }
         else {
            if( dt.score() > scr ) {
               curr->sc_ = dt.score();
               conf->push_back( curr );
            }
            //~ else { printf("score too low: %.2f < %.2f\n",dt.score(), scr); }
         }
      }
      else if( this->perturbTested_ != allTested ) {
         break;
      }
      else {
         --i;
      }
      curr = NULL;
   }

   //~ conf->print();
   printf("done bMC...%d | %d\n",conf->num(),conf->cap());
   bool done = false;
   if( !conf->empty() ) {
      conf->sort();
      done = true;
   }
   
   for( int i=0; i < num; ++i ) {
      delete dtstock[i];
      dtstock[i] = NULL;
   }
   delete [] dtstock;
   dtstock = NULL;
   
   pep = NULL;
   conf = NULL;
   sor = NULL;
   //~ conf->print();
   //~ printf("%.2f\t%.2f\n", (*conf)[0].sc_, (*conf)[1].sc_);
   //~ printf("%u\t%u\n%u\n%u\n\n", &(*conf)[0], &(*conf)[1], conf->beg(),conf->end());

   return done;// > 0 ? true : false;
}

void bMC::setupPerturb( const int min[3], const int max[3] ) {
   int num = this->orgPep_->pp_->numPnts_;

   // Calculate the range & the total number of options
   int range[3] = { max[0] - min[0], max[1] - min[1], max[2] - min[2] };
   ++range[0]; ++range[1]; ++range[2];
   this->numOpt_ = range[0]; this->numOpt_ *= range[1]; this->numOpt_ *= range[2]; --this->numOpt_;
   
   // Setup the perturbation matrix
   this->perturbSchema_ = new int* [this->numOpt_];
   int x = 0;
   for( int i=min[0]; i <= max[0]; ++i ) {
      for( int k=min[1]; k <= max[1]; ++k ) {
         for( int m=min[2]; m <= max[2]; ++m ) {
            if( i == 0 && k == 0 && m == 0 ) continue;
            this->perturbSchema_[x] = new int [3];
            this->perturbSchema_[x][0] = i;
            this->perturbSchema_[x][1] = k;
            this->perturbSchema_[x][2] = m;
            ++x;
         }
      }
   }
   
   // Setup the perturbation checklist
   this->perturbDone_ = new bHex* [num];
   for( int i=0; i < num; ++i ) {
      this->perturbDone_[i] = new bHex( this->numOpt_ );
   }
   this->currPnt_ = 0;
   this->currPtb_ = 0;
   
   // Save the variables
   this->perturbRangeMin_[0] = min[0]; this->perturbRangeMin_[1] = min[1]; this->perturbRangeMin_[2] = min[2];
   this->perturbRangeMax_[0] = max[0]; this->perturbRangeMax_[1] = max[1]; this->perturbRangeMax_[2] = max[2];
   this->numOpt_ = x;
   this->havePerturbSchema_ = true;
   return;
}

int bMC::perturb( bConf &mod ) {
   if( !this->havePerturbSchema_ ) { this->setupPerturb( this->perturbRangeMin_, this->perturbRangeMax_); }
   
   
   //~ printf("\nPertub:\n");
   //~ for( uint i=0; i < perturbDone_.size(); ++i ) {
      //~ printf("\t[%u]\t",i);
      //~ if( perturbDone_[i].empty() ) { printf("\n"); continue; }
      //~ for( list<ushort>::iterator it = perturbDone_[i].begin(); it != perturbDone_[i].end(); ++it ) {
         //~ printf("%u\t",*it);
      //~ }
      //~ printf("\n");
   //~ }

   
   float ptrb[3] = { 0.0, 0.0, 0.0 };
   bPoints* mpt = mod.pp_;
   int whichPnt = 0;
   ushort whichPtb = 0;
   bool valid = false;
   int numAttempts = 0;
   //~ printf("perturb [%lu] ", this->perturbDone_.size());
   //~ printf("[%lu]", this->perturbDone_[currPnt_].size());
   //~ printf("num     [%d] [%d]\n", mpt->numPnts_, this->numOpt_);
   do {
      // Randomly choose the point and perturbation
      //~ printf("num: %d\n", mpt->numPnts_);
      whichPnt = bCalc::getRandomNumber( mpt->numPnts_ );
      whichPtb = bCalc::getRandomNumber( this->numOpt_ );
      
      // Deterministic
      //~ whichPnt = this->currPnt_;
      //~ whichPtb = this->currPtb_;
      //~ printf("[%u ~ %u]\t[%u ~ %u]\n", whichPnt, currPnt_, whichPtb, currPtb_);
      
      //~ mpt->print();
      //~ printf("size of perturb: %d, [ %d ]\n", this->perturbDone_.size(), whichPnt);

      // Check that the perturbation hasn't been done before
      if( !(*(this->perturbDone_[ whichPnt ]) & whichPtb) ) {
         //~ printf("No elements [%u]\t",whichPnt);
         valid = true;
         *(this->perturbDone_[ whichPnt ]) |= whichPtb;
      }
      //~ else {
         //~ printf("Elements   [%u, %lu]\t",whichPnt,this->perturbDone_[whichPnt].size());
         //~ it = this->perturbDone_[ whichPnt ].begin();
         //~ int cnt = 0;
         
         //~ while( *it < whichPtb && it != this->perturbDone_[ whichPnt ].end() ) { ++it; }//printf("%u > %u\n",*it, whichPtb); if( ++cnt > 10 ) exit(1); }
         //~ if( *it != whichPtb ) {
            //~ valid = true;
            //~ this->perturbDone_[ whichPnt ].insert( it, whichPtb );
         //~ }
      //~ }
      
      if( valid ) {
         //~ printf("Doing      : [%u, %u]\n", whichPnt, whichPtb);
         int wpX = whichPnt; wpX *= 3;
         ptrb[0] = mpt->pnts_[wpX]; ++wpX;
         ptrb[1] = mpt->pnts_[wpX]; ++wpX;
         ptrb[2] = mpt->pnts_[wpX];
      
         ptrb[0] += this->perturbSchema_[ whichPtb ][0];
         ptrb[1] += this->perturbSchema_[ whichPtb ][1];
         ptrb[2] += this->perturbSchema_[ whichPtb ][2];
         
         valid = this->orgGrd_->isaGridPoint( ptrb );
      }
      //~ else {
         //~ printf("Already done: [%u, %u]\n", whichPnt, whichPtb);
      //~ }

      // Deterministic
      //~ if( !valid ) {
         //~ ++this->currPtb_;
         //~ if( this->currPtb_ == this->numOpt_ ) { ++this->currPnt_; this->currPtb_ = 0;}
         //~ if( this->currPnt_ == mpt->numPnts_ ) break;
      //~ }
      if( ++numAttempts > 600 ) { break; }
   } while( !valid );

   if( valid ) {
      int wpX = whichPnt; wpX *= 3;
      mpt->pnts_[wpX] = ptrb[0]; ++wpX;
      mpt->pnts_[wpX] = ptrb[1]; ++wpX;
      mpt->pnts_[wpX] = ptrb[2];
      //~ whichPnt /= 3;
   }
   else {
      this->perturbTested_ |= ( 1 >> whichPnt );
      whichPnt = -1;
   }
   //~ printf("perturb: %d | %d [ %d, %d, %d ]\n",whichPnt, whichPtb, perturbSchema_[ whichPtb ][0], perturbSchema_[ whichPtb ][1], perturbSchema_[ whichPtb ][2] );
   //~ printf("num: %d\n",mod.numPnts_);

   // Deterministic
   //~ ++this->currPtb_;
   //~ if( this->currPtb_ == this->numOpt_ ) { ++this->currPnt_; this->currPtb_ = 0; }

   //~ printf("\n");
   return whichPnt;
}

//~ void bMC::pymol( FILE* op, char name[], char color[] ) {
   //~ printf("here\n");
   //~ uint size = this->conf_->size() < 5 ? this->conf_->size() : 5;
   //~ uint len = strlen( name ); len += 5;
   
   //~ // Setup color gradient
   //~ char* gradient[size];
   //~ float red = 1.0;
   //~ float grn = 0.0;
   //~ float inc = 1.0 / size;
   //~ for( uint i=0; i < size; ++i ) {
      //~ gradient[i] = new char[ len ];
      //~ memset( gradient[i], '\0', len );
      //~ sprintf( gradient[i], "%s%04u", name, i );
      //~ fprintf( op, "cmd.set_color(\"%s\", [1.00, %.4f, 0.0])\n",gradient[i],grn );
      //~ grn += inc;
   //~ }
   
   //~ char nameFull[32];
   //~ char scr[7];
   //~ memset( nameFull, '\0', 32 );
   //~ bConf* it = this->conf_->beg();
   //~ bConf* bk = this->conf_->end();
   //~ for( uint i=0; it != bk && i < size; ++i, it=it->next() ) {
      //~ sprintf( nameFull, "%s_conf%02u", name, i );
      //~ sprintf( scr, "%.2f", it->sc_ );
      //~ it->pp_->setToNormalSpace();
      //~ it->pp_->pymolConnectedPseudoatoms( op, nameFull, gradient[i], 1, scr );
   //~ }
//~ }

/******************************************************* bMC_MetTable */

bMC_MetTable::bMC_MetTable() :
   delt_(NULL), prob_(NULL)
{}

bMC_MetTable::~bMC_MetTable() {
   if( this->delt_ != NULL ) { delete [] this->delt_; this->delt_ = NULL; }
   if( this->prob_ != NULL ) { delete [] this->prob_; this->prob_ = NULL; }
}

float bMC_MetTable::operator[]( float find ) {
   if( this->delt_ == NULL || this->prob_ == NULL ) { throw "[bMC::MT] Please identify a Metropolis table"; }
   
   ushort end = this->size_;
   ushort med = end; med >>= 1;
   ushort beg = 0;
   
   ushort at = end;
   bool found = false;
   
   for( uint i=0; i < 4 && (beg != med || end != med); ++i ) {
      if( this->delt_[med] > find ) {
         beg = med; med >>= 1; med += beg;
      }
      else if( this->delt_[med] < find ) {
         end = med; med >>= 1;
      }
      else { at = med; i = 2; found = true;}
   }
   
   for( ushort i=beg; i < end && !found; ++i ) { if( find <= this->delt_[i] ) { at = i; found = true; } }
   //~ if( at > 0 ) { --at; }
   return this->prob_[at];
}

void bMC_MetTable::read( char file[] ) {
   // open the file
   ifstream ip;
   ip.open(file, ifstream::in);
   if(!ip) { return; }

   // Read in data
   deque<string> data;
   string bffr;
   uint i = 0;
   for( i = 0; getline( ip, bffr ); ++i ) { data.push_back(bffr); }
   ip.close();

   // Resize data structures
   if( this->delt_ != NULL ) { delete [] this->delt_; this->delt_ = NULL; }
   if( this->prob_ != NULL ) { delete [] this->prob_; this->prob_ = NULL; }
   this->delt_ = new float[i];
   this->prob_ = new float[i];

   // loop through each line and save the coordinates
   for( i = 0; i < data.size(); ++i ) {
      istringstream ss(data[i]);
      ss >> this->delt_[i];
      ss >> this->prob_[i];
   }
   return;
}