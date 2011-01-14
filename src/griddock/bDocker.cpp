
#include <omp.h>
double omp_get_wtime(void);

#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//~ #include <ctime>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "bDocker.h"
#include "bAA.h"
#include "bGrid.h"
#include "bPoints.h"
#include "bDelTess.h"
#include "bCentroid.h"
#include "bSort.h"
#include "bSys.h"
#include "bCalc.h"
#include "bMC.h"
#include "bList.h"

using namespace std;
using namespace bStd;

#define PI 3.14159265

uint bDocker::_mcitr_ = 0;
uint bDocker::_mctot_ = 100;

const char* bDocker::usage() {
   string usage = "[Usage] griddock.e -p [pdb code] [options]\n";
   usage += "\t-d\t(data) pdb code, i.e., 1AWQ_alone or 1IFR\n";
   usage += "\t-p\t(path) default: ../../grid/\n\t-r\t(resolution) grid spacing; default: 0.5\u212B\n";
   usage += "\t-f\t(fit) grid proximity to protein; default: 4\u212B\n";
   usage += "\t-t\t(thickness) grid depth; default: 6\u212B\n";
   return usage.c_str();
}

/* */
/* Constructors and Destructors */
bDocker::bDocker() : makeMovie_(false) {
   // Set Flags
   haveGrd_ = false;
   havePrt_ = false;
   havePep_ = false;
   haveDT_ = false;
   havePSeq_ = false;
   havePScr_ = false;
   
   this->pepScr_ = NULL;
   this->pepOrd_ = NULL;

   //~ rndNum_.seed((unsigned int)time(NULL)); // seeds the Mersenne Twister RNG
   //~ pep_.clear();
   this->numConf_ = 0;
   this->numPepIter_ = 10;
}
bDocker::~bDocker() {
   
   //~ this->pep_.clear();
   //~ this->pepNear_.clear();
   //~ this->oPepSeq_.clear();
   //~ subs = omp_get_wtime();
   //~ this->pepTess_.clear();
   //~ sube = omp_get_wtime();
   //~ printf("~bDelTess
   //~ this->pepGrd_.clear();
   
   for( uint i=0; i < this->oPep_.size(); ++i ) {
      delete this->oPep_[i];
      this->oPep_[i] = NULL;
      
      delete this->oPepNear_[i];
      this->oPepNear_[i] = NULL;
   }
   this->oPep_.clear();
   this->oPepNear_.clear();
   
   for( uint i=0; i < this->pep_.size(); ++i ) {
      delete this->pep_[i];
      this->pep_[i] = NULL;
      
      //~ delete this->pepNear_[i];
      //~ this->pepNear_[i] = NULL;
   }
   this->pep_.clear();
   //~ this->pepNear_.clear();
   
   delete [] this->pepScr_;
   this->pepScr_ = NULL;
   
   delete [] this->pepOrd_;
   this->pepOrd_ = NULL;
   
   delete [] this->pepRMSD_;
   this->pepRMSD_ = NULL;
   
}

/* */
/* Handle Command Line Arguements */
bool bDocker::cla(int numArg, char** argArray) {
   
   // loop through all given parameters (i=1 to skip command)
   char* src = new char [64]; memset( src, '\0', 64 );
   char* pth = new char [64]; memset( pth, '\0', 64 );
   float param[3];
   for(int i=1; i<numArg; ++i) {
      // catches a filename...not the best
      // but this allows for no flags at all

      // check for a flag
      if(argArray[i][0] == '-') {
         switch(argArray[i][1]) {
            case 'd':
               strcpy( src, argArray[++i] );
               break;
            case 'p':
               strcpy( pth, argArray[++i] );
               break;
            case 'r':
               param[0] = atof(argArray[++i]);
               break;
            case 'f':
               param[1] = atof(argArray[++i]);
               break;
            case 't':
               param[2] = atof(argArray[++i]);
               break;
            case 'i':
               this->numPepIter_ = atoi( argArray[++i] );
               break;
            default:
               (void)printf("\tUnknown parameter: '%s'\n\n",argArray[i]);
               throw usage();
               break;
         }
      }
      else if( i == 1 ) {
         strcpy( src, argArray[i] );
      }
      else {
         throw usage();
         //~ (void)printf("Usage: ./main.e [-i] <prot file> [<tet file>] [options]\n");
         //~ (void)printf("\tr: resolution [1.0,0.5]\n");
         //~ (void)printf("\tf: fit; defines exclusion radius\n");
         //~ (void)printf("\tt: thickness; defines thickness of grid\n\n");
         exit(1);
      }
   }

   try{ 
      if( pth != NULL ) { this->prt_.setSrc( src, pth, true ); }
      else { this->prt_.setSrc( src ); }
      this->grd_.setParam( param );
   }
   catch( const char* e ) { printf("%s\n",e); }
   
   delete [] src; src = NULL;
   delete [] pth; pth = NULL;

   return 1;
}


/* */
/* Do It All -- temporary Driver */
void bDocker::doItAll() {
   double start, end;
   start = omp_get_wtime();

/* Grid */  this->grd_.setPointObject( this->prt_ );

/* Prot */  printf("[bDocker] Reading protein...\n");
            if( !(this->prt_.haveFile()) ) {
               bCentroid::pdb2Centroids( this->prt_.pntBase_, this->prt_.pntPath_ );
            }
            this->havePrt_ = prt_.readPoints();

/* Tess */  printf("[bDocker] Tessellating...\n");
            bPoints** src = new bPoints*[1];
            src[0] = & this->prt_;
            this->haveDT_ = this->prtDT_.tessellate_full( src, 1 );
            //~ catch( const char* e ) { printf("[bDelTess] %s.\n", e); }

/* Prot */  printf("[bDocker] Getting tetrahedral centroids...\n");
            bDelTess partial;
            partial.addSrc( this->prt_ );
            partial.prep();
            partial.tessellate();
            partial.verify();
            partial.removeExcess();
            partial.trim( 10.0 );
            partial.onion();
            partial.removeExcess();
            partial.findEdges(); 
            if( bCentroid::dt2Centroids( partial, this->prt_.tets_ ) == 0 ) {
               this->prt_.setTets();
               printf("tet file: %s\n", this->prt_.tets_->pntBase_);
               if( this->prt_.tets_->haveFile() ) {
                  printf("havefile\n");
                  bCentroid::dt2Centroids( partial, this->prt_.pntBase_, this->prt_.pntPath_ );
               }
               this->prt_.haveTets_ = this->prt_.tets_->readPoints();
            }
            else { this->prt_.haveTets_ = true; }
            this->prt_.prep();

/* Grid */  this->prt_.setToDockingSpace(); // used to be w/ prot
            this->grd_.initializeGrid();
            this->haveGrd_ = this->grd_.stampPoints();

/* Seq */   (void)printf("[bDocker] Reading in peptide sequences...\n");
            this->havePSeq_ = this->findOrgPep();
//~ int z = 0;

/* Pep */   (void)printf("[bDocker] Generating random peptides...\n");
            //~ printf("[%d]\n",z); ++z; // 0
            float coordRestr[3] = { 0, 120, 0 };
            int howMany = this->numPepIter_;
            int howLong = strlen( this->oPep_[0]->aaSeq_ );
            //~ printf("how long: %d\n",howLong);
            this->pepScr_ = new float[ howMany ];
            this->pepOrd_ = new int[ howMany ];
            this->pepRMSD_ = new double[ howMany ];
            int count = 0;
            int tenPercent = howMany / 10;
            int percentDone = tenPercent;
            //~ printf("[%d]\n",z); ++z; // 0
            for( int i=0; i < howMany; ++i ) {
               ++count;
               bPoints* pep = new bPoints;
               bPoints* pepN = new bPoints;
               bDelTess* pepT = new bDelTess;
               if( i == percentDone ) {
                  printf("|"); fflush(stdout);
                  percentDone += tenPercent;
               }
               
               // initialize peptide
               pep->setSeq( this->oPep_[0]->aaSeq_ );
               pep->sizeAndSpace( this->prt_, howLong );

               // Generate random peptide
               int startPt = 0;
               do{ startPt = bCalc::getRandomNumber(howLong - 1) + 1; } while( startPt >= howLong );// random start point (non-index)
               bool haveChain = generateRandomChain( this->grd_, *pep, coordRestr, startPt );

               // Find Nearby Points
               bool haveNearby = pep->findNearbyPoints( this->prt_, *pepN );

               // Tessellate
               if( haveChain && haveNearby ) {
                  this->pep_.push_back( pep );
                  //~ this->pepNear_.push_back( pepN );

                  bPoints* bp[2] = { pepN, pep };
                  pepT->addSrc( bp, 2 );
                  pepT->prep();
                  bool haveDT = false;
                  haveDT = pepT->tessellate();
                  haveDT &= pepT->verify();
                  if( !haveDT ) { pepT->retry(); }
                  
                  if( haveDT ) {
                     pepT->focus( 2 );
                     pepT->trim( 11 );
                  }
                  //~ else { printf("still no good!\n"); }

                  if( haveDT && haveChain && (pepT->score() > 5) ) {

                     ++this->numConf_;
                     this->pepScr_[i] = pepT->score_;
                     this->pepOrd_[i] = i;
                     this->pepRMSD_[i] = bPoints::rmsd( *pep, *(this->oPep_[0]) );
                     //~ pepT->removeExcess();
                     //~ pepT->findEdges();
                     //~ pepT->findTypes();
                     //~ this->pepTess_.push_back( pepT );
                     //~ this->pepTess_.back().removeExcess();
                     //~ this->pepTess_.back().findEdges();
                     //~ this->pepTess_.back().findTypes();
                  }
                  else {
                     delete this->pep_.back();
                     this->pep_.pop_back();
                     //~ this->pepNear_.pop_back();
                     if( ! haveDT || ! haveChain ) { printf("*"); } //printf("== bad %d ==\n",i);
                     --i;
                  }
               }
               else { --i; }
               
               delete pepT;
               pepT = NULL;
               
               delete pepN;
               pepN = NULL;
            }
            this->havePep_ = !(this->pep_.empty());
            printf(" iterations: %d\n",count);
            this->havePScr_ = true;
            bSort::qsort( this->pepScr_, this->pepOrd_, howMany, true );
            //~ for( int i=0; i < howMany; ++i ) { printf("[%d] %4d %.2f\n", i, this->pepOrd_[i], this->pepScr_[this->pepOrd_[i]]); }

/* MMC *   printf("[bDocker] Performing Metropolis Monte Carlo...\n");
            _mcitr_ = 0;
            uint numMC = this->pep_.size();
            numMC = 5;
            //~ uint count = 0;
            bList pep;
            for( uint k=0; k < numMC; ++k ) {
               bConf* conf = new bConf( *(this->pep_[this->pepOrd_[k] ]), this->pepScr_[this->pepOrd_[k]] );
               //~ bList* list = new bList;
               
               //~ bMC mc;
               //~ mc.setup( this->grd_, this->prt_, *conf );
               //~ bool haveMet = false;
               //~ haveMet = mc.metropolis();
            
               //~ if( haveMet ) {
                  pep.push_back( conf );
               //~ }
            }
            
            printf("pep: [ %u, %u ]\n", pep.size(), pep.cap());
            pep.print();
            char name[128] = "plist";
            FILE *op = fopen("d2.py", "w");
            bList::pymolGradient( op );
            
            runMC( &pep, op, name );
            
            fclose(op);
            
            //~ char name[32] = "plist";
            //~ FILE *op = fopen("d2.py", "w");
            //~ pep.pymolGradient( op );
            //~ pep.pymol( op, name );
            //~ fclose(op);
            //~ printf("size: %d | %d | %d\n", dtstock[i]->src_[0]->numPnts_,dtstock[i]->src_[1]->numPnts_,dtstock[i]->src_[2]->numPnts_);

   /* Done */

   end = omp_get_wtime();
   this->printPyMol();
   printf("Work took %f sec. time.\n", end-start);

}

bool bDocker::findOrgPep() {
   char path[64];
   strcpy( path, this->prt_.pntPath_ );
   strcat( path, "orgPep/");

   // Open directory
   DIR* dp;
   struct dirent* dir;
   if( (dp = opendir( path )) == NULL ) {
      throw "[bDocker] Unable to open pepSeq directory.";
      return false;
   }

   // Get file bases
   deque<string> files;
   char prev[64]; memset( prev, '\0', 64 );
   bool areEqual = false;
   while( (dir = readdir(dp)) != NULL ) {
      if( dir->d_name[0] == '.' ) { continue; }
      
      // Determine if it's the same data set
      areEqual = false;
      uint size = strlen( dir->d_name );
      if( size == strlen( prev ) ) {
         if( memcmp( dir->d_name, prev, size - 3 ) == 0 ) {
            areEqual = true;
         }
      }
      
      // Save the data set
      if( !areEqual ) {
         memset( prev, '\0', 64 );
         memmove( prev, dir->d_name, size );
         memset( dir->d_name + size - 4, '\0', 4 );
         files.push_back( string(dir->d_name) );
      }
   }
   closedir(dp);

   // Opening files
   char fname[96];
   for( uint i=0; i < files.size(); ++i ) {
      memset( fname, '\0', 96 );
      memmove( fname, files[i].c_str(), files[i].size() );
      
      // Check for PDB file
      memmove( fname + files[i].size(), ".pdb", 4 );

      bool exists = bSys::fileExists( fname, path );
      if( exists ) { exists = bCentroid::pdb2Centroids( fname, path ); }
      
      if( exists ) {
         bPoints* pep = new bPoints;
         pep->setSrc( files[i].c_str(), path );
         pep->readPoints();
         pep->sizeAndSpace( this->prt_ );
         this->oPep_.push_back( pep );
         //~ this->oPepSeq_.push_back( pep->aaSeq_ );
         
         bPoints* pepN = new bPoints;
         if( pep->findNearbyPoints( this->prt_, *pepN ) ) { this->oPepNear_.push_back( pepN ); }

         bDelTess tess;
         bPoints* bp[2] = { pepN, pep };
         tess.addSrc( bp, 2 );
         tess.prep();
         tess.tessellate();
         tess.verify();
         tess.removeExcess();
         tess.focus( bp[1] );
         tess.trim( 12 );
         tess.removeExcess();
         tess.findEdges();
         tess.findTypes();
         tess.score();
         this->oPepTess_.push_back( tess );
         
         printf("original size:  %lu\n", tess.simplex_.size() );
         printf("original score: %.2f\n", tess.score() );
      }
   }
   return true;
}

/* Set Peptide Size and Space */
void bDocker::initializePeptides( int nC, int nP ) {
   if( this->pep_.size() < (unsigned int)nC) { this->pep_.resize(nC); }
   for( int i=0; i < nC; ++i ) { this->pep_[i]->sizeAndSpace( this->prt_, nP ); }
   return;
}

/* generate chain */
bool bDocker::generateRandomChain(bGrid &grd, bPoints &pep, float sRestr[], int startPt, int size) {
   if(size == -1) { size = pep.capPnts_; } // default
   
   bool validPeptide = false;
   int count = 0;
   bGrid pepGrd;
   pepGrd.overlay( grd );
   //~ pepGrd.fit_ = 3;
   //~ pepGrd.initializeStamps();
   pepGrd.setPointObject( pep );
   //~ this->pepGrd_.setPointObject(p);
   //~ this->pepGrd_.clearGrid();
   while(!validPeptide && count < 5000) {
      validPeptide = true; // set...this will change if the points fail
      pep.clear(); // clear out the chain each time
      pepGrd.clear();
      bool yup = false;
      int i=0;
      do {
         yup = this->generateRandomGridPoint( grd, pepGrd, pep, startPt ); // set initial random point
      } while( !yup );
      // these two loops cover the entire chain
      // -- if point generation ever fails (i.e. not on grid w/in ~numIter)
      //    then we start over again
      for(int i=startPt+1; i<=size && validPeptide; ++i) {
         //~ printf("pt: %d\n",i);
         validPeptide = this->generateDisplacedRandomGridPoint( grd, pepGrd, pep, i, sRestr );
      }
      for(int i=startPt-1; i>0 && validPeptide; --i) {
         //~ printf("pt: %d\n",i);
         validPeptide = this->generateDisplacedRandomGridPoint( grd, pepGrd, pep, i, sRestr );
      }
      ++count;
   }
   //~ this->tempGrd_.push_back( this->pepGrd_ );

   //~ if( validPeptide ) { this->pepGrd_.push_back( pepGrd ); }
   return validPeptide;
}

/* Generate a random point on the grid */
bool bDocker::generateRandomGridPoint( bGrid &grd, bGrid &pgrd, bPoints &pep, int whichPoint ) {
   bool isValidPoint = false; // flag

   // randomly generate point until we get a hit
   int cnt = 0;
   float newPoint[3] = { 0, 0, 0 };
   while(!isValidPoint && cnt < 5000) {
      for(int i=0; i<3; ++i) { // randomly seed new coordinates
         newPoint[i] = bCalc::getRandomNumber( grd.max_[i] - grd.buffer_ );
         newPoint[i] += grd.buffer_;
      }

      if( grd.isaGridPoint( newPoint ) ) { // check for validity
         isValidPoint = true;
         pep.addPointAtPos( newPoint, whichPoint );
         pgrd._stampPoint( pgrd.fStmp_, newPoint );
      }
      ++cnt;
   }

   return isValidPoint;
}
/* Generate a Directed Random Point on the Grid
   -- note that sCoord represents [r,phi,psi]
   -- if only r!=0, assume any value of phi and psi work; use one previous point
   -- if only r!=0 && phi!=0, we assume any psi works; use two previous points
   -- not too sure what to do if given phi
*/
bool bDocker::generateDisplacedRandomGridPoint( bGrid &grd, bGrid &pgrd, bPoints &pep, int whichPoint, float sCoord[] ) {

   int capPoints = pep.capPnts_; // count number of points
   int newX = (whichPoint-1)*3; // get the point index

   // identify which direction previous points are from
   int direction = 1;

   if(whichPoint == capPoints) {
      direction = -1; // last point -- go backwards
   }
   else if(whichPoint == 1) {
      direction = 1; // first point - go forwards
   }
   else{
      int tmpIndex = newX - 3;
      int sum = pep.pnts_[ tmpIndex ];
      sum += pep.pnts_[ ++tmpIndex ];
      sum += pep.pnts_[ ++tmpIndex ];
      if( sum != 0.0 ) {
         direction = -1; // points behind -- go backwards
      }
      else {
         direction = 1; // points ahead -- go forwards
      }
      //~ else {
         //~ return generateRandomGridPoint(grd,pgrd,pnt,whichPoint); // no points...should not be here
         //~ // note: this assumes that we won't ever have a point at the origin
      //~ }
   }

   int orgX = newX + direction*3; // declare adjacent point indices
   int mrrX = newX + direction*6;

   bool have2ndPnt = false; // check for second adjacent point
   if( mrrX >= 0 && mrrX < (capPoints * 3) ) {
      float sum = pep.pnts_[mrrX] + pep.pnts_[mrrX + 1] + pep.pnts_[mrrX + 2];
      if( sum > 0.0 ) { have2ndPnt = true; }
   }

   float orgPnt[3] = { pep.pnts_[orgX], pep.pnts_[orgX + 1], pep.pnts_[orgX + 2] };
   float mrrPnt[3] = { 0, 0, 0 };
   float newPnt[3] = { 0, 0, 0 };

   if( have2ndPnt ) {
      mrrPnt[0] = pep.pnts_[mrrX]; // save point
      mrrPnt[1] = pep.pnts_[mrrX + 1];
      mrrPnt[2] = pep.pnts_[mrrX + 2];

      mrrPnt[0] -= orgPnt[0]; // translate to the origin
      mrrPnt[1] -= orgPnt[1];
      mrrPnt[2] -= orgPnt[2];
      bPoints::c2s( mrrPnt );
   }

   // set r restrictions
   if(!sCoord[0]) { sCoord[0] = ( grd.fit_ * 1.2); }

   // generate random angles until we get a valid point
   bool isValidPoint = false;
   int counter = 0;
   while(!isValidPoint && counter < 500) {

      //~ printf("cnt: %d\n", counter);
      // set r
      newPnt[0] = sCoord[0];

      // set theta [-90,90]
      // -- we assume the adjacent point is at the origin
      // -- a second point is directly below (if we have it)
      // -- assumes the restriction is obtuse
      if(have2ndPnt && sCoord[1] != 0) { // restrict theta
         newPnt[1] = sCoord[1] - 90;
      }
      else { newPnt[1] = (float)(90 - bCalc::getRandomNumber(180)); }
         // we don't have a second point, so randomly generate the spherical coordinates

      // set phi [-180,180]
      if(have2ndPnt && sCoord[2] != 0) { newPnt[2] = sCoord[2]; } // restrict phi
      else { newPnt[2] = (float)(180 - bCalc::getRandomNumber(360)); }
         // we don't have a second point, so randomly generate the spherical coordinates

      bPoints::s2c( newPnt ); // convert spherical to cartesian

      if(have2ndPnt) { // rotate the coordinates if we need to (i.e. if we have a second point)
         bPoints::rotateTheta( newPnt, (90 + mrrPnt[1]) );
         bPoints::rotatePhi( newPnt, (-mrrPnt[2]) );
      }

      newPnt[0] += orgPnt[0];
      newPnt[1] += orgPnt[1];
      newPnt[2] += orgPnt[2];

      newPnt[0] = (int)newPnt[0];
      newPnt[1] = (int)newPnt[1];
      newPnt[2] = (int)newPnt[2];

      if( grd.isaGridPoint( newPnt ) && !pgrd.isaGridPoint(newPnt) ){ // check if the point is on the grid
         isValidPoint = true;
         pep.addPointAtPos( newPnt, whichPoint );
         pgrd._stampPoint( pgrd.fStmp_, newPnt );
      }
      ++counter; // add one to bad point count
      //~ printf("counter: [%d] %4d\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\n",whichPoint,counter, newPnt[0], newPnt[1], newPnt[2], mrrPnt[1], mrrPnt[2]);
   }

   return isValidPoint;
}


/* Uses the Mersenne Twister algorithm to get a random number
   -- uses the nextInt method from Java to constrain w/in [0-range]
*/



void bDocker::runMC( bList* init, FILE* op, char name[] ) {
   if( _mcitr_ > _mctot_ ) { return; }
   else { ++_mcitr_; }

   deque< bList* > todo;
   todo.push_back( init );

   //~ for( uint i=0; i < 100 && ! todo.empty(); ++i ) {
   uint k=0;
   while( ++k < 10 && ! todo.empty() ) {
      bMC mc;
      bList* set = todo[0];
      todo.pop_front();
      printf("[%u]\n",k);
      set->print();
      for( uint i=0; i < set->size(); ++i ) {
         //~ printf("[i, %u]\n",i);(*set)[i].print();
         mc.setup( this->grd_, this->prt_, (*set)[i] );
         if( mc.metropolis() ) {
            sprintf(name,"%s.%u",name,i);
            //~ (*set)[i].chld()->print();
            (*set)[i].chld()->trim();
            //~ (*set)[i].chld()->print();
            (*set)[i].chld()->pymol( op, name );
            todo.push_back( (*set)[i].chld() );
            //~ runMC( (*set)[i].chld(), op, name );
         }
      }
      set = NULL;
      printf("todo: %u, %lu\n",todo.empty(), todo.size());
   }
   return;
}

/* */
/* Create PyMol Python Script */
void bDocker::printPyMol() {

   // Reset to normal space
   bPoints::setToNormalSpace( this->prt_ );
   bPoints::setToNormalSpace( this->pep_, this->numConf_ );
   //~ bPoints::setToNormalSpace( this->pepNear_, this->pepNear_.size() );
   bPoints::setToNormalSpace( this->oPep_, this->oPep_.size() );
   bPoints::setToNormalSpace( this->oPepNear_, this->oPepNear_.size() );

   // colors
   char colorRed[4] = "red";
   char colorPnk[5] = "pink";
   char colorRub[5] = "ruby";
   char colorMar[7] = "marine";
   //~ char colorBlu[5] = "blue";
   //~ char colorYOr[13] = "yelloworange";
   //~ char colorPYw[11] = "paleyellow";
   //~ char colorPGr[10] = "palegreen";
   //~ char colorG60[7] = "gray60";

   char colorGre[6] = "green";
   char colorSky[8] = "skyblue";
   char colorPur[7] = "purple";
   //~ char colorHPk[8] = "hotpink";
   char colorOrg[7] = "orange";
   char colorG10[7] = "gray10";

   // create the filename
   char *spdb = new char[64];
   char *sprt = new char[64];
   strcpy( spdb, this->prt_.pntPath_ );
   strcpy( sprt, this->prt_.pntPath_ );
   strcat( spdb, this->prt_.pntBase_ );
   strcat( sprt, this->prt_.pntBase_ );
   strcat( spdb, ".pdb" );
   strcat( sprt, "_grid.py" );
   (void)printf("[bDocker] Writing to: '%s'\n", sprt);
   // open the file and check
   FILE *op = fopen(sprt, "w");
   if(!op) {
      (void)printf("%s did not open!\n", sprt);
   }

   // Write header
   pymolHeader(op,1);

   // Write pdb
   pymolPdb( op, this->prt_.pntBase_, this->prt_.pntPath_ );

   // Write origin
   //~ this->pymolOrigin(op);
   //~ this->pymolGridLines(op);

   /* Variables */   char *name = new char[64];
                     float size = 0.2;

   /* Protein */ if( this->havePrt_ && this->prt_.havePnts_ ) {
      strcpy(name,"proteins");
      bPoints::pymolPoints( op, prt_.pnts_, prt_.numPnts_, name, colorRed, size );
   }

   /* Protein Centroids */ if( this->havePrt_ && this->prt_.haveTets_ ) {
      //~ printf("have em\n");
      size = 0.1;
      strcpy(name,"centroids");
      bPoints::pymolPoints( op, prt_.tets_->pnts_, prt_.tets_->numPnts_, name, colorPnk, size );//,color,size);
   }

   /* Tessellation */ if( this->haveDT_ ) {
      strcpy( name, "tess" );
      this->prtDT_.pymol( op, name, colorRub );
      //~ for( uint i=0; i < tmpDT_.size(); ++i ) {
         //~ sprintf( name, "tmp_%u", i );
         //~ this->tmpDT_[i].pymol( op, name, colorMar );
      //~ }
   }

   /* Grid */ if( this->haveGrd_ ) {
      strcpy(name,"grid");
      this->grd_.pymol3dGrid( op, name, colorMar);
      //~ strcpy(name,"grid2");
      //~ this->grd_.pymol3dGrid( op, name, colorBlu, true );
      (void)fprintf(op,"cmd.disable(\"%s\")\n",name);
      
      //~ grd_.removeInternal( op );
      //~ bPoints test;
      //~ test.sizeAndSpace( prt_ );
      //~ test.isInRes_ = true;
      //~ test.isInPos_ = true;
      //~ for( uint z=0; z < 10; ++z ) {
         //~ uint pnt1[3];
         //~ float fpnt1[3];
         //~ do {
            //~ pnt1[0] = bCalc::getRandomNumber( grd_.length_ ); pnt1[1] = bCalc::getRandomNumber( grd_.height_ ); pnt1[2] = bCalc::getRandomNumber( grd_.depth_ );
            //~ fpnt1[0] = (float)pnt1[0]; fpnt1[1] = (float)pnt1[1]; fpnt1[2] = (float)pnt1[2];
         //~ } while( ! grd_.isaGridPoint(fpnt1) );
         //~ test.addPoint( fpnt1 );
         //~ grd_.countNearby( pnt1,op );
      //~ }
      //~ strcpy( name, "nearby" );
      //~ test.setToNormalSpace();
      //~ test.pymolPseudoatoms( op, name, colorHPk );
   }

   /* OrgPep */ if( this->oPep_.size() > 0 ) {
      for( uint i=0; i < this->oPep_.size(); ++i ) {
         sprintf( name, "pOrg_%02d", i );
         char scr[7];
         sprintf( scr, "%.2f", this->oPepTess_[i].score_ );
         this->oPep_[i]->pymolConnectedPseudoatoms( op, name, colorOrg, 1, scr );
         sprintf( name, "pOrg_%02d_dt", i );
         this->oPepTess_[i].pymol( op, name, colorOrg );
      }
      
   }
   int num2Print = this->pep_.size() > 200 ? 200 : this->pep_.size();
   /* Scores */ if( this->havePScr_ ) {
      this->pymolScores( op );
   }

   /* Peptides */ if( this->havePep_ ) {
      for( int i=0; i < num2Print; ++i ) {

         // Color
         float percent = (float)i / (float)this->pep_.size();
         char *color;
         if( percent < .05 ) { color = colorGre; }
         else if( percent < .10 ) { color = colorSky; }
         else if( percent < .20 ) { color = colorPur; }
         else if( percent < .30 ) { color = colorG10; }
         else { color = colorOrg; }

         char scr[7];
         sprintf( scr, "%.2f", this->pepScr_[ this->pepOrd_[i] ] );
         (void)sprintf(name,"p%03d", i); //this->pepOrd_[i]);
         this->pep_[ this->pepOrd_[i] ]->pymolConnectedPseudoatoms( op, name, color, 1, scr );
         //~ if( showPepDetail && i < 15 ) {
            //~ try{
            //~ (void)sprintf(name,"p%03u_dt", i); //this->pepOrd_[i]);
            //~ this->pepTess_[ this->pepOrd_[i] ].pymol( op, name, color );
            //~ (void)fprintf(op,"cmd.disable(\"%s\")\n",name);
            //~ }
            //~ catch( const char* e ) { printf("%s\n", e); }
         //~ }
      }
      fprintf(op,"cmd.show(\"sticks\",\"p*\")\n");
      fprintf(op,"cmd.hide(\"sticks\",\"*_dt\")\n");
      fprintf(op,"cmd.show(\"sticks\",\"top*\")\n");
   }

   /* MonteCarlo */ if( ! this->mc_.empty() ) {
      for( uint i=0; i < this->mc_.size(); ++i ) {
         //~ printf("hi\n");
         //~ sprintf(name,"mcTop%d",i);
         //~ char scr[7];
         //~ sprintf( scr, "%.2f", this->mc_[i].conf_.front().dt_->score() );
         //~ this->mc_[i].conf_.front().pp_->pymolConnectedPseudoatoms( op, name, colorOrg, 1, scr );

         //~ sprintf(name,"mc%03u",i);
         //~ this->mc_[i].pymol( op, name, colorBlu );

      }
   }

   // Movie -- Rotate!
   /* PyMol Movie */ if( this->makeMovie_ ) {
      this->pymolMovie(op);
   }

   // Close the file handle
   //~ (void)fprintf(op,"cmd.clip(\"slab\",20)\n");
   //~ (void)fprintf(op,"cmd.hide(\"everything\",\"%s\")\n", this->prt_.pntBase_);
   (void)fprintf( op, "cmd.orient(\"%s\")\n", this->prt_.pntBase_ );
   //~ (void)fprintf(op,"cmd.hide(\"labels\")\n");

   fclose(op);
   
   delete [] spdb;
   delete [] sprt;
   delete [] name;
   
   // Reset to Docking space
   //~ bPoints::setToDockingSpace( this->prt_ );
   //~ bPoints::setToDockingSpace( this->pep_, this->numConf_ );
   //~ bPoints::setToDockingSpace( this->oPep_, this->oPep_.size() );
   //~ bPoints::setToDockingSpace( this->oPepNear_, this->oPepNear_.size() );

   // Run the command
   // ('cuz windows sucks...)
   //system("pymolwin -r \"D:\\My Dropbox\\grid\\1AWQ_alone_grid.py\"");
   return;
}

/* Write PyMol Header */
void bDocker::pymolHeader(FILE* op,int white) {
   (void)fprintf(op,"from pymol import cmd\n");
   (void)fprintf(op,"from pymol.cgo import *\n");
   (void)fprintf(op,"\n");
   if(white) {
      (void)fprintf(op,"cmd.bg_color(\"white\")\n");
   }
   (void)fprintf(op,"\n");
   
   // set parameters
   (void)fprintf(op,"cmd.set(\"auto_show_spheres\", \"on\")\n");
   (void)fprintf(op,"cmd.set(\"sphere_scale\", .25)\n");
   (void)fprintf(op,"cmd.set(\"auto_show_lines\", \"on\")\n");
   (void)fprintf(op,"cmd.set(\"label_position\",(1.5,1.5,1.5))\n");
   (void)fprintf(op,"cmd.set(\"label_size\",8)\n");
   (void)fprintf(op,"\n");

   return;
}

/* Write PyMol PDB */
void bDocker::pymolPdb( FILE* op, char* pdb, char* path ) {
   char* pdbfile;
   int size = 0;
   if( path != NULL ) {
      size = strlen(pdb) + strlen(path) + 5;
      pdbfile = new char[ size ];
      memset( pdbfile, '\0', size );
      memmove( pdbfile, path, strlen(path) );
      memmove( pdbfile + strlen(path), pdb, strlen(pdb) );
   }
   else {
      size = strlen(pdb) + 5;
      pdbfile = new char[ size ];
      memset( pdbfile, '\0', size );
      memmove( pdbfile, pdb, strlen(pdb) );
   }
   if( ! bSys::fileExists( pdbfile ) ) { memmove( pdbfile + strlen(pdbfile), ".pdb", 4 ); }
   
   (void)fprintf(op,"cmd.load(\"%s\")\n", pdbfile);
   (void)fprintf(op,"cmd.color(\"firebrick\", \"%s\")\n", pdb );
   (void)fprintf(op,"cmd.hide(\"everything\", \"%s\")\n", pdb );
   (void)fprintf(op,"cmd.show(\"cartoon\", \"%s\")\n", pdb );
   
   //1AWQ only
   (void)fprintf(op,"cmd.fetch(\"1AWQ\")\n");
   (void)fprintf(op,"cmd.hide(\"everything\", \"1AWQ\")\n");
   (void)fprintf(op,"cmd.show(\"cartoon\", \"1AWQ\")\n");
   (void)fprintf(op,"cmd.color(\"firebrick\", \"1AWQ//A//\")\n");
   (void)fprintf(op,"cmd.color(\"yellow\", \"resi 1054\")\n");
   (void)fprintf(op,"cmd.color(\"yellow\", \"resi 1055\")\n");
   (void)fprintf(op,"cmd.color(\"yellow\", \"resi 1063\")\n");
   (void)fprintf(op,"cmd.color(\"yellow\", \"resi 1072\")\n");
   (void)fprintf(op,"cmd.color(\"yellow\", \"resi 1102\")\n");
   (void)fprintf(op,"cmd.color(\"yellow\", \"resi 1121\")\n");
   (void)fprintf(op,"cmd.color(\"yellow\", \"resi 1126\")\n");

   delete [] pdbfile;
   pdbfile = NULL;
   
   return;
}

/* Write PyMol Origin */
void bDocker::pymolOrigin(FILE* op) {
   float mi[3];
   float ma[3];
   for(int i=0; i<3; ++i) {
      mi[i] = grd_.min_[i]*grd_.res_ + prt_.planeDisplacement_[i];
      ma[i] = grd_.max_[i]*grd_.res_ + prt_.planeDisplacement_[i];
   }
   (void)fprintf(op,"orgn = [\n");
   (void)fprintf(op,"\tLINEWIDTH, 2.0,\n");
   (void)fprintf(op,"\tBEGIN, LINES,\n");
   (void)fprintf(op,"\tCOLOR, 0.2, 0.2, 0.2,\n");
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",mi[0],mi[1],mi[2]);
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",mi[0],mi[1],ma[2]);
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",mi[0],mi[1],mi[2]);
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",mi[0],ma[1],mi[2]);
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",mi[0],mi[1],mi[2]);
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",ma[0],mi[1],mi[2]);
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",ma[0],mi[1],mi[2]);
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",ma[0],ma[1],mi[2]);
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",mi[0],ma[1],mi[2]);
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",ma[0],ma[1],mi[2]);
   (void)fprintf(op,"\tEND\n");
   (void)fprintf(op,"\t]\n");
   (void)fprintf(op,"cmd.load_cgo(orgn,'orgn')\n");
   (void)fprintf(op,"cmd.disable(\"orgn\")\n");

   return;
}
/* Write PyMol Origin */
void bDocker::pymolTrueOrigin(FILE* op, int max) {

   (void)fprintf(op,"orgn = [\n");
   (void)fprintf(op,"\tLINEWIDTH, 2.0,\n");
   (void)fprintf(op,"\tBEGIN, LINES,\n");
   (void)fprintf(op,"\tCOLOR, 0.0, 0.0, 0.0,\n");
      (void)fprintf(op,"\tVERTEX, %d, %d, %d, \n",0,0,0);
      (void)fprintf(op,"\tVERTEX, %d, %d, %d, \n",0,0,max);

   (void)fprintf(op,"\tCOLOR, 0.4, 0.4, 0.4,\n");
      (void)fprintf(op,"\tVERTEX, %d, %d, %d, \n",0,0,0);
      (void)fprintf(op,"\tVERTEX, %d, %d, %d, \n",0,max,0);

      (void)fprintf(op,"\tVERTEX, %d, %d, %d, \n",0,0,0);
      (void)fprintf(op,"\tVERTEX, %d, %d, %d, \n",max,0,0);
   (void)fprintf(op,"\tEND\n");
   (void)fprintf(op,"\t]\n");
   (void)fprintf(op,"cmd.load_cgo(orgn,'orgn')\n");

   return;
}
/* Write PyMol Grid */
void bDocker::pymolGridLines(FILE* op) {
   int mi[3];
   int ma[3];
   for(int i=0; i<3; ++i) {
      mi[i] = grd_.min_[i]*grd_.res_ + prt_.planeDisplacement_[i];
      ma[i] = grd_.max_[i]*grd_.res_ + prt_.planeDisplacement_[i];
   }
   (void)fprintf(op,"gridLines = [\n");
   (void)fprintf(op,"\tLINEWIDTH, 1.5,\n");
   (void)fprintf(op,"\tBEGIN, LINES,\n");
   (void)fprintf(op,"\tCOLOR, 0.8, 0.8, 0.8,\n");
   int increment = 10;//*res_;
   for(int k=mi[1]; k<=ma[1];k=k+increment) {//ma[1]; ++k) {
      for(int i=mi[0]; i<=ma[0];i=i+increment) {//ma[0]; ++i) {
         for(int z=mi[2]; z<=ma[2];z=z+increment) {//ma[2]; ++z) {

            // extend in the z direction
            (void)fprintf(op,"\tVERTEX, %d, %d, %d, \n",i,k,mi[2]);
            (void)fprintf(op,"\tVERTEX, %d, %d, %d, \n",i,k,ma[2]);

            // extend in the x direction
            (void)fprintf(op,"\tVERTEX, %d, %d, %d, \n",mi[0],k,z);
            (void)fprintf(op,"\tVERTEX, %d, %d, %d, \n",ma[0],k,z);

            // extend in the y direction
            (void)fprintf(op,"\tVERTEX, %d, %d, %d, \n",i,mi[1],z);
            (void)fprintf(op,"\tVERTEX, %d, %d, %d, \n",i,ma[1],z);

         }
      }
   }
   //~ (void)fprintf(op,"\tCOLOR, 0.8, 0.8, 0.8,\n");
      //~ (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",mi[0],mi[1],ma[2]);
      //~ (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",mi[0],mi[1],mi[2]);
      //~ (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",mi[0],ma[1],mi[2]);
      //~ (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",mi[0],mi[1],mi[2]);
      //~ (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",ma[0],mi[1],mi[2]);
      //~ (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",ma[0],mi[1],mi[2]);
      //~ (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",ma[0],ma[1],mi[2]);
      //~ (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",mi[0],ma[1],mi[2]);
      //~ (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",ma[0],ma[1],mi[2]);
   (void)fprintf(op,"\tEND\n");
   (void)fprintf(op,"\t]\n");
   (void)fprintf(op,"cmd.load_cgo(gridLines,'gridLines')\n");
   (void)fprintf(op,"cmd.disable(\"gridLines\")\n");

   return;

}

/* Write PyMol Movie */
void bDocker::pymolMovie(FILE* op) {
   // -- set up the frames
   (void)fprintf(op,"cmd.mclear()\n");

   // Rotate and zoom
   // Correct zoom depending on the grid size -- should always go just inside the grid
   (void)fprintf(op,"cmd.mset(\"1 x360\")\n");
   float zoom = ((grd_.max_[2] - grd_.min_[2])/2*11);
   zoom /= 400;

   // start out rotating
   (void)fprintf(op,"for i in range(0,60):\n");
   (void)fprintf(op,"\tcmd.mdo(i,'turn x,1; turn y,1; turn z,0')\n");

   // continue rotating; zoom in
   (void)fprintf(op,"for i in range(60,180):\n");
   (void)fprintf(op,"\tcmd.mdo(i,'turn x,1; turn y,1; turn z,0; move z,%.2f')\n",zoom);

   // continue rotating; zoom out
   (void)fprintf(op,"for i in range(180,300):\n");
   (void)fprintf(op,"\tcmd.mdo(i,'turn x,1; turn y,1; turn z,0; move z,-%.2f')\n",zoom);

   // continue rotating; no zoom
   (void)fprintf(op,"for i in range(300,360):\n");
   (void)fprintf(op,"\tcmd.mdo(i,'turn x,1; turn y,1; turn z,0')\n");
   
   // Rotate Y-axis
   //~ (void)fprintf(op,"cmd.mset(\"1 x180\")\n");
   //~ (void)fprintf(op,"for i in range(0,180):\n");
   //~ (void)fprintf(op,"\tcmd.mdo(i,'turn y,2')\n");

   return;
}

/* Miscellaneous Options */
void bDocker::pymolOptions(FILE *op) {
   //~ (void)fprintf(op,"cmd.clip(\"slab\",100)\n"); 
   return;
}

void bDocker::pymolScores( FILE* op ) {
   char** color = new char*[5];
   color[0] = new char[6]; memmove( color[0], "green\0", 6);
   color[1] = new char[8]; memmove( color[1], "skyblue\0", 8);
   color[2] = new char[7]; memmove( color[2], "purple\0", 7);
   color[3] = new char[8]; memmove( color[3], "hotpink\0", 8);
   color[4] = new char[7]; memmove( color[4], "orange\0", 7);
   char name[8];

   memmove( name, "top10\0", 6 );
   int topN[5] = { 0, 10, 25, 50, 100 };
   register uint cnt = 0;
   register int size = 0;
   register int pos = 0;
   for( int k=1; k < 5; ++k ) {
      size = topN[k] - topN[ k - 1 ];
      sprintf( name, "top%d", topN[k] );
      for( int i=0; i < size; ++i ) {
         if( cnt >= this->pep_.size() ) { i = size; k = 5; continue; }
         pos = this->pepOrd_[cnt];
         if( k < 3 ) printf("[%d] %.2f, %.2f\n", cnt, this->pepScr_[pos], this->pepRMSD_[pos]);

         char label[16];
         sprintf( label, "%.2f, %.2f", this->pepScr_[pos], this->pepRMSD_[pos] );
         this->pep_[pos]->pymolConnectedPseudoatoms( op, name, color[k-1], i, label );
         ++cnt;
      }
   }

   delete [] color[0]; color[0] = NULL;
   delete [] color[1]; color[1] = NULL;
   delete [] color[2]; color[2] = NULL;
   delete [] color[3]; color[3] = NULL;
   delete [] color[4]; color[4] = NULL;
   delete [] color;
   color = NULL;
   return;
}



void bDocker::pymolBox(FILE *op, float* min, float *max, char name[], char color[]) {
   (void)fprintf(op,"%s = [\n", name);
   (void)fprintf(op,"\tLINEWIDTH, 1.0,\n");
   (void)fprintf(op,"\tBEGIN, LINES,\n");

      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",min[0],min[1],min[2]); // A - E
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",min[0],min[1],max[2]);

      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",max[0],min[1],min[2]); // B - F
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",max[0],min[1],max[2]);

      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",max[0],max[1],min[2]); // C - G
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",max[0],max[1],max[2]);

      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",min[0],max[1],min[2]); // D - H
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",min[0],max[1],max[2]);


      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",min[0],min[1],min[2]); // A - B
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",max[0],min[1],min[2]);

      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",max[0],min[1],min[2]); // B - C
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",max[0],max[1],min[2]);

      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",max[0],max[1],min[2]); // C - D
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",min[0],max[1],min[2]);
      
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",min[0],max[1],min[2]); // D - A
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",min[0],min[1],min[2]);


      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",min[0],min[1],max[2]); // E - F
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",max[0],min[1],max[2]);

      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",max[0],min[1],max[2]); // F - G
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",max[0],max[1],max[2]);

      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",max[0],max[1],max[2]); // G - H
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",min[0],max[1],max[2]);

      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",min[0],max[1],max[2]); // H - E
      (void)fprintf(op,"\tVERTEX, %.2f, %.2f, %.2f, \n",min[0],min[1],max[2]);

   (void)fprintf(op,"\tEND\n");
   (void)fprintf(op,"\t]\n");
   (void)fprintf(op,"cmd.load_cgo(%s,'%s')\n", name, name);
   (void)fprintf(op,"cmd.color(\"%s\", \"%s\")\n", color, name);
   (void)fprintf(op,"cmd.disable(\"orgn\")\n");
   return;
}
