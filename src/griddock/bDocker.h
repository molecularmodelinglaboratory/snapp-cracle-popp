#ifndef BDOCKER_H
#define BDOCKER_H

//~ #include <valarray>
#include <vector>
#include <deque>
#include "bGrid.h"
#include "bPoints.h"
#include "bDelTess.h"
#include "bMC.h"
#include "MersenneTwister.h"
#include "bList.h"

namespace bStd { class bDocker; };

class bStd::bDocker {
   typedef unsigned int uint;

public:
   bDocker();
   ~bDocker();
   
   // Command Line Arguements
   const char* usage();
   bool cla( int, char** );

   // Peptides
   bool findOrgPep();
   bool findPepSeq();
   bool findPepConf();
   void initializePeptides( int, int );
   bool scorePeptides();

   // Random Point generators

   bool generateRandomGridPoint( bGrid&, bGrid&, bPoints&, int );
   bool generateDisplacedRandomGridPoint( bGrid&, bGrid&, bPoints&, int, float[] );
   bool generateRandomChain( bGrid&, bPoints&, float[], int=1, int=-1 );
   bool generateChains( bGrid&, std::vector<bPoints>&, float[], int=-1, int=-1 );

   void runMC( bList*, FILE*, char[] );

   // PyMol
   void printPyMol();
   static void pymolHeader( FILE*, int=0 );
   static void pymolPdb( FILE*, char*, char* =(NULL) );
   void pymolOrigin( FILE* );
   void pymolTrueOrigin( FILE*, int );
   void pymolGridLines( FILE* );
   void pymolMovie( FILE* );
   void pymolOptions( FILE* );
   void pymolBox( FILE*, float[], float[], char[], char[] );
   void pymolScores( FILE* );

   // Temporary
   void doItAll();
      
protected:
   // Counters
   int numConf_;
   static uint _mcitr_;
   static uint _mctot_;

   // Sources
   bGrid grd_;
   bPoints prt_;
   bDelTess prtDT_;
   //~ std::deque<bDelTess> tmpDT_;
   std::deque<bPoints*> pep_;
   //~ std::deque<bPoints*> pepNear_;
   //~ std::deque<bDelTess> pepTess_;
   //~ std::deque<bGrid> pepGrd_;
   //~ std::deque<char*> oPepSeq_;
   std::deque<bPoints*> oPep_;
   std::deque<bPoints*> oPepNear_;
   std::deque<bDelTess> oPepTess_;
   //~ std::deque<float> pepScr_;
   float*  pepScr_;
   int*    pepOrd_;
   double* pepRMSD_;
   bPoints* test_;
   std::deque<bMC> mc_;
   ushort numPepIter_;

   // Flags
   bool makeMovie_;
   bool haveGrd_;
   bool havePrt_;
   bool havePep_;
   bool haveDT_;
   bool havePSeq_;
   bool havePScr_;
   bool haveMC_;

   //~ bStamp ext_;

private:

};




#endif
