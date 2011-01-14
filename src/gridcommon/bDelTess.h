#ifndef BDELTESS_H
#define BDELTESS_H

#include <deque>
#include <stdio.h>

#include "bPoints.h"
#include "bSimplex.h"
#include "bHex.h"
#include "MersenneTwister.h"

namespace bStd { struct ptData; class bDelTess; };

struct bStd::ptData {
   unsigned short pos_;
   unsigned short chn_;

   ptData() {
      pos_ = 0;
      chn_ = 0;
   }
   ptData( unsigned short p, unsigned short c ) {
      pos_ = p;
      chn_ = c;
   }
   ptData( const ptData &p ) {
      pos_ = p.pos_;
      chn_ = p.chn_;
   }
};

class bStd::bDelTess : virtual bStd::bVBase {
   typedef std::deque<float> def;
   typedef std::deque<short> des;
   typedef std::deque<int> dei;
   typedef std::deque<bHex> deH;
   typedef std::deque<ptData> deP;
   typedef std::deque<bSimplex> deX;

   friend class bCentroid;
   friend class bDocker;
   friend class bMC;

   private:

      // Self
      deH edges_;
      deP vrtx_;
      deX simplex_;
      deH simplID_;
      dei simplTy_;
      def simplSc_;
      int max_;
      int min_;
      int size_;
      float score_;
   
      // Debugging IDs
      deH simTrm_;
      deH simOni_;
      deH simInf_;
      deH simInv_;
      deH simRem_;
      deP pntInv_;

      // Sources
      bPoints** src_;
      int       numChains_;
      deP       toAdd_;

      // Information (constantly changing)
      double* currPt_;
      int    currPtId_;
      bHex   sick_;
      bHex   well_;
      deH    newID_;
      dei    newList_;
      dei    hullID_;
      dei    open_;

      // Flags
      bool haveType_;
      bool haveEdge_;
      bool isTrim_;
      bool isSlim_;
      bool isVerified_;
      bool isReady_;
      bool haveScore_;

      // PyMol
      bool printEdg_;
      bool printVal_;
      bool printInv_;
      bool printTrm_;
      bool printRem_;
      bool printOni_;
      bool printInf_;

      // Extras
      MTRand rndNum_;
      static int numInst_;

   protected:

   public:

      // Constructors
      bDelTess();
      bDelTess( const bDelTess&, const bool =(true) );
      ~bDelTess();
      bDelTess& operator=( const bDelTess& );
      void erase();
      void clear();
      void reset();
      bool prep();

      // Setup
      void addSrc( bPoints& );
      void addSrc( bPoints**, int );
      bool relink( bPoints*, bPoints * );
      void randomizePts();
      int  getRandomNumber( const int );
      void defineAA3to1();
   
      // Simplex Manipulation
      int addSimplex( const int [] );
      int addSimplex( bHex & );
      void delSimplex( const int );

      // Tessellation methods
      bool tessellate_full( bPoints&, const char* =(NULL), const char* =(NULL) );
      bool tessellate_full( bPoints**, int, const char* =(NULL), const char* =(NULL), float =(12.0), bool =(false) );
      bool tessellate_std( bPoints**, int, int =(0) );
      bool tessellate_l1o( bPoints**, int, int, int);
      bool tessellate_l1o_fin( bPoints*, int, int =(0) );

      // Tessellation
      bool tessellate();
      bool checkPoint( int );
      void checkSimplex( bSimplex& );
      void outbreak( bSimplex*, bSimplex* );
         void infectedHull( bSimplex* );
         void noInfection( bSimplex*, bSimplex* );
      bool treat();
         bool addNewSimplexes();
         void identifyNeighbors();
         void removeSimplexes();
      bool cleanIteration();
      void immunize();
         void reconcile( dei&, bHex& );
      void resetHandlers();
      bool isNew( bSimplex* );
      
      void dumbCheck();
      void printNghbr( int i );
      void chkNghbr();
      
      // Finalization
      bool verify();
      void retry();
      void onion();
      void findEdges();
      void removeExcess();
      void trim( float =( (float)10.0 ) );
      void trimTet( int, float );
      bool skip( int ) const;
      bool nullNgh( int );
      void focus( int );
      void focus( bPoints* );
      void intersect();
      float score();

      // Output
      void print( FILE* =(stdout) ) const;
      void printTetRes( FILE* =(stdout) );
      void printTetPos( FILE* =(stdout) ) const;
      void printTetPnt( FILE* =(stdout) ) const;
      void printTetCat( FILE* =(stdout) ) const;
      
      // PyMol
      void pymol( FILE*, char[], char[] =("ruby"), char[] =(NULL) );
      void pymolPseudoatoms( FILE*, char[], char[] );
      void pymolEdges( FILE*, char[], char[], char[] =(NULL) );
      void pymolSimplexes( FILE*, char[], char[], deH&, bool =(false) );

      // Not Fully Tested
      void findTypes();
      int  findType( int );
};

#endif
