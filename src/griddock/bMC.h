#ifndef BMC_H
#define BMC_H

#include <deque>
#include <list>

#include "bHex.h"
#include "bDelTess.h"
#include "bGrid.h"
#include "bPoints.h"
#include "bList.h"

namespace bStd { class bMC; class bMC_MetTable; };

typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned short ushort;


class bStd::bMC {
   //~ typedef std::deque< uint > dei;
   //~ typedef std::deque< float > def;
   typedef std::deque< ulong > del;
   typedef std::deque< bHex > deH;
   typedef std::deque< bPoints > deP;
   typedef std::deque< bDelTess > deT;
   typedef std::list< bConf > liC;
   typedef std::list< ushort > lis;
   
   friend class bDocker;
   
private:

   //~ bList*  conf_;
   //~ bool    myList_;

   const bGrid*   orgGrd_;
   const bPoints* orgSrc_;
         bConf*   orgPep_;

   bool   havePerturbSchema_;
   int**  perturbSchema_;
   bHex** perturbDone_;
   ushort numOpt_;
   int    perturbRangeMin_[3];
   int    perturbRangeMax_[3];
   ushort perturbTested_;

   float  threshold_;
   ushort numIter_;

   
   uint currPnt_;
   uint currPtb_;
   
   static bool haveGradient_;


public: 
   bMC( bGrid* =(NULL), bPoints* =(NULL), bConf* =(NULL) );
   bMC( const bMC& );
   ~bMC();

   void clear();
   void setup( const bGrid&, const bPoints&, bConf& );

   bool metropolis( const bPoints&, const bPoints& );
   bool metropolis();
   ulong savePepId( const bPoints* );
   
   void setupPerturb( const int[3], const int[3] );
   int  perturb( bConf& );

   void pymol( FILE*, char[], char[] );
};

/****************************** MetTable */
class bStd::bMC_MetTable{
   
private:
   float* delt_;
   float* prob_;
   ushort size_;

public:
   bMC_MetTable();
   ~bMC_MetTable();

   float operator[]( float );

   void read( char[] );

   
};

#endif
