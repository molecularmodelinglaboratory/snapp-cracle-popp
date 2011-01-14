#ifndef BCALC_H
#define BCALC_H

#include "MersenneTwister.h"

namespace bStd { class bCalc; };

class bStd::bCalc {
private:
   static MTRand rndNum_;
   static bool   isSeeded_;

public:
   bCalc();
   ~bCalc();

   static int rng( int );
   static int getRandomNumber( int );
};

#endif
