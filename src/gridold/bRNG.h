#ifndef BCALC_H
#define BCALC_H

#include "MersenneTwister.h"

namespace bStd { class bCalc; };

class bCalc {
   private:
      MTRand rndNum_;

   public:
      int rng( int );
      int getRandomNumber( int );
};

#endif
