
#include <ctime>

#include "bCalc.h"
#include "MersenneTwister.h"

using namespace std;
using namespace bStd;

MTRand bCalc::rndNum_;
bool bCalc::isSeeded_ = false;

bCalc::bCalc() {}
bCalc::~bCalc() {
}

/* Mersenne Twister Interface */
MTRand rndNum_;
int bCalc::rng( int range ) { return bCalc::getRandomNumber(range); }
int bCalc::getRandomNumber( int range ) {
   if( !isSeeded_ ) { rndNum_.seed((unsigned int)time(NULL)); isSeeded_ = true; }
   if(range <= 0) { return -1; } // check range...we could let it be lower than zero; just adjust

   if((range & -range) == range) { // re[ar]range
      return (int)((range*(long)rndNum_.randInt()) >> 31);
   }

   int bits = 0; // put in range...longer implementation
   int val = 0;
   do {
      bits = rndNum_.randInt();
      val = bits % range;
   } while(bits - val + (range-1) < 0);

   return val;
}
