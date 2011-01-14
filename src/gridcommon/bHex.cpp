#include <stdio.h>
#include <string.h>
#include "bHex.h"
using namespace std;
using namespace bStd;

unsigned int bHex::hamming_[65536]; // 16-bit hamming table
bool bHex::isHamRdy_;

/* Constructors */
/* default */
bHex::bHex() : size_( 0 ), numBits_(0), numFlip_(0), hex_(NULL), active_(NULL) {
   this->resize(0);
   this->isVerified_ = false;
   this->haveActive_ = false;
   this->haveNumFlip_ = false;
}

/* somewhat default */
bHex::bHex(int nB) : size_( 0 ), numBits_(0), numFlip_(0), hex_(NULL), active_(NULL) {
   this->resize(nB);
   this->isVerified_ = false;
   this->haveActive_ = false;
   this->haveNumFlip_ = false;
}


/* deep copy constructor */
bHex::bHex( const bHex &h ) : size_( 0 ), numBits_(0), numFlip_(0), hex_(NULL), active_(NULL) {
   this->size_ = h.size_;
   this->numBits_ = h.numBits_;
   this->numFlip_ = h.numFlip_;
   
   this->hex_ = (h.size_ == 0) ? NULL : new ulong[h.size_];
   
   for( int i=0; i < h.size_; ++i ) {
      this->hex_[i] = h.hex_[i];
   }
   
   this->isVerified_ = false;
   this->haveActive_ = false;
   this->haveNumFlip_ = false;
   this->getActive();
}

/* Destructor */
bHex::~bHex() {
   if( this->hex_ != NULL ) {
      delete [] this->hex_;
      this->hex_ = NULL;
   }
   if( this->active_ != NULL ) {
      delete [] this->active_;
      this->active_ = NULL;
   }
}

/* Assignment operators [Copy] */
bHex& bHex::operator=( const bHex &rhs ) {
   if( this == &rhs ) { return *this; }

   delete [] this->hex_;
   this->hex_ = NULL;

   delete [] this->active_;
   this->active_ = NULL;

   this->size_ = rhs.size_;
   this->numBits_ = rhs.numBits_;
   this->numFlip_ = rhs.numFlip_;

   this->hex_ = rhs.size_ == 0 ? NULL : new ulong[rhs.size_];

   for( int i=0; i < rhs.size_; ++i ) {
      this->hex_[i] = rhs.hex_[i];
   }

   this->isVerified_ = false;
   this->haveActive_ = false;
   this->haveNumFlip_ = false;
   this->getActive(); 
   return *this;
}
void bHex::erase() {
   for( int i=0; i < this->size_; ++i ) { this->hex_[i] &= 0x0; }
   if( this->active_ != NULL ) { delete [] this->active_; this->active_ = NULL; }
   this->numFlip_ = 0;
   this->haveActive_ = false;
   this->haveNumFlip_ = false;
   this->isVerified_ = false;
   return;
}

bHex& bHex::operator=( const int rhs ) {
   this->erase();
   if( rhs != 0 ) {
      int n = rhs;
      int section = this->adjust( n );
      ulong translator = 1;
      translator <<= n; // *
      this->hex_[section] |= translator;
   }
   
   this->isVerified_ = false;
   return *this;
}
/* OR ASSIGNMENT operator 
   Flips a single bit at location 'n'. Allows access to higher order
   bits without a full function
*/
bHex& bHex::operator|=( const int rhs ) {
   int n = rhs;
   int section = this->adjust( n );
   
   ulong translator = 1;
   translator <<= n; // *
   this->hex_[section] |= translator;
   
   this->isVerified_ = false;
   return *this;
}

bHex& bHex::operator|=( const bHex& rhs ) {
   if( this == &rhs ) { return *this; }
   int size = this->size_ > rhs.size_ ? rhs.size_ : this->size_;
   for( int i=0; i < size; ++i )
      this->hex_[i] |= rhs.hex_[i];
   
   this->isVerified_ = false;
   return *this;
}


/* AND ASSIGNMENT */
bHex& bHex::operator&=( const int rhs ) {
   int n = rhs;
   int section = this->adjust( n );
   
   ulong translator = 1;
   translator <<= n; // *
   this->hex_[section] &= translator;
   
   this->isVerified_ = false;
   return *this;
}

bHex& bHex::operator&=( const bHex& rhs ) {
   if( this == &rhs ) { return *this; }
   int size = this->size_ > rhs.size_ ? rhs.size_ : this->size_;
   for( int i=0; i < size; ++i ) {
      this->hex_[i] &= rhs.hex_[i];
   }
   this->isVerified_ = false;
   return *this;
}

/* XOR ASSIGNMENT */
bHex& bHex::operator^=( const int rhs ) {
   int n = rhs;
   int section = this->adjust( n );
   
   ulong translator = 1;
   translator <<= n; // *
   this->hex_[section] ^= translator;
   this->isVerified_ = false;
   return *this;
}

bHex& bHex::operator^=( const bHex& rhs ) {
   if( this == &rhs ) { return *this; }
   int size = this->size_ > rhs.size_ ? rhs.size_ : this->size_;
   for( int i=0; i < size; ++i ) {
      this->hex_[i] ^= rhs.hex_[i];
   }
   this->isVerified_ = false;
   return *this;
}

/* NOT ASSIGNMENT */
bHex& bHex::flipBits() {
   //~ if( this == &rhs ) { return *this; }
   for( int i=0; i < this->size_; ++i ) {
      this->hex_[i] = ~this->hex_[i];
   }
   this->isVerified_ = false;
   return *this;
}

/* ADJUST */
int bHex::adjust( int &adj ) {
   int section = 0;
   while( adj > 63 ) {
      adj -= 64;
      ++section;
   }
   return section;
}



/* Comparison operator */
bool bHex::operator==( const bHex &h ) const {
   bool equal = true;
   if( this->size_ != h.size_ ) { return false; }
   for( int i=0; i < this->size_ && equal; ++i ) {
      equal = ( this->hex_[i] == h.hex_[i] );
   }
   return equal;
}
bool bHex::operator!=( const bHex &h ) const {
   return !( *this == h );
}

bool bHex::operator==( const ulong rhs ) const {
   return (this->hex_[0] == rhs);
}
bool bHex::operator!=( const ulong rhs ) const {
   return !( *this == rhs );
}


bool bHex::operator<( const bHex &rhs ) const {
   bool lt = true;
   int size = this->size_;
   if( size != rhs.size_ ) {
      if( size > rhs.size_ ) { return lt; }
      else { return false; }
   }
   for( int i=(size - 1); i >= 0; --i ) {
      if( this->hex_[i] > rhs.hex_[i] ) { break; }
      else if( this->hex_[i] < rhs.hex_[i] ) { lt = false; break; }
      else {}
   }
   return lt;
}

bool bHex::operator<=( const bHex &rhs  ) const {
   bool le = true;
   if( this->size_ != rhs.size_ ) {
      if( this->size_ > rhs.size_ ) { return le; }
      else { return false; }
   }
   for( int i=this->size_ - 1; i >= 0; --i ) {
      if( this->hex_[i] <= rhs.hex_[i] ) { break; }
      else if( this->hex_[i] > rhs.hex_[i] ) { le = false; break; }
      else {}
   }
   return le;
}
bool bHex::operator>( const bHex &rhs ) const { return !( *this <= rhs ); }
bool bHex::operator>=( const bHex &rhs ) const { return !( *this < rhs ); }

const bHex bHex::operator&( const bHex &rhs ) const { return bHex( *this ) &= rhs; }
const bHex bHex::operator|( const bHex &rhs ) const { return bHex( *this ) |= rhs; }
const bHex bHex::operator^( const bHex &rhs ) const { return bHex( *this ) ^= rhs; }

bool bHex::operator&( const ulong rhs ) const {
   return this->hex_[0] & rhs;
}

bool bHex::operator&( const int num ) const {
   int section = 0;
   int n = num;
   while( n > 63 ) {
      n -= 64;
      ++section;
   }
   
   
   ulong translator = 1;
   translator <<= n; // *
   return this->hex_[section] & translator;
}




int bHex::size() const {
   return this->numBits_; // tad confusing; meant to work with resize
}


/* OR
   -- Performs the expected | function by assigning multiple bits
   This implementation allows for higher order access
*/
void bHex::assign( const ulong num, int depth ) {
   if( num == 0 ) { return; }

   int section = adjust( depth );
   //~ int d = depth; // - 64 * section; // in adjust
   ulong n1 = num;
   ulong n2 = num;

   n1 <<= depth;
   this->hex_[section] |= n1;

   //~ this->print();
   //~ this->_printHex( num );
   //~ printf("\tdepth: %d\n", depth);

   if( depth != 0 ) {
      int rev = 64 - depth;
      n2 >>= rev;
      ++section;
      //~ printf("size: %d\tsection: %d\t%#lx\n", size_, section, this->hex_[section] );
      if( (n2 > 0) && (this->size_ > section) ) { this->hex_[section] |= n2; }
   }
   this->isVerified_ = false;
   return;
}

void bHex::filter( const ulong num, int depth ) {
   int section = adjust( depth );
   //~ int d = depth; // - 64 * section;
   ulong n1 = num;
   ulong n2 = num;
   
   int pre = 0;
   while( pre < section ) {
      this->hex_[pre] = 0;
      ++pre;
   }
   
   n1 <<= depth;
   this->hex_[section] &= n1;
   
   if( depth != 0 ) {
      int rev = 64 - depth;
      n2 >>= rev;
      this->hex_[++section] &= n2;
   }
   
   while( ++section < this->size_ ) {
      this->hex_[section] = 0;
   }
   this->isVerified_ = false;
   return;
}

void bHex::remove( const ulong num, int depth ) {
   int section = adjust( depth );
   //~ int d = depth - 64 * section;
   ulong n1 = num;
   ulong n2 = num;
   
   n1 <<= depth;
   this->hex_[section] |= n1;
   this->hex_[section] ^= n1;
   
   if( depth != 0 ) {
      int rev = 64 - depth;
      n2 >>= rev;
      this->hex_[++section] |= n2;
      this->hex_[++section] ^= n2;
   }
   this->isVerified_ = false;
   return;
}

/* Resize */
void bHex::resize( int needBits ) {
   if( this->numBits_ >= needBits ) { return; } // don't resize if we don't need to

   int numSec = 1;
   --needBits;
   while( needBits > 63 ) {
      needBits -= 64;
      ++numSec;
   }

   if( this->hex_ == NULL ) {
      this->hex_ = new ulong[numSec];
      for(int i=0; i < numSec; ++i) { this->hex_[i] = 0x0; }
   }
   else {
      bHex t(*this);
      delete [] this->hex_; // fixed this with a new copy constructor
      delete [] this->active_;
      this->hex_ = NULL;
      this->active_ = NULL;

      this->hex_ = new ulong[numSec];
      for( int i=0; i < t.size_; ++i ) { this->hex_[i] = t.hex_[i]; }

      //~ this->findActive();
      //~ this->active_ = new int[this->numFlip_];
      //~ for( int i=0; i < t.numFlip_; ++i ) {
         //~ this->active_[i] = t.active_[i];
      //~ }
      
   }

   this->numBits_ = 64 * numSec; // do we need a last used?
   --this->numBits_; // last index
   this->size_ = numSec;

   this->isVerified_ = false;
   return;
}


/* Empty */
bool bHex::empty() const {
   for( int i=0; i < this->size_; ++i ) {
      if( this->hex_[i] != 0x0 ) {
         return false;
      }
   }
   return true;
}

/* Print */
void bHex::print( FILE* op ) const {
   for( int i=0; i < this->size_; ++i ) {
      if( this->hex_[i] == 0 ) { continue; }
      fprintf(op, "[%d] ",i);
      this->printBit( this->hex_[i], op);
      fprintf(op, " ");
   }
   printf("\n");
   
   return;
}
void bHex::print_full( FILE* op ) {
   for( int i=0; i < this->size_; ++i ) {
      if( this->hex_[i] == 0 ) { continue; }
      fprintf(op,  "%d: %#18lx\t", i, this->hex_[i] );
      this->printBit( this->hex_[i], op);
      fprintf(op, "\n");
   }
   
   return;
}

void bHex::printBit( ulong t, FILE* op ) {
   if(!isHamRdy_) { setupHammingTable(); }
   int amt = 64;
   while(t && amt > 0) {
      int cnt = bHex::getHammingWeight( t );
      --t;
      cnt -= bHex::getHammingWeight( t );
      cnt *= -1;
      cnt += 2; // +2 to account for the space needed for '1'
      fprintf(op, "%0*d",cnt,1);
      t >>= cnt;
      amt -= cnt;
   }
   
   if( amt > 0 ) {
      fprintf(op, "%0*d",amt,0);
   }
   
   return;
}

void bHex::printHex( FILE* op ) {
   printf("0x");
   for( int i=(this->size_ - 1); i >= 0; --i ) {
      printHex( this->hex_[i], op, false );
   }
   printf("\n");
}

void bHex::printHex( ulong h, FILE* op, bool standalone ) {
   ulong hex[4];
   register int i = 0;
   for( i=0; h; ++i ) { hex[i] = h & 0xFFFF; h >>= 16; }
   for( int j=i; j < 4; ++j ) { hex[j] = 0; }
   if( standalone ) { fprintf(op, "0x"); }
   for( int j=3; j >= 0; --j ) { fprintf(op, "%04lx ", hex[j]); }
   return;
}


/* Count Number of Flipped Bits */
int bHex::getNumFlip() {
   if( this->isVerified_ && this->haveNumFlip_ ) { return this->numFlip_; }
   int flipped = 0;
   for( int i=0; i < this->size_; ++i) {
      flipped += this->hex_[i] == 0 ? 0 : getHammingWeight( this->hex_[i] );
   }
   this->numFlip_ = flipped;
   
   this->haveNumFlip_ = true;
   this->isVerified_ = true;
   return this->numFlip_;
}

bHex::ulong bHex::getRow( int i ) {
   return this->hex_[i];
}

/* List Active
>> Description
   Returns a reference to an array listing the active indices
*/
int* bHex::getActive() {
   if( this->isVerified_ && this->haveActive_ ) { return this->active_; }
   if(!isHamRdy_) { setupHammingTable(); }
   if( this->active_ != NULL ) {
      delete [] this->active_;
      this->active_ = NULL;
   }

   if( this->getNumFlip() == 0 ) { return NULL; }
   this->active_ = new int[ this->numFlip_ ];
   
   // find all positions
   int pos = 0;
   for( int i=0; i < this->size_; ++i ) {
      ulong t = this->hex_[i];
      int correction = i*64;
      while(t) {
         int cnt = bHex::getHammingWeight( t ); // num of ones
         ulong init = t; // save old value of t
         --t;        // flip bits up to and incl lowest
         cnt -= bHex::getHammingWeight( t ); // recount
         cnt *= -1;  // get positive
         ++cnt;      // increment for position
         t &= init;  // erase previous
         
         this->active_[pos] = cnt + correction;
         ++pos;
      }
   }
   
   this->isVerified_ = true;
   this->haveActive_ = true;
   return this->active_;
}

void bHex::getActive( int active[] ) const {
   if(!isHamRdy_) { setupHammingTable(); }
   
   if( this->active_ != NULL ) {
      active[0] = this->active_[0];
      active[1] = this->active_[1];
      active[2] = this->active_[2];
      active[3] = this->active_[3];
   }
   else {
      int pos = 0;
      for( int i=0; i < this->size_; ++i ) {
         ulong t = this->hex_[i];
         int correction = i*64;
         while(t) {
            int cnt = bHex::getHammingWeight( t ); // num of ones
            ulong init = t; // save old value of t
            --t;        // flip bits up to and incl lowest
            cnt -= bHex::getHammingWeight( t ); // recount
            cnt *= -1;  // get positive
            ++cnt;      // increment for position
            t &= init;  // erase previous
            
            active[pos] = cnt + correction;
            ++pos;
         }
      }
   }
   return;
}

void bHex::getActive( ulong ul, int active[] ) {
   if(!isHamRdy_) { setupHammingTable(); }
   int pos = 0;
   while( ul ) {
      int cnt = bHex::getHammingWeight( ul ); // num of ones
      ulong init = ul; // save old value of t
      --ul;        // flip bits up to and incl lowest
      cnt -= bHex::getHammingWeight( ul ); // recount
      cnt *= -1;  // get positive
      ++cnt;      // increment for position
      ul &= init;  // erase previous
      
      active[pos] = cnt;
      ++pos;
   }
   return;
}

/* Setup Hamming Table
>> Description
   If we can store a lookup table of the hamming function of every 16 bit integer, we can do the following to compute the Hamming weight of every 32 bit integer.
      static unsigned char wordbits[65536] = { bitcounts of ints between 0 and 65535 };
      static int popcount( unsigned int i )
      {
          return( wordbits[i&0xFFFF] + wordbits[i>>16] );
      }
*/
void bHex::setupHammingTable() {
   register int a = 65536;
   for(int i =0; i<a; ++i) {
      bHex::hamming_[i] = calculateHammingWeight(i);
   }
   bHex::isHamRdy_ = true;
   return;
}

/*	Calculate Hamming Weight */
int bHex::calculateHammingWeight(int x) {
   int i=0;
   while(x) {
      x &= x-1;
      ++i;
   }
   return i;
}

/*	Get Hamming Weight for 64 bit numbers */
int bHex::getHammingWeight(unsigned long x) {
   int weight = 0;
   for(int i=0; i<4 && x > 0; ++i) {
      weight += bHex::hamming_[int(x & 0xFFFF)];
      x >>= 16;
   }
   return weight;
}

int bHex::getHammingWeight() {
   this->numFlip_ = 0;
   for( int i=0; i < this->size_; ++i) {
      this->numFlip_ += getHammingWeight( this->hex_[i] );
   }
   return this->numFlip_;
}




/* EOF */


