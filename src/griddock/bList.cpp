
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bPoints.h>
#include "bList.h"


using namespace std;
using namespace bStd;

float bList::_resizeAmt_ = 1.2;
int   bList::_step_ = 20;

/****** bCONF ******/

bConf::bConf()
   : sc_(0.0), next_( NULL ), prev_( NULL ), pare_(NULL), chld_(NULL)
{
   this->pp_ = NULL;
}

bConf::bConf( const bConf &rhs )
   : sc_(0.0), next_( NULL ), prev_( NULL ), pare_(NULL), chld_(NULL)
{
   this->pp_ = new bPoints( *(rhs.pp_) );
   this->sc_ = rhs.sc_;
   this->next_ = rhs.next_;
   this->prev_ = rhs.prev_;
   this->pare_ = rhs.pare_;
   this->chld_ = new bList( *(rhs.chld_) );
}

bConf::bConf( const bPoints &rhs )
   : sc_(0.0), next_( NULL ), prev_( NULL ), pare_(NULL), chld_(NULL)
{
   this->pp_ = new bPoints( rhs );
}

bConf::bConf( bPoints &rhs, float s )
   : sc_(s), next_( NULL ), prev_( NULL ), pare_(NULL), chld_(NULL)
{
   this->pp_ = new bPoints( rhs );
}

bConf::~bConf() {
   this->sc_ = 0.0;
   this->pare_ = NULL;
   delete this->pp_; this->pp_ = NULL;
   delete this->chld_; this->chld_ = NULL;
}

void bConf::clear() {
   if( this->pp_ != NULL ) { delete this->pp_; this->pp_ = NULL; }
   this->sc_ = 0.0;
   this->next_ = NULL;
   this->prev_ = NULL;
   this->pare_ = NULL;
   if( this->chld_ != NULL ) { delete this->chld_; this->chld_ = NULL; }
   return;
}

/****** Assignment */
bConf& bConf::operator=( const bConf &rhs ) {
   if( this != &rhs ) {
      this->clear();
      this->pp_ = new bPoints( *(rhs.pp_) );
      this->sc_ = rhs.sc_;
      this->next_ = rhs.next_;
      this->prev_ = rhs.prev_;
      this->pare_ = rhs.pare_;
      this->chld_ = new bList( *(rhs.chld_) );
   }
   return *this;
}

/****** Conditional */
bool bConf::operator< ( const bConf &rhs ) { return rhs.sc_ <  this->sc_; }
bool bConf::operator> ( const bConf &rhs ) { return rhs.sc_ >  this->sc_; }
bool bConf::operator==( const bConf &rhs ) { return rhs.sc_ == this->sc_; }

/****** Get */
bConf* bConf::next() { return this->next_; }
bConf* bConf::prev() { return this->prev_; }
bConf* bConf::pare() { return this->pare_->pare_; }
bList* bConf::home() { return this->pare_; }
bList* bConf::chld() { return this->chld_; }
bList* bConf::new_chld() { 
   if( this->chld_ != NULL ) {
      this->chld_->clear();
      delete this->chld_;
      this->chld_ = NULL;
   }
   this->chld_ = new bList;
   this->chld_->pare( this );
   return this->chld_;
}
float bConf::score() { return this->sc_; }
bPoints* bConf::pts() { return this->pp_; }

/****** Set */
void bConf::next( bConf* n ) { this->next_ = n; return; }
void bConf::prev( bConf* n ) { this->prev_ = n; return; }
void bConf::pare( bList* l ) { this->pare_ = l; return; }
void bConf::chld( bList* l ) {
   if( this->chld_ != NULL ) {
      delete this->chld_;
      this->chld_ = NULL;
   }
   this->chld_ = new bList( *l );
   this->chld_->pare( this );
   return;
}
void bConf::chld( bConf* n ) {
   if( this->chld_ == NULL ) {
      this->chld_ = new bList( n );
      this->chld_->pare( this );
   }
   else { this->chld_->push_back( n ); }
   return;
}
void  bConf::score( float s ) { this->sc_ = s; return; }


void bConf::pymol( FILE* op, char name[], char color[] ) {
   char label[8];
   sprintf( label, "%.2f", this->sc_ );
   //~ pp_->print();
   
   //~ this->pp_->pymolConnectedPseudoatoms( op, name, color, 1, label );
   this->pp_->setToNormalSpace();
   //~ sprintf(name,"%s-NS",name);
   this->pp_->pymolConnectedPseudoatoms( op, name, color, 1, label );
   if( this->chld_ != NULL ) {
      this->chld_->pymol( op, name );
   }
   this->pp_->setToDockingSpace();
   //~ sprintf(name,"%s-DS",name);
   //~ this->pp_->pymolConnectedPseudoatoms( op, name, color, 1, label );
   //~ pp_->print();
}

void bConf::print() {
   printf("score: %.2f\n", this->sc_);
   this->pp_->print();
}


/************************** bLIST ********************************************/
/************************** bLIST ********************************************/
/************************** bLIST ********************************************/
bList::bList( uint c )
   : cap_(c), num_(0), pare_(NULL), beg_(NULL), end_(NULL), list_(NULL)
{
   this->list_ = new bConf*[ this->cap_ ];
}

bList::bList( const bList &rhs, uint c )
   : cap_(c), num_(0), pare_(NULL), beg_(NULL), end_(NULL), list_(NULL)
{
   this->cap_ = rhs.num_ < rhs.cap_ ? rhs.cap_ : rhs.num_ * _resizeAmt_;
   this->num_ = rhs.num_;
   this->pare_ = rhs.pare_;

   this->list_ = new bConf*[ this->cap_ ];
   this->list_[0] = new bConf( *(rhs.list_[0]) );
   uint i=1, k=0;
   this->beg_ = this->list_[k];
   for( i=1; i < this->num_; ++i, ++k ) {
      this->list_[i] = new bConf( *(rhs.list_[i]) );
      this->list_[i]->prev_ = this->list_[k];
      this->list_[k]->next_ = this->list_[i];
      this->list_[i]->pare_ = this;
   }
   this->end_ = this->list_[i];
   for( i=i; i < this->cap_; ++i ) {
      this->list_[i] = NULL;
   }

   this->list_[0]->prev_ = NULL;
   this->list_[0]->pare_ = this;
   this->list_[i]->next_ = NULL;
}

bList::bList( const bConf &rhs, uint c )
   : cap_(c), num_(1), pare_(NULL), beg_(NULL), end_(NULL), list_(NULL)
{
   this->list_ = new bConf*[ this->cap_ ];
   this->list_[0] = new bConf( rhs );
   this->list_[0]->pare( this );
   this->list_[0]->next( NULL );
   this->list_[0]->prev( NULL );
   this->beg_ = this->list_[0];
   this->end_ = this->list_[0];
   ++this->num_;
}

bList::bList( bConf* rhs, uint c )
   : cap_(c), num_(1), pare_(NULL), beg_(rhs), end_(rhs), list_(NULL)
{
   this->list_ = new bConf*[ this->cap_ ];
   this->list_[0] = rhs;
   this->list_[0]->pare( this );
   this->list_[0]->next( NULL );
   this->list_[0]->prev( NULL );
   this->beg_ = this->list_[0];
   this->end_ = this->list_[0];
   ++this->num_;
}

bList::~bList() {
   bConf* curr = this->beg_;
   bConf* next = NULL;
   for( uint i=0; i < this->num_; ++i ) {
      next = curr->next();
      delete curr;
      curr = next;
   }
   delete [] this->list_;
   this->list_ = NULL;
   
   this->cap_ = 0;
   this->num_ = 0;
   this->pare_ = NULL;
   this->beg_ = NULL;
   this->end_ = NULL;
}

/***** DATA MANAGEMENT ******/
void bList::clear() {
   bConf* next;
   bConf* curr = this->beg_;
   while( curr != this->end_ ) {
      next = curr->next_;
      curr->clear();
      delete curr;
      curr = next;
   }
}

/****** OVERLOADED ******/
bConf& bList::operator[]( uint x ) {
   bConf* access = NULL;
   if( x < this->cap_ ) { access = this->list_[x]; }
   else {
      uint pos = this->cap_; --pos;
      access = this->list_[pos]->next();
      while( ++pos != x ) {
         access = access->next();
      }
   }
   return *access;
}

/****** BOOL ******/
bool bList::empty() { return (this->num_ == 0) ? true : false; }

/*** GET */
bConf* bList::pare() { return this->pare_; }
bConf* bList::beg()  { return this->beg_; }
bConf* bList::end()  { return this->end_; }
uint   bList::cap()  { return this->cap_; }
uint   bList::num()  { return this->num_; }
uint   bList::size() { return this->num_; }
void bList::resize( uint rs ) {
   if( this->list_ != NULL ) { delete [] this->list_; this->list_ = NULL; }

   // Make sure it's big enough
   if( rs < this->num_ ) { rs = this->num_; rs *= _resizeAmt_; }
   this->list_ = new bConf*[rs];
   this->cap_ = rs;

   // Repopulate the list
   bConf* access = this->beg_;
   for( uint i=0; i < this->num_; ++i ) {
      this->list_[i] = access;
      access = access->next();
   }

   // Ensure null-ness
   for( uint i=this->num_; i < rs; ++i ) {
      this->list_[i] = NULL;
   }
   return;
}


/*** SET */
void bList::pare( bConf* n ) { this->pare_ = n; return; }
void bList::cap( uint n ) { this->resize( n ); return; }

/*** MODIFY */
void bList::push_back( bConf &c ) {
   bConf* n = new bConf( c );
   this->push_back( n );
   return;
}

void bList::push_back( bConf* n ) {
   // Resize if needed
   uint tooBig = this->cap_; tooBig *= _resizeAmt_;
   if( this->num_ > tooBig ) { tooBig *= _resizeAmt_; this->resize( tooBig ); }

   // Save in array if room
   if( this->num_ < this->cap_ ) {
      this->list_[ this->num_ ] = n;
   }

   // Maintain integrity
   if( this->num_ == 0 ) {
      this->end_ = n;
      this->beg_ = n;
   }
   else {
      this->end_->next( n );
      n->prev( this->end_ );
      this->end_ = n;
   }

   ++this->num_;
   return;
}

void bList::del( uint pos ) {
   // Remove from pointers
   bConf* access = &((*this)[ pos ]);
   access->next()->prev( access->prev() );
   access->prev()->next( access->next() );

   // Check for extremeties
   if( pos == 0 ) { this->beg_ = access->next(); }
   else if ( pos == this->num_ - 1 ) { this->end_ = access->prev(); }
   else {}

   // Maintain the list
   bConf* toDel = access;
   while( pos < this->cap_ && pos < this->num_ ) {
      access = access->next();
      this->list_[pos] = access;
      ++pos;
   }

   // Delete element
   delete toDel;
   toDel = NULL;
   access = NULL;
   return;
}

void bList::insert( bConf &n ) {
   return;
}

/****** Sort */
void bList::sort( bool rev ) {
   int right = this->num_; --right;
   if( this->num_ > this->cap_ ) { this->resize( this->num_ * _resizeAmt_ ); }
   rev ? this->_sortRev( 0, right ) : this->_sort( 0, right );

   // Pointer Maintenance
   uint i = 0;
   this->beg_ = this->list_[0];
   bConf* prev = this->beg_;
   this->list_[0]->prev( NULL );
   for( i=1; i < this->num_; ++i ) {
      prev->next( this->list_[i] );
      this->list_[i]->prev( prev );
      prev = this->list_[i];
   }
   this->list_[--i]->next( NULL ); // strange --i...
   this->end_ = this->list_[i];

   return;
}

void bList::_sortRev( int left, int right ) {
   if( left == right ) { return; } // return if we're done
   else if( (left + 1) == right ) {
      if( list_[left]->sc_ < list_[right]->sc_ ) { swap( left, right ); }
      return;
   }
   else {}

   int pivot = (right + left) >> 1; // get the pivot (i.e. the middle point)
   if( left > right ) exit(1);
   float pivotVal = this->list_[pivot]->sc_; // save the value

   swap( pivot, right ); // move the pivot element to the end  (temporarily)
   pivot = right;

   int store = left; // start at the left 
   for(int i = left; i < right; ++i) {
      if( this->list_[i]->sc_ == pivotVal) {
         swap( i--, --right );
      }
      else if( this->list_[i]->sc_ > pivotVal) {    // if less than pivot
         swap( i, store );  // swap points
         ++store;  // move the store pointer
      }
      else {}
   }

   int newStore = store;
   if( store < right ) {
      while( right < pivot ) {
         swap( store, right );
         ++store;
         ++right;
      }
   }

   swap( pivot, store ); // put the pivot back in place
   _sortRev( left, newStore ); // recursively iterate through both sides
   _sortRev( store, right );

   return;
}

void bList::_sort( int left, int right ) {
   if( left == right ) { return; } // return if we're done
   else if( (left + 1) == right ) { return; }
   else {}

   int pivot = (right + left) >> 1; // get the pivot (i.e. the middle point)
   float pivotVal = this->list_[pivot]->sc_; // save the value

   swap( pivot, right ); // move the pivot element to the end  (temporarily)
   pivot = right;

   int store = left; // start at the left 
   for(int i = left; i < right; ++i) {
      if( this->list_[i]->sc_ == pivotVal) {
         swap( i--, --right );
      }
      else if( this->list_[i]->sc_ < pivotVal) {    // if less than pivot
         swap( i, store );  // swap points
         ++store;  // move the store pointer
      }
      else {}
   }

   int newStore = store;
   if( store < right ) {
      while( right < pivot ) {
         swap( store, right );
         ++store;
         ++right;
      }
   }

   swap( pivot, store ); // put the pivot back in place
   _sort( left, newStore ); // recursively iterate through both sides
   _sort( store, right );

   return;
}

void bList::swap( int a, int b ) {
   if( a == b ) { return; }
   else {}
   bConf* t = this->list_[a];
   this->list_[a] = this->list_[b];
   this->list_[b] = t;
   t = NULL;
   return;
}

/****** UTILITY ******/

void bList::trim( uint t ) {
   //~ printf("TRIM\n");
   bConf* c = NULL;
   bConf* p = NULL;
   if( t >= this->num_ ) { return; }

   // Destroy old list and recreate
   delete [] this->list_;
   this->list_ = NULL;
   this->list_ = new bConf*[ t ];

   // Save new list
   uint i = 1;
   this->list_[0] = this->beg_;
   c = this->list_[0];
   for( ; i < t; ++i ) {
      this->list_[i] = c->next();
      c = this->list_[i];
   }
   this->end_ = c;

   // Erase old list
   p = c->next();
   while( p != NULL ) {
      c = p->next();
      p->clear();
      delete p;
      if( c != NULL ) { p = c->next(); }
      else { p = NULL; }
   }

   this->num_ = t;
   this->cap_ = t;
   this->end_->next( NULL );
   return;
}

/****** PYMOL ******/
void bList::pymol( FILE* op, char  name[] ) {
   int tc = 0;
   int gc = 0;
   int gr = 0;
   int grStep = this->num_ < 100 ? this->num_ : 100; grStep /= (_step_);
   char color[8];
   memset( color, '\0', 8 );
   memmove( color, "grad0", 5 );
   
   char n[32];
   memset( n, '\0', 32 );
   
   //~ printf("PYMOL:\n");
   //~ print();
   bConf* it = this->beg_;
   while( it != NULL ) {
      if( tc > 100 ) break;
      sprintf( n, "%s.%d", name, tc );
      it->pymol( op, n, color );
      
      
      it = it->next();
      ++tc;
      if( ++gc > grStep ) {
         gc = 0; ++gr;
         sprintf( color, "grad%d", gr );
      }
   }
}

void bList::pymolGradient( FILE *op, int step ) {
   // Setup color gradient
   _step_ = step;
   char gradient[8];
   float red = 1.0;
   float grn = 0.0;
   float blu = 0.0;
   float inc = 1.0 / step;
   for( int i=0; i < step; ++i ) {
      memset( gradient, '\0', 8 );
      sprintf( gradient, "grad%d", i );
      fprintf( op, "cmd.set_color(\"%s\", [1.00, %.4f, 0.0])\n",gradient,grn );
      grn += inc;
      red += 0.0;
      blu += 0.0;
   }
}

void bList::print() {
   printf("beg: %p\t[%u, %u]\n", this->beg_, this->num_, this->cap_);
   bConf* it = this->beg_;
   uint i= 0;
   for( ; i < 20; ++i ) {
      if( i >= this->num_ ) { break; }
      printf("[%3d] %p [%p]\tnpt: %d\tscr: %.2f\n", i,it, it, it->pp_->size(), it->sc_);
      it = it->next();
   }
   for( ; i < this->num_ - 10; ++i ) {
      if( i >= this->num_ ) { break; }
      it = it->next();
   }
   if( i == this->num_ - 10 ) { printf("...\n"); }
   for( ; i < this->num_; ++i ) {
      if( i >= this->num_ ) { break; }
      printf("[%3d]: %p [%p]\tnpt: %d\tscr: %.2f\n", i,it, it, it->pp_->size(), it->sc_);
      it = it->next();
   }
   printf("end: %p\n\n", this->end_);
}
