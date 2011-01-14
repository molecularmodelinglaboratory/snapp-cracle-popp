
#include <stdio.h>

#include "bPoints.h"
#include "bDelTess.h"
#include "bHex.h"
#include "bConf.h"

using namespace std;
using namespace bStd;

uint bConf::cnt_ = 0;

bConf::bConf() {
   //~ this->id_ = new bHex;
   this->pp_ = new bPoints;
   this->nb_ = new bPoints;
   //~ this->dt_ = new bDelTess;
   this->sc_ = 0.0;
   this->pare_ = NULL;
   this->chld_ = NULL;
   this->numChld_ = 0;
}
bConf::bConf( const bConf &rhs ) {
   //~ this->id_ = new bHex( *rhs.id_ );
   this->pp_ = new bPoints( *rhs.pp_ );
   this->nb_ = new bPoints( *rhs.nb_ );
   //~ this->dt_ = new bDelTess( *rhs.dt_ );
   this->sc_ = rhs.sc_;
   this->pare_ = rhs.pare_;
   
   this->numChld_ = rhs.numChld_;
   this->chld_ = new bConf*[ this->numChld_ ];
   for( uint i=0; i < this->numChld_; ++i ) {
      this->chld_[i] = rhs.chld_[i];
   }
}
bConf::~bConf() {
   //~ delete this->id_; this->id_ = NULL;
   delete this->pp_; this->pp_ = NULL;
   delete this->nb_; this->nb_ = NULL;
   //~ delete this->dt_; this->dt_ = NULL;
   this->sc_ = 0.0;
   this->pare_ = NULL;
   for( uint i=0; i < this->numChld_; ++i ) { this->chld_[i] = NULL; }
   delete [] this->chld_; this->chld_ = NULL;
}

bConf& bConf::operator=( const bConf &rhs ) {
   //~ delete this->id_; this->id_ = NULL;
   delete this->pp_; this->pp_ = NULL;
   delete this->nb_; this->nb_ = NULL;
   this->pare_ = NULL;
   for( uint i=0; i < this->numChld_; ++i ) { this->chld_[i] = NULL; }
   delete [] this->chld_; this->chld_ = NULL;
   //~ delete this->dt_; this->dt_ = NULL;
   //~ this->id_ = new bHex( *rhs.id_ );
   this->pp_ = new bPoints( *rhs.pp_ );
   this->nb_ = new bPoints( *rhs.nb_ );
   this->sc_ = rhs.sc_;
   this->pare_ = rhs.pare_;
   this->numChld_ = rhs.numChld_;
   this->chld_ = new bConf*[ this->numChld_ ];
   for( uint i=0; i < this->numChld_; ++i ) {
      this->chld_[i] = rhs.chld_[i];
   }
   //~ this->dt_ = new bDelTess( *rhs.dt_ );
   return *this;
}

bool bConf::operator<( const bConf &rhs ) {
   //~ ++cnt_;
   //~ printf("%3u] %.2f < %.2f = %s\n", cnt_, rhs.sc_, this->sc_, rhs.sc_ < this->sc_ ? "true" : "false");
   return rhs.sc_ < this->sc_;
}