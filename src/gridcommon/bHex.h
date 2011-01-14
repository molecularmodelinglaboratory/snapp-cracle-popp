#ifndef BHEX_H
#define BHEX_H

#include <stdio.h>
#include "bVBase.h"

namespace bStd { class bHex; };

class bStd::bHex : virtual bVBase {
   public:
      typedef ulong* uhex;
      static uint hamming_[65536];
      static bool isHamRdy_;
   
      bHex();
      bHex( int );
      bHex( const bHex& );
      ~bHex();
   
      bHex& operator=( const bHex& );
      bHex& operator=( const int );
      bHex& operator|=( const bHex& );
      bHex& operator|=( const int );
      bHex& operator&=( const bHex& );
      bHex& operator&=( const int );
      bHex& operator^=( const bHex& );
      bHex& operator^=( const int );
      bHex& flipBits();
      int adjust( int& );
   
      bool operator==( const bHex& ) const ;
      bool operator!=( const bHex& ) const ;
      bool operator<( const bHex& ) const ;
      bool operator>( const bHex& ) const ;
      bool operator<=( const bHex& ) const ;
      bool operator>=( const bHex& ) const ;
      const bHex operator&( const bHex& ) const ;
      const bHex operator|( const bHex& ) const ;
      const bHex operator^( const bHex& ) const ;
      
      bool operator==( const ulong ) const ;
      bool operator!=( const ulong ) const ;
      bool operator&( const int ) const ;
      
      bool operator&( const ulong ) const ;

      void assign( const ulong, int =(0) );
      void filter( const ulong, int =(0) );
      void remove( const ulong, int =(0) );
      bool empty() const;
      void erase();

      int size() const;
      void resize();
      void resize( int );

      void print( FILE* =(stdout) ) const;
      void print_full( FILE* =(stdout) );
      void printHex( FILE* =(stdout) );
      static void printBit( ulong, FILE* =(stdout) );
      static void printHex( ulong, FILE* =(stdout), bool =(true) );
      
      ulong getRow( int );
      int getNumFlip();
      int* getActive();
      void getActive( int [] ) const;
      static void getActive( ulong, int [] );
      
      static void setupHammingTable();
      static int calculateHammingWeight( int );
      static int getHammingWeight( ulong );
      int getHammingWeight();
   
   protected:
      
   private:
      int size_;
      int numBits_;
      int numFlip_;

      uhex hex_;
      int* active_;

      // flags
      bool haveActive_;
      bool haveNumFlip_;
      bool isVerified_;
};



#endif
