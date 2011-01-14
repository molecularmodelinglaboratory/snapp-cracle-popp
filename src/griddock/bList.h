#ifndef BLIST_H
#define BLIST_H

#include <bPoints.h>
//~ #include "bPoints.h"

namespace bStd { class bConf; class bList; };

class bStd::bConf {
   typedef unsigned int uint;

   friend class bDocker;
   friend class bMC;
   friend class bList;

   private:
   bPoints* pp_;
   float sc_;
   bConf* next_;
   bConf* prev_;
   bList* pare_;
   bList* chld_;
   static uint cnt_;

   public:
   bConf();
   bConf( const bConf& );
   bConf( const bPoints& );
   bConf( bPoints&, float =(0.0) );
   ~bConf();
   
   void clear();

   bConf& operator=( const bConf& );
   bool   operator<( const bConf& );
   bool   operator>( const bConf& );
   bool   operator==( const bConf& );


   bConf* next();
   bConf* prev();
   bConf* pare();
   bList* home();
   bList* chld();
   bList* new_chld();
   float  score();
   bPoints* pts();
   
   void next( bConf* );
   void prev( bConf* );
   void pare( bList* );
   void chld( bList* );
   void chld( bConf* );
   void score( float );
   
   void pymol( FILE*, char[], char[] );
   void print();
};

class bStd::bList {
   typedef unsigned int uint;

   friend class bConf;

   private:
   uint cap_;
   uint num_;
   bConf* pare_;
   bConf* beg_;
   bConf* end_;
   bConf** list_;
   static float _resizeAmt_;
   
   // pymol
   static int _step_;

   public:
   bList( uint =(10) );
   bList( const bList&, uint=(10) );
   bList( const bConf&, uint=(10) );
   bList( bConf*, uint=(10) );
   ~bList();

   void clear();
   void trim( uint =(10) );

   bConf& operator[]( uint );
   bool   empty();

   bConf* pare();
   bConf* beg();
   bConf* end();
   uint   cap();
   uint   num();
   uint   size();
   
   void resize( uint );

   void pare( bConf* );
   void cap ( uint );

   
   void push_back( bConf& );
   void push_back( bConf* );
   void del( uint );
   
   void insert( bConf& );
   void insert( bConf* );
   
   void sort( bool =(true) );
   void _sort( int, int );
   void _sortRev( int, int );
   void swap( int, int );
   
   void pymol( FILE*, char[] );
   static void pymolGradient( FILE*, int =(20) );
   
   void print();
};



#endif

