#ifndef BSTAMP_H
#define BSTAMP_H

namespace bStd { class bStamp; };

class bStd::bStamp {
   typedef unsigned long ulong;

   // declare friends
   friend class bGrid;
   friend class bDocker;

   public:
   bStamp();
   bStamp( const bStamp& );
   ~bStamp();

   bStamp& operator=( const bStamp& );

   bool initializeStamp( const int, const int );
   void print();
   static float pointDistance( const float*, const float* );

   protected:
   ulong* stamp_;
   ulong** slice_;
   int radius_;
   int size_;

   private:
};


#endif

