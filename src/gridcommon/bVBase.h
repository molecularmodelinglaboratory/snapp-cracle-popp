
#ifndef BVBASE_H
#define BVBASE_H

#include <stdio.h>
#include <stdlib.h>

   extern const char _HOME_[];
   extern const char _SLASH_;
   extern const short _MAXBIT_;
   extern const double _PI_;

#ifdef _OPENMP
#include <omp.h>
#endif

namespace bStd { class bVBase; };

class bStd::bVBase {
   protected:
      
   typedef unsigned short ushort;
   typedef unsigned int   uint;
   typedef unsigned long  ulong;
   
   public:

   virtual void print( FILE* =(stdout) ) const =0;
   static  void pymolHeader( FILE*, int=0 );
};


#endif
