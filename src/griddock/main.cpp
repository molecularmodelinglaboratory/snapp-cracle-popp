
//~ #include <cstdio>
//~ #include <vector>
//~ #include <claVararray>
//~ #include <iostream>
//~ #include <fstream>
//~ #include <string>

#include <stdio.h>
//~ #include <string.h>
//~ #include <sys/stat.h>

#include "bDocker.h"
//~ #include "bSimplex.h"
//~ #include "bAA.h"
//~ #include "bSys.h"
#include "bSort.h"
//~ #include "bTree.h"
//~ #include "bCentroid.h"
#include "bHex.h"

using namespace std;
using namespace bStd;


/*******************************************

   1) Read in protein
      -- recenter: origin in center of protein?
      -- find max (min) for each axis

      _width_i = _length_i + (2*f)/r + (2*t)r

   2) Create exlusion matrix

      _width_em = 1 + (2*f)/r


*******************************************/

int main( int argc, char **argv ) {

	//~ bTree test;
	//~ return;
   
   //~ unsigned long test = 0;
   
   //~ for( int i=0; i < 64; ++i ) { test = i; bHex::printBit( test ); printf(": %lu\n",test); }
    //~ printf("\n");
   
   //~ bHex n(64);
   //~ for( int i=0; i < 64; ++i ) { n.remove(0xFFF); n.assign(i); n.print(); }
   
   //~ test = 2000000;
   //~ printf("here\n");
   //~ for( int i=0; i < 64; ++i ) { ++test; bHex::printHex( test ); printf(": %lu\n",test); }
    //~ printf("\n");
   
   //~ return 0;

   //~ int listarr[2] = {4,1};
   //~ int* list = listarr;
   //~ for( int i=0; i < 2; ++i ) { printf("%2d ", list[i]); } printf("\n");
   //~ bSort::qsort( list, 2 );
   //~ for( int i=0; i < 2; ++i ) { printf("%2d ", list[i]); } printf("\n");
   //~ return 1;
   /* PPI Docker */
   try{ 
      bDocker d;
      d.cla(argc,argv);
      d.doItAll();
      printf("\n\tdone.\n\n");
   }
   catch( const char* e ) { printf("%s\n",e); }

   return 0;
}
