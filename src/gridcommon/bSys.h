#ifndef BSYS_H
#define BSYS_H

#include "bVBase.h"

namespace bStd { class bSys; }

class bStd::bSys : virtual bVBase {

public:

   bSys();
   ~bSys();

   static uint cla( int, char**, char**&, char**& );


   static bool fileExists( char*, char* );
   static bool fileExists( char* );

   static bool fileExists( const char*, const char* );
   static bool fileExists( const char* );

   static bool  dirExists( const char*, const char* );
   static bool  dirExists( const char* );

   static uint fileList( const char*, char**& );

   static bool removeExt( char* );
   static bool multicat( char*&, const char*, const char* =(NULL), const char* =(NULL), const char* =(NULL) );

};





#endif
