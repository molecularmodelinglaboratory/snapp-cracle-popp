#ifndef BSNAPP_H
#define BSNAPP_H


namespace bStd { class bSnapp; };

class bStd::bSnapp {
   typedef unsigned int uint;

   private:
      static bool isSet_;
      static char file_[256];
      static char defdir_[64];
      static char deffile_[64];
      static bool haveFile_;
      static bool findFile( char[] );

   public:
      static float snapp_[26][26][26][26][5];

      bSnapp();
      bSnapp( const bSnapp & );
      ~bSnapp();

      static void readSnapp( const char[] =(file_) );
      static float score( char[], int );
      static float score( int[], int );
   
      static int res2num( char );
};




#endif
