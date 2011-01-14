#ifndef BAA_H
#define BAA_H

namespace bStd { class bAA; }

class bStd::bAA {
   typedef unsigned short ushort;
   
   private:
      static char res3_[26][4];
      static char res1_[26];
      static char resNS_tb[492][4];
      static ushort resNS_lu[492];
      static ushort resNS_dx[28];

      static bool have3_;
      static bool have1_;
      static bool haveNS_;
      static bool haveAll_;

   public:
      bAA();
      ~bAA();
   
      static void set();
      static void set31();
      static void setNS();
      
      static int   aa2num( char );
      static char* aa1to3( char );
      static char  aa3to1( char* );
      static char  aa3to1( const char* );
   
      static void  print( FILE* =(stdout) );
};





#endif
