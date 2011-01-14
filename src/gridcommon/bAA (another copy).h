#ifndef BAA_H
#define BAA_H

namespace bStd { class bAA; }

class bStd::bAA {
   private:
      static char res3_[20][4];
      static char res1_[21];
   
      static bool have3_;
      static bool have1_;
   
   public:
      bAA();
      ~bAA();
   
      static void set1();
      static void set3();
      
      static int   aa2num( char );
      static char* aa1to3( char );
      static char  aa3to1( char* );
      static char  aa3to1( const char* );
   
      static void print();
};





#endif
