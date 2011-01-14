
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include "bAA.h"
using namespace std;
using namespace bStd;

char bAA::res3_[20][4];
char bAA::res1_[21];

bool bAA::have3_ = false;
bool bAA::have1_ = false;

bAA::bAA() {
   
}

bAA::~bAA() {}

void bAA::set3() {
   memmove( res3_[0x00], "ALA\0", 4 ); memmove( res3_[0x01], "CYS\0", 4 );
   memmove( res3_[0x02], "ASP\0", 4 ); memmove( res3_[0x03], "GLU\0", 4 );
   memmove( res3_[0x04], "PHE\0", 4 ); memmove( res3_[0x05], "GLY\0", 4 );
   memmove( res3_[0x06], "HIS\0", 4 ); memmove( res3_[0x07], "ILE\0", 4 );
   memmove( res3_[0x08], "LYS\0", 4 ); memmove( res3_[0x09], "LEU\0", 4 );
   memmove( res3_[0x0A], "MET\0", 4 ); memmove( res3_[0x0B], "ASN\0", 4 );
   memmove( res3_[0x0C], "PRO\0", 4 ); memmove( res3_[0x0D], "GLN\0", 4 );
   memmove( res3_[0x0E], "ARG\0", 4 ); memmove( res3_[0x0F], "SER\0", 4 );
   memmove( res3_[0x10], "THR\0", 4 ); memmove( res3_[0x11], "VAL\0", 4 );
   memmove( res3_[0x12], "TRP\0", 4 ); memmove( res3_[0x13], "TYR\0", 4 );

   have3_ = true;
   return;
}

void bAA::set1() {
   res1_[0x00] = 'A'; res1_[0x01] = 'C';
   res1_[0x02] = 'D'; res1_[0x03] = 'E';
   res1_[0x04] = 'F'; res1_[0x05] = 'G';
   res1_[0x06] = 'H'; res1_[0x07] = 'I';
   res1_[0x08] = 'K'; res1_[0x09] = 'L';
   res1_[0x0A] = 'M'; res1_[0x0B] = 'N';
   res1_[0x0C] = 'P'; res1_[0x0D] = 'Q';
   res1_[0x0E] = 'R'; res1_[0x0F] = 'S';
   res1_[0x10] = 'T'; res1_[0x11] = 'V';
   res1_[0x12] = 'W'; res1_[0x13] = 'Y';
   res1_[0x14] = '\0';
   have1_ = true;
   return;
}

char bAA::aa3to1( char* r ) { const char s[4] = { r[0], r[1], r[2], r[3] }; return aa3to1( s ); }
char bAA::aa3to1( const char* r ) {
   char one = '\0';
   char res[4];
   res[0] = ( r[0] < 65 ) ? toupper( r[0] ) : r[0];
   res[1] = ( r[1] < 65 ) ? toupper( r[1] ) : r[1];
   res[2] = ( r[2] < 65 ) ? toupper( r[2] ) : r[2];
   res[3] = '\0';
   switch( res[0] ) {
      case 'A': 
         switch( res[1] ) {
            case 'L': /* Alanine */
               one = 'A'; break;
            case 'R': /* Arginine */
               one = 'R'; break;
            case 'S':
               switch( res[2] ) {
                  case 'N': /* Asparagine */
                     one = 'N'; break;
                  case 'P': /* Aspartic Acid */
                     one = 'D'; break;
                  default: break;
               }
            default: break;
         }
         break;
      case 'C': /* Cysteine */
         one = 'C'; break;
      case 'G':
         switch( res[2] ) {
            case 'N': /* Glutamine */
               one = 'Q'; break;
            case 'U': /* Glutamic Acid */
               one = 'E'; break;
            case 'Y': /* Glycine */
               one = 'G'; break;
            default: break;
         }
         break;
      case 'H': /* Histidine */
         one = 'H'; break;
      case 'I': /* Isoleucine */
         one = 'I'; break;
      case 'L':
         switch( res[1] ) {
            case 'E': /* Leucine */
               one = 'L'; break;
            case 'Y': /* Lysine */
               one = 'K'; break;
            default: break;
         }
         break;
      case 'M': /* Methionine */
         one = 'M'; break;
      case 'P':
         switch( res[1] ) {
            case 'H': /* Phenylalanine */
               one = 'F'; break;
            case 'R': /* Proline */
               one = 'P'; break;
            default: break;
         }
         break;
      case 'S': /* Serine */
         one = 'S'; break;
      case 'T':
         switch( res[1] ) {
            case 'H': /* Threonine */
               one = 'T'; break;
            case 'R': /* Tryptophan */
               one = 'W'; break;
            case 'Y': /* Tyrosine */
               one = 'Y'; break;
            default: break;
         }
         break;
      case 'V': /* Valine */
         one = 'V'; break;
      default: break;
      
   }
   return one;
}

char* bAA::aa1to3( char res ) {
   if( !have3_ ) { set3(); }
   int ref = aa2num( res );
   return res3_[ref];

   if( res < 65 ) { res = toupper( res ); }
   int num = res;
        if( num > 'X' ) { num -= 5; }
   else if( num > 'U' ) { num -= 4; }
   else if( num > 'O' ) { num -= 3; }
   else if( num > 'J' ) { num -= 2; }
   else if( num > 'B' ) { --num; }
   else {}
   num -= 66;
   return res3_[num];
}

int bAA::aa2num( char res ) {
   if( !have1_ ) { set1(); }
   char* pc = (char*) memchr( res1_, res, 20 );
   return pc - res1_;
}

void bAA::print() {
   printf("[bAA]\n\t");
   for( int i=0; i < 26; ++i ) {
      printf("[%2d] %c:%s   ", i, res1_[i], res3_[i]);
      if( (i % 5) == 0 ) { printf("\n\t"); }
   }
   printf("\n");
   return;
}

