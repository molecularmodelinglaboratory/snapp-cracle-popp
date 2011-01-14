
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include "bAA.h"
using namespace std;
using namespace bStd;

char bAA::res3_[26][4];
char bAA::res1_[26];

char bAA::resNS_tb[492][4];
bAA::ushort bAA::resNS_lu[492];
bAA::ushort bAA::resNS_dx[28];

bool bAA::have3_ = false;
bool bAA::have1_ = false;
bool bAA::haveNS_ = false;
bool bAA::haveAll_ = false;

bAA::bAA()
{
}

bAA::~bAA()
{
}

void bAA::set() {
   if( !haveAll_ ) {
      if( !have3_ || !have1_ ) set31();
      if( !haveNS_) setNS();
   }
   haveAll_ = true;
   return;
}

void bAA::set31() {
   memmove( res3_[0x00], "ALA\0", 4 );   res1_[0x00] = 'A';
   memmove( res3_[0x01], "ASX\0", 4 );   res1_[0x01] = 'B'; // *
   memmove( res3_[0x02], "CYS\0", 4 );   res1_[0x02] = 'C';
   memmove( res3_[0x03], "ASP\0", 4 );   res1_[0x03] = 'D';
   memmove( res3_[0x04], "GLU\0", 4 );   res1_[0x04] = 'E';
   memmove( res3_[0x05], "PHE\0", 4 );   res1_[0x05] = 'F';
   memmove( res3_[0x06], "GLY\0", 4 );   res1_[0x06] = 'G';
   memmove( res3_[0x07], "HIS\0", 4 );   res1_[0x07] = 'H';
   memmove( res3_[0x08], "ILE\0", 4 );   res1_[0x08] = 'I';
   memmove( res3_[0x09], "XLE\0", 4 );   res1_[0x09] = 'J'; // *
   memmove( res3_[0x0A], "LYS\0", 4 );   res1_[0x0A] = 'K';
   memmove( res3_[0x0B], "LEU\0", 4 );   res1_[0x0B] = 'L';
   memmove( res3_[0x0C], "MET\0", 4 );   res1_[0x0C] = 'M';
   memmove( res3_[0x0D], "ASN\0", 4 );   res1_[0x0D] = 'N';
   memmove( res3_[0x0E], "PYL\0", 4 );   res1_[0x0E] = 'O'; // *
   memmove( res3_[0x0F], "PRO\0", 4 );   res1_[0x0F] = 'P';
   memmove( res3_[0x10], "GLN\0", 4 );   res1_[0x10] = 'Q';
   memmove( res3_[0x11], "ARG\0", 4 );   res1_[0x11] = 'R';
   memmove( res3_[0x12], "SER\0", 4 );   res1_[0x12] = 'S';
   memmove( res3_[0x13], "THR\0", 4 );   res1_[0x13] = 'T';
   memmove( res3_[0x14], "SEC\0", 4 );   res1_[0x14] = 'U'; // *
   memmove( res3_[0x15], "VAL\0", 4 );   res1_[0x15] = 'V';
   memmove( res3_[0x16], "TRP\0", 4 );   res1_[0x16] = 'W';
   memmove( res3_[0x17], "XAA\0", 4 );   res1_[0x17] = 'X'; // *
   memmove( res3_[0x18], "TYR\0", 4 );   res1_[0x18] = 'Y';
   memmove( res3_[0x19], "GLX\0", 4 );   res1_[0x19] = 'Z'; // *

   have3_ = true;
   have1_ = true;
   return;
}

char bAA::aa3to1( char* r ) {
   const char s[4] = { r[0], r[1], r[2], r[3] }; return aa3to1( s );
}
char bAA::aa3to1( const char* r ) {
   if( !haveAll_ ) { set(); }
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
                  case 'X': /* One of the prev two */
                     one = 'B';
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
            case 'X': /* One of the first two */
               one = 'Z'; break;
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
            case 'Y': /* ?? */
               one = 'O'; break;
            default: break;
         }
         break;
      case 'S': /* Serine */
         switch( res[2] ) {
            case 'R': /* Serine */
               one = 'S'; break;
            case 'C': /* Selenocysteine? */
               one = 'U'; break;
            default: break;
         }
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
      case 'X':
         switch( res[1] ) {
            case 'L': /* ILE or LEU */
               one = 'J'; break;
            case 'A': /* Any */
               one = 'X'; break;
            default: break;
         }
         break;
      default: break;
   }

   if( one == '\0' ) {
      int beg = aa2num( res[0] );
      int end = 0;
      if( beg < 0 ) {
         beg = 0;
         end = resNS_dx[ 0 ];
      }
      else {
         end = beg + 1;
         beg = resNS_dx[ beg ];
         end = resNS_dx[ end ];
      }
      for( int i=beg; i < end; ++i ) {
         if( memcmp( res, resNS_tb[i], 3 ) == 0 ) {
            one = res1_[ resNS_lu[i] ];
            break;
         }
      }
      if( one == '\0' ) {
         beg = resNS_dx[26];
         end = resNS_dx[27];
         for( int i=beg; i < end; ++i ) {
            if( memcmp( res, resNS_tb[i], 3 ) == 0 ) {
               one = res1_[ resNS_lu[i] ];
               break;
            }
         }
      }
      if( one == '\0' ) { printf("**************************** Unknown amino acid: %s\n", res); }
   }

   return one;
}

char* bAA::aa1to3( char res ) {
  if( !have3_ ) { set(); }
   int ref = aa2num( res );
   return res3_[ref];

   if( res < 65 ) { res = toupper( res ); }
   int num = res;
   num -= 66;
  return res3_[num];
}

int bAA::aa2num( char res ) {
   if( !have1_ ) { set(); }
   char* pc = (char*) memchr( res1_, res, 26 );
   if( pc == NULL ) { *pc = '1'; }
   return pc - res1_;
}

void bAA::setNS() {
   memmove( resNS_tb[0], "033\0", 4 ); resNS_lu[0] = 0x15;
   memmove( resNS_tb[1], "0A0\0", 4 ); resNS_lu[1] = 0x03;
   memmove( resNS_tb[2], "0A1\0", 4 ); resNS_lu[2] = 0x18;
   memmove( resNS_tb[3], "0A2\0", 4 ); resNS_lu[3] = 0x0A;
   memmove( resNS_tb[4], "0A5\0", 4 ); resNS_lu[4] = 0x0D;
   memmove( resNS_tb[5], "0A8\0", 4 ); resNS_lu[5] = 0x02;
   memmove( resNS_tb[6], "0A9\0", 4 ); resNS_lu[6] = 0x05;
   memmove( resNS_tb[7], "0AA\0", 4 ); resNS_lu[7] = 0x15;
   memmove( resNS_tb[8], "0AB\0", 4 ); resNS_lu[8] = 0x15;
   memmove( resNS_tb[9], "0AC\0", 4 ); resNS_lu[9] = 0x06;
   memmove( resNS_tb[10], "0AF\0", 4 ); resNS_lu[10] = 0x16;
   memmove( resNS_tb[11], "0AG\0", 4 ); resNS_lu[11] = 0x0B;
   memmove( resNS_tb[12], "0AH\0", 4 ); resNS_lu[12] = 0x12;
   memmove( resNS_tb[13], "0AK\0", 4 ); resNS_lu[13] = 0x03;
   memmove( resNS_tb[14], "0AY\0", 4 ); resNS_lu[14] = 0x0A;
   memmove( resNS_tb[15], "0AZ\0", 4 ); resNS_lu[15] = 0x0F;
   memmove( resNS_tb[16], "0CS\0", 4 ); resNS_lu[16] = 0x00;
   memmove( resNS_tb[17], "143\0", 4 ); resNS_lu[17] = 0x02;
   memmove( resNS_tb[18], "1LU\0", 4 ); resNS_lu[18] = 0x0B;
   memmove( resNS_tb[19], "1PA\0", 4 ); resNS_lu[19] = 0x05;
   memmove( resNS_tb[20], "1TQ\0", 4 ); resNS_lu[20] = 0x16;
   memmove( resNS_tb[21], "1TY\0", 4 ); resNS_lu[21] = 0x18;
   memmove( resNS_tb[22], "200\0", 4 ); resNS_lu[22] = 0x05;
   memmove( resNS_tb[23], "23F\0", 4 ); resNS_lu[23] = 0x05;
   memmove( resNS_tb[24], "2AS\0", 4 ); resNS_lu[24] = 0x03;
   memmove( resNS_tb[25], "2FM\0", 4 ); resNS_lu[25] = 0x0C;
   memmove( resNS_tb[26], "2LU\0", 4 ); resNS_lu[26] = 0x0B;
   memmove( resNS_tb[27], "2ML\0", 4 ); resNS_lu[27] = 0x0B;
   memmove( resNS_tb[28], "2MR\0", 4 ); resNS_lu[28] = 0x11;
   memmove( resNS_tb[29], "2MT\0", 4 ); resNS_lu[29] = 0x0F;
   memmove( resNS_tb[30], "2TY\0", 4 ); resNS_lu[30] = 0x18;
   memmove( resNS_tb[31], "2VA\0", 4 ); resNS_lu[31] = 0x15;
   memmove( resNS_tb[32], "3AH\0", 4 ); resNS_lu[32] = 0x07;
   memmove( resNS_tb[33], "3MD\0", 4 ); resNS_lu[33] = 0x03;
   memmove( resNS_tb[34], "4BF\0", 4 ); resNS_lu[34] = 0x18;
   memmove( resNS_tb[35], "4DP\0", 4 ); resNS_lu[35] = 0x16;
   memmove( resNS_tb[36], "4FB\0", 4 ); resNS_lu[36] = 0x0F;
   memmove( resNS_tb[37], "4FW\0", 4 ); resNS_lu[37] = 0x16;
   memmove( resNS_tb[38], "4HT\0", 4 ); resNS_lu[38] = 0x16;
   memmove( resNS_tb[39], "5CS\0", 4 ); resNS_lu[39] = 0x02;
   memmove( resNS_tb[40], "5HP\0", 4 ); resNS_lu[40] = 0x04;
   memmove( resNS_tb[41], "6CL\0", 4 ); resNS_lu[41] = 0x0A;
   memmove( resNS_tb[42], "6CW\0", 4 ); resNS_lu[42] = 0x16;
   resNS_dx[0x00] = 43;
   memmove( resNS_tb[43], "AA3\0", 4 ); resNS_lu[43] = 0x00;
   memmove( resNS_tb[44], "AA4\0", 4 ); resNS_lu[44] = 0x00;
   memmove( resNS_tb[45], "AAR\0", 4 ); resNS_lu[45] = 0x11;
   memmove( resNS_tb[46], "ABA\0", 4 ); resNS_lu[46] = 0x00;
   memmove( resNS_tb[47], "ACB\0", 4 ); resNS_lu[47] = 0x03;
   memmove( resNS_tb[48], "ACL\0", 4 ); resNS_lu[48] = 0x11;
   memmove( resNS_tb[49], "AEI\0", 4 ); resNS_lu[49] = 0x03;
   memmove( resNS_tb[50], "AFA\0", 4 ); resNS_lu[50] = 0x0D;
   memmove( resNS_tb[51], "AGM\0", 4 ); resNS_lu[51] = 0x11;
   memmove( resNS_tb[52], "AHB\0", 4 ); resNS_lu[52] = 0x0D;
   memmove( resNS_tb[53], "AHL\0", 4 ); resNS_lu[53] = 0x11;
   memmove( resNS_tb[54], "AHO\0", 4 ); resNS_lu[54] = 0x00;
   memmove( resNS_tb[55], "AHP\0", 4 ); resNS_lu[55] = 0x00;
   memmove( resNS_tb[56], "AIB\0", 4 ); resNS_lu[56] = 0x00;
   memmove( resNS_tb[57], "AKL\0", 4 ); resNS_lu[57] = 0x03;
   memmove( resNS_tb[58], "ALC\0", 4 ); resNS_lu[58] = 0x00;
   memmove( resNS_tb[59], "ALG\0", 4 ); resNS_lu[59] = 0x11;
   memmove( resNS_tb[60], "ALM\0", 4 ); resNS_lu[60] = 0x00;
   memmove( resNS_tb[61], "ALN\0", 4 ); resNS_lu[61] = 0x00;
   memmove( resNS_tb[62], "ALO\0", 4 ); resNS_lu[62] = 0x13;
   memmove( resNS_tb[63], "ALS\0", 4 ); resNS_lu[63] = 0x00;
   memmove( resNS_tb[64], "ALT\0", 4 ); resNS_lu[64] = 0x00;
   memmove( resNS_tb[65], "ALY\0", 4 ); resNS_lu[65] = 0x0A;
   memmove( resNS_tb[66], "APH\0", 4 ); resNS_lu[66] = 0x00;
   memmove( resNS_tb[67], "API\0", 4 ); resNS_lu[67] = 0x0A;
   memmove( resNS_tb[68], "APK\0", 4 ); resNS_lu[68] = 0x0A;
   memmove( resNS_tb[69], "AR4\0", 4 ); resNS_lu[69] = 0x04;
   memmove( resNS_tb[70], "ARM\0", 4 ); resNS_lu[70] = 0x11;
   memmove( resNS_tb[71], "ARO\0", 4 ); resNS_lu[71] = 0x11;
   memmove( resNS_tb[72], "AS2\0", 4 ); resNS_lu[72] = 0x03;
   memmove( resNS_tb[73], "ASA\0", 4 ); resNS_lu[73] = 0x03;
   memmove( resNS_tb[74], "ASB\0", 4 ); resNS_lu[74] = 0x03;
   memmove( resNS_tb[75], "ASI\0", 4 ); resNS_lu[75] = 0x03;
   memmove( resNS_tb[76], "ASK\0", 4 ); resNS_lu[76] = 0x03;
   memmove( resNS_tb[77], "ASL\0", 4 ); resNS_lu[77] = 0x03;
   memmove( resNS_tb[78], "ASQ\0", 4 ); resNS_lu[78] = 0x03;
   memmove( resNS_tb[79], "AYA\0", 4 ); resNS_lu[79] = 0x00;
   memmove( resNS_tb[80], "AZK\0", 4 ); resNS_lu[80] = 0x0A;
   resNS_dx[0x01] = 81;
   memmove( resNS_tb[81], "B1F\0", 4 ); resNS_lu[81] = 0x05;
   memmove( resNS_tb[82], "B2A\0", 4 ); resNS_lu[82] = 0x00;
   memmove( resNS_tb[83], "B2F\0", 4 ); resNS_lu[83] = 0x05;
   memmove( resNS_tb[84], "B2I\0", 4 ); resNS_lu[84] = 0x08;
   memmove( resNS_tb[85], "B2V\0", 4 ); resNS_lu[85] = 0x15;
   memmove( resNS_tb[86], "B3A\0", 4 ); resNS_lu[86] = 0x00;
   memmove( resNS_tb[87], "B3D\0", 4 ); resNS_lu[87] = 0x03;
   memmove( resNS_tb[88], "B3E\0", 4 ); resNS_lu[88] = 0x04;
   memmove( resNS_tb[89], "B3K\0", 4 ); resNS_lu[89] = 0x0A;
   memmove( resNS_tb[90], "B3S\0", 4 ); resNS_lu[90] = 0x12;
   memmove( resNS_tb[91], "B3X\0", 4 ); resNS_lu[91] = 0x0D;
   memmove( resNS_tb[92], "B3Y\0", 4 ); resNS_lu[92] = 0x18;
   memmove( resNS_tb[93], "BAL\0", 4 ); resNS_lu[93] = 0x00;
   memmove( resNS_tb[94], "BBC\0", 4 ); resNS_lu[94] = 0x02;
   memmove( resNS_tb[95], "BCS\0", 4 ); resNS_lu[95] = 0x02;
   memmove( resNS_tb[96], "BCX\0", 4 ); resNS_lu[96] = 0x02;
   memmove( resNS_tb[97], "BFD\0", 4 ); resNS_lu[97] = 0x03;
   memmove( resNS_tb[98], "BHD\0", 4 ); resNS_lu[98] = 0x03;
   memmove( resNS_tb[99], "BIF\0", 4 ); resNS_lu[99] = 0x05;
   memmove( resNS_tb[100], "BLE\0", 4 ); resNS_lu[100] = 0x0B;
   memmove( resNS_tb[101], "BLY\0", 4 ); resNS_lu[101] = 0x0A;
   memmove( resNS_tb[102], "BMT\0", 4 ); resNS_lu[102] = 0x13;
   memmove( resNS_tb[103], "BNN\0", 4 ); resNS_lu[103] = 0x00;
   memmove( resNS_tb[104], "BOR\0", 4 ); resNS_lu[104] = 0x11;
   memmove( resNS_tb[105], "BPE\0", 4 ); resNS_lu[105] = 0x02;
   memmove( resNS_tb[106], "BSE\0", 4 ); resNS_lu[106] = 0x12;
   memmove( resNS_tb[107], "BTA\0", 4 ); resNS_lu[107] = 0x0B;
   memmove( resNS_tb[108], "BTC\0", 4 ); resNS_lu[108] = 0x02;
   memmove( resNS_tb[109], "BTR\0", 4 ); resNS_lu[109] = 0x16;
   memmove( resNS_tb[110], "BUC\0", 4 ); resNS_lu[110] = 0x02;
   memmove( resNS_tb[111], "BUG\0", 4 ); resNS_lu[111] = 0x0B;
   resNS_dx[0x02] = 112;
   memmove( resNS_tb[112], "C1X\0", 4 ); resNS_lu[112] = 0x0A;
   memmove( resNS_tb[113], "C3Y\0", 4 ); resNS_lu[113] = 0x02;
   memmove( resNS_tb[114], "C5C\0", 4 ); resNS_lu[114] = 0x02;
   memmove( resNS_tb[115], "C6C\0", 4 ); resNS_lu[115] = 0x02;
   memmove( resNS_tb[116], "CAB\0", 4 ); resNS_lu[116] = 0x00;
   memmove( resNS_tb[117], "CAF\0", 4 ); resNS_lu[117] = 0x02;
   memmove( resNS_tb[118], "CAS\0", 4 ); resNS_lu[118] = 0x02;
   memmove( resNS_tb[119], "CAY\0", 4 ); resNS_lu[119] = 0x02;
   memmove( resNS_tb[120], "CCL\0", 4 ); resNS_lu[120] = 0x0A;
   memmove( resNS_tb[121], "CCS\0", 4 ); resNS_lu[121] = 0x02;
   memmove( resNS_tb[122], "CEA\0", 4 ); resNS_lu[122] = 0x02;
   memmove( resNS_tb[123], "CGU\0", 4 ); resNS_lu[123] = 0x04;
   memmove( resNS_tb[124], "CHG\0", 4 ); resNS_lu[124] = 0x00;
   memmove( resNS_tb[125], "CHP\0", 4 ); resNS_lu[125] = 0x06;
   memmove( resNS_tb[126], "CIR\0", 4 ); resNS_lu[126] = 0x11;
   memmove( resNS_tb[127], "CLB\0", 4 ); resNS_lu[127] = 0x00;
   memmove( resNS_tb[128], "CLD\0", 4 ); resNS_lu[128] = 0x00;
   memmove( resNS_tb[129], "CLE\0", 4 ); resNS_lu[129] = 0x0B;
   memmove( resNS_tb[130], "CLG\0", 4 ); resNS_lu[130] = 0x0A;
   memmove( resNS_tb[131], "CLH\0", 4 ); resNS_lu[131] = 0x0A;
   memmove( resNS_tb[132], "CME\0", 4 ); resNS_lu[132] = 0x02;
   memmove( resNS_tb[133], "CMH\0", 4 ); resNS_lu[133] = 0x02;
   memmove( resNS_tb[134], "CML\0", 4 ); resNS_lu[134] = 0x02;
   memmove( resNS_tb[135], "CMT\0", 4 ); resNS_lu[135] = 0x02;
   memmove( resNS_tb[136], "CR5\0", 4 ); resNS_lu[136] = 0x06;
   memmove( resNS_tb[137], "CRU\0", 4 ); resNS_lu[137] = 0x04;
   memmove( resNS_tb[138], "CS1\0", 4 ); resNS_lu[138] = 0x02;
   memmove( resNS_tb[139], "CS3\0", 4 ); resNS_lu[139] = 0x02;
   memmove( resNS_tb[140], "CS4\0", 4 ); resNS_lu[140] = 0x02;
   memmove( resNS_tb[141], "CSA\0", 4 ); resNS_lu[141] = 0x02;
   memmove( resNS_tb[142], "CSB\0", 4 ); resNS_lu[142] = 0x02;
   memmove( resNS_tb[143], "CSD\0", 4 ); resNS_lu[143] = 0x02;
   memmove( resNS_tb[144], "CSE\0", 4 ); resNS_lu[144] = 0x02;
   memmove( resNS_tb[145], "CSI\0", 4 ); resNS_lu[145] = 0x06;
   memmove( resNS_tb[146], "CSO\0", 4 ); resNS_lu[146] = 0x02;
   memmove( resNS_tb[147], "CSP\0", 4 ); resNS_lu[147] = 0x02;
   memmove( resNS_tb[148], "CSR\0", 4 ); resNS_lu[148] = 0x02;
   memmove( resNS_tb[149], "CSS\0", 4 ); resNS_lu[149] = 0x02;
   memmove( resNS_tb[150], "CSU\0", 4 ); resNS_lu[150] = 0x02;
   memmove( resNS_tb[151], "CSW\0", 4 ); resNS_lu[151] = 0x02;
   memmove( resNS_tb[152], "CSX\0", 4 ); resNS_lu[152] = 0x02;
   memmove( resNS_tb[153], "CSZ\0", 4 ); resNS_lu[153] = 0x02;
   memmove( resNS_tb[154], "CTH\0", 4 ); resNS_lu[154] = 0x13;
   memmove( resNS_tb[155], "CWR\0", 4 ); resNS_lu[155] = 0x12;
   memmove( resNS_tb[156], "CXM\0", 4 ); resNS_lu[156] = 0x0C;
   memmove( resNS_tb[157], "CY0\0", 4 ); resNS_lu[157] = 0x02;
   memmove( resNS_tb[158], "CY1\0", 4 ); resNS_lu[158] = 0x02;
   memmove( resNS_tb[159], "CY3\0", 4 ); resNS_lu[159] = 0x02;
   memmove( resNS_tb[160], "CY4\0", 4 ); resNS_lu[160] = 0x02;
   memmove( resNS_tb[161], "CYA\0", 4 ); resNS_lu[161] = 0x02;
   memmove( resNS_tb[162], "CYD\0", 4 ); resNS_lu[162] = 0x02;
   memmove( resNS_tb[163], "CYF\0", 4 ); resNS_lu[163] = 0x02;
   memmove( resNS_tb[164], "CYG\0", 4 ); resNS_lu[164] = 0x02;
   memmove( resNS_tb[165], "CYM\0", 4 ); resNS_lu[165] = 0x02;
   memmove( resNS_tb[166], "CYQ\0", 4 ); resNS_lu[166] = 0x02;
   memmove( resNS_tb[167], "CYR\0", 4 ); resNS_lu[167] = 0x02;
   memmove( resNS_tb[168], "CZ2\0", 4 ); resNS_lu[168] = 0x02;
   memmove( resNS_tb[169], "CZZ\0", 4 ); resNS_lu[169] = 0x02;
   resNS_dx[0x03] = 170;
   memmove( resNS_tb[170], "DAB\0", 4 ); resNS_lu[170] = 0x00;
   memmove( resNS_tb[171], "DAH\0", 4 ); resNS_lu[171] = 0x05;
   memmove( resNS_tb[172], "DAL\0", 4 ); resNS_lu[172] = 0x00;
   memmove( resNS_tb[173], "DAR\0", 4 ); resNS_lu[173] = 0x11;
   memmove( resNS_tb[174], "DAS\0", 4 ); resNS_lu[174] = 0x03;
   memmove( resNS_tb[175], "DBS\0", 4 ); resNS_lu[175] = 0x12;
   memmove( resNS_tb[176], "DBU\0", 4 ); resNS_lu[176] = 0x00;
   memmove( resNS_tb[177], "DBY\0", 4 ); resNS_lu[177] = 0x18;
   memmove( resNS_tb[178], "DBZ\0", 4 ); resNS_lu[178] = 0x00;
   memmove( resNS_tb[179], "DCY\0", 4 ); resNS_lu[179] = 0x02;
   memmove( resNS_tb[180], "DDE\0", 4 ); resNS_lu[180] = 0x07;
   memmove( resNS_tb[181], "DGL\0", 4 ); resNS_lu[181] = 0x04;
   memmove( resNS_tb[182], "DGN\0", 4 ); resNS_lu[182] = 0x10;
   memmove( resNS_tb[183], "DHA\0", 4 ); resNS_lu[183] = 0x00;
   memmove( resNS_tb[184], "DHI\0", 4 ); resNS_lu[184] = 0x07;
   memmove( resNS_tb[185], "DHN\0", 4 ); resNS_lu[185] = 0x15;
   memmove( resNS_tb[186], "DIL\0", 4 ); resNS_lu[186] = 0x08;
   memmove( resNS_tb[187], "DIR\0", 4 ); resNS_lu[187] = 0x11;
   memmove( resNS_tb[188], "DIV\0", 4 ); resNS_lu[188] = 0x15;
   memmove( resNS_tb[189], "DLE\0", 4 ); resNS_lu[189] = 0x0B;
   memmove( resNS_tb[190], "DLS\0", 4 ); resNS_lu[190] = 0x0A;
   memmove( resNS_tb[191], "DLY\0", 4 ); resNS_lu[191] = 0x0A;
   memmove( resNS_tb[192], "DMH\0", 4 ); resNS_lu[192] = 0x0D;
   memmove( resNS_tb[193], "DMK\0", 4 ); resNS_lu[193] = 0x03;
   memmove( resNS_tb[194], "DNE\0", 4 ); resNS_lu[194] = 0x0B;
   memmove( resNS_tb[195], "DNG\0", 4 ); resNS_lu[195] = 0x0B;
   memmove( resNS_tb[196], "DNL\0", 4 ); resNS_lu[196] = 0x0A;
   memmove( resNS_tb[197], "DNM\0", 4 ); resNS_lu[197] = 0x0B;
   memmove( resNS_tb[198], "DNP\0", 4 ); resNS_lu[198] = 0x00;
   memmove( resNS_tb[199], "DNS\0", 4 ); resNS_lu[199] = 0x0A;
   memmove( resNS_tb[200], "DOH\0", 4 ); resNS_lu[200] = 0x03;
   memmove( resNS_tb[201], "DON\0", 4 ); resNS_lu[201] = 0x0B;
   memmove( resNS_tb[202], "DP1\0", 4 ); resNS_lu[202] = 0x11;
   memmove( resNS_tb[203], "DPH\0", 4 ); resNS_lu[203] = 0x05;
   memmove( resNS_tb[204], "DPL\0", 4 ); resNS_lu[204] = 0x0F;
   memmove( resNS_tb[205], "DPN\0", 4 ); resNS_lu[205] = 0x05;
   memmove( resNS_tb[206], "DPP\0", 4 ); resNS_lu[206] = 0x00;
   memmove( resNS_tb[207], "DPR\0", 4 ); resNS_lu[207] = 0x0F;
   memmove( resNS_tb[208], "DSE\0", 4 ); resNS_lu[208] = 0x12;
   memmove( resNS_tb[209], "DSG\0", 4 ); resNS_lu[209] = 0x0D;
   memmove( resNS_tb[210], "DSN\0", 4 ); resNS_lu[210] = 0x12;
   memmove( resNS_tb[211], "DSP\0", 4 ); resNS_lu[211] = 0x03;
   memmove( resNS_tb[212], "DTH\0", 4 ); resNS_lu[212] = 0x13;
   memmove( resNS_tb[213], "DTR\0", 4 ); resNS_lu[213] = 0x16;
   memmove( resNS_tb[214], "DTY\0", 4 ); resNS_lu[214] = 0x18;
   memmove( resNS_tb[215], "DVA\0", 4 ); resNS_lu[215] = 0x15;
   resNS_dx[0x04] = 216;
   memmove( resNS_tb[216], "EFC\0", 4 ); resNS_lu[216] = 0x02;
   memmove( resNS_tb[217], "EHP\0", 4 ); resNS_lu[217] = 0x05;
   memmove( resNS_tb[218], "EPM\0", 4 ); resNS_lu[218] = 0x0C;
   memmove( resNS_tb[219], "ESC\0", 4 ); resNS_lu[219] = 0x0C;
   resNS_dx[0x05] = 220;
   memmove( resNS_tb[220], "FCL\0", 4 ); resNS_lu[220] = 0x05;
   memmove( resNS_tb[221], "FGL\0", 4 ); resNS_lu[221] = 0x06;
   memmove( resNS_tb[222], "FGP\0", 4 ); resNS_lu[222] = 0x12;
   memmove( resNS_tb[223], "FLA\0", 4 ); resNS_lu[223] = 0x00;
   memmove( resNS_tb[224], "FLE\0", 4 ); resNS_lu[224] = 0x0B;
   memmove( resNS_tb[225], "FLT\0", 4 ); resNS_lu[225] = 0x18;
   memmove( resNS_tb[226], "FME\0", 4 ); resNS_lu[226] = 0x0C;
   memmove( resNS_tb[227], "FOE\0", 4 ); resNS_lu[227] = 0x02;
   memmove( resNS_tb[228], "FOG\0", 4 ); resNS_lu[228] = 0x05;
   memmove( resNS_tb[229], "FPA\0", 4 ); resNS_lu[229] = 0x05;
   memmove( resNS_tb[230], "FT6\0", 4 ); resNS_lu[230] = 0x16;
   memmove( resNS_tb[231], "FTR\0", 4 ); resNS_lu[231] = 0x16;
   memmove( resNS_tb[232], "FTY\0", 4 ); resNS_lu[232] = 0x18;
   resNS_dx[0x06] = 233;
   memmove( resNS_tb[233], "G3A\0", 4 ); resNS_lu[233] = 0x00;
   memmove( resNS_tb[234], "GAU\0", 4 ); resNS_lu[234] = 0x04;
   memmove( resNS_tb[235], "GGL\0", 4 ); resNS_lu[235] = 0x04;
   memmove( resNS_tb[236], "GHG\0", 4 ); resNS_lu[236] = 0x10;
   memmove( resNS_tb[237], "GHP\0", 4 ); resNS_lu[237] = 0x06;
   memmove( resNS_tb[238], "GL3\0", 4 ); resNS_lu[238] = 0x06;
   memmove( resNS_tb[239], "GLH\0", 4 ); resNS_lu[239] = 0x10;
   memmove( resNS_tb[240], "GLQ\0", 4 ); resNS_lu[240] = 0x04;
   memmove( resNS_tb[241], "GLZ\0", 4 ); resNS_lu[241] = 0x06;
   memmove( resNS_tb[242], "GMA\0", 4 ); resNS_lu[242] = 0x04;
   memmove( resNS_tb[243], "GPL\0", 4 ); resNS_lu[243] = 0x0A;
   memmove( resNS_tb[244], "GSC\0", 4 ); resNS_lu[244] = 0x06;
   memmove( resNS_tb[245], "GSU\0", 4 ); resNS_lu[245] = 0x04;
   memmove( resNS_tb[246], "GT9\0", 4 ); resNS_lu[246] = 0x02;
   resNS_dx[0x07] = 247;
   memmove( resNS_tb[247], "H1D\0", 4 ); resNS_lu[247] = 0x0C;
   memmove( resNS_tb[248], "H5M\0", 4 ); resNS_lu[248] = 0x0F;
   memmove( resNS_tb[249], "HAC\0", 4 ); resNS_lu[249] = 0x00;
   memmove( resNS_tb[250], "HAR\0", 4 ); resNS_lu[250] = 0x11;
   memmove( resNS_tb[251], "HBN\0", 4 ); resNS_lu[251] = 0x07;
   memmove( resNS_tb[252], "HIA\0", 4 ); resNS_lu[252] = 0x07;
   memmove( resNS_tb[253], "HIC\0", 4 ); resNS_lu[253] = 0x07;
   memmove( resNS_tb[254], "HIP\0", 4 ); resNS_lu[254] = 0x07;
   memmove( resNS_tb[255], "HIQ\0", 4 ); resNS_lu[255] = 0x07;
   memmove( resNS_tb[256], "HLU\0", 4 ); resNS_lu[256] = 0x0B;
   memmove( resNS_tb[257], "HMF\0", 4 ); resNS_lu[257] = 0x00;
   memmove( resNS_tb[258], "HMR\0", 4 ); resNS_lu[258] = 0x11;
   memmove( resNS_tb[259], "HPC\0", 4 ); resNS_lu[259] = 0x05;
   memmove( resNS_tb[260], "HPE\0", 4 ); resNS_lu[260] = 0x05;
   memmove( resNS_tb[261], "HPQ\0", 4 ); resNS_lu[261] = 0x05;
   memmove( resNS_tb[262], "HRG\0", 4 ); resNS_lu[262] = 0x11;
   memmove( resNS_tb[263], "HRP\0", 4 ); resNS_lu[263] = 0x16;
   memmove( resNS_tb[264], "HSE\0", 4 ); resNS_lu[264] = 0x12;
   memmove( resNS_tb[265], "HSL\0", 4 ); resNS_lu[265] = 0x12;
   memmove( resNS_tb[266], "HSO\0", 4 ); resNS_lu[266] = 0x07;
   memmove( resNS_tb[267], "HTI\0", 4 ); resNS_lu[267] = 0x02;
   memmove( resNS_tb[268], "HTR\0", 4 ); resNS_lu[268] = 0x16;
   memmove( resNS_tb[269], "HV5\0", 4 ); resNS_lu[269] = 0x00;
   memmove( resNS_tb[270], "HY3\0", 4 ); resNS_lu[270] = 0x0F;
   memmove( resNS_tb[271], "HYI\0", 4 ); resNS_lu[271] = 0x0C;
   memmove( resNS_tb[272], "HYP\0", 4 ); resNS_lu[272] = 0x0F;
   resNS_dx[0x08] = 273;
   memmove( resNS_tb[273], "I58\0", 4 ); resNS_lu[273] = 0x0A;
   memmove( resNS_tb[274], "IAM\0", 4 ); resNS_lu[274] = 0x00;
   memmove( resNS_tb[275], "IAS\0", 4 ); resNS_lu[275] = 0x03;
   memmove( resNS_tb[276], "IGL\0", 4 ); resNS_lu[276] = 0x06;
   memmove( resNS_tb[277], "IIL\0", 4 ); resNS_lu[277] = 0x08;
   memmove( resNS_tb[278], "ILG\0", 4 ); resNS_lu[278] = 0x04;
   memmove( resNS_tb[279], "ILX\0", 4 ); resNS_lu[279] = 0x08;
   memmove( resNS_tb[280], "IML\0", 4 ); resNS_lu[280] = 0x08;
   memmove( resNS_tb[281], "IOY\0", 4 ); resNS_lu[281] = 0x05;
   memmove( resNS_tb[282], "IPG\0", 4 ); resNS_lu[282] = 0x06;
   memmove( resNS_tb[283], "IYR\0", 4 ); resNS_lu[283] = 0x18;
   memmove( resNS_tb[284], "IYT\0", 4 ); resNS_lu[284] = 0x13;
   resNS_dx[0x09]  = 285;
   resNS_dx[0x0A] = 285;
   memmove( resNS_tb[285], "K1R\0", 4 ); resNS_lu[285] = 0x02;
   memmove( resNS_tb[286], "KCX\0", 4 ); resNS_lu[286] = 0x0A;
   memmove( resNS_tb[287], "KGC\0", 4 ); resNS_lu[287] = 0x0A;
   memmove( resNS_tb[288], "KOR\0", 4 ); resNS_lu[288] = 0x0C;
   memmove( resNS_tb[289], "KST\0", 4 ); resNS_lu[289] = 0x0A;
   memmove( resNS_tb[290], "KYN\0", 4 ); resNS_lu[290] = 0x00;
   resNS_dx[0x0B] = 292;
   memmove( resNS_tb[291], "LAL\0", 4 ); resNS_lu[291] = 0x00;
   memmove( resNS_tb[292], "LCX\0", 4 ); resNS_lu[292] = 0x0A;
   memmove( resNS_tb[293], "LDH\0", 4 ); resNS_lu[293] = 0x0A;
   memmove( resNS_tb[294], "LED\0", 4 ); resNS_lu[294] = 0x0B;
   memmove( resNS_tb[295], "LEF\0", 4 ); resNS_lu[295] = 0x0B;
   memmove( resNS_tb[296], "LEN\0", 4 ); resNS_lu[296] = 0x0B;
   memmove( resNS_tb[297], "LLP\0", 4 ); resNS_lu[297] = 0x0A;
   memmove( resNS_tb[298], "LLY\0", 4 ); resNS_lu[298] = 0x0A;
   memmove( resNS_tb[299], "LME\0", 4 ); resNS_lu[299] = 0x04;
   memmove( resNS_tb[300], "LPD\0", 4 ); resNS_lu[300] = 0x0F;
   memmove( resNS_tb[301], "LPG\0", 4 ); resNS_lu[301] = 0x06;
   memmove( resNS_tb[302], "LPS\0", 4 ); resNS_lu[302] = 0x12;
   memmove( resNS_tb[303], "LTR\0", 4 ); resNS_lu[303] = 0x16;
   memmove( resNS_tb[304], "LVG\0", 4 ); resNS_lu[304] = 0x06;
   memmove( resNS_tb[305], "LYM\0", 4 ); resNS_lu[305] = 0x0A;
   memmove( resNS_tb[306], "LYN\0", 4 ); resNS_lu[306] = 0x0A;
   memmove( resNS_tb[307], "LYR\0", 4 ); resNS_lu[307] = 0x0A;
   memmove( resNS_tb[308], "LYX\0", 4 ); resNS_lu[308] = 0x0A;
   memmove( resNS_tb[309], "LYZ\0", 4 ); resNS_lu[309] = 0x0A;
   resNS_dx[0x0C] = 310;
   memmove( resNS_tb[310], "M0H\0", 4 ); resNS_lu[310] = 0x02;
   memmove( resNS_tb[311], "M3L\0", 4 ); resNS_lu[311] = 0x0A;
   memmove( resNS_tb[312], "MAA\0", 4 ); resNS_lu[312] = 0x00;
   memmove( resNS_tb[313], "MAI\0", 4 ); resNS_lu[313] = 0x11;
   memmove( resNS_tb[314], "MBQ\0", 4 ); resNS_lu[314] = 0x18;
   memmove( resNS_tb[315], "MC1\0", 4 ); resNS_lu[315] = 0x12;
   memmove( resNS_tb[316], "MCL\0", 4 ); resNS_lu[316] = 0x0A;
   memmove( resNS_tb[317], "MCS\0", 4 ); resNS_lu[317] = 0x02;
   memmove( resNS_tb[318], "MEA\0", 4 ); resNS_lu[318] = 0x05;
   memmove( resNS_tb[319], "MEG\0", 4 ); resNS_lu[319] = 0x04;
   memmove( resNS_tb[320], "MEN\0", 4 ); resNS_lu[320] = 0x0D;
   memmove( resNS_tb[321], "MEQ\0", 4 ); resNS_lu[321] = 0x10;
   memmove( resNS_tb[322], "MEU\0", 4 ); resNS_lu[322] = 0x06;
   memmove( resNS_tb[323], "MFN\0", 4 ); resNS_lu[323] = 0x04;
   memmove( resNS_tb[324], "MGG\0", 4 ); resNS_lu[324] = 0x11;
   memmove( resNS_tb[325], "MGN\0", 4 ); resNS_lu[325] = 0x10;
   memmove( resNS_tb[326], "MGY\0", 4 ); resNS_lu[326] = 0x06;
   memmove( resNS_tb[327], "MHL\0", 4 ); resNS_lu[327] = 0x0B;
   memmove( resNS_tb[328], "MHO\0", 4 ); resNS_lu[328] = 0x0C;
   memmove( resNS_tb[329], "MHS\0", 4 ); resNS_lu[329] = 0x07;
   memmove( resNS_tb[330], "MIS\0", 4 ); resNS_lu[330] = 0x12;
   memmove( resNS_tb[331], "MLE\0", 4 ); resNS_lu[331] = 0x0B;
   memmove( resNS_tb[332], "MLL\0", 4 ); resNS_lu[332] = 0x0B;
   memmove( resNS_tb[333], "MLY\0", 4 ); resNS_lu[333] = 0x0A;
   memmove( resNS_tb[334], "MLZ\0", 4 ); resNS_lu[334] = 0x0A;
   memmove( resNS_tb[335], "MME\0", 4 ); resNS_lu[335] = 0x0C;
   memmove( resNS_tb[336], "MNL\0", 4 ); resNS_lu[336] = 0x0B;
   memmove( resNS_tb[337], "MNV\0", 4 ); resNS_lu[337] = 0x15;
   memmove( resNS_tb[338], "MPQ\0", 4 ); resNS_lu[338] = 0x06;
   memmove( resNS_tb[339], "MSA\0", 4 ); resNS_lu[339] = 0x06;
   memmove( resNS_tb[340], "MSE\0", 4 ); resNS_lu[340] = 0x0C;
   memmove( resNS_tb[341], "MSL\0", 4 ); resNS_lu[341] = 0x0C;
   memmove( resNS_tb[342], "MSO\0", 4 ); resNS_lu[342] = 0x0C;
   memmove( resNS_tb[343], "MTY\0", 4 ); resNS_lu[343] = 0x18;
   memmove( resNS_tb[344], "MVA\0", 4 ); resNS_lu[344] = 0x15;
   resNS_dx[0x0D] = 345;
   memmove( resNS_tb[345], "N10\0", 4 ); resNS_lu[345] = 0x12;
   memmove( resNS_tb[346], "N7P\0", 4 ); resNS_lu[346] = 0x0F;
   memmove( resNS_tb[347], "NAL\0", 4 ); resNS_lu[347] = 0x00;
   memmove( resNS_tb[348], "NAM\0", 4 ); resNS_lu[348] = 0x00;
   memmove( resNS_tb[349], "NBQ\0", 4 ); resNS_lu[349] = 0x18;
   memmove( resNS_tb[350], "NC1\0", 4 ); resNS_lu[350] = 0x12;
   memmove( resNS_tb[351], "NCB\0", 4 ); resNS_lu[351] = 0x00;
   memmove( resNS_tb[352], "NDF\0", 4 ); resNS_lu[352] = 0x05;
   memmove( resNS_tb[353], "NEM\0", 4 ); resNS_lu[353] = 0x07;
   memmove( resNS_tb[354], "NEP\0", 4 ); resNS_lu[354] = 0x07;
   memmove( resNS_tb[355], "NFA\0", 4 ); resNS_lu[355] = 0x05;
   memmove( resNS_tb[356], "NHL\0", 4 ); resNS_lu[356] = 0x04;
   memmove( resNS_tb[357], "NIY\0", 4 ); resNS_lu[357] = 0x18;
   memmove( resNS_tb[358], "NLE\0", 4 ); resNS_lu[358] = 0x0B;
   memmove( resNS_tb[359], "NLN\0", 4 ); resNS_lu[359] = 0x0B;
   memmove( resNS_tb[360], "NLO\0", 4 ); resNS_lu[360] = 0x0B;
   memmove( resNS_tb[361], "NLP\0", 4 ); resNS_lu[361] = 0x0B;
   memmove( resNS_tb[362], "NLQ\0", 4 ); resNS_lu[362] = 0x10;
   memmove( resNS_tb[363], "NMC\0", 4 ); resNS_lu[363] = 0x06;
   memmove( resNS_tb[364], "NNH\0", 4 ); resNS_lu[364] = 0x11;
   memmove( resNS_tb[365], "NPH\0", 4 ); resNS_lu[365] = 0x02;
   memmove( resNS_tb[366], "NTY\0", 4 ); resNS_lu[366] = 0x18;
   memmove( resNS_tb[367], "NVA\0", 4 ); resNS_lu[367] = 0x15;
   memmove( resNS_tb[368], "NZH\0", 4 ); resNS_lu[368] = 0x07;
   resNS_dx[0x0E] = 369;
   memmove( resNS_tb[369], "OAS\0", 4 ); resNS_lu[369] = 0x12;
   memmove( resNS_tb[370], "OCS\0", 4 ); resNS_lu[370] = 0x02;
   memmove( resNS_tb[371], "OCY\0", 4 ); resNS_lu[371] = 0x02;
   memmove( resNS_tb[372], "OHS\0", 4 ); resNS_lu[372] = 0x03;
   memmove( resNS_tb[373], "OLT\0", 4 ); resNS_lu[373] = 0x13;
   memmove( resNS_tb[374], "OMT\0", 4 ); resNS_lu[374] = 0x0C;
   memmove( resNS_tb[375], "OPR\0", 4 ); resNS_lu[375] = 0x11;
   memmove( resNS_tb[376], "ORN\0", 4 ); resNS_lu[376] = 0x00;
   memmove( resNS_tb[377], "ORQ\0", 4 ); resNS_lu[377] = 0x11;
   memmove( resNS_tb[378], "OSE\0", 4 ); resNS_lu[378] = 0x12;
   memmove( resNS_tb[379], "OTY\0", 4 ); resNS_lu[379] = 0x18;
   memmove( resNS_tb[380], "OXX\0", 4 ); resNS_lu[380] = 0x03;
   resNS_dx[0x0F] = 381;
   memmove( resNS_tb[381], "P1L\0", 4 ); resNS_lu[381] = 0x02;
   memmove( resNS_tb[382], "P2Y\0", 4 ); resNS_lu[382] = 0x0F;
   memmove( resNS_tb[383], "PAQ\0", 4 ); resNS_lu[383] = 0x18;
   memmove( resNS_tb[384], "PAS\0", 4 ); resNS_lu[384] = 0x03;
   memmove( resNS_tb[385], "PAT\0", 4 ); resNS_lu[385] = 0x16;
   memmove( resNS_tb[386], "PAU\0", 4 ); resNS_lu[386] = 0x00;
   memmove( resNS_tb[387], "PBB\0", 4 ); resNS_lu[387] = 0x02;
   memmove( resNS_tb[388], "PBF\0", 4 ); resNS_lu[388] = 0x05;
   memmove( resNS_tb[389], "PCA\0", 4 ); resNS_lu[389] = 0x04;
   memmove( resNS_tb[390], "PCC\0", 4 ); resNS_lu[390] = 0x0F;
   memmove( resNS_tb[391], "PCS\0", 4 ); resNS_lu[391] = 0x05;
   memmove( resNS_tb[392], "PEC\0", 4 ); resNS_lu[392] = 0x02;
   memmove( resNS_tb[393], "PF5\0", 4 ); resNS_lu[393] = 0x05;
   memmove( resNS_tb[394], "PFF\0", 4 ); resNS_lu[394] = 0x05;
   memmove( resNS_tb[395], "PG1\0", 4 ); resNS_lu[395] = 0x12;
   memmove( resNS_tb[396], "PG9\0", 4 ); resNS_lu[396] = 0x06;
   memmove( resNS_tb[397], "PGY\0", 4 ); resNS_lu[397] = 0x06;
   memmove( resNS_tb[398], "PHA\0", 4 ); resNS_lu[398] = 0x05;
   memmove( resNS_tb[399], "PHD\0", 4 ); resNS_lu[399] = 0x03;
   memmove( resNS_tb[400], "PHI\0", 4 ); resNS_lu[400] = 0x05;
   memmove( resNS_tb[401], "PHL\0", 4 ); resNS_lu[401] = 0x05;
   memmove( resNS_tb[402], "PHM\0", 4 ); resNS_lu[402] = 0x05;
   memmove( resNS_tb[403], "PLE\0", 4 ); resNS_lu[403] = 0x0B;
   memmove( resNS_tb[404], "PM3\0", 4 ); resNS_lu[404] = 0x05;
   memmove( resNS_tb[405], "POM\0", 4 ); resNS_lu[405] = 0x0F;
   memmove( resNS_tb[406], "PPH\0", 4 ); resNS_lu[406] = 0x0B;
   memmove( resNS_tb[407], "PPN\0", 4 ); resNS_lu[407] = 0x05;
   memmove( resNS_tb[408], "PR3\0", 4 ); resNS_lu[408] = 0x02;
   memmove( resNS_tb[409], "PRR\0", 4 ); resNS_lu[409] = 0x00;
   memmove( resNS_tb[410], "PRS\0", 4 ); resNS_lu[410] = 0x0F;
   memmove( resNS_tb[411], "PSA\0", 4 ); resNS_lu[411] = 0x05;
   memmove( resNS_tb[412], "PSH\0", 4 ); resNS_lu[412] = 0x07;
   memmove( resNS_tb[413], "PTH\0", 4 ); resNS_lu[413] = 0x18;
   memmove( resNS_tb[414], "PTM\0", 4 ); resNS_lu[414] = 0x18;
   memmove( resNS_tb[415], "PTR\0", 4 ); resNS_lu[415] = 0x18;
   memmove( resNS_tb[416], "PVH\0", 4 ); resNS_lu[416] = 0x07;
   memmove( resNS_tb[417], "PYA\0", 4 ); resNS_lu[417] = 0x00;
   memmove( resNS_tb[418], "PYX\0", 4 ); resNS_lu[418] = 0x02;
   resNS_dx[0x10] = 419;
   resNS_dx[0x11] = 419;
   memmove( resNS_tb[419], "R1A\0", 4 ); resNS_lu[419] = 0x02;
   memmove( resNS_tb[420], "R1B\0", 4 ); resNS_lu[420] = 0x02;
   memmove( resNS_tb[421], "R1F\0", 4 ); resNS_lu[421] = 0x02;
   memmove( resNS_tb[422], "R7A\0", 4 ); resNS_lu[422] = 0x02;
   memmove( resNS_tb[423], "RCY\0", 4 ); resNS_lu[423] = 0x02;
   resNS_dx[0x12] = 424;
   memmove( resNS_tb[424], "S1H\0", 4 ); resNS_lu[424] = 0x12;
   memmove( resNS_tb[425], "SAC\0", 4 ); resNS_lu[425] = 0x12;
   memmove( resNS_tb[426], "SAH\0", 4 ); resNS_lu[426] = 0x02;
   memmove( resNS_tb[427], "SAR\0", 4 ); resNS_lu[427] = 0x06;
   memmove( resNS_tb[428], "SBD\0", 4 ); resNS_lu[428] = 0x12;
   memmove( resNS_tb[429], "SBL\0", 4 ); resNS_lu[429] = 0x12;
   memmove( resNS_tb[430], "SCH\0", 4 ); resNS_lu[430] = 0x02;
   memmove( resNS_tb[431], "SCS\0", 4 ); resNS_lu[431] = 0x02;
   memmove( resNS_tb[432], "SCY\0", 4 ); resNS_lu[432] = 0x02;
   memmove( resNS_tb[433], "SDP\0", 4 ); resNS_lu[433] = 0x12;
   memmove( resNS_tb[434], "SEB\0", 4 ); resNS_lu[434] = 0x12;
   memmove( resNS_tb[435], "SEC\0", 4 ); resNS_lu[435] = 0x00;
   memmove( resNS_tb[436], "SEG\0", 4 ); resNS_lu[436] = 0x00;
   memmove( resNS_tb[437], "SEL\0", 4 ); resNS_lu[437] = 0x12;
   memmove( resNS_tb[438], "SEP\0", 4 ); resNS_lu[438] = 0x12;
   memmove( resNS_tb[439], "SET\0", 4 ); resNS_lu[439] = 0x12;
   memmove( resNS_tb[440], "SHC\0", 4 ); resNS_lu[440] = 0x02;
   memmove( resNS_tb[441], "SHP\0", 4 ); resNS_lu[441] = 0x06;
   memmove( resNS_tb[442], "SHR\0", 4 ); resNS_lu[442] = 0x0A;
   memmove( resNS_tb[443], "SIB\0", 4 ); resNS_lu[443] = 0x02;
   memmove( resNS_tb[444], "SLZ\0", 4 ); resNS_lu[444] = 0x0A;
   memmove( resNS_tb[445], "SMC\0", 4 ); resNS_lu[445] = 0x02;
   memmove( resNS_tb[446], "SME\0", 4 ); resNS_lu[446] = 0x0C;
   memmove( resNS_tb[447], "SMF\0", 4 ); resNS_lu[447] = 0x05;
   memmove( resNS_tb[448], "SNC\0", 4 ); resNS_lu[448] = 0x02;
   memmove( resNS_tb[449], "SOC\0", 4 ); resNS_lu[449] = 0x02;
   memmove( resNS_tb[450], "SOY\0", 4 ); resNS_lu[450] = 0x12;
   memmove( resNS_tb[451], "STY\0", 4 ); resNS_lu[451] = 0x18;
   memmove( resNS_tb[452], "SVA\0", 4 ); resNS_lu[452] = 0x12;
   resNS_dx[0x13] = 453;
   memmove( resNS_tb[453], "TAV\0", 4 ); resNS_lu[453] = 0x03;
   memmove( resNS_tb[454], "TBG\0", 4 ); resNS_lu[454] = 0x06;
   memmove( resNS_tb[455], "TBM\0", 4 ); resNS_lu[455] = 0x13;
   memmove( resNS_tb[456], "THC\0", 4 ); resNS_lu[456] = 0x13;
   memmove( resNS_tb[457], "TIH\0", 4 ); resNS_lu[457] = 0x00;
   memmove( resNS_tb[458], "TMB\0", 4 ); resNS_lu[458] = 0x13;
   memmove( resNS_tb[459], "TMD\0", 4 ); resNS_lu[459] = 0x13;
   memmove( resNS_tb[460], "TNB\0", 4 ); resNS_lu[460] = 0x02;
   memmove( resNS_tb[461], "TNR\0", 4 ); resNS_lu[461] = 0x12;
   memmove( resNS_tb[462], "TOX\0", 4 ); resNS_lu[462] = 0x16;
   memmove( resNS_tb[463], "TPL\0", 4 ); resNS_lu[463] = 0x16;
   memmove( resNS_tb[464], "TPO\0", 4 ); resNS_lu[464] = 0x13;
   memmove( resNS_tb[465], "TQQ\0", 4 ); resNS_lu[465] = 0x16;
   memmove( resNS_tb[466], "TRF\0", 4 ); resNS_lu[466] = 0x16;
   memmove( resNS_tb[467], "TRN\0", 4 ); resNS_lu[467] = 0x16;
   memmove( resNS_tb[468], "TRO\0", 4 ); resNS_lu[468] = 0x16;
   memmove( resNS_tb[469], "TRQ\0", 4 ); resNS_lu[469] = 0x16;
   memmove( resNS_tb[470], "TRW\0", 4 ); resNS_lu[470] = 0x16;
   memmove( resNS_tb[471], "TRX\0", 4 ); resNS_lu[471] = 0x16;
   memmove( resNS_tb[472], "TTQ\0", 4 ); resNS_lu[472] = 0x16;
   memmove( resNS_tb[473], "TTS\0", 4 ); resNS_lu[473] = 0x18;
   memmove( resNS_tb[474], "TY2\0", 4 ); resNS_lu[474] = 0x18;
   memmove( resNS_tb[475], "TY3\0", 4 ); resNS_lu[475] = 0x18;
   memmove( resNS_tb[476], "TYB\0", 4 ); resNS_lu[476] = 0x18;
   memmove( resNS_tb[477], "TYI\0", 4 ); resNS_lu[477] = 0x18;
   memmove( resNS_tb[478], "TYN\0", 4 ); resNS_lu[478] = 0x18;
   memmove( resNS_tb[479], "TYO\0", 4 ); resNS_lu[479] = 0x18;
   memmove( resNS_tb[480], "TYQ\0", 4 ); resNS_lu[480] = 0x18;
   memmove( resNS_tb[481], "TYS\0", 4 ); resNS_lu[481] = 0x18;
   memmove( resNS_tb[482], "TYT\0", 4 ); resNS_lu[482] = 0x18;
   memmove( resNS_tb[483], "TYY\0", 4 ); resNS_lu[483] = 0x18;
   resNS_dx[0x14] = 484;
   memmove( resNS_tb[484], "UMA\0", 4 ); resNS_lu[484] = 0x00;
   resNS_dx[0x15] = 485;
   memmove( resNS_tb[485], "VAD\0", 4 ); resNS_lu[485] = 0x15;
   memmove( resNS_tb[486], "VAF\0", 4 ); resNS_lu[486] = 0x15;
   resNS_dx[0x16] = 487;
   resNS_dx[0x17] = 487;
   memmove( resNS_tb[487], "XX1\0", 4 ); resNS_lu[487] = 0x0A;
   resNS_dx[0x18] = 488;
   memmove( resNS_tb[488], "YCM\0", 4 ); resNS_lu[488] = 0x02;
   memmove( resNS_tb[489], "YOF\0", 4 ); resNS_lu[489] = 0x18;
   resNS_dx[0x19] = 490;
   resNS_dx[0x1A] = 490;
   memmove( resNS_tb[490], "SAM\0", 4 ); resNS_lu[490] = 0x0C;
   memmove( resNS_tb[491], "UNK\0", 4 ); resNS_lu[491] = 0x17;
   resNS_dx[0x1B] = 492;

   return;
}

void bAA::print( FILE* op ) {
   fprintf(op, "[bAA]\n\t");
   for( int i=0; i < 26; ++i ) {
      fprintf(op, "[%2d] %c:%s   ", i, res1_[i], res3_[i]);
      if( (i % 5) == 0 ) { fprintf(op, "\n\t"); }
   }
   fprintf(op, "\n");
   return;
}
