
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "bVBase.h"
#include "bSnapp.h"
#include "bAA.h"
#include "bSys.h"

using namespace std;
using namespace bStd;

/* STATIC */
bool bSnapp::isSet_ = false;
char bSnapp::file_[256] = "\0";
char bSnapp::defdir_[64] = "~/Dropbox/";
char bSnapp::deffile_[64] = "lib/SNAPP_VALUES_BALA.txt";
bool bSnapp::haveFile_ = findFile( bSnapp::file_ );

float bSnapp::snapp_[26][26][26][26][5];

bool bSnapp::findFile( char f[] ) {
   char* idir = getenv("GRID_DIR");
   if( idir == NULL ) {
      idir = defdir_;
   }
   int dsize = strlen( idir );
   int fsize = strlen( deffile_ );
   memset( f, '\0', 256 );
   memmove( f, idir, dsize );
   memmove( f+dsize, deffile_, fsize );
   printf("IDIR: %s\n", getenv("GRID_DIR"));
   printf("IDIR: %s\n", idir);
   printf("SNAPP file: %s\n", f);
   return bSys::fileExists(f);
}

/****** CONSTRUCTOR & DESTRUCTOR */
bSnapp::bSnapp()
{
}

bSnapp::bSnapp( const bSnapp &rhs )
{
}

bSnapp::~bSnapp()
{
}

/****** READ SNAPP */

void bSnapp::readSnapp( const char file[] ) {
   if( !haveFile_ ) {
      printf("no snapp file: %s\n", file);
      FILE* op = fopen( file, "r");
      if( !op ) { printf("nope...no file\n"); }
   }
   if( isSet_ || !haveFile_ ) return;

   // open the file
   ifstream ip;
   ip.open(file, ifstream::in);
   if(!ip) { printf(" No snapp file: %s\n",file_); return; }

   // read in the file... not the most efficient, but we need the size
   // -- we read in w/o doing anything to get a cound of the residues
   //~ deque<string> data;
   string bffr;
   char res;
   for( uint i = 0; getline( ip, bffr ); ++i ) {
      istringstream ss(bffr);
      int v[4];
      ss >> res; v[0] = res2num( res );
      ss >> res; v[1] = res2num( res );
      ss >> res; v[2] = res2num( res );
      ss >> res; v[3] = res2num( res );
      //~ = { res2num( ss >> res ), res2num( ss >> res ),
                   //~ res2num( ss >> res ), res2num( ss >> res ) };
      ss >> snapp_[ v[0] ][ v[1] ][ v[2] ][ v[3] ][0];
      ss >> snapp_[ v[0] ][ v[1] ][ v[2] ][ v[3] ][1];
      ss >> snapp_[ v[0] ][ v[1] ][ v[2] ][ v[3] ][2];
      ss >> snapp_[ v[0] ][ v[1] ][ v[2] ][ v[3] ][3];
      ss >> snapp_[ v[0] ][ v[1] ][ v[2] ][ v[3] ][4];

      //~ printf("[%-4d,", v[0]);
      //~ printf("%-4d,", v[1]);
      //~ printf("%-4d,", v[2]);
      //~ printf("%-4d,]", v[3]);
      //~ printf("{%-6.2f,", snapp_[ v[0] ][ v[1] ][ v[2] ][ v[3] ][0]);
      //~ printf("%-6.2f,", snapp_[ v[0] ][ v[1] ][ v[2] ][ v[3] ][1]);
      //~ printf("%-6.2f,", snapp_[ v[0] ][ v[1] ][ v[2] ][ v[3] ][2]);
      //~ printf("%-6.2f }\n", snapp_[ v[0] ][ v[1] ][ v[2] ][ v[3] ][3]);
   }
   ip.close();
//~ exit(1);
   isSet_ = true;
   return;
}



float bSnapp::score( char res[4], int type ) {
   int v[4] = { res2num( res[0] ), res2num( res[1] ),
                res2num( res[2] ), res2num( res[3] ) };
   return score( v, type );
}

float bSnapp::score( int res[4], int type ) {
   if( !isSet_ ) readSnapp();
   return snapp_[ res[0] ][ res[1] ][ res[2] ][ res[3] ][ type ];
}

int bSnapp::res2num( char res ) {
   res = toupper( res );
   int num = res;
   num -= 'A';
   //~ switch( res ) {
      //~ case 'A': num =  0; break; case 'B': num =  1; break;
      //~ case 'C': num =  2; break; case 'D': num =  3; break;
      //~ case 'E': num =  4; break; case 'F': num =  5; break;
      //~ case 'G': num =  6; break; case 'H': num =  7; break;
      //~ case 'I': num =  8; break; case 'J': num =  9; break;
      //~ case 'K': num = 10; break; case 'L': num = 11; break;
      //~ case 'M': num = 12; break; case 'N': num = 13; break;
      //~ case 'O': num = 14; break; case 'P': num = 15; break;
      //~ case 'Q': num = 16; break; case 'R': num = 17; break;
      //~ case 'S': num = 18; break; case 'T': num = 19; break;
      //~ case 'U': num = 20; break; case 'V': num = 21; break;
      //~ case 'W': num = 22; break; case 'X': num = 23; break;
      //~ case 'Y': num = 24; break; case 'Z': num = 25; break;
      //~ default:  num = -1; break;
   //~ };
   return num;
}

void print( FILE* op ) {
   return;
}