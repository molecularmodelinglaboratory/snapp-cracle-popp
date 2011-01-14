#include <deque>
#include <string>

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "bVBase.h"
#include "bSys.h"

using namespace std;
using namespace bStd;

/****** Constructors */
bSys::bSys()
{
}

bSys::~bSys()
{
}


/****** Command Line Arguements */
uint bSys::cla( int numArg, char** arrArg, char** &flgSaved, char** &valSaved ) {
   // Description: Takes the command line arguements & sorts them into
   //    two corresponding lists -- one w/ flags, the other w/ values.
   //    WRONG -- It will allow for values without flags (only at the beginning)
   //    WRONG -- or for flags w/o or w/ many values.
   //    If there is no value, the value is set as an empty string (not NULL).

   if( numArg == 1 ) { return 0; }

   deque<string> flg;
   deque<string> val;
   string tmpFlg;
   string tmpVal;
   for( int i=1; i < numArg; ++i ) {
      printf("[%d] %s :: ", i, arrArg[i] );
      if( arrArg[i][0] == '-' ) {
         printf(" flg :: ");
         ++arrArg[i];
         tmpFlg = arrArg[i];
         //~ bool hasVar = true;
         //~ while( hasVar && (i + 1) < numArg ) {
         if( arrArg[(i + 1)][0] != '-' ) {
            //~ printf("(%d) ",i);
            printf(" has var :: ");
            tmpVal = arrArg[ ++i ];
            flg.push_back( tmpFlg );
            val.push_back( tmpVal );
         }
         //~ }
         else {
            tmpVal = "";
            printf(" no  var :: ");
            flg.push_back( tmpFlg );
            val.push_back( tmpVal );
            break;
         }
      }
      else {
         printf(" var :: ");
         tmpFlg = "";
         tmpVal = arrArg[i];
         flg.push_back( tmpFlg );
         val.push_back( tmpVal );
      }
      printf("\n");
   }

   flgSaved = (flg.size() == 0) ? NULL : new char*[ flg.size() ];
   valSaved = (val.size() == 0) ? NULL : new char*[ val.size() ];

   for( uint i=0; i < flg.size(); ++i ) {
      flgSaved[i] = new char[ flg[i].size() + 1 ];
      memset( flgSaved[i], '\0', flg[i].size() + 1 );
      memmove( flgSaved[i], flg[i].c_str(), flg[i].size() );

      valSaved[i] = new char[ val[i].size() + 1 ];
      memset( valSaved[i], '\0', val[i].size() + 1 );
      memmove( valSaved[i], val[i].c_str(), val[i].size() );

      printf("found: [%lu] %s => %s :: %s => %s\n", flg[i].size(), flg[i].c_str(), val[i].c_str(), flgSaved[i], valSaved[i]);
   }

   return flg.size();
}

/****** File Checking */

bool bSys::fileExists( char* file, char* path) { return fileExists( (const char*) file, (const char*) path ); }
bool bSys::fileExists( const char* file, const char* path) {
   int ps = strlen( path );
   int fs = strlen( file );
   char fullname[ ps + fs + 1 ];
   memset( fullname, '\0', ps + fs + 1 );
   memmove( fullname, path, ps );
   memmove( fullname + ps, file, fs );
   return fileExists( fullname );
}

bool bSys::fileExists( char* file ) { return fileExists( (const char*) file ); }
bool bSys::fileExists( const char* file ) {
   // modified from <http://www.techbytes.ca/techbyte103.html>
   struct stat info;
   bool exists = false;
   int status;

   /* Get attrib */ status = stat( file, &info );
   /* Does exist */ if( status == 0 && !(info.st_mode & S_IFDIR) ) { exists = true; }
   /* May not    */ else {}
  return exists;
}

bool bSys::dirExists( const char* dir, const char* path) {
   int ps = strlen( path );
   int fs = strlen( dir );
   char fullname[ ps + fs + 1 ];
   memset( fullname, '\0', ps + fs + 1 );
   memmove( fullname, path, ps );
   memmove( fullname + ps, dir, fs );
   return dirExists( fullname );
}
bool bSys::dirExists( const char* dir ) {
   struct stat info;
   bool exists = false;
   int status;

   /* Get attrib */ status = stat( dir, &info );
   /* Does exist */ if( status == 0 && (info.st_mode & S_IFDIR) ) { exists = true; }
   /* May not    */ else {}
  return exists;
}


uint bSys::fileList( const char* dir2check, char** &savelist ) {
   if( savelist != NULL ) { throw "[bSys::fileList] Need empty list."; }

   // Open directory
   DIR* dp;
   struct dirent* dir;
   if( (dp = opendir( dir2check )) == NULL ) {
      char e[96];
      sprintf( e, "[bSys] Unable to open directory <%s>", dir2check );
      throw e;
   }

   // Get file bases
   deque<string> files;
   char prev[64]; memset( prev, '\0', 64 );
   while( (dir = readdir(dp)) != NULL ) {
      //~ printf("[bSys] file: %s", dir->d_name);
      if( dir->d_name[0] == '.' ) { continue; }
      else if( dirExists( dir->d_name, dir2check ) ) { continue; }
      else if( fileExists( dir->d_name, dir2check ) ) {
         files.push_back( string(dir->d_name) );
      }
      else {}
   }
   closedir(dp);

   // Save the file list
   savelist = new char*[ files.size() ];
   for( uint i=0; i < files.size(); ++i ) {
      savelist[i] = new char[ files[i].size() + 1 ];
      memset( savelist[i], '\0', files[i].size() + 1 );
      memmove( savelist[i], files[i].c_str(), files[i].size() );
   }

   return files.size();
}
/****** File Handling */
bool bSys::removeExt( char* file ) {
   int size = strlen( file );
   int rpos = strcspn( file, "." );
   int cnt = size - rpos;
   
   //~ --rpos;
   memset( file + rpos, '\0', cnt );
   return true;

   while( file[ size - 1 ] != '.' ) { --size; ++cnt; }
   --size; ++cnt;
   memset( file + size, '\0', cnt );
   return true;
}

bool bSys::multicat( char* &save, const char* c1, const char* c2, const char* c3, const char* c4 ) {
   char* temp = NULL;
   uint s0 = 0;
   uint s1 = (c1 == NULL) ? 0 : strlen( c1 );
   uint s2 = (c2 == NULL) ? 0 : strlen( c2 );
   uint s3 = (c3 == NULL) ? 0 : strlen( c3 );
   uint s4 = (c4 == NULL) ? 0 : strlen( c4 );
   uint ss = s1; ss += s2; ss += s3; ss += s4; ++ss;
   if( save != NULL ) {
      s0 = strlen( save );
      ss += s0;
      temp = new char[ s0 + 1 ];
      memmove( temp, save, s0 ); temp[ s0 ] = '\0';
      delete [] save; save = NULL;
   }
   save = new char[ ss ];
   memset( save, '\0', ss );
   if( s0 ) { memmove( save, temp, s0 );}
   if( s1 ) { memmove( save + s0, c1, s1 ); s1 += s0; }
   if( s2 ) { memmove( save + s1, c2, s2 ); s2 += s1; }
   if( s3 ) { memmove( save + s2, c3, s3 ); s3 += s2; }
   if( s4 ) { memmove( save + s3, c4, s4 ); }
   
   if( temp != NULL ) {
      delete [] temp;
      temp = NULL;
   }
   return (ss > 0) ? true : false;
}