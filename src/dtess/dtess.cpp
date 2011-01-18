// Filename: dtess.cpp
//Author:   Stephen J. Bush <sjbush at unc dot edu>
// Date:     11.16.10
//
// Delaunay Tessellation

#include <omp.h>
double omp_get_wtime(void);

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "bDelTess.h"
#include "bPoints.h"
#include "bCentroid.h"

using namespace std;
using namespace bStd;

void usage() {
   printf("Usage: ./dtess -p <file> <options>\n");
   printf("\t-p\tpdbfile (will calculate side chain centroids) \t[can handle multiple chains]\n");
   printf("\t  \tOR point coordinates \t\t\t\t[will act as one chain]\n");
   printf("\t-t\ttrim length (float)\n");
   printf("\t-w\twhite background for pymol file\n");
   printf("\t-h\tthis menu\n");
   exit(0);
}

int main( int argc, char **argv ) {

   char* file = new char[128]; memset( file, '\0', 128 );
   float trim = 10.0;
   bool white = false;
   bool isPDB = true;

   // Get input options
   if( argc < 2 ) { usage(); }
   for( int i=1; i < argc; ++i ) {
         if(argv[i][0] == '-') {
         switch(argv[i][1]) {
            case 'p':
               strcpy( file, argv[++i] );
               break;
            case 't':
               trim = atof(argv[++i]);
               break;
            case 'w':
               white = true;
               break;
            case 'h':
               usage();
            default:
               (void)printf("\nUnknown parameter: '%s'\n\n",argv[i]);
               usage();
               break;
         }
      }
   }


   // Check for valid file
   FILE* op = fopen( file, "r" );
   if( !op ) { printf("Unable to open '%s'\n",file); exit(1); }
   fclose(op);

   double start, end;
   start = omp_get_wtime();

   // Get file parts
   char* endpath = strrchr( file, '/' ); ++endpath;
   char* endname = strrchr( file, '.' );
   int basesize = endname - endpath;
   char* base = new char[ ++basesize ];
   memset( base, '\0', basesize );
   memmove( base, endpath, --basesize );

   int pathsize = endpath - file;
   char* path = new char[ ++pathsize ];
   memset( path, '\0', pathsize );
   memmove( path, file, --pathsize );

   printf("file: %s\t", file);
   printf("path: %s\t", path);
   printf("base: %s\n", base);
   //~ printf("\n");

   // Check for PDB file
   if( memcmp( endname, ".pdb", 4 ) != 0 ) { isPDB = false; }

   // Read in the protein and calculate the centroids
   bPoints** protchains = NULL;
   int numchains = 0;
   if( isPDB ) {
      printf("here\n");
      bCentroid::pdb2Centroids( protchains, numchains, base, path );
      //~ printf("numchains: %d\n", numchains);
      //~ for( int i=0; i < numchains; ++i ) {
         //~ protchains[i]->print();
      //~ }
   }
   else {
      numchains = 1;
      protchains = new bPoints*[numchains];
      protchains[0] = new bPoints;
      protchains[0]->readPoints( file );
      //~ protchains[0]->print();
   }


   for( int k=0; k < numchains; ++k ) {
      printf("%d\n",protchains[k]->size());
      for( int i=0; i < protchains[k]->size(); ++i ) {
         printf("\t\t%u\n",protchains[k]->pos( i ));
      }
   }
   //~ exit(1);

   // Tessellate
   bDelTess dt;
   bool valid = dt.tessellate_full( protchains, numchains, base, path, trim, true );
   if( !valid ) { printf("Unable to tessellate '%s'\n",file); exit(1); }

   // Print out visual
   int nsize = strlen(base) + 4;
   char name[ nsize ];
   memset( name, '\0', nsize );
   memmove( name, base, nsize - 4 );
   memmove( name + nsize - 4, "_dt", 3 );
   int osize = pathsize + nsize + 4;
   char* outf = new char[ osize ];
   memset( outf, '\0', osize );
   memmove( outf, path, pathsize );
   strcat( outf, base );
   strcat( outf, ".py" );
   printf("wrote: %s\n", outf);
   op = fopen(outf, "w");
   (void)fprintf(op,"from pymol import cmd\n");
   (void)fprintf(op,"from pymol.cgo import *\n");
   (void)fprintf(op,"\n");
   if(white) { (void)fprintf(op,"cmd.bg_color(\"white\")\n"); }
   (void)fprintf(op,"cmd.set(\"auto_show_spheres\", \"on\")\n");
   (void)fprintf(op,"cmd.set(\"sphere_scale\", .25)\n");
   (void)fprintf(op,"cmd.set(\"auto_show_lines\", \"on\")\n");
   (void)fprintf(op,"cmd.set(\"label_position\",(1.5,1.5,1.5))\n");
   (void)fprintf(op,"cmd.set(\"label_size\",8)\n");
   (void)fprintf(op,"\n");
   dt.pymol( op, name );
   (void)fprintf( op, "cmd.orient(\"%s\")\n", name );
   fclose(op);

   // Free memory
   for( int i=0; i < numchains; ++i ) { delete protchains[0]; protchains[0] = NULL; }
   delete [] protchains; protchains = NULL;
   delete [] file; file = NULL;
   delete [] base; base = NULL;
   delete [] path; path = NULL;
   delete [] outf; outf = NULL;

   end = omp_get_wtime();
   printf("WRuntime: %f sec. time.\n\n", end-start);


   return valid ? 0 : 1;
}
