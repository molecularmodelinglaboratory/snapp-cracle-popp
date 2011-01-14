
#include <stdlib.h>
#include <stdio.h>
#include "bVBase.h"

#ifndef _SLASH
#define _SLASH
#if defined(_WIN32) || defined(_WIN64)
# error WIN not supported yet.
   const char _HOME_[13] = "%UserProfile%";
   const char _SLASH_    = '\\';
#else
   const char _HOME_[6] = "$HOME";
   const char _SLASH_    = '/';
#endif
#endif

#ifndef _MAXBIT
#define _MAXBIT
#if defined(_M_X64) || defined(__amd64__) || defined(_LP64)
   const short _MAXBIT_ = 64;
#else
# error 32-bit not supported yet.
   const short _MAXBIT_ = 32;
#endif
#endif

#ifndef PI
#define PI
   const double _PI_ = 3.14159265358979323846264338327950288;
#endif


using namespace std;
using namespace bStd;

/* Write PyMol Header */
void bVBase::pymolHeader(FILE* op,int white) {
   (void)fprintf(op,"from pymol import cmd\n");
   (void)fprintf(op,"from pymol.cgo import *\n");
   (void)fprintf(op,"\n");
   if(white) {
      (void)fprintf(op,"cmd.bg_color(\"white\")\n");
   }
   (void)fprintf(op,"\n");
   
   // set parameters
   (void)fprintf(op,"cmd.set(\"auto_show_spheres\", \"on\")\n");
   (void)fprintf(op,"cmd.set(\"sphere_scale\", .25)\n");
   (void)fprintf(op,"cmd.set(\"auto_show_lines\", \"on\")\n");
   (void)fprintf(op,"cmd.set(\"label_position\",(1.5,1.5,1.5))\n");
   (void)fprintf(op,"cmd.set(\"label_size\",8)\n");
   (void)fprintf(op,"\n");

   return;
}