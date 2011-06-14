#ifndef BPOINTS_H
#define BPOINTS_H

#include <vector>
#include <deque>
#include <fstream>

#include "bVBase.h"



namespace bStd { class bPoints; class bPVar; };

class bStd::bPoints : virtual bStd::bVBase {
   typedef float* fp;
   //~ typedef unsigned int uint;

   friend class bGrid;
   friend class bDocker;
   friend class bDelTess;
   friend class bSort; // extraneous?
   friend class bSimplex;
   friend class bMC;
   friend class bCentroid;
   friend class bConf;
   
public:
   /* Constructors */
   bPoints( int =(0) );
   bPoints( const bPoints& );
   ~bPoints( );

   /* Overloaded Operators */
   bPoints& operator=( const bPoints& );

   /* Object Manipulation */
   void clear();
   int  size( int =(0) );
   void resize( int );
   void resizeCopy( int );

   /* Data Handling */
   bool haveFile();
   bool cla( int, char** );
   bool setSrc( const char[], const char[] =("../../grid/"), bool =(false) );
   bool setTets();
   void setGridParam( float, float, float, float, float );
   static bool fileExists( char*, char* );
   static bool fileExists( char* );
   static bool dirExists( char* );
      void _saveFilename( char[], char[] );
      void _saveFilenameNoExt( char[], char[] );

   /* Data Input */
   bool readPoints( );
   //~ bool readPoints( char* );
   //~ bool readPoints( char*, char* );
   int readPoints( char* );
   static int _readPoints( char*, float* & );
   void addPoint( float[], int =(-1) );
   void addPointAtPos( float[], int );
   void addPoints( float[], int, bool =(false) );
   void getSeq();
   void setSeq( char[] );
   void prep();

   /* Spatial Considerations */
          bool isInDockingSpace();
   void sizeAndSpace( const bPoints&, int =(0) );
          bool setToDockingSpace();
   static bool setToDockingSpace( bPoints& );
   static bool setToDockingSpace( std::deque<bPoints*>&, int );
          bool setToNormalSpace();
   static bool setToNormalSpace( bPoints& );
   static bool setToNormalSpace( std::deque<bPoints*>&, int );
   bool _translatePointPlane( );
   bool _changePointSpace( );
   //~ void set2SameSpace( const bPoints& );

   /* Measurement */
   void _findMinMax();
   bool findNearbyPoints( const bPoints&, bPoints&, bool =(false) );
   static float pointDistance( float*, float* );
   double rmsd( const bPoints& );
   static double rmsd( const bPoints&, const bPoints & );

   /* Rotation Functions */
   static void rotate( float[], const double[] );
   static bool rotateTheta( float[], int );
   static bool rotatePhi( float[], int );
   static bool rotateThetaPhi( float[], int, int );
   static bool rotatePhiTheta( float[], int, int );

   /* Conversion Functions */
   static bool c2s( float* );
   static bool s2c( float* );

   /* PyMol */
   void pymol( FILE*, char[], char[] =("ruby"), int =(1), char[] =(NULL) ) const;
   void pymol_dc( FILE*, char[], char[] =("ruby"), int =(1), char[] =(NULL) ) const;
   static void pymolPseudoatoms( FILE*, bPoints**, int, char[], char[], char** =(NULL) );
          void pymolPseudoatoms( FILE*, char[], char[], int =(1), char[] =(NULL) ) const ;
   static void pymolPseudoatoms( FILE*, float*, int, char[], char[], char[] =(NULL), int =(1), char[] =(NULL) );

   static void pymolConnectedPseudoatoms( FILE*, bPoints**, int, char[], char[], char** =(NULL) );
          void pymolConnectedPseudoatoms( FILE*, char[], char[], int =(1), char[] =(NULL) ) const;
   static void pymolConnectedPseudoatoms( FILE*, float*, int, char[], char[], char[] =(NULL), int =(1), char[] =(NULL) );

   static void pymolPoints( FILE*, float*, int, char[] =( "points" ) );//, float*, float );
   static void pymolPoints( FILE*, float*, int, char[], char[], float );//, float*, float );

   static void pymolConnectedPoints( FILE*, float*, int, char[] =( "points" ) );//, float*, float );
   static void pymolConnectedPoints( FILE*, float*, int, char[], char[], float );//, float*, float );

   // Print
   void print( FILE* =(stdout) ) const;
   static void printPoint( const int* );
   static void printPoint( const float* );
   static void printPoints( const float*, int, FILE* =(stdout) );

protected:
   // Flags
   bool havePnts_;
   bool haveTets_;
   bool isInRes_;
   bool isInPos_;
   bool haveMM_;
   bool haveSeq_;

   // Counters
   int numPnts_;
   int capPnts_;

   // Data
   char* pntBase_;
   char* pntPath_;
   float min_[3];
   float max_[3];
   float planeDisplacement_[3]; // amount required to translate to (+) plane
   fp pnts_;
   bPoints* tets_;
   uint* pos_;
   char* aaSeq_;

   // Parameters
   float res_; // grid resolution (in angstroms)
   float fit_; // exclusion radius
   float thk_; // inclusion radius
   float offset_;
   float buffer_;

private:

};

/********************************************************************
   bPVar
********************************************************************/

class bStd::bPVar {
   friend class bPoints;

   private:

   uint pos_;
   char var_;
   char res_;
   double pnt_[];

   public:
   bPVar();
   bPVar(  double[], uint =(0), char =(NULL), char =(NULL) );
   ~bPVar();

   void set(  double[], uint =(0), char =(NULL), char =(NULL) );
   void setPos( uint );
   void setVar( char );
   void setRes( char );
   void setPnt( double[] );

   uint pos();
   char var();
   char res();
   void pnt( double[] );
};




#endif
