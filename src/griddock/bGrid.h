#ifndef BGRID_H
#define BGRID_H

#include <bPoints.h>
#include "bStamp.h"
#include "bHex.h"

namespace bStd { class bGrid; };

class bStd::bGrid {
   typedef unsigned int uint;
   friend class bPoints;
   friend class bDocker;
   
public:
   // Construct | Destruct
   bGrid();
   bGrid( const bGrid& );
   ~bGrid();
   void clear();

   // Assignment
   void copyBasicParam( const bGrid& );
   bGrid& operator=( const bGrid& );
   bGrid& operator-=( const bGrid& );
   bGrid& operator+=( const bGrid& );
   bGrid& operator/=( const bGrid& );

   // Validation
   bool sameDim( const bGrid& );
   void overlay( const bGrid& );

   // Command Line Arguements
   bool cla(int,char**);
   void setParam( float[] );

   // Protein
   void setPointObject( bPoints& );

   // Grid and Stamps
   void initializeGrid();
   bool initializeStamps();

   // Create the Grid
   bool stampPoints();
      bool _stampPoints( bStamp&, float*, int );
      void _stampPoint( bStamp&, float[] );
      void transposeIndex( int&, int& );

   // PyMol Output
   void pymol3dGrid( FILE*, const char[] =("grid"), const char[] =("blue"), bool =(false) );
   void pymol3dGrid( FILE*, bHex*, int, int, const char[], const char[], float );
   
   // Debugging
   void printFStmp();
   void printTStmp();
   void printProtein();
   void print2dProtein();
   void print2dGrid( bGrid*, int);

   // Toolkit
   uint countNearby( uint[], FILE* );
   void removeInternal( FILE* );
   bool isaGridPoint( float[] ) const;

protected:

   // Flags
   bool haveFStmp_;
   bool haveTStmp_;
   bool isStamped_; // have we stamped the grid?

   // Dimensions
   int size_;
   int length_;
   int height_;
   int depth_;
   int min_[3]; // min index of the grid (should always be 0 in resolution)
   int max_[3]; // max index of the grid (note: index, not count)

   // Data
   bPoints* prt_; // list of protein coordinates
   bHex* grd_;
   bStamp fStmp_;
   bStamp tStmp_;

   // Parameters
   float res_; // grid resolution (in angstroms)
   float fit_; // exclusion radius
   float thk_; // inclusion radius
   float offset_;
   float buffer_;

private:

};


#endif

