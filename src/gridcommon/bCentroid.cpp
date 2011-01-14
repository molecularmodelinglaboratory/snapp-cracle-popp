#include <deque>
#include <fstream>
#include <sstream>
#include <string>

#include <string.h>
#include <stdio.h>

#include "bCentroid.h"
#include "bDelTess.h"
#include "bAA.h"
#include "bSys.h"

using namespace std;
using namespace bStd;

std::deque<float>  bCentroid::centroids;
std::deque<string> bCentroid::residues;
std::deque<bCentroid::uint>   bCentroid::residuesPos;
std::deque<bCentroid::uint>   bCentroid::chainPos;
std::deque<bCentroid::uint>   bCentroid::chainSize;

bCentroid::bCentroid() {
   
}
bCentroid::bCentroid( char base[], char path[] ) {
   strcpy( this->fileBase_, base );
   strcpy( this->filePath_, path );
}
bCentroid::bCentroid( const bCentroid &rhs ) {
   
}

bCentroid::~bCentroid() {
   
}

void bCentroid::setSrc( char base[], char path[] ) {
   strcpy( this->fileBase_, base );
   strcpy( this->filePath_, path );
   return;
}

bool bCentroid::pdb2Centroids() { return pdb2Centroids( this->fileBase_, this->filePath_ ); }
bool bCentroid::pdb2Centroids( char fileBase[], char filePath[], char pathOut[] ) {


   // Open the file
   ifstream ip;
   char* file = new char[96];
   strcpy( file, filePath );
   strcat( file, fileBase );
   //~ if( memcmp( file + strlen(file) - 4, ".pdb", 4 ) != 0 ) strcat( file, ".pdb" );
   //~ ip.open(file, ifstream::in);
   //~ if(!ip) { throw "[bCentroid] Bad file."; return false; }

   // Print centroids
   bool haveCent = false;
   try{ haveCent = _pdb2Centroids( fileBase, filePath ); }
   catch( const char* e ) { printf("%s\n",e); }

   if( haveCent ) {
      //~ printf("have centroids\n");

      if( pathOut != NULL ) {
         strcpy( file, pathOut );
         strcat( file, fileBase );
      }

      int size = strlen(file);
      if( memcmp( file + size - 4, ".pdb", 4 ) == 0 ) {
         //~ strcat( file, ".pdb" );
         size -= 4;
         memset( file + size, '\0', 4 );
         //~ --size;
      }

      memmove( file + size, ".out", 4 );
      //~ printf("file: %s\n", file);
      FILE* op = fopen( file, "w" );
      for( uint i=0; i < centroids.size(); ++i ) {
         fprintf( op, "%.3f ", centroids[i] );
         ++i; fprintf( op, "%.3f ", centroids[i] );
         ++i; fprintf( op, "%.3f\n", centroids[i] );
      }
      fclose( op );
      
      // Print sequence
      memmove( file + size, ".seq", 4 );
      //~ printf("file: %s\n", file);
      op = fopen( file, "w" );
      for( uint i=0; i < residues.size(); ++i ) {
         fprintf( op, "%s %d\n", residues[i].c_str(), residuesPos[i] );
      }
      fclose( op );
   }
   //~ else { printf(" no centroids\n"); }
   
   delete [] file;
   file = NULL;
   return haveCent;
}

bool bCentroid::pdb2Centroids( bPoints** &chainPtr, int &numChain, char fileBase[], char filePath[] ) {
   bool haveCent = _pdb2Centroids( fileBase, filePath );
   if( haveCent ) {
      numChain = chainSize.size();
      chainPtr = new bPoints*[ numChain ];
      uint cX = 0;
      uint pX = 0;
      uint rX = 0;
      for( int i=0; i < numChain; ++i ) {
         chainPtr[i] = new bPoints;
         //~ printf("chain pos: %d\n", chainPos[i]);
         char  seqSave[ chainSize[i] + 1 ];
         float pntSave[ 3 * chainSize[i] ];
         memset( seqSave, '\0', chainSize[i] + 1 );
         pX = 0;
         rX = chainPos[i];
         for( uint k=0; k < chainSize[i]; ++k ) {
            //~ printf("pX|cX: %d | %d [%d, %lu] \t%c\t%.2f, %.2f, %.2f\n", pX, cX, 3*chainSize[i], centroids.size(), bAA::aa3to1( residues[k].c_str() ), centroids[cX], centroids[cX+1], centroids[cX+2]);
            pntSave[pX] = centroids[cX];
            pntSave[++pX] = centroids[++cX];
            pntSave[++pX] = centroids[++cX];
            seqSave[k] = bAA::aa3to1( residues[rX].c_str() );
            //~ printf("%d: %s\n", rX, residues[rX].c_str() );
            ++pX; ++cX; ++rX;
            //~ memmove( ugh, residues[k].c_str
         }
         //~ bPoints::printPoints( pntSave, chainSize[i] );
         //~ printf("seq: %s\n", seqSave);
         chainPtr[i]->addPoints( pntSave, chainSize[i] );
         chainPtr[i]->setSeq( seqSave );
         chainPtr[i]->_findMinMax();
         //~ chainPtr[i]->print();
      }
   }
//~ printf("done\n");
   return haveCent;
}

bool bCentroid::_pdb2Centroids( char fileBase[], char filePath[] ) {
   centroids.clear();
   residues.clear();
   residuesPos.clear();
   chainPos.clear();
   chainSize.clear();

   // Make sure the path works
   char* file = new char[96];
   strcpy( file, filePath );
   strcat( file, fileBase );
   if( ! bSys::fileExists( file ) ) {
      strcat( file, ".pdb" );
      //~ printf("[bC] file: %s\n", file);
      if( ! bSys::fileExists( file ) ) { printf("not found...\n"); throw "[bC::pdb2C] Invalid File."; return false; }
   }
   //~ if( memcmp( file + strlen(file) - 4, ".pdb", 4 ) != 0 ) strcat( file, ".pdb" );

   // Open the file
   ifstream ip;
   ip.open(file, ifstream::in);
   if(!ip) { throw "[bCentroid] Bad file."; return false; }
   
   // delete memory
   delete [] file;
   file = NULL;

   // Storage
   //~ deque<float>  centroids;
   //~ deque<string> residues;
   //~ deque<uint>   residuesPos;
   //~ deque<uint>   chainPos;
   //~ deque<uint>   chainSize;
   string        bffr;

   // Read (in order)
   string type;
   int    posAtom;
   string atom;
   string res;
   char   chain;
   int    posRes;
   float  coord[3] = { 0.0, 0.0, 0.0 };

   // Data handling
   float  centr[3] = { 0.0, 0.0, 0.0 };
   int    prevPos = 0;
   bool   prevCa = false;
   bool   prevCCarbonyl = false;
   int    cntAtom = 0;
   char   currChain = '\0';
          //~ numChain = 0;
   uint   numRes = 0;

   // Line by line
   chainPos.push_back(0);
   for( uint i = 0; getline( ip, bffr ); ++i ) { 
      istringstream ss( bffr );

      // Skip unless it's an atom
      ss >> type;
      if( type.compare("ATOM") ) { continue; }

      // Read up to coordinates
      ss >> posAtom;
      ss >> atom;
      ss >> res;
      ss >> chain;
      ss >> posRes;
      //~ printf("%s\n",bffr.c_str());
      //~ if( bffr[25] != ' ') { printf("%s\n",bffr.c_str()); }
      if( bffr[26] != ' ') { continue; }
      //~ if( bffr[27] != ' ') { printf("%s\n",bffr.c_str()); }

      // Check for backbone N (new residue, handle at Ca)
      if( !atom.compare("N") && posRes != prevPos ) {
         prevPos = posRes;
         //~ printf("\t\tNITROGEN!\n");
         continue;
      }

      // Check for C-alpha (new residue!)
      else if( !atom.compare("Ca") || !atom.compare("CA") ) {
         //~ printf("\t\tALPHA CARBON!");

         // Check for new chain
         if( chain != currChain ) {
            if( currChain != '\0' ) {
               chainPos.push_back( residues.size() );
               chainSize.push_back( numRes );
               numRes = 0;
               //~ ++chain;
            }
            currChain = chain;
            //~ printf("new chain!\n");
         }

         // Save the residue
         //~ printf("%d: %s\n", numRes, res.c_str());
         ++numRes;
         residues.push_back( res );
         residuesPos.push_back( posRes );
         prevCa = true;

         if( cntAtom != 0 ) {
            // Calculate the centroid, save, & reset
            centr[0] /= cntAtom;
            centr[1] /= cntAtom;
            centr[2] /= cntAtom;
            centroids.push_back( centr[0] );
            centroids.push_back( centr[1] );
            centroids.push_back( centr[2] );
            centr[0] = 0.0;
            centr[1] = 0.0;
            centr[2] = 0.0;
            cntAtom = 0;
         }
      }

      // Check for backbone C
      else if( !atom.compare("C") && prevCa ) {
         prevCa = false;
         prevCCarbonyl = true;
         //~ printf("\t\tCARBONYL!\n");
         continue;
      }

      // Check for backbone C
      else if( !atom.compare("O") && prevCCarbonyl ) {
         prevCCarbonyl = false;
         //~ printf("\t\tCARBONYL!\n");
         continue;
      }

      // Ignore extra carboxylate atom
      else if( !atom.compare("OXT") ) {
         continue;
      }

      else {}
      //~ printf("\n");

      // Weight
      float weight = 1.0;

      // Save the coordinates
      cntAtom += weight;
      ss >> coord[0];
      ss >> coord[1];
      ss >> coord[2];
      coord[0] *= weight;
      coord[1] *= weight;
      coord[2] *= weight;

      centr[0] += coord[0];
      centr[1] += coord[1];
      centr[2] += coord[2];
   }
   ip.close();

   // Save last centroid
   if( numRes != 0 ) {
      chainPos.push_back( residues.size() );
      chainSize.push_back( numRes );
   }


   if( cntAtom != 0 ) {
      // Calculate the centroid, save, & reset
      centr[0] /= cntAtom;
      centr[1] /= cntAtom;
      centr[2] /= cntAtom;
      centroids.push_back( centr[0] );
      centroids.push_back( centr[1] );
      centroids.push_back( centr[2] );
   }

   return (centroids.size() > 0) ? true : false;
}


bCentroid::uint bCentroid::dt2Centroids( bDelTess &dt ) {

   //~ deque<float> centroids;
   centroids.clear();
   float centroid[3] = { 0.0, 0.0, 0.0 };
   for( uint i=0; i < dt.simplID_.size(); ++i ) {
      if( dt.skip( i ) ) { continue; }
      
      int *a = dt.simplID_[i].getActive(); // get points
      ushort chn[4] = { dt.vrtx_[a[0]].chn_, dt.vrtx_[a[1]].chn_, dt.vrtx_[a[2]].chn_, dt.vrtx_[a[3]].chn_ };
      ushort pos[4] = { dt.vrtx_[a[0]].pos_, dt.vrtx_[a[1]].pos_, dt.vrtx_[a[2]].pos_, dt.vrtx_[a[3]].pos_ };
      pos[0] *= 3; pos[1] *= 3; pos[2] *= 3; pos[3] *= 3;


      centroid[0] = dt.src_[chn[0]]->pnts_[pos[0]];
      centroid[0] += dt.src_[chn[1]]->pnts_[pos[1]];
      centroid[0] += dt.src_[chn[2]]->pnts_[pos[2]];
      centroid[0] += dt.src_[chn[3]]->pnts_[pos[3]];
         
      centroid[1] = dt.src_[chn[0]]->pnts_[++pos[0]];
      centroid[1] += dt.src_[chn[1]]->pnts_[++pos[1]];
      centroid[1] += dt.src_[chn[2]]->pnts_[++pos[2]];
      centroid[1] += dt.src_[chn[3]]->pnts_[++pos[3]];
         
      centroid[2] = dt.src_[chn[0]]->pnts_[++pos[0]];
      centroid[2] += dt.src_[chn[1]]->pnts_[++pos[1]];
      centroid[2] += dt.src_[chn[2]]->pnts_[++pos[2]];
      centroid[2] += dt.src_[chn[3]]->pnts_[++pos[3]];

      centroid[0] /= 4;
      centroid[1] /= 4;
      centroid[2] /= 4;

      centroids.push_back( centroid[0] );
      centroids.push_back( centroid[1] );
      centroids.push_back( centroid[2] );
   }

   return centroids.size() / 3;
}

bCentroid::uint bCentroid::dt2Centroids( bDelTess &dt, bPoints* &tet ) {
   uint numCentroids = dt2Centroids( dt );

   // Save centroids
   if( tet == NULL ) { tet = new bPoints; }
   else { tet->clear(); }
   float pntSave[ 3 * numCentroids ];
   uint pX = 0;
   for( uint k=0; k < numCentroids; ++k ) {
      pntSave[pX] = centroids[pX]; ++pX;
      pntSave[pX] = centroids[pX]; ++pX;
      pntSave[pX] = centroids[pX]; ++pX;
   }
   tet->addPoints( pntSave, numCentroids );

   return numCentroids;
}

bCentroid::uint bCentroid::dt2Centroids( bDelTess &dt, char base[], char path[] ) {
   uint numCentroids = dt2Centroids( dt );

   // Print centroids
   char file[64];
   strcpy( file, path );
   strcat( file, base );
   strcat( file, "_CENTROIDS.out" );
   FILE* op = fopen( file, "w" );
   for( uint i=0; i < centroids.size(); ++i ) {
      fprintf( op, "%.4f ", centroids[i] );
      ++i; fprintf( op, "%.4f ", centroids[i] );
      ++i; fprintf( op, "%.4f\n", centroids[i] );
   }
   fclose( op );

   return numCentroids;
}


