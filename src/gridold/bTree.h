#ifndef BTREE_H
#define BTREE_H

/* Documentation

   bTree is a hash / tree data structure. It may be traversed in the typical
   node-to-node fashion or by levels of the tree. Using the bracket operators
   will return either the first Child on the indicated level (int) or the
   Child of the appropriate node from the very first level (bHex).
   
   Child are structures which point to the first Child of a given parent,
   the next and previous Child on the level (circular, with a bool
   indicating the first instance), and the Child of the first node on
   the current level. Using the bracket operators will return the node
   indicated (int). 
   Notes: [] may potentially loop over Child on the level
   
   Nodes correspond to typical nodes and point to their parent, next and
   previous Child (circular, with a bool indicating the first instance),
   and their Child.

*/

#include "bHex.h"

namespace bStd { class bTree; class bChild; class bHexNode; };
typedef unsigned short ushort;

/* bTree */
template <class T>
class bTree {

private:
   bChild* root_;
   ushort lvl_;
   ushort chl_;
   ushort ttl_;


public:
   bTree();
   bTree( const bTree& );
   ~bTree();
   
   bChild& operator[]( const int );
   bChild& operator[]( const T& );
   
   ushort insert( const T& );
};

/* bChild */
class bChild {
private:
   bNode** list_;
   bChild* lf_;
   bChild* rt_;
   bChild* up_;
   bChild* dn_;
   
   ushort cap_;
   ushort num_;
   
   bool first_;

public:
   bChild();
   bChild( const bChild& );
   ~bChild();
   
   bHexNode& operator[]( const int );
};

/* bNode */
class bNode {
private:
   bNode* lf_;
   bNode* rt_;
   bNode* up_;
   bChild* dn_;

   bool first_;

public:
   bNode();
   bNode( const bNode& );
   ~bNode();

   bNode& prev();
   bNode& next();
   bNode& prnt();
   bChild& chld();
}

/* bCharNode */

/* bHexNode */
class bHexNode : public bNode {
private:
   bHex key_;

public:
   bHexNode();
   bHexNode( const bHexNode& );
   ~bHex();


}



#endif
