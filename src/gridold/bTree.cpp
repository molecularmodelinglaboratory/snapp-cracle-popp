
#include bTree.h

using namespace std;
using namespace bStd;




/*********************** bChild */
bChild::bChild() :
   list_( new bNode*[2] ), this->cap_( 2 ),
   lf_(NULL), rt_(NULL), up_(NULL), dn_(NULL)
{
}

bChild::bChild( const bChild& rhs ) :
   list_(NULL), this->cap_(0),
   lf_(NULL), rt_(NULL), up_(NULL), dn_(NULL)
{
}

bChild::~bChild() {
   if( this->list_ != NULL ) {
      for( int i=0; i < this->cap_; ++i ) {
         delete this->list_[i];
         this->list_[i] = NULL;
      }
      delete [] this->list_;
      this->list_ = NULL;
   }
}

/*********************** bNode */
bNode::bNode() :
   lf_(NULL), rt_(NULL), up_(NULL), dn_(NULL)
{
}
bNode::bNode( const bNode &rhs )
   lf_(NULL), rt_(NULL), up_(NULL), dn_(NULL)
{
}
bNode::~bNode() {
   this->lf_ = NULL;
   this->rt_ = NULL;
   this->up_ = NULL;
   this->dn_ = NULL;
}

/*********************** bHexNode */
bHexNode::bHexNode() :
   lf_(NULL), rt_(NULL), up_(NULL), dn_(NULL)
{
}
bHexNode::bHexNode( const bHexNode &rhs ) : 
   bNode(rhs), lf_(NULL), rt_(NULL), up_(NULL), dn_(NULL)
{
}
bHexNode::~bHexNode() {
   this->lf_ = NULL;
   this->rt_ = NULL;
   this->up_ = NULL;
   this->dn_ = NULL;
}
