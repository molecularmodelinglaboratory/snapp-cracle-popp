#ifndef BLIST_H
#define BLIST_H

namespace bStd { class bNode, bList; }

template <class T>
class bStd::bNode {
   typedef bList<T> L;
   
   private:
   T* self_;
   T* next_;
   T* back_;
   T* pare_;
   L* chld_;
   L* list_;
   
   public:
   bNode();
   bNode( const bNode<T>& );
   bNode( const T& );
   ~bNode();
   
   bool operator==( const bNode<T>& );
   bool operator!=( const bNode<T>& );
   bool operator< ( const bNode<T>& );
   bool operator<=( const bNode<T>& );
   bool operator> ( const bNode<T>& );
   bool operator>=( const bNode<T>& );
   
   bNode<T>& back();
   bNode<T>& next();
   bNode<T>& pare();
   L& chld();
   
   void add_back( T* );
   void add_next( T* );
   void add_pare( T* );
   void add_chld( T* );
   void add_chld( L* );
};


/***** bNode *********************************************/
template <class T>
bNode<T>::bNode() :
   self_( NULL ), next_(NULL), back_(NULL), pare_(NULL), child_(NULL)
{
};

template <class T>
bNode<T>::bNode( const bNode<T>& rhs ) :
   self_( NULL ), next_(NULL), back_(NULL), pare_(NULL), child_(NULL)
{
   this->self_ = new T( *(rhs.self_) );
   this->next_ = rhs.next_;
   this->back_ = rhs.back_;
   this->pare_ = rhs.pare_;
   this->chld_ = new L( *(ths.chld_) );
};

template <class T>
bNode<T>::bNode( const T& rhs ) :
   self_( NULL ), next_(NULL), back_(NULL), pare_(NULL), child_(NULL)
{
   this->self_ = new T( rhs );
};

template <class T>
bNode<T>::~bNode() {
   if( this->self_ != NULL ) { delete this->self_; this->self_ = NULL; }
   this->next_ = NULL;
   this->back_ = NULL;
   this->pare_ = NULL;
   if( this->chld_ != NULL ) { delete this->chld_; this->chld_ = NULL; }
}

template <class T> bool bNode<T>::operator==( const bNode<T>& rhs ) { return this->self_ == rhs.self_; }
template <class T> bool bNode<T>::operator!=( const bNode<T>& rhs ) { return !(this->self_ == rhs.self_;) }
template <class T> bool bNode<T>::operator< ( const bNode<T>& rhs ) { return this->self_ <  rhs.self_; }
template <class T> bool bNode<T>::operator<=( const bNode<T>& rhs ) { return !(this->self_ >  rhs.self_;) }
template <class T> bool bNode<T>::operator> ( const bNode<T>& rhs ) { return this->self_ >  rhs.self_; }
template <class T> bool bNode<T>::operator>=( const bNode<T>& rhs ) { return !(this->self_ <  rhs.self_;) }

template <class T> bNode<T>& bNode<T>::next() { return *(this->next_); }
template <class T> bNode<T>& bNode<T>::back() { return *(this->back_); }
template <class T> bNode<T>& bNode<T>::pare() { return *(this->pare_); }
template <class T> bNode<T>& bNode<T>::chld() { return *(this->chld_); }

template <class T> void bNode<T>::add_next( T* rhs ) { this->next_ = rhs; return; }
template <class T> void bNode<T>::add_back( T* rhs ) { this->back_ = rhs; return; }
template <class T> void bNode<T>::add_pare( T* rhs ) { this->pare_ = rhs; return; }
template <class T> void bNode<T>::add_chld( L* rhs ) { this->chld_ = rhs; return; }
template <class T> void bNode<T>::add_chld( T* rhs ) {
   if( this->chld_ == NULL ) { this->chld_ = new L; }
   this->chld_->add( rhs );
   return;
}

/***** bList *********************************************/
template <class T>
class bStd::bList {
   typedef bNode<T> N;
   
   bNode<T>*  root_;
   bNode<T>** list_;
   
   size_t cap_;
   size_t num_;
   
   bList();
   bList( const bList& );
   bList( const T& );
   ~bList();
   
   T* operator[]( uint );
   
   bNode<T>* begin();
   bNode<T>* end();
   T* front();
   T* back();
   
   void add( T& );
   void push_back( T& );
   void push_front( T& );

   void add( T* );
   void push_back( T* );
   void push_front( T* );

   void del( uint );
   void pop_back();
   void pop_front();
   
   size_t size();
   void resize( size_t );
};

template<T>
bList<T>::bList() :
   root_(NULL), list_(NULL), cap_(0), num_(0)
{
}

template<T>
bList<T>::bList( const bList<T>& rhs ) :
   root_(NULL), list_(NULL), cap_(0), num_(0)
{
   this->cap_ = rhs.cap_;
   this->num_ = rhs.num_;
   this->list_ = new bNode<T>*[ this->cap_ ];
   
   for( size_t i=0; i < rhs.num_; ++i ) {
      this->push_back( rhs[i] );
   }
   
   for( size_t i=this->nuM_; i < this->cap_; ++i ) {
      this->list_[i] = NULL;
   }

   this->root_ = this->list_[0];
}

template<T>
bList<T>::bList( const T& rhs ) :
   root_(NULL), list_(NULL), cap_(5), num_(0)
{
   this->list_ = new bNode<T>*[ this->cap_ ];
   this->push_back( &rhs );
   this->root_ = this->list_[0];
}

template<T>
bList<T>::~bList() {
   for( size_t i=0; i < this->num_; ++i ) {
      delete this->list_[i];
      this->list_[i] = NULL;
   }
   delete [] this->list_;
   this->list_ = NULL;

   this->root_ = NULL;
   this->cap_ = 0;
   this->num_ = 0;
}

template<T> T* bList<T>::operator[]( uint pos ) { return this->list_[pos]->self_; }

template<T> bNode<T>* bList<T>::begin() { return this->list_[0]; }
template<T> bNode<T>* bList<T>::end  () { return this->list_[(this->num_ - 1)]->self_; }

template<T> T* bList<T>::front() { return this->list_[0]->self_; }
template<T> T* bList<T>::back () { return this->list_[(this->num_ - 1)]->self_; }

template<T>
void bList<T>::add( T& rhs ) {
   if( this->num_ == this->cap_ ) { this->resize( this->cap_ + 5 ); }
   
   size_t num = this->num_;
   this->list_[num] = new bNode<T>( T );
   bNode<T>* curr = this->list_[num];
   
   for( size_t i=0; i < num; ++i ) {
      if( *(this->list_[i]) > rhs ) {
         // need to check for sorting
         // need to figure out a good way to insert sorted...
      }
   }
   
   return;
}


#endif