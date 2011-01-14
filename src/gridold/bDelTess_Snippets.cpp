      /* IDENTIFY NEIGHBORS -- from function tessellate */
      for( int i=0; i < newList.size(); ++i ) {
         
         int nghbrCnt = 0;
         bHex test( this->simplID_[ newList[i] ] );
         //~ printf("neighboring %d\n",newList[i]);
         
         // Check for neighbors amongst other new simplexes
         for( int k=(i + 1); k < newList.size() && nghbrCnt < 4; ++k ) {
            test &= this->simplID_[ newList[k] ];
            if( test.getFlipped() == 3 ) {
               while( this->simplex_[ newList[i] ].nghbr_[ nghbrCnt ] != NULL ) { ++nghbrCnt; }
               this->simplex_[ newList[i] ].nghbr_[ nghbrCnt ] = &( this->simplex_[ newList[k] ] );
               ++nghbrCnt;
               for( int m=0; m < 4; ++m ) {
                  if( this->simplex_[ newList[k] ].nghbr_[m] == NULL ) {
                     this->simplex_[ newList[k] ].nghbr_[m] = &(this->simplex_[ newList[i] ]);
                     m = 4;
                  }
               }
            }
            test = this->simplID_[ newList[i] ];
         }
         
         // Check for neighbors amongst existing simplexes
         for( int k=0; k < wellFlip && nghbrCnt < 4; ++k ) {
            test &= this->simplID_[ well[k] ];
            if( test.getFlipped() == 3 ) {
               while( this->simplex_[ newList[i] ].nghbr_[ nghbrCnt ] != NULL ) { ++nghbrCnt; }
               this->simplex_[ newList[i] ].nghbr_[ nghbrCnt ] = &( this->simplex_[ well[k] ] );
               ++nghbrCnt;
               for( int m=0; m < 4; ++m ) {
                  if( this->simplex_[ well[k] ].nghbr_[m] == NULL ) {
                     this->simplex_[ well[k] ].nghbr_[m] = &(this->simplex_[ newList[i] ]);
                     m = 4;
                  }
               }
            }
            test = this->simplID_[ newList[i] ];
         }
         //~ status( newList[i] );
      } // end neighbors

            /*
      deque<int> newList;
      //~ printf("-- NEW --\n");
      for( int i=0; i < this->newID_.size(); ++i ) {
         this->newID_[i] |= this->vrtx_.size() - 1; // Add the latest point
         
         try {
            int pos = this->addSimplex( this->newID_[i] );
            newList.push_back( pos );
         }
         catch( int e ) {
            isCoPlanar = true;
            for( int k=0; k < newList.size(); ++k ) {
               this->delSimplex( newList[k] );
            }
            break;
         }
         //~ printf("p|i: %d | %d\t", pos, this->simplex_[pos].id_ );
         //~ this->newID_[i].print();
         //~ printf("new simplex: %d\n", newList.back());
      }
      
      // If co-planar, reset
      if( isCoPlanar ) {
         this->toAdd_.push_back( this->vrtx_.back() );
         this->vrtx_.pop_back();
         continue;
      }
      */