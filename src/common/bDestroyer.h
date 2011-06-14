#ifndef BDESTROYER_H
#define BDESTROYER_H

/// This class was gleaned from
/// http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt

/// This is a Singleton destruction class
/// Two issues:
///   1) don't try to delete before the end of the program.
///   2) hope and pray the doomed do not depend on each other.

namespace griddock { template <class DOOMED> class bDestroyer; };

template <class DOOMED>
class griddock::bDestroyer
{
   public:
      bDestroyer(DOOMED* = 0);
      ~bDestroyer();

      void setDoomed(DOOMED*);

   private:
      // Prevent users from making copies of a
      // bDestroyer to avoid double deletion:
      bDestroyer(const bDestroyer<DOOMED>&);
      void operator=(const bDestroyer<DOOMED>&);

   private:
      DOOMED* doomed_;
};

template <class DOOMED>
griddock::bDestroyer<DOOMED>::bDestroyer (DOOMED* d)
{
   doomed_ = d;
}

template <class DOOMED>
griddock::bDestroyer<DOOMED>::~bDestroyer ()
{
   delete doomed_;
}

template <class DOOMED>
void griddock::bDestroyer<DOOMED>::setDoomed (DOOMED* d)
{
   doomed_ = d;
}

#endif
