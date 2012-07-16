#ifndef ROD_UTIL_H
#define ROD_UTIL_H

#include "swak/Swak.h"

namespace RodUtil
{
  //template <typename Rod, int64 (*end_func_ptr)(const Rod *) = DefaultGetEnd>

  template <typename Rod,
            typename LessType,
            typename NextType,
            int (*start_func)(const Rod&),
            int (*end_func)(const Rod&)>
  class RodWalker
  {
    typedef typename set<Rod, LessType>::const_iterator RodSetIter;

  public:
    RodWalker(NextType &next_obj_);

    void JumpTo(int new_pos);

    int GetPos();

    /* 
     * GetOverlaps() return value:
     *   Returns false if this function will do nothing because there's no overlaps to check to clear out
     *   and there's no more rods we could add.  Usually this occurs after calling get overlaps after
     *   the last entry on a chrom is encountered
     */
    
    /* Assumes overlaps was cleared before first call. The overlaps set is used to maintain state
     * between calls so do not manipulate the set between calls and dont mix and match the set version
     * with the vector version of this function. We deem this the "set version" of the function. */
    bool GetOverlaps(set<Rod, LessType> &overlaps);
    
    bool HasNext() { return has_next; }

    const Rod &NextRod() { return next_rod; }


  private:
    int Start(const Rod &rod);

    int End(const Rod &rod);

    bool GetOverlapsImpl(set<Rod, LessType> &overlaps);

    int pos;
    bool has_next;
    Rod next_rod;
    NextType next_obj;
  };

};

#include "RodWalkerImpl.h"

#endif
