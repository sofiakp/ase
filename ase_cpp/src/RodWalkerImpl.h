#ifndef ROD_WALKER_IMPL_H
#define ROD_WALKER_IMPL_H

namespace RodUtil
{
  //#define rod_walker_template template <typename Rod, typename RodIterType, int64 (*start_func_ptr)(const Rod *), int64 (*end_func_ptr)(const Rod *)>
  //#define rod_walker_scope RodWalker<Rod, RodIterType, start_func_ptr, end_func_ptr>

  //#define rod_walker_template template <typename Rod, typename LessType, typename RodSet>
  //#define rod_walker_scope RodWalker<Rod, LessType, RodSet>

  #define rod_walker_template template <typename Rod, typename LessType, typename NextType, int (*start_func)(const Rod&), int (*end_func)(const Rod&)>
  #define rod_walker_scope RodWalker<Rod, LessType, NextType, start_func, end_func>

  // vvvvv  Constructors  vvvvv

  rod_walker_template
  rod_walker_scope::RodWalker(NextType &next_obj_) : next_obj(next_obj_)
  {
    pos = -1;
  }

  // ^^^^^  Constructors  ^^^^^


  // vvvvv  Utility fxns  vvvvv

  rod_walker_template
  void rod_walker_scope::JumpTo(int new_pos)
  {
    if (new_pos < pos || new_pos < 0)
      throw runtime_error("In RodWalker::JumpTo(new_pos), new_pos must be >= than current pos and non-negative! " + sprint(pos) + sprint(new_pos));

    // Initialize if this is the first call to jumpto
    if (pos < 0)
      has_next = next_obj.GetNext(next_rod);
    pos = new_pos;
  }

  rod_walker_template
  int rod_walker_scope::GetPos()
  {
    return pos;
  }

  rod_walker_template
  int rod_walker_scope::Start(const Rod &rod)
  {
    return (*start_func)(rod);
  }

  rod_walker_template
  int rod_walker_scope::End(const Rod &rod)
  {
    return (*end_func)(rod);
  }


  // ^^^^^  Utility fxns  ^^^^^

  // vvvvv  Long fxns  vvvvv

  rod_walker_template
  bool rod_walker_scope::GetOverlapsImpl(set<Rod, LessType> &overlaps)
  {
    AssertMsg(pos >= 0, "Must call JumpTo() before GetOverlaps()");

    if (overlaps.empty() && !has_next)
      return false;

    // Remove things from overlaps if we hae moved past them
    //cerr << "checking if we should remove any..." << endl;
    for (RodSetIter oiter = overlaps.begin(); oiter != overlaps.end(); )
    {
      int end = End(*oiter);
      //cerr << "Pos: " << pos << " End: " << end << endl;

      if (end <= pos)
      {
        //cerr << "Removing "  << endl;
        overlaps.erase(oiter++);
      }
      else
        break;
    }

    // Add anything new that we have encountered
    //cerr << "checking if we should add any..." << endl;
    while (has_next)
    {
      int start = Start(next_rod);
      int end = End(next_rod);

      if (start > pos)
      {
        //cerr << "Too early: (won't add any)" << endl;
        //cerr << "@" << pos << ": Pausing bc " << new_rod_ptr->Name << ". Starts at " << start << " (" << new_rod_ptr->Position << ")" << endl;
        break;
      }

      //cerr << "Pos: " << pos << " Start: " << start << " End: " << end << endl;
      if (pos >= start && pos < end)
      {
        //cerr << "Adding" << endl;
        //cerr << "@" << pos << ": Adding " << new_rod_ptr->Name << ". Starts at " << start << " (" << new_rod_ptr->Position << ")" << endl;
        overlaps.insert(next_rod);
      }

      has_next = next_obj.GetNext(next_rod);
    }

    //cerr << "@" << pos << ": Overlaps size now: " << overlaps.size() << endl;

    return true; // At the least we at least did something (checked whether to remove overlaps or to add new things)
  }

  rod_walker_template
  bool rod_walker_scope::GetOverlaps(set<Rod, LessType> &overlaps)
  {
    return GetOverlapsImpl(overlaps);
  }

};

#endif
