#include "swak/PolyExe.h"

namespace PolyExe
{
  map<string, int (*)(const vector<string> &)> __main_map;
  vector<string> __main_vec;
  vector<string> __main_descs;

  bool HaveMain(const string &exe)
  {
    return __main_map.find(exe) != __main_map.end();
  }

  int NumMains()
  {
    return __main_vec.size();
  }

  string MainNames(int i)
  {
    if (i >= __main_vec.size() || i < 0)
    {
      cerr << "Illegal main index" << i << " when there are only " << NumMains() << endl;
      exit(1);
    }
    return __main_vec[i];
  }

  string MainDescs(int i)
  {
    if (i >= __main_descs.size() || i < 0)
    {
      cerr << "Illegal main index" << i << " when there are only " << NumMains() << endl;
      exit(1);
    }
    return __main_descs[i];
  }

  // Assumes argv[0] is the name of the parent executable and argv[1] is name of sub-executable (main)
  int CallMain(int argc, char ** argv)
  {
    if (argc < 2)
    {
      cerr << "Expected argc to be >= 2" << endl;
      return 1;
    }

    if (!HaveMain(argv[1]))
    {
      cerr << "No such main exists" << endl;
      return 1;
    }

    vector<string> args;
    args.assign(argv, argv + argc);
    return (*__main_map[argv[1]])(args);
  }

  VecS MainNamesVec()
  {
    return __main_vec;
  }

  VecS MainDescsVec()
  {
    return __main_descs;
  }
};
