#ifndef SWAK_POLY_EXE_H
#define SWAK_POLY_EXE_H

#include "Swak.h"

#define DeclarePolyExe(prefix, x, desc) do { extern int prefix##_##x(const vector<string> &); \
  PolyExe::__main_map[#x] = &prefix##_##x; PolyExe::__main_vec.push_back(#x); PolyExe::__main_descs.push_back(desc); } while(false)

namespace PolyExe
{
  extern map<string, int (*)(const vector<string> &)> __main_map;
  extern vector<string> __main_vec;
  extern vector<string> __main_descs;

  bool HaveMain(const string &exe);

  int NumMains();

  string MainNames(int i);
  string MainDescs(int i);

  VecS MainNamesVec();
  VecS MainDescsVec();

  // Assumes argv[0] is the name of the parent executable and argv[1] is name of sub-executable (main)
  int CallMain(int argc, char ** argv);
};

#endif
