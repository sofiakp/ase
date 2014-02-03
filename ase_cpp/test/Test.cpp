#include "swak/Swak.h"
#include "swak/PolyExe.h"
#include "swak/System.h"
#include "swak/Helpers.h"
#include "swak/StlUtil.h"

int main(int argc, char ** argv)
{
  DeclarePolyExe(main, test_rod_util, "Test ROD utils");

  if (argc <= 1 || string(argv[1]) == string("-h") || !PolyExe::HaveMain(argv[1]))
  {
    string usage_line = Basename(argv[0]) + " <command> [options] <args>";
    string desc = "";
    PrintUsage(usage_line, desc, PolyExe::MainNamesVec(), 100, "Commands");
    return 0;
  }

  return PolyExe::CallMain(argc, argv);
}
