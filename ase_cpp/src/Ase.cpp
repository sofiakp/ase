#include "swak/Swak.h"
#include "swak/PolyExe.h"
#include "swak/System.h"
#include "swak/Helpers.h"
#include "swak/StlUtil.h"

int main(int argc, char ** argv)
{
  DeclarePolyExe(main, reconcile, "Chooses the best alignment between mappings to parallel genomes (eg. b6 and cast)");
  DeclarePolyExe(main, asequantmultirg, "Counts alleles for one BAM file with multiple RGs");
  DeclarePolyExe(main, asequantmultibam, "Counts alleles for multiple BAM files");
  DeclarePolyExe(main, aseregion, "Counts reads for each region in a BED file");

  if (argc <= 1 || string(argv[1]) == string("-h") || !PolyExe::HaveMain(argv[1]))
  {
    string usage_line = Basename(argv[0]) + " <command> [options] <args>";
    string desc = "";
    PrintUsage(usage_line, desc, PolyExe::MainNamesVec(), 100, "Commands");
    return 0;
  }

  return PolyExe::CallMain(argc, argv);
}
