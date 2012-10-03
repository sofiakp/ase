#include "swak/Swak.h"
#include "swak/PolyExe.h"
#include "swak/System.h"
#include "swak/Helpers.h"
#include "swak/StlUtil.h"

int main(int argc, char ** argv)
{
    DeclarePolyExe(main, separatedbSNP, "Merely separate dbSNP into sub-datasets according to chromesomes to facilitate further analysis.");
    DeclarePolyExe(main, LDDataPreProcess, "Preprocess LD data from HapMap.");
    DeclarePolyExe(main, splitLD, "Split LD datasets to prepare for randomization.");
    DeclarePolyExe(main, mapIDAlleles, "Map the information for the AS-SNPs");
    DeclarePolyExe(main, annotateeQTLs, "Map the information for the eQTLs");
    DeclarePolyExe(main, processGWASCatalog, "Parse the GWASCatalog and pre-filter redundant records.");
    DeclarePolyExe(main, generateFairBinAssociationSNPs, "Generate a fair GWAS set or eQTL set according to LD structure and bin this set for randomized matching.");
    DeclarePolyExe(main, computeOverlapping, "Compute the number of overlapping SNPs without filtering.");
    DeclarePolyExe(main, generateBinFilterNULLSet, "Generate the matched randomized null SNP set.");
    DeclarePolyExe(main, filterASSNPs, "Filter the allele-specific SNPs and compute the size of the filtered set.");
    DeclarePolyExe(main, computeStatistics, "Compute enrichments compared to null set.");
    if (argc <= 1 || string(argv[1]) == string("-h") || !PolyExe::HaveMain(argv[1]))
    {
        string usage_line = Basename(argv[0]) + " <command> [options] <args>";
        string desc = "";
        PrintUsage(usage_line, desc, PolyExe::MainNamesVec(), 100, "Commands");
        return 0;
    }
    
    return PolyExe::CallMain(argc, argv);
}
