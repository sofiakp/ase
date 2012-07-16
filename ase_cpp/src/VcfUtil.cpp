#include "VcfUtil.h"

namespace VcfUtil
{
  void GetNext(vcf_file &vfile, uint entry_index, vcf_entry &entry)
  {
    Assert (entry_index < vfile.N_entries);

    string line;
    vfile.get_vcf_entry(entry_index, line);
    entry.reset(line);

    // For the locus
    bool PARSE_ALT=true;
    bool PARSE_FILTER=true;
    bool PARSE_INFO=true;
    entry.parse_basic_entry(PARSE_ALT, PARSE_FILTER, PARSE_INFO);

    // For the individuals
    bool PARSE_FORMAT=true;
    entry.parse_full_entry(PARSE_FORMAT);
  }
};
