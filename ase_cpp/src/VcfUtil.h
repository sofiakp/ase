#ifndef VCF_UTIL_H
#define VCF_UTIL_H

#include "vcf_file.h"
#include "swak/Swak.h"

namespace VcfUtil
{
  void GetNext(vcf_file &vfile, uint entry_index, vcf_entry &entry);
};

#endif
