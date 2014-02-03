#ifndef SWAK_FASTQ_READER_H
#define SWAK_FASTQ_READER_H

#include "Swak.h"

namespace Swak
{
  class Fastq
  {
    public:
      string header, seq, plusline, qual;
  };

  class FastqReader
  {
    public:

    void Init();

    void Init(istream * in_stream_ptr_);

    FastqReader();

    bool GetNext(Fastq &fq);

    private:
    istream * in_stream_ptr;
  };
};

ostream& operator<< (ostream &os, const Swak::Fastq &fq);

#endif
