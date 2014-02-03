#include "swak/FastqReader.h"

namespace Swak
{
  void FastqReader::Init()
  {
    in_stream_ptr = NULL;
  }

  void FastqReader::Init(istream * in_stream_ptr_)
  {
    Init();
    in_stream_ptr = in_stream_ptr_;
  }

  FastqReader::FastqReader()
  {
    Init();
  }

  bool FastqReader::GetNext(Fastq &r)
  {
    if (!getline(*in_stream_ptr, r.header))
      return false;
    if (!getline(*in_stream_ptr, r.seq))
      return false;
    if (!getline(*in_stream_ptr, r.plusline))
      return false;
    if (!getline(*in_stream_ptr, r.qual))
      return false;
    return true;
  }
};

ostream& operator<< (ostream &os, const Swak::Fastq &fq)
{
  os << fq.header << '\n' << fq.seq << '\n' << fq.plusline
     << '\n' << fq.qual << endl;
  return os;
}
