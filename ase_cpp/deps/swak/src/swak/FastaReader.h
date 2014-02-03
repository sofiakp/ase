#ifndef FASTA_READER_H
#define FASTA_READER_H

#include "Swak.h"

const int MAX_FASTA_LINE_LEN = 500000;

class FastaReader
{
public:
  typedef enum {First, Eof, Header, Other} LineState;

  FastaReader();

  void Open(istream &is);

  LineState ReadNextLine();

  bool GetNext();

  string header;
  string seq;

private:
  istream * in_stream;

  string next_line;
  LineState next_state;

  bool HasNext();

  void ReadNext();
};

#endif
