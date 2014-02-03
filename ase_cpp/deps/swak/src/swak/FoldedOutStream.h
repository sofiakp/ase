#ifndef FOLDED_OUT_STREAM_H
#define FOLDED_OUT_STREAM_H

#include "Swak.h"

class FoldedOutStream
{
public:

  FoldedOutStream(int64 line_width_=50);

  FoldedOutStream(ostream &os, int64 line_width_=50);

  void Init(int64 line_width_=50);

  void Init(ostream &os, int64 line_width_=50);

  void Write(const string &val);

private:
  int64 line_width; // Max number of chars per line (not including newline)
  int64 line_chars_left;  // How many chars we have left to print on this line
  ostream * out_stream; // Stream to write to
  string buf;

};
#endif
