#ifndef SWAK_FIELD_READER_H
#define SWAK_FIELD_READER_H

#include "Swak.h"

class FieldReader
{
public:
  FieldReader();
  FieldReader(const string &pathname);
  FieldReader(istream &in);
  ~FieldReader();

  void Init(const string &pathname);
  void Init(istream &in);
  void Close();
  void Reset();
  bool IsOpen() const;
  bool GetNext(vector<string> &fields);
  void SetMergeAdjacent(bool merge_adjacent_delims)
  {
    merge_adj_delims = merge_adjacent_delims;
  }

  void SetDelim(const string &new_delim)
  {
    delim = new_delim;
  }

  string LastLine()
  {
    return line;
  }

  int64 LastLineNum()
  {
    return line_num;
  }

private:
  string delim;
  string line;
  int64 line_num;
  istream *in_stream;
  bool merge_adj_delims;
  enum { Closed, File, Stream } ftype;
};

#endif
