#include "swak/FieldReader.h"

FieldReader::FieldReader()
{
  ftype = Closed;
  delim = "\t";
  line_num = 0;
  merge_adj_delims = false;
}

FieldReader::~FieldReader()
{
  Close();
}

FieldReader::FieldReader(const string &pathname)
{
  delim = "\t";
  ftype = Closed;
  Init(pathname);
}

FieldReader::FieldReader(istream &in)
{
  delim = "\t";
  ftype = Closed;
  Init(in);
}

void FieldReader::Init(const string &pathname)
{
  Close();

  in_stream = (istream *)new ifstream(pathname.c_str(), ifstream::in);
  if (!in_stream || !((ifstream *)in_stream)->is_open())
    throw runtime_error("Unable to open file \"" + pathname + "\"");
  ftype = File;
  line_num = 0;
}

void FieldReader::Init(istream &in)
{
  in_stream = &in;
  ftype = Stream;
  line_num = 0;
}

void FieldReader::Close()
{
  switch (ftype)
  {
  case Closed:
  case Stream:
    break;

  case File:
    ((ifstream *)in_stream)->close();
    delete in_stream;
    break;
  }
  ftype = Closed;
}

void FieldReader::Reset()
{
  if (ftype != Closed)
  {
    in_stream->clear();
    switch (ftype)
    {
    case File:
      ((ifstream *)in_stream)->seekg(0, ifstream::beg);
      break;

    case Stream:
      throw runtime_error("Can't reset stream");
    case Closed: // To suppress warning
      break;
    }
  }
}

bool FieldReader::IsOpen() const
{
  return ftype != Closed;
}

bool FieldReader::GetNext(vector<string> &fields)
{
  while(getline(*in_stream, line))
  {
    ++line_num;

    if (line.empty() || line[0] == '#')
      continue;

    if (merge_adj_delims)
      SplitString2(line, delim, fields);
    else
      SplitString(line, delim, fields);

    return true;
  }

  return false;
}
