#include "swak/FoldedOutStream.h"

FoldedOutStream::FoldedOutStream(int64 line_width_)
{
  Init(line_width_);
}

FoldedOutStream::FoldedOutStream(ostream &os, int64 line_width_)
{
  Init(os, line_width_);
}

void FoldedOutStream::Init(int64 line_width_)
{
  if (line_width_ <= 0)
    throw runtime_error("Illegal line width <= 0 in FoldedOutStream");
  else
    line_width = line_width_;
  line_chars_left = line_width;
}

void FoldedOutStream::Init(ostream &os, int64 line_width_)
{
  Init(line_width_);
  out_stream = &os;
}

void FoldedOutStream::Write(const string &val)
{
  if (out_stream == NULL)
    throw runtime_error("FoldedOutStream's stream is NULL");
  if (!out_stream->good())
    throw runtime_error("Can't write to FoldedOutStream's stream");

  int64 pos = 0;

  while (pos < val.size())
  {
    if (line_chars_left == 0)
    {
      line_chars_left = line_width;
      (*out_stream) << '\n';

      // If we needed a newline because of line width but the next char in the file is a newline, we can skip the one in the file
      if (val[pos] == '\n')
      {
        ++pos;
        continue;
      }
    }

    int64 chars_to_write = min(int64(val.size() - pos), line_chars_left);

    buf = val.substr(pos, chars_to_write);

    // If we are printing a newline, update line_chars_left and only print up to and including the newline
    int64 newline_pos = buf.find('\n');
    if (newline_pos != string::npos)
    {
      chars_to_write = newline_pos + 1;
      buf.resize(chars_to_write);
      line_chars_left = line_width; // We know we just printed a newline so reset line_chars_left
    }
    else
      line_chars_left -= chars_to_write;

    (*out_stream) << buf;

    pos += chars_to_write;
  }
}
