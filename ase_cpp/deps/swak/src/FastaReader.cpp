#include "swak/FastaReader.h"

FastaReader::FastaReader()
{
}

void FastaReader::Open(istream &is)
{
  if (!is.good())
    throw runtime_error("Bad fasta stream");

  in_stream = &is;
  next_state = First;
}

FastaReader::LineState FastaReader::ReadNextLine()
{
  if (next_state == Eof)
    return Eof;

  do 
  {
    char buffer[MAX_FASTA_LINE_LEN];

    if (!in_stream->getline(buffer, MAX_FASTA_LINE_LEN))
    {
      if (in_stream->gcount() == 0)
        return (next_state = Eof);
      throw runtime_error("Fasta line too long");
    }
    next_line = buffer;

    int end = next_line.length();
    while ((end > 0) && isspace(next_line[end - 1]))
      --end;
    next_line.resize(end);

  } while (next_line.empty());

  if (next_line[0] == '>')
    return (next_state = Header);

  return (next_state = Other);
}

bool FastaReader::GetNext()
{ 
  if (!HasNext())
    return false;
  ReadNext();
  return true;
}

bool FastaReader::HasNext()
{
  while (next_state == First)
    ReadNextLine();

  if (next_state == Eof)
    return false;

  if (next_state != Header)
    throw runtime_error("Invalid fasta file");
  return true;
}

void FastaReader::ReadNext()
{
  if (!HasNext())
    throw runtime_error("Can't ReadNext when !HasNext");

  header = next_line.substr(1, next_line.length() - 1);
  seq.clear();

  if (ReadNextLine() == Eof)
    throw runtime_error("Truncated fasta file");
  do
  {
    seq.append(next_line);
  } while (ReadNextLine() == Other);
}
