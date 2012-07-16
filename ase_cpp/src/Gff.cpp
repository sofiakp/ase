#include "Gff.h"
#include "swak/Interval.h"

GffEntry::GffEntry()
{
  start = -1;
  stop = -1;
  score = -1;
  strand = '?';
  frame = -1;
}

void GffEntry::GetAttributes(map<string, string> &attrib_map) const
{
  attrib_map.clear();
  VecS fields = SplitString2(StripString(group), ';');

  for (int i = 0; i < fields.size(); ++i)
  {
    string &attrib = fields[i];
    if (attrib.empty())
      continue;
    if (attrib[0] == ' ') // Trim leading space
      attrib = attrib.substr(1);

    size_t attrib_name_end = attrib.find(' ');
    if (attrib_name_end == string::npos)
      throw runtime_error(string("Gff error: No space in attribute"));

    string name = attrib.substr(0, attrib_name_end);

    if (attrib_name_end >= attrib.size() - 1)
      throw runtime_error("Gff error: Attribute had no value!");

    string value;
    StripString(attrib.substr(attrib_name_end + 1), value, '"');

    attrib_map[name] = value;
  }
}

MapSS GffEntry::GetAttributes() const
{
  MapSS attribs;
  GetAttributes(attribs);
  return attribs;
}

void GffEntry::SetAttributes(const MapSS &attr)
{
  stringstream ss;

  int i = 0;
  StlForMapConst(string, string, iter, attr)
  {
    if (i != 0)
      ss << ' ';

    const string &key = iter->first;
    const string &value = iter->second;
    ss << key;

    bool needs_quotes = StringContains(key, ' ');

    if (needs_quotes)
      ss << " \"" << value << "\";";
    else
      ss << ' ' << value << ';';

    ++i;
  }

  group = ss.str();
}

// Doesn't write newline
ostream& operator<< (ostream &os, const GffEntry &g)
{
  os << g.chrom << '\t'
     << (g.source.empty() ? "." : g.source) << '\t'
     << g.feature << '\t'
     << (g.start + 1) << '\t' // Stored as 0-based so we +1 since text file should be 1-based
     << (g.stop + 1) << '\t'; // Stored as 0-based so we +1 since text file should be 1-based
  if (g.score < 0)
    os << ".\t";
  else
    os << g.score << '\t';
  os << g.strand << '\t';
  if (g.frame < 0)
    os << '.' << '\t';
  else
    os << g.frame << '\t';
  os << g.group;
  return os;
}

bool GffEntry::operator<(const GffEntry &other) const
{
  if (chrom == other.chrom)
  {
    Interval my_int(start, stop + 1);
    Interval other_int(other.start, other.stop + 1);

    return my_int < other_int;
  }
  return chrom < other.chrom;
}

bool GffEntry::Overlaps(int64 other_start, int64 other_end) const
{
  return start < other_end && other_start < stop + 1;
}

bool GffEntry::Overlaps(const GffEntry &other) const
{
  return chrom == other.chrom && Overlaps(other.start, other.stop + 1);
}

string GffEntry::GetAttribute(const string &key) const
{
  MapSS attribs = GetAttributes();
  return attribs[key];
}

bool GffEntry::RemoveAttribute(const string &key)
{
  MapSS attr = GetAttributes();
  if (attr.find(key) != attr.end())
  {
    attr.erase(key);
    SetAttributes(attr);
    return true;
  }
  return false;
}
