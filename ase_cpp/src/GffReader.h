#ifndef GFF_READER_H
#define GFF_READER_H

#include "swak/System.h"
#include "Gff.h"
#include "swak/FieldReader.h"
#include "swak/StringUtil.h"

class GffReader
{
public:
  GffReader()
  {
  }

  GffReader(istream &is)
  {
    Init(is);
  }

  GffReader(const string &filename)
  {
    reader.Init(filename);
  }

  void Init(istream &is)
  {
    reader.Init(is);
  }

  void Init(const string &filename)
  {
    reader.Init(filename);
  }

  bool GetNext(GffEntry &gff_entry)
  {
    bool got_next = reader.GetNext(fields);

    if (!got_next)
      return false;

    if (fields.size() != 9)
      throw runtime_error(string("Unexepected number of fields in a GFF @ line ") + ToStr(reader.LastLineNum()) + ": " + reader.LastLine()); 

    gff_entry.chrom = fields[0];
    gff_entry.source = fields[1];
    gff_entry.feature = fields[2];
    gff_entry.start = StrTo<int64>(fields[3]) - 1; // We convert to 0-based from 1-based
    gff_entry.stop = StrTo<int64>(fields[4]) - 1;  // We convert to 0-based from 1-based
    gff_entry.score = StrTo<int64>(fields[5]);

    if (fields[6].size() != 1)
      throw runtime_error(string("Illegal GFF strand: ") + fields[6]);
    gff_entry.strand = fields[6][0];

    if (fields[7] == ".")
      gff_entry.frame = -1;
    else
    {
      gff_entry.frame = StrTo<int>(fields[7]);
      if (gff_entry.frame < 0 || gff_entry.frame > 2)
        throw runtime_error(string("Illegal GFF frame: ") + fields[7]);
    }

    gff_entry.group = fields[8];

    return true;
  }

  vector<GffEntry> GetAll()
  {
    vector<GffEntry> vec;
    GffEntry gff;
    while (GetNext(gff))
      vec.push_back(gff);

    return vec;
  }

  vector<GffEntry> &GetAll(vector<GffEntry> &vec)
  {
    vec.clear();

    GffEntry gff;
    while (GetNext(gff))
      vec.push_back(gff);

    return vec;
  }

private:
  FieldReader reader;
  vector<string> fields;
};

#endif
