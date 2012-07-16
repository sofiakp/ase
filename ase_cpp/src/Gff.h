#ifndef GFF_H
#define GFF_H

#include "swak/Swak.h"

class GffEntry
{
public:
  string chrom;
  string source;
  string feature;
  int64 start; // 0-based
  int64 stop;  // 0-based, inclusive
  float score;
  char strand;
  int frame;
  string group; // Eg gene name, whatever

  /*
  int64 start() { return start; }
  int64 stop() { return stop; }
  */
  int64 end() const { return stop + 1; } 
  int64 size() const { return stop - start + 1; }

  GffEntry();

  //friend ostream& operator <<(ostream &os, const GffEntry &g);

  // For GTF format
  void GetAttributes(MapSS &attrib_map) const;
  MapSS GetAttributes() const;
  string GetAttribute(const string &key) const;
  void SetAttributes(const MapSS &attr);

  template <typename T>
  void SetAttribute(const string &key, const T &value)
  {
    MapSS attr = GetAttributes();
    attr[key] = ToStr(value);
    SetAttributes(attr);
  }

  bool RemoveAttribute(const string &key);

  bool operator<(const GffEntry &other) const;
  bool Overlaps(int64 other_start, int64 other_end) const;
  bool Overlaps(const GffEntry &other) const;
};

ostream& operator<< (ostream &os, const GffEntry &g);


#endif
