#ifndef GENE_H
#define GENE_H

#include "swak/Swak.h"
#include "swak/Interval.h"
#include "Gff.h"

typedef Interval Exon;
typedef Interval CDS;

// Forward declarations
class Gene; 

// BEGIN:Splice
class Splice
{
public:
  // Relationship to the fake splices chrom -- must be set when that chrom is created!
  string chrom;     // Name of the fake "splices" chrom
  int64 start;      // What are the beginning and end of the junction within the "splices" chrom (0-based)
  int64 end;        // (0-based)
  char strand;
  int64 junction_pos; // The actual splice occurs between splice_pos and splice_pos + 1 on the "splices" chrom.
                      // Should be between coord.start and coord.end (0-based)

  // Attributes of the transcript/gene
  string gene_chrom;
  string gene_id;
  vector<string> tx_ids;
  Exon partial_exon1;
  Exon partial_exon2;

  void Init();

  Splice();

  Splice(const string &gene_id_, const VecS &tx_ids_, const Exon &exon1_, const Exon &exon2_);

  Splice(const GffEntry &gff);

  void Init(const GffEntry &gff);

  // Doesn't write newline
  GffEntry ToGff() const;

  // Unique string describing this junction
  string JunctionStr() const;
};

vector<Splice> GetSplicesFromFile(const string &filename);



class Transcript
{
public:
  string id;
  vector<Exon> exons;
  vector<CDS> cds_regions;

  string chrom;
  string parent_gene_id;
  char strand;

  int64 start() const; // 0-based
  int64 end() const;   // 0-based

  // Assumes exons are sorted
  // Order by leftmost exon's left then by the lastmost exon's right
  bool operator<(const Transcript &other) const;

  // NOTE: Doesn't set chrom, start, end of the junction
  // extension: how far on either side of the splice junction you want
  vector<Splice> GetSplices(int extension) const;

  string GetSequence(const string &ref);

  vector<GffEntry> GetExonGffs() const;
};

ostream& operator<< (ostream &os, const Transcript &t);


class Gene
{
public:
  string chrom; 
  string id;  // gene_id
  string name; // optional: gene_name
  char strand;
  vector<Transcript> transcripts;

  int64 start() const; // 0-based
  int64 end() const;   // 0-based

  // Assumes transcripts are in order
  bool operator<(const Gene &other) const;

  // NOTE: Doesn't set chrom, start, end of the junction
  vector<Splice> GetSplices(int extension) const;
};

ostream& operator<< (ostream &os, const Gene &g);

class TxExon // An exon on the transcriptome
{
public:
  string chrom;
  int64 start;
  int64 end;
  string gene_id;
  string tx_id;
  int exon_index;
  char strand;

  TxExon();

  TxExon(const GffEntry &gff);

  void Init();

  void Init(const GffEntry &gff);

  // Assumes on same chrom;
  bool operator<(const TxExon &other) const;
};


namespace GeneUtils
{
  void AddTxEntries(const string &chrom_seq, const Transcript &tx, const GffEntry &gff_template, string &cur_tx_seq, vector<GffEntry> &gffs);
};

#endif
