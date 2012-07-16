#include "Gene.h"
#include "GffReader.h"
#include "swak/StlUtil.h"

// BEGIN:Splice

void Splice::Init()
{
  chrom = "chrSp";
  start = 1;
  end = 2;
  junction_pos = 1;
  strand = '+';
}

Splice::Splice()
{
  Init();
}

Splice::Splice(const string &gene_id_, const VecS &tx_ids_, const Exon &exon1_, const Exon &exon2_)
{
  Init();
  gene_id = gene_id_;
  tx_ids = tx_ids_;
  partial_exon1 = exon1_;
  partial_exon2 = exon2_;
}

Splice::Splice(const GffEntry &gff)
{
  Init(gff);
}

void Splice::Init(const GffEntry &gff)
{
  Init();
  MapSS attr;
  chrom = gff.chrom;
  start = gff.start;
  end = gff.stop + 1;
  strand = gff.strand;

  gff.GetAttributes(attr);
  gene_id = attr["gene_id"];
  tx_ids = SplitString2(attr["tx_ids"], ",");


  // These are all 1-based in the file, so - 1 as we read them in
  partial_exon1.SetStart(StrTo<int64>(attr["pexon1_start"]) - 1);
  partial_exon1.SetEnd(StrTo<int64>(attr["pexon1_stop"]) - 1 + 1);
  partial_exon2.SetStart(StrTo<int64>(attr["pexon2_start"]) - 1);
  partial_exon2.SetEnd(StrTo<int64>(attr["pexon2_stop"]) - 1 + 1);
  junction_pos = StrTo<int64>(attr["jxn_pos"]) - 1;
  gene_chrom = attr["gene_chrom"];

}

GffEntry Splice::ToGff() const
{
  GffEntry gff;
  gff.chrom = chrom;
  gff.start = start;
  gff.stop = end - 1;
  gff.feature = "splice";
  gff.strand = '+';

  MapSS attr;
  attr["gene_id"] = gene_id;
  attr["tx_ids"]  = JoinFields(tx_ids, ",");

  // These are all 1-based in the text file, so + 1 as we output them
  attr["pexon1_start"] = ToStr(partial_exon1.start() + 1);
  attr["pexon1_stop"]  = ToStr(partial_exon1.end() - 1 + 1);
  attr["pexon2_start"] = ToStr(partial_exon2.start() + 1);

  attr["pexon2_stop"]  = ToStr(partial_exon2.end() - 1 + 1);
  attr["jxn_pos"]      = ToStr(junction_pos + 1);
  attr["gene_chrom"]   = gene_chrom;

  gff.SetAttributes(attr);

  return gff;
}

vector<Splice> GetSplicesFromFile(const string &filename)
{
  GffReader reader(filename);

  vector<Splice> splices;
  GffEntry gff;
  Splice splice;

  while(reader.GetNext(gff))
  {
    if (gff.feature == "splice")
    {
      splice.Init(gff);
      splices.push_back(splice);
    }
  }

  return splices;
}

string Splice::JunctionStr() const
{
  stringstream ss;
  ss << gene_id << "@" << gene_chrom << ":" << partial_exon1.start() << "-" << partial_exon1.end() << "x" << partial_exon2.start() << "-" << partial_exon2.end();
  return ss.str();
}

// END:Splice
// BEGIN:Transcript

int64 Transcript::start() const
{
  return exons[0].start();
}

int64 Transcript::end() const
{
  return exons[exons.size() - 1].end();
}

// Assumes exons are sorted
// Order by leftmost exon's left then by the lastmost exon's right
bool Transcript::operator<(const Transcript &other) const
{
  if (start() == other.start())
    return end() < other.end();

  return start() < other.start();
}

// NOTE: Doesn't set chrom, start, end of the junction
// extension: how far on either side of the splice junction you want
vector<Splice> Transcript::GetSplices(int extension) const
{
  vector<Splice> junctions;

  junctions.resize(exons.size() - 1);

  for(int i = 1; i < exons.size(); ++i)
  {
    Splice &sj = junctions[i - 1];
    sj.gene_chrom = chrom;
    sj.gene_id = parent_gene_id;
    sj.tx_ids.push_back(id);

    // We must clamp the extension in case the exon is too short.
    int64 start_for_exon1 = max(exons[i-1].end() - extension, exons[i-1].start());
    int64 end_for_exon2 = min(exons[i].start() + extension, exons[i].end());

    sj.partial_exon1.SetStart(start_for_exon1);
    sj.partial_exon1.SetEnd(exons[i-1].end());

    sj.partial_exon2.SetStart(exons[i].start());
    sj.partial_exon2.SetEnd(end_for_exon2);
  }

  return junctions;
}

ostream& operator<< (ostream &os, const Transcript &t)
{
  vector<GffEntry> exon_gffs = t.GetExonGffs();

  StlFor(i, exon_gffs)
    os << exon_gffs[i] << endl;

  StlFor(i, t.cds_regions)
  {
    const CDS &e = t.cds_regions[i];

    GffEntry gff;
    gff.chrom = t.chrom;
    gff.feature = "CDS";
    gff.start = e.start();
    gff.stop = e.end() - 1;
    gff.strand = t.strand;

    MapSS attr;
    attr["gene_id"] = t.parent_gene_id;
    attr["transcript_id"] = t.id;
    gff.SetAttributes(attr);

    os << gff << endl;
  }
  return os;
}

vector<GffEntry> Transcript::GetExonGffs() const
{
  vector<GffEntry> gffs;

  StlFor(i, exons)
  {
    const Exon &e = exons[i];

    GffEntry gff;
    gff.chrom = chrom;
    gff.feature = "exon";
    gff.start = e.start();
    gff.stop = e.end() - 1;
    gff.strand = strand;

    MapSS attr;
    attr["gene_id"] = parent_gene_id;
    attr["transcript_id"] = id;
    gff.SetAttributes(attr);

    gffs.push_back(gff);
  }

  return gffs;
}

// END:Transcript
// BEGIN:Gene

int64 Gene::start() const
{
  return transcripts[0].start();
}

int64 Gene::end() const
{
  int64 rightmost = -1;
  StlFor(i, transcripts)
    rightmost = max(rightmost, transcripts[i].end());
  return rightmost;
}

// Assumes transcripts are in order
bool Gene::operator<(const Gene &other) const
{
  if (start() == other.start())
    return end() < other.end();
  return start() < other.start();
}


// NOTE: Doesn't set chrom, start, end of the junction
vector<Splice> Gene::GetSplices(int extension) const
{
  map<string, Splice> jxn_hash; // So we can avoid adding a junction twice

  // Add each splice, skipping duplicates (adding each redundant transcript to tx_ids)
  StlFor(i, transcripts)
  {
    vector<Splice> jxns = transcripts[i].GetSplices(extension);

    StlFor(j, jxns)
    {
      string jxn_str = jxns[j].JunctionStr();

      if (!Map::Contains(jxn_hash, jxn_str))
        jxn_hash[jxn_str] = jxns[j];
      else
      {
        Assert(jxns[j].tx_ids.size() == 1);
        jxn_hash[jxn_str].tx_ids.push_back(jxns[j].tx_ids[0]);
      }
    }
  }

  // Copy the hash values into a vector
  vector<Splice> junctions;
  StlForMap(string, Splice, iter, jxn_hash)
    junctions.push_back(iter->second);
  return junctions;
}

ostream& operator<< (ostream &os, const Gene &g)
{
  StlFor(i, g.transcripts)
    os << g.transcripts[i];
  return os;
}

// END:GENE
// BEGIN:TxExon
TxExon::TxExon()
{
  Init();
}

TxExon::TxExon(const GffEntry &gff)
{
  Init(gff);
}

void TxExon::Init()
{
  chrom = "chr?";
  start = -1;
  end = -1;
  exon_index = 0;
  strand = '+';
}

void TxExon::Init(const GffEntry &gff)
{
  chrom = gff.chrom;
  start = gff.start;
  end = gff.end();
  gene_id = gff.GetAttribute("gene_id");
  tx_id = gff.GetAttribute("transcript_id");
  exon_index = StrTo<int>(gff.GetAttribute("exon_index"));
  strand = gff.strand;
}

bool TxExon::operator<(const TxExon &other) const
{
  return end <= other.start;
}


// End:TxExon

// BEGIN:GENEUTILS
namespace GeneUtils
{
  // Only uses exons, not CDS etc
  void AddTxEntries(const string &chrom_seq, const Transcript &tx, const GffEntry &gff_template, string &cur_tx_seq, vector<GffEntry> &gffs)
  {
    int64 pos = cur_tx_seq.size();

    vector<Exon> exons = tx.exons;
    sort(exons.begin(), exons.end());

    StlFor(e, exons)
    {
      // Assumed to have chrom, feature, strand
      GffEntry gff = gff_template;
      gff.SetAttribute("gene_id", tx.parent_gene_id);
      gff.SetAttribute("transcript_id", tx.id);
      gff.SetAttribute("exon_index", e);
      gff.start = pos;
      gff.stop = pos + exons[e].size() - 1;

      cur_tx_seq += chrom_seq.substr(exons[e].start(), exons[e].size());
      gffs.push_back(gff);
      pos += exons[e].size();
    }
  }
};

// END:GENEUTILS
