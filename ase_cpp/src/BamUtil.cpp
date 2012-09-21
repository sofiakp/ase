#include "BamUtil.h"

namespace BamUtil
{
  using namespace BamTools;

  vector<CigarOp> ParseCigarStr(const string &cigar)
  {
    vector<CigarOp> vec;
    CigarOp op;
    string len;

    StlFor(i, cigar)
    {
      char c = cigar[i];
      if (c >= '0' && c <= '9')
      {
        len += c;
      }
      else
      {
        op.Type = cigar[i];
        op.Length = StrTo<int>(len);

        vec.push_back(op);
        len.clear();
      }
    }

    return vec;
  }

  string CigarToStr(const vector<CigarOp> &cigar, const string &spacer)
  {
    string s;

    StlFor(i, cigar)
    {
      s += ToStr(cigar[i].Length);
      s += cigar[i].Type;
      if (i != cigar.size() - 1)
        s += spacer;
    }

    return s;
  }

  Hit::Hit() {}

  void Hit::Init(const string &chrom_, char strand_, int pos_, const vector<CigarOp> &cigar_, int edit_dist_)
  {
    chrom = chrom_;
    strand = strand_;
    pos = pos_;
    cigar = cigar_;
    edit_dist = edit_dist_;
  }

  Hit::Hit(const string &chrom_, char strand_, int pos_, const vector<CigarOp> &cigar_, int edit_dist_)
  {
    Init(chrom_, strand_, pos_, cigar_, edit_dist_);
  }

  void Hit::Init(const string &alt_hit_str)
  {
    //chrTx,+7371,100M,0
    VecS fields;
    SplitString(alt_hit_str, ',', fields);
    Assert(fields[0].size() >= 1);
    Assert(fields[1].size() >= 2);
    Assert(fields[2].size() >= 2);
    Assert(fields[3].size() >= 1);

    int p = StrTo<int>(fields[1].substr(1)) - 1;
    Init(fields[0], fields[1][0], p, ParseCigarStr(fields[2]), StrTo<int>(fields[3]));
  }

  Hit::Hit(const string &alt_hit_str)
  {
    Init(alt_hit_str);
  }

  string Hit::ToStr() const
  {
    stringstream ss;
    ss << chrom << ',' << strand << pos << ',' << CigarToStr(cigar) << ',' << edit_dist;
    return ss.str();
  }


  void GetAltHits(const BamAlignment &bam, vector<Hit> &alt_hits)
  {
    alt_hits.clear();

    string val;
    bam.GetTag("XA:Z", val);

    VecS alt_hit_strings;
    SplitString(val.substr(0, val.size() - 1), ';', alt_hit_strings);
    StlFor(i, alt_hit_strings)
      alt_hits.push_back(Hit(alt_hit_strings[i]));
  }

  int AlignSize(const BamAlignment &bam)
  {
    int aligned = 0;
    StlFor(i, bam.CigarData)
    {
      if (bam.CigarData[i].Type == 'M' || bam.CigarData[i].Type == 'D' || bam.CigarData[i].Type == 'N')
        aligned += bam.CigarData[i].Length;
    }

    return aligned;
  }

  int AlignSize(const Hit &hit)
  {
    int aligned = 0;
    StlFor(i, hit.cigar)
    {
      if (hit.cigar[i].Type == 'M' || hit.cigar[i].Type == 'D' || hit.cigar[i].Type == 'N')
        aligned += hit.cigar[i].Length;
    }

    return aligned;
  }

  int AlignStart(const BamAlignment &bam)
  {
    return bam.Position;
  }

  int64 AlignStart(const BamAlignment * bam)
  {
    return bam->Position;
  }

  int AlignEnd(const BamAlignment &bam)
  {
    return bam.Position + AlignSize(bam);
  }

  int64 AlignEnd(const BamAlignment * bam)
  {
    return AlignEnd(*bam);
  }

  int AlignEnd(const Hit &hit)
  {
    return hit.pos + AlignSize(hit);
  }

  string BamToCoreStr(const BamAlignment &bam, const vector<RefData> &ref_vec)
  {
    stringstream ss;
    ss << ref_vec[bam.RefID].RefName << ": " << (bam.Position + 1);
    return ss.str();
  }

  void OpenBam(BamReader &bam_reader, const string &bam_file)
  {
    AssertMsg(bam_reader.Open(bam_file), "Couldn't open bam file for reading: " + bam_file);
  }

  void OpenBam(BamWriter &bam_writer, const string &bam_file, const SamHeader &header, const vector<RefData> &chroms)
  {
    AssertMsg(bam_writer.Open(bam_file, header, chroms), "Couldn't open bam file for writing: " + bam_file);
  }

  void MakeAlignmentRepetitive(BamAlignment &bam)
  {
    //cout << "Making alignment repetitive" << endl;

    // Store edit dist if it had it
    int edit_dist = 0;
    bool had_edit_dist = bam.GetTag("NM", edit_dist);
    //bam.TagData.clear();
    ClearBwaTags(bam);
    if (had_edit_dist)
      Assert(bam.EditTag("NM", "i", edit_dist)); // Store edit dist

    // Mark it as repeptive
    // TODO: UNDO when "A" tags are fixed
    //Assert(bam.EditTag("XT", "A", 'R'));
    bam.MapQuality = 0;
    bam.SetIsMapped(true);
  }

  void MakeAlignmentUnmapped(BamAlignment &bam)
  {
    //cout << "Making alignment unmapped" << endl;
    bam.MapQuality = 0;
    bam.Position = -1;
    bam.RefID = -1;
    bam.CigarData.clear();
    bam.SetIsMapped(false);
    //bam.TagData.clear();
    ClearBwaTags(bam);
  }

  void FlipAlignmentStrand(BamAlignment &bam)
  {
    // (1) flip the flag
    bam.SetIsReverseStrand(!bam.IsReverseStrand());

    // (2) revcomp the seq
    reverse(bam.QueryBases.begin(), bam.QueryBases.end());
    StlFor(i, bam.QueryBases)
    {
      char old = bam.QueryBases[i];
      if (old == 'A')
        bam.QueryBases[i] = 'T';
      else if (old == 'C')
        bam.QueryBases[i] = 'G';
      else if (old == 'G')
        bam.QueryBases[i] = 'C';
      else if (old == 'T')
        bam.QueryBases[i] = 'A';
      else if (old == 'a')
        bam.QueryBases[i] = 't';
      else if (old == 'c')
        bam.QueryBases[i] = 'g';
      else if (old == 'g')
        bam.QueryBases[i] = 'c';
      else if (old == 't')
        bam.QueryBases[i] = 'a';
    }

    // (3) reverse the qual
    reverse(bam.Qualities.begin(), bam.Qualities.end());

    // (4) reverse the cigar
    reverse(bam.CigarData.begin(), bam.CigarData.end());
  }

  string GetReadGroup(const BamAlignment &bam)
  {
    string rg;
    bam.GetTag("RG", rg);
    return rg;
  }

  void ClearBwaTags(BamAlignment &bam)
  {
/*
      │NM  │ Edit distance                                  │
      │MD  │ Mismatching positions/bases                    │
      │AS  │ Alignment score                                │
      │BC  │ Barcode sequence                               │
      ├────┼────────────────────────────────────────────────┤
      │X0  │ Number of best hits                            │
      │X1  │ Number of suboptimal hits found by BWA         │
      │XN  │ Number of ambiguous bases in the referenece    │
      │XM  │ Number of mismatches in the alignment          │
      │XO  │ Number of gap opens                            │
      │XG  │ Number of gap extentions                       │
      │XT  │ Type: Unique/Repeat/N/Mate-sw                  │
      │XA  │ Alternative hits; format: (chr,pos,CIGAR,NM;)* │
      ├────┼────────────────────────────────────────────────┤
      │XS  │ Suboptimal alignment score                     │
      │XF  │ Support from forward/reverse alignment         │
      │XE  │ Number of supporting seeds    
      */

    int num_tags = 4 + 8 + 3;
    const char * tags[] = {"NM", "MD", "AS", "BC", "X0", "X1", "XN", "XM", "XO", "XG", "XT", "XA", "XS", "XF", "XE"};
    for (int i = 0; i < num_tags; ++i)
    {
      if (bam.HasTag(tags[i]))
        AssertMsg(bam.RemoveTag(tags[i]), bam.Name);
    }
  }

  string BamChrom(const BamAlignment &bam, const BamReader &bam_reader)
  {
    const vector<RefData> &refs = bam_reader.GetReferenceData();

    if (bam.IsMapped() && bam.RefID >= 0 && bam.RefID < refs.size())
      return refs[bam.RefID].RefName;

    string chrom;
    return chrom;
  }
};
