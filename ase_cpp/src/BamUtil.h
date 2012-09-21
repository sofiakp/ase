#ifndef BAM_UTIL_H
#define BAM_UTIL_H

#include "BamAux.h"
#include "BamReader.h"
#include "BamWriter.h"
#include "swak/Swak.h"

namespace BamUtil
{
  using namespace BamTools;

  vector<CigarOp> ParseCigarStr(const string &cigar);

  string CigarToStr(const vector<CigarOp> &cigar, const string &spacer=" ");

  inline bool IsAligned(char cigar_char)
  {
    return cigar_char == 'M' || cigar_char == 'D';
  }

  inline bool IsAligned(const CigarOp &cigar)
  {
    return IsAligned(cigar.Type);
  }

  inline bool IsFromRead(char cigar_char)
  {
    return cigar_char != 'D';
  }

  inline bool IsFromRead(const CigarOp &cigar)
  {
    return IsFromRead(cigar.Type);
  }

  inline bool IsInsertionLike(const CigarOp &op)
  {
    return op.Type == 'I' || op.Type == 'H' || op.Type == 'S';
  }


  struct Hit
  {
    string chrom;
    char strand;
    int pos; // 0-based
    vector<CigarOp> cigar;
    int edit_dist;

    Hit();

    void Init(const string &chrom_, char strand_, int pos_, const vector<CigarOp> &cigar_, int edit_dist_);

    Hit(const string &chrom_, char strand_, int pos_, const vector<CigarOp> &cigar_, int edit_dist_);

    void Init(const string &alt_hit_str);

    Hit(const string &alt_hit_str);
    
    string ToStr() const;
  };


  void GetAltHits(const BamAlignment &bam, vector<Hit> &alt_hits);

  int AlignSize(const BamAlignment &bam);

  int AlignSize(const Hit &hit);

  int AlignStart(const BamAlignment &bam);

  int AlignEnd(const BamAlignment &bam);

  int AlignEnd(const Hit &hit);

  inline void AlignedString(const BamAlignment &bam, string &out_str)
  {
    out_str.clear();

    const string &bases = bam.QueryBases;

    int offset = 0; // Offset into bases

    // Assume bam has a valid cigar string
    for (int i = 0; i < bam.CigarData.size(); ++i)
    {
      char type = bam.CigarData[i].Type;
      int len = bam.CigarData[i].Length;

      // Only a match contributes letters from the read
      if (type == 'M')
        out_str += bases.substr(offset, len);
      else if (type == 'D' || type == 'N')
        out_str += string(len, '-');

      // Anything but a deletion consumes letter from the read
      if (type != 'D' && type != 'N')
        offset += len;

      // Ignore any inserted/clipped letters entirely
      // NOTE: does NOT support 'padding'
    }
  }

  inline string AlignedString(const BamAlignment &bam)
  {
    string s;
    AlignedString(bam, s);
    return s;
  }

  // Uses mapping quality for deletion qualities
  inline void AlignedQualities(const BamAlignment &bam, string &out_str, int ascii_offset=33)
  {
    out_str.clear();

    char del_qual = char(ascii_offset) + bam.MapQuality;

    int offset = 0; // Offset into qualities

    // Assume bam has a valid cigar string
    for (int i = 0; i < bam.CigarData.size(); ++i)
    {
      char type = bam.CigarData[i].Type;
      int len = bam.CigarData[i].Length;

      // Only a match contributes letters from the read
      if (type == 'M')
        out_str += bam.Qualities.substr(offset, len);
      else if (type == 'D' || type == 'N')
        out_str += string(len, del_qual);

      // Anything but a deletion consumes letter from the read
      if (type != 'D' && type != 'N')
        offset += len;

      // Ignore any inserted/clipped letters entirely
      // NOTE: does NOT support 'padding'
    }
  }

  inline string AlignedQualities(const BamAlignment &bam, int ascii_offset=33)
  {
    string s;
    AlignedQualities(bam, s, ascii_offset);
    return s;
  }

  string BamToCoreStr(const BamAlignment &bam, const vector<RefData> &ref_vec);

  void OpenBam(BamReader &bam_reader, const string &bam_file);

  void OpenBam(BamWriter &bam_writer, const string &bam_file, const SamHeader &header, const vector<RefData> &chroms);

  void MakeAlignmentRepetitive(BamAlignment &bam);

  void MakeAlignmentUnmapped(BamAlignment &bam);

  void FlipAlignmentStrand(BamAlignment &bam);

  string GetReadGroup(const BamAlignment &bam);

  void ClearBwaTags(BamAlignment &bam);

  string BamChrom(const BamAlignment &bam, const BamReader &bam_reader);

  // Consumes letters in a cigar string that are aligned to the genome and any inserted letters in the process
  inline void ConsumeAlignedLetters(vector<CigarOp> &ops, int num_letters, vector<CigarOp> &out_ops)
  {
    out_ops.clear();
    
    // To make this faster, we will reverse the ops and work from the end forward
    reverse(ops.begin(), ops.end());

    while (num_letters > 0)
    {
      AssertMsg(ops.size() > 0, "Ran out of cigar operations while consuming aligned letters");
      CigarOp &op = ops[ops.size() - 1];

      if (IsInsertionLike(op))
      {
        out_ops.push_back(op);
        ops.pop_back();
        // No need to update num_letters since we didnt remove any aligned letters
      }
      else if (op.Length <= num_letters)
      {
        // Consume entire op
        out_ops.push_back(op);
        ops.pop_back();
        num_letters -= op.Length;
      }
      else // op.Length > num_letters (we are done)
      {
        // Consume partial op
        CigarOp op_part(op.Type, num_letters);
        out_ops.push_back(op_part);

        ops[ops.size() - 1].Length -= num_letters; 

        num_letters = 0;
      }
    }

    AssertMsg(num_letters == 0, "Last cigar op wasnt big enough when consuming aligned letters");

    // Undo previous reverse
    reverse(ops.begin(), ops.end());
  }

  // Checks to see if the cigar string has the right number of letters used
  inline bool CigarIsValid(const BamAlignment &bam, const vector<CigarOp> &cigar)
  {
    int read_bases = bam.QueryBases.size();

    for (int i = 0; i < cigar.size(); ++i)
    {
      char type = cigar[i].Type;

      if (type == 'M' || type == 'I' || type == 'S' || type == 'H')
        read_bases -= cigar[i].Length;
    }

    return read_bases == 0;
  }

  inline bool CigarIsValid(const BamAlignment &bam)
  {
    return CigarIsValid(bam, bam.CigarData);
  }

  class BamLess
  {
  public:
    bool operator() (const BamAlignment &lhs, const BamAlignment &rhs) const
    {
      int l_end = AlignEnd(lhs);
      int r_end = AlignEnd(rhs);

      if (l_end == r_end)
        return lhs.Name < rhs.Name;
      else
        return l_end < r_end;
    }
  };

  class BamStreamer
  {
  public:

    BamStreamer(BamReader &bam_reader_) : bam_reader(bam_reader_)
    {
      JumpTo(0);
      SetMinQual(0);
    }

    void JumpTo(int chrom)
    {
      cur_chrom = chrom;
      bam_reader.Jump(cur_chrom);
    }

    void SetMinQual(int new_min_qual)
    {
      min_qual = new_min_qual;
    }

    bool GetNext(BamAlignment &bam)
    {
      while ((bam_reader.GetNextAlignment(bam)) && bam.RefID == cur_chrom) 
      {
        if (bam.MapQuality >= min_qual)
          return true;
      }

      return false;
    }

  private:

    BamReader &bam_reader;
    int cur_chrom;
    int min_qual;
  };

  
};

#endif
