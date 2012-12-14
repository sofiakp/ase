#include "swak/Swak.h"
#include "swak/System.h"
#include "swak/Helpers.h"
#include "RodUtil.h"
#include "BamUtil.h"
#include "VcfUtil.h"
#include "swak/FastaReader.h"
#include "swak/FieldReader.h"
#include "swak/BinaryIO.h"
#include <vector>

using namespace BamTools;
using namespace BamUtil;
using namespace RodUtil;
using namespace VcfUtil;

namespace AseQuantMultiRG
{
    struct Counts
    {
        int num_ref;
        int num_alt;
        int num_other;
        
        Counts()
        {
            num_ref = num_alt = num_other = 0;
        }
        
        int NumTotal()
        {
            return num_other + num_alt + num_ref;
        }
    };
    
    struct Snp
    {
        int pos; // 0-based
        char ref;
        char alt;
        map<string, Counts> fwd, rev;
        
        Snp()
        {
            pos = 0;
            alt = ref = '?';
        }
    };
    
    // User args
    string bam_file;
    string snp_file;
    int min_map_qual = 10;
    int min_base_qual = 20;
    int min_ref_reads = 0;
    int min_alt_reads = 0;
    int ascii_offset = 33;
    int min_cover = 1;
    
    // Program vars
    BamReader bam_reader;
    map<string, vector<Snp> > snps_by_chrom;
    
    void ReadVcf()
    {
        cerr << "* Opening vcf file " << endl;
        vcf_file vfile(snp_file);
        cerr << "* Done Opening vcf file " << endl;
        vcf_entry entry(vfile.N_indv);
        
        cerr << "* Reading vcf file " << endl;
        for (uint i = 0; i < vfile.N_entries; ++i)
        {
            VcfUtil::GetNext(vfile, i, entry);
            
            // AssertMsg(entry.get_REF().size() == 1, "Variant at " + ToStr(entry.get_POS()) + " has a non-snp ref allele");
            // AssertMsg(entry.get_ALT().size() == 1, "Variant at " + ToStr(entry.get_POS()) + " has a non-snp alt allele");
            
            Snp snp;
            
            snp.pos = entry.get_POS() - 1; // Make 0-based
            snp.ref = entry.get_REF()[0];
            snp.alt = entry.get_ALT()[0];
            
            string chrom = entry.get_CHROM();
            
            // AssertMsg(snps_by_chrom[chrom].size() == 0 || snp.pos > snps_by_chrom[chrom].rbegin()->pos, "Expected all SNPs to be in order on each chrom (and no duplicates) " + ToStr(snp.pos) );
            
            if (snps_by_chrom[chrom].size() == 0)
                snps_by_chrom[chrom].reserve(1000000);
            
            snps_by_chrom[chrom].push_back(snp);
        }
        cerr << "* Done reading vcf file " << endl;
    }
    
    void Init(const vector<string> &all_args)
    {
        VecS args;
        
        string prog_desc = "Counts alleles at SNPs for each Readgroup. Outputs VCF-like file with counts of SNPs. The counts are seperated by ',' and the order for each column is REF_FORWART_STRAND_COUNTS, ALT_FORWART_STRAND_COUNTS, OTHER_FORWART_STRAND_COUNTS, REF_REVERSE_STRAND_COUNTS, ALT_REVERSE_STRAND_COUNTS,OTHER_REVERSE_STRAND_COUNTS.\nOutput is directed to stdout.";
        Swak::OptionParser op("<bam> <known_snps.vcf>", prog_desc);
        
        op.AddOpt(min_map_qual, 'm', "min-map-qual", "INT", "Min mapping qual to use [" + ToStr(min_map_qual) + "]");
        op.AddOpt(min_base_qual, 'b', "min-base-qual", "INT", "Min base qual to vote or be counted [" + ToStr(min_base_qual) + "]");
        //op.AddOpt(min_ref_reads, 'r', "min-ref-reads", "INT", "Min reads required for ref allele [" + ToStr(min_ref_reads) + "]");
        //op.AddOpt(min_alt_reads, 'a', "min-alt-reads", "INT", "Min reads required for alt allele [" + ToStr(min_alt_reads) + "]");
        op.AddOpt(ascii_offset, 'a', "ascii-offset", "INT", "Ascii offset of the BAM's read qualities (note: bam transforms qualities to ascii33) [" + ToStr(ascii_offset) + "]");
        //op.AddOpt(min_cover, 'c', "min-cover", "INT", "Minimum number of reads overlapping a SNP to print the entry [" + ToStr(min_cover) + "]");
        
        if (!op.Parse(all_args, args, 2) || args.size() < 2)
        {
            op.PrintUsage();
            exit(1);
        }
        
        bam_file = args[0];
        snp_file = args[1];
        
        ReadVcf();
        
        StlForMap(string, vector<Snp>, iter, snps_by_chrom)
        cerr << "* Found " << iter->second.size() << " snps on " << iter->first << endl;
    }
};

using namespace AseQuantMultiRG;

int main_asequantmultirg(const vector<string> &all_args)
{
    Init(all_args);
    
    cerr << "* Reading bam file " << endl;
    OpenBam(bam_reader, bam_file);
    bam_reader.OpenIndex(bam_file + ".bai");
    
    vector<string> readGroupVector;
    SamHeader header = bam_reader.GetHeader();
    SamReadGroupDictionary headerRG = header.ReadGroups;
    for (SamReadGroupIterator it = headerRG.Begin(); it != headerRG.End(); it ++)
    {
        readGroupVector.push_back(it -> ID);
    }
    
    
    vector<RefData> chroms = bam_reader.GetReferenceData();
    
    cout << "#CHROM" << "\t" << "POS" << "\t" << "REF" << "\t" << "ALT";
    for (vector<string>::iterator it = readGroupVector.begin(); it != readGroupVector.end(); it ++)
    {
        cout << "\t" << *it;
    }
    cout << endl;
    
    StlFor(chrom_idx, chroms)
    {
        string &chrom = chroms[chrom_idx].RefName;
        vector<Snp> snps = snps_by_chrom[chrom];
        
        int s = 0; // Index into snp array
        
        BamAlignment bam;
        bam_reader.Jump(chrom_idx);
        
        string align;
        string qualities;
        
        cerr << "* On chrom " << chrom << endl;

        while (bam_reader.GetNextAlignment(bam) && bam.RefID == chrom_idx) 
        {
	  if (bam.MapQuality < min_map_qual || !bam.IsMapped())
                continue;
       
            string currentRG;
            Assert(bam.GetReadGroup(currentRG));
            
            int start = AlignStart(bam);
            int end = AlignEnd(bam);
            
            // Move the current SNP pointer so that it is ahead of the read's start (since bam alignments are in sorted order)
            while (s < snps.size() && snps[s].pos < start)
                ++s;
            
            // Stop everything if we have visited all SNPs on this chrom
            if (s >= snps.size())
                break;
            
            // Find any/all SNPs that are within the bam alignment
            int n = 0; // Number of SNPs overlapped
            while ((s + n) < snps.size() && snps[s + n].pos < end) // Then it overlaps!
                ++n;
            
            // Now, look at each SNP and see which way it votes
            AlignedString(bam, align);
            AlignedQualities(bam, qualities);
            Assert(align.size() == qualities.size());

            // Now, tally votes
            for (int i = 0; i < n; ++i)
            {
                Snp &snp = snps[s + i];
                char base = align[snp.pos - start]; // Base from the read
                int qual = int(qualities[snp.pos - start]) - ascii_offset; // Base from the read
                
                //AssertMsg(qual >= 0 && qual <= 100, ToStr(qual) + "\n" + bam.Name + "\n" + CigarToStr(bam.CigarData) + "\n" + bam.QueryBases + "\n" + bam.Qualities);
                
                if (base == '-' || qual < min_base_qual)
                    continue;
                
                map<string, Counts> &RG_counts = bam.IsReverseStrand() ? snp.rev : snp.fwd;
                
                map<string, Counts>::iterator searchIt = RG_counts.find(currentRG);
                
                if (searchIt == RG_counts.end())
                {
                    if (base == snp.ref)
                    {
                        RG_counts[currentRG].num_ref = 1;
                        RG_counts[currentRG].num_alt = 0;
                        RG_counts[currentRG].num_other = 0;
                    }
                    else if (base == snp.alt)
                    {
                        RG_counts[currentRG].num_ref = 0;
                        RG_counts[currentRG].num_alt = 1;
                        RG_counts[currentRG].num_other = 0;
                    }
                    else
                    {
                        RG_counts[currentRG].num_ref = 0;
                        RG_counts[currentRG].num_alt = 0;
                        RG_counts[currentRG].num_other = 1;
                    }
                }
                else
                {
                    if (base == snp.ref)
                    {
                        searchIt -> second.num_ref += 1;
                    }
                    else if (base == snp.alt)
                    {
                        searchIt -> second.num_alt += 1;
                    }
                    else
                    {
                        searchIt -> second.num_other += 1;
                    }
                }
            }
        }
        
        // Output counts
        for (int s = 0; s < snps.size(); ++s)
        {
            cout << chrom << "\t" << snps[s].pos + 1 << "\t" << snps[s].ref << "\t" << snps[s].alt;
            for (vector<string>::iterator it = readGroupVector.begin(); it != readGroupVector.end(); it ++)
            {
                map<string, Counts>::iterator searchIt = snps[s].fwd.find(*it);
                if (searchIt != snps[s].fwd.end())
                {
                    cout << "\t" << searchIt -> second.num_ref << "," << searchIt -> second.num_alt << "," << searchIt -> second.num_other << ",";
                }
                else
                {
                    cout << "\t" << "0,0,0,";
                }
                searchIt = snps[s].rev.find(*it);
                if (searchIt != snps[s].rev.end())
                {
                    cout << searchIt -> second.num_ref << "," << searchIt -> second.num_alt << "," << searchIt -> second.num_other;
                }
                else
                {
                    cout << "0,0,0";
                }
            }
            cout << endl;
        }
    }
    
    cerr << "* Done" << endl;
    
    return 0;
}

