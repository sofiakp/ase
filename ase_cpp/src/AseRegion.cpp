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

namespace AseRegion
{
    struct GenomicRegion
    {
        unsigned long start; //0-based
        unsigned long end; // 0-based
        map<string, int> fwdCounts;
        map<string, int> revCounts;
    };
    
    // User args
    string bam_file;
    string bed_file;
    int min_map_qual = 10;
    int min_base_qual = 20;
    //int min_ref_reads = 0;
    //int min_alt_reads = 0;
    int ascii_offset = 33;
    //int min_cover = 1;
    
    // Program vars
    BamReader bam_reader;
    map<string, vector<GenomicRegion> > chrom_genomicRegions;
    
    void ReadBed()
    {
        ifstream infile;
        infile.open(bed_file.c_str());
        cerr << "* Reading BED file " << endl;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            vector<string> info;
            stringstream convert(currentLine);
            string tmp;
            while (getline(convert, tmp, '\t'))
            {
                info.push_back(tmp);
            }
            GenomicRegion tmpGenomicRegion;
            tmpGenomicRegion.start = strtoul(info[1].c_str(), NULL, 0);
            tmpGenomicRegion.end = strtoul(info[2].c_str(), NULL, 0);
            if (chrom_genomicRegions.find(info[0].c_str()) == chrom_genomicRegions.end())
            {
                vector<GenomicRegion> tmpVectorGenomicRegion;
                tmpVectorGenomicRegion.push_back(tmpGenomicRegion);
                chrom_genomicRegions[info[0]] = tmpVectorGenomicRegion;
            }
            else
            {
                chrom_genomicRegions[info[0]].push_back(tmpGenomicRegion);
            }
        }
        cerr << "* Done reading BED file " << endl;
    }
    
    void Init(const vector<string> &all_args)
    {
        VecS args;
        
        string prog_desc = "Count reads that overlapping each regions from a BED file for each Readgroup in a BAM file. Output into stdout BED-like file with extra columns indicating counts for forward strand and reverse strand in order grouped by different readgroups. Both BED file and bam file should be sorted by position.";
        Swak::OptionParser op("<bam> <bed>", prog_desc);
        
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
        bed_file = args[1];
        
        ReadBed();
    }
};

using namespace AseRegion;

int main_aseregion(const vector<string> &all_args)
{
    Init(all_args);
    
    cerr << "* Reading bam file " << endl;
    OpenBam(bam_reader, bam_file);
    bam_reader.OpenIndex(bam_file + ".bai");
    
    vector<string> readGroupVector; //Obtain all the readgroups.
    SamHeader header = bam_reader.GetHeader();
    SamReadGroupDictionary headerRG = header.ReadGroups;
    for (SamReadGroupIterator it = headerRG.Begin(); it != headerRG.End(); it ++)
    {
        readGroupVector.push_back(it -> ID);
    }
    
    cout << "#CHROM" << "\t" << "StartPos" << "\t" << "EndPos";
    for (vector<string>::iterator it = readGroupVector.begin(); it != readGroupVector.end(); it ++)
    {
        cout << "\t" << *it;
    }
    cout << endl;
    
    vector<RefData> chroms = bam_reader.GetReferenceData();
    
    StlFor(chrom_idx, chroms)
    {
        string &chrom = chroms[chrom_idx].RefName;
        cerr << "* On chrom " << chrom << endl;
        
        map<string, vector<GenomicRegion> >::iterator searchIt = chrom_genomicRegions.find(chrom);
        
        BamAlignment startPointer; // This pointer will point to the region immediately before the start of current regions under inspection.
        bam_reader.Jump(chrom_idx);
        if (!bam_reader.GetNextAlignment(startPointer))
            break;
        
        int count = 0;
        // For each region, walk through all the reads correspoinding to this region and count the reads.
        for (vector<GenomicRegion>::iterator it = searchIt -> second.begin(); it != searchIt -> second.end(); ++it)
        {
            bam_reader.Jump(chrom_idx, startPointer.Position); // Fix the reading pointer.
            if (!bam_reader.GetNextAlignment(startPointer))
                break;
            int flag = 0;
            while (true)
            {
                int startEnd = startPointer.GetEndPosition();
                if (startEnd < it -> start)
                {
                    if (!bam_reader.GetNextAlignment(startPointer))
                    {
                        flag = 1;
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
            
            if (flag == 1)
            {
                break;
            }
            // Now startPointer assumes its rightful position.
            BamAlignment nextPointer = startPointer; //This pointer traverse through all reads that align to the current genomic region in bed file and the iteration ends when this pointer pass through the end of the region.
            
            while (true)
            {
                int nextStart = nextPointer.Position;
                if (nextStart > it -> end)
                {
                    break; // This iteration is done.
                }
                
                if (nextPointer.MapQuality < min_map_qual)
                {
                    if (!bam_reader.GetNextAlignment(nextPointer))
                    {
                        break;
                    }
                    continue;
                }
                
                string currentRG;
                Assert(nextPointer.GetReadGroup(currentRG));
                
                map<string, int> &RG_counts = nextPointer.IsReverseStrand() ? it -> revCounts : it -> fwdCounts;
                map<string, int>::iterator searchItForRG = RG_counts.find(currentRG);
                if (searchItForRG == RG_counts.end())
                {
                    RG_counts[currentRG] = 1;
                }
                else
                {
                    ++ RG_counts[currentRG];
                }
                if (!bam_reader.GetNextAlignment(nextPointer))
                {
                    break;
                }
            }
            count ++;
            if (count % 1000 == 0)
                cerr << "Processed" << "\t" << count << endl;
        }
        
        // Output the counts
        for (vector<GenomicRegion>::iterator it = searchIt -> second.begin(); it != searchIt -> second.end(); ++it)
        {
            cout << chrom << "\t" << it -> start << "\t" << it -> end;
            for (vector<string>::iterator subIt = readGroupVector.begin(); subIt != readGroupVector.end(); ++subIt)
            {
                map<string, int>::iterator searchItForRG = it -> fwdCounts.find(*subIt);
                if (searchItForRG != it -> fwdCounts.end())
                {
                    cout << "\t" << searchItForRG -> second << ",";
                }
                else
                {
                    cout << "\t" << "0,";
                }
                searchItForRG = it -> revCounts.find(*subIt);
                if (searchItForRG != it -> revCounts.end())
                {
                    cout << searchItForRG -> second;
                }
                else
                {
                    cout << "0";
                }
            }
            cout << endl;
        }
    }
    
    cerr << "* Done" << endl;
    
    return 0;
}

