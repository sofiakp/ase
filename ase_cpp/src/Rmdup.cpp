#include "swak/Swak.h"
#include "swak/System.h"
#include "swak/Helpers.h"
#include "BamReader.h"
#include "BamWriter.h"
#include "BamUtil.h"

using namespace BamTools;
using namespace BamUtil;

namespace Reconcile
{
    // User args
    string bam_file;
    
    // Program vars
    vector<string> args;
    YAML::Node options;
    
    BamReader bam_reader;
  BamWriter bam_writer;
    void ProcessCmdLine(const vector<string> &all_args)
    {
        int num_args = 1;
        int num_prog_names = 2;
        
        if (!ProcessInput(all_args, num_prog_names, num_args, options, args))
        {
            string usage = Basename(all_args[0]) + " rmdup [options] <in.bam>";
            string desc = "Assumes in.bam is sorted by coordinate";
            
            VecS opt_lines;
            
            PrintUsage(usage, desc, opt_lines); 
            exit(0);
        }
        
        bam_file = args[0];
    }
    
    void Init()
    {
        OpenBam(bam_reader, bam_file);
             
        SamHeader header = bam_reader.GetHeader();
             
        OpenBam(bam_writer, '/dev/stdout', header, bam_reader.GetReferenceData());
        
        // TODO: Make deterministic
        srand48(0xCAFEBABA);
    }
    
};


using namespace Reconcile;

int main_reconcile(const vector<string> &all_args)
{
    ProcessCmdLine(all_args);
    Init();
    
    BamAlignment bam;
    BamAlignment bam_prev;
    
    int64 num_reads = 0;
    int64 num_bam1 = 0;
    int64 num_bam2 = 0;
    int64 num_amb = 0;
    
    while (bam_reader.GetNextAlignment(bam))
    {
        num_reads ++;
        
        int edit_dist1_f = INT_MAX;
        int edit_dist1_l = INT_MAX;
        int edit_dist2_f = INT_MAX;
        int edit_dist2_l = INT_MAX;
        
        BamAlignment * bam_ptr_f = NULL;
        BamAlignment * bam_ptr_l = NULL;
        
        if (!bam.IsMapped())
        {
	  //output alignment; continue
	}
	if(){
	  // the previous is not empty and this pos is the same as the previous
	  //choose best
	}else if(){
	  // the previous is not empty
	  //output the previous
	  //previous = current
	    if (lrand48() % 2 == 0)
	    {
	      bam_ptr_f = &bam1_f;
	      bam_ptr_l = &bam1_l;
	    }
	    else
	    {
	      bam_ptr_f = &bam2_f;
	      bam_ptr_l = &bam2_l;
	    }
	    bam_ptr_f->RemoveTag("XC");
            bam_ptr_f->AddTag("XC", "Z", "both_unmap");
            bam_ptr_l->RemoveTag("XC");
            bam_ptr_l->AddTag("XC", "Z", "both_unmap");
            bam_ptr_f->EditTag("RG", "Z", "ambiguous");
            bam_ptr_l->EditTag("RG", "Z", "ambiguous");
            num_amb ++;
        }
        else if ((bam1_f.IsMapped() && bam1_l.IsMapped()) && (!bam2_f.IsMapped() || !bam2_l.IsMapped()))
        {
            bam_ptr_f = &bam1_f;
            bam_ptr_l = &bam1_l;
            bam_ptr_f->RemoveTag("XC");
            bam_ptr_f->AddTag("XC", "Z", "one_unmap");
            bam_ptr_l->RemoveTag("XC");
            bam_ptr_l->AddTag("XC", "Z", "one_unmap");
            bam_ptr_f->EditTag("RG", "Z", readgroup_1);
            bam_ptr_l->EditTag("RG", "Z", readgroup_1);
            num_bam1 ++;
        }
        else if ((!bam1_f.IsMapped() || !bam1_l.IsMapped()) && (bam2_f.IsMapped() && bam2_l.IsMapped()))
        {
            bam_ptr_f = &bam2_f;
            bam_ptr_l = &bam2_l;
            
            bam_ptr_f->RemoveTag("XC");
            bam_ptr_f->AddTag("XC", "Z", "one_unmap");
            bam_ptr_l->RemoveTag("XC");
            bam_ptr_l->AddTag("XC", "Z", "one_unmap");
            bam_ptr_f->EditTag("RG", "Z", readgroup_2);
            bam_ptr_l->EditTag("RG", "Z", readgroup_2);
            num_bam2 ++;
        }
        else
        {
            Assert(bam1_f.GetTag("NM", edit_dist1_f));
            Assert(bam1_l.GetTag("NM", edit_dist1_l));
            Assert(bam2_f.GetTag("NM", edit_dist2_f));
            Assert(bam2_l.GetTag("NM", edit_dist2_l));
            if ((edit_dist1_f < edit_dist2_f && edit_dist1_l <= edit_dist2_l) || (edit_dist1_f <= edit_dist2_f && edit_dist1_l < edit_dist2_l))
            {
                bam_ptr_f = &bam1_f;
                bam_ptr_l = &bam1_l;
                bam_ptr_f->EditTag("RG", "Z", readgroup_1);
                bam_ptr_l->EditTag("RG", "Z", readgroup_1);
                num_bam1 ++;
            }
            else if ((edit_dist1_f > edit_dist2_f && edit_dist1_l >= edit_dist2_l) ||(edit_dist1_f >= edit_dist2_f && edit_dist1_l > edit_dist2_l))
            {
                bam_ptr_f = &bam2_f;
                bam_ptr_l = &bam2_l;
                bam_ptr_f->EditTag("RG", "Z", readgroup_2);
                bam_ptr_l->EditTag("RG", "Z", readgroup_2);
                num_bam2 ++;
            }
            else
            {
                if (lrand48() % 2 == 0)
                {
                    bam_ptr_f = &bam1_f;
                    bam_ptr_l = &bam1_l;
                }
                else
                {
                    bam_ptr_f = &bam2_f;
                    bam_ptr_l = &bam2_l;
                }
                bam_ptr_f->EditTag("RG", "Z", "ambiguous");
                bam_ptr_l->EditTag("RG", "Z", "ambiguous");
                num_amb ++;
            }
            if (bam1_f.Position == bam2_f.Position && bam1_l.Position == bam2_l.Position)
            {
                bam_ptr_f->RemoveTag("XC");
                bam_ptr_f->AddTag("XC", "Z", "same_pos");
                bam_ptr_l->RemoveTag("XC");
                bam_ptr_l->AddTag("XC", "Z", "same_pos");
            }
            else
            {
                bam_ptr_f->RemoveTag("XC");
                bam_ptr_f->AddTag("XC", "Z", "diff_pos");
                bam_ptr_l->RemoveTag("XC");
                bam_ptr_l->AddTag("XC", "Z", "diff_pos");
                string RGtag;
                bam_ptr_f->GetTag("RG", RGtag);
                if (RGtag == "ambiguous")
                {
                    bam_ptr_f->MapQuality = 0;
                    bam_ptr_l->MapQuality = 0;
                }
            }
        }
        bam_writer.SaveAlignment(*bam_ptr_f);
        bam_writer.SaveAlignment(*bam_ptr_l);
        
        if (num_reads % 1000 == 0)
        {
            cerr << "* Alignments processed so far: " << num_reads << "\r";
            cerr.flush();
        }
    }
    cerr << endl;
    Assert(!bam_reader1.GetNextAlignment(bam1_f) && !bam_reader2.GetNextAlignment(bam2_f));
    cout << "* Total alignments processed: \t" << num_reads << endl;
    
    cout << "* Num in readgroup " << ToStrL(readgroup_1 + ":", 11) << "\t" << ToStrL(num_bam1, 10) << "\t" << ToPercent(num_bam1, num_reads) << endl;
    cout << "* Num in readgroup " << ToStrL(readgroup_2 + ":", 11) << "\t" << ToStrL(num_bam2, 10) << "\t" << ToPercent(num_bam2, num_reads) << endl;
    cout << "* Num in readgroup amb:       \t" << ToStrL(num_amb, 10) << "\t" << ToPercent(num_amb, num_reads) << endl;
    
    return 0;
}

