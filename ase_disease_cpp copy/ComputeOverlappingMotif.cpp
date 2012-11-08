#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include "boost/algorithm/string.hpp"
#include "boost/program_options.hpp"
namespace po = boost::program_options;

using namespace std;
using namespace boost;

#define MOTIFDIRECTORY  "/home/yulingl/ENCODE_ChIP-seq_motif/motif/"
#define CHIPSEQDIRECTORY "/home/yulingl/ENCODE_ChIP-seq_motif/Gm12878_ChIP-seq_preprocessed/"
#define LDDIRECTORY "/home/yulingl/ase_diseases/hapmap_ld/"

namespace computeOverlappingMotif
{
    //User arguments
    vector<string> inputFileNames, inputDataSetNames, outputDirectory;
    string population;
    bool ChIPseqOn, motifOn;
    
    void init(const vector<string> &all_args)
    {
        po::options_description desc("Compute the number of overlapping hits of AS-SNPs with motifs and ChIP-seq peaks. You can specify multiple datasets, but the subjects should be of the same population. Output is directed to stdout.");
        desc.add_options()
        ("help,h", "display help message")
        ("inputs,i", po::value< vector<string> >(&inputFileNames), "specify all the AS-SNP files separated by space")
        ("datasets,d", po::value< vector<string> >(&inputDataSetNames), "specify all the name of the datasets in the same order as the inputs")
        ("outputs,o", po::value<string>(&population), "specify the output directory")
        ("c", po::value<bool>(&ChIPseqOn), "if overlapping with ChIP-seq is needed, please specify 1. Otherwise, please specify 0")
        ("m", po::value<bool>(&motifOn), "if overlapping with motif is needed, please specify 1. Otherwise, please specify 0")
        ("population,p", po::value<string>(&population), "specify the population");
        po::variables_map vm;
        if (all_args.size() <= 2)
        {
            cout << desc << "\n";
            exit(0);
        }
        
        const char* av[all_args.size()];
        for (int i = 0; i < all_args.size(); ++ i)
        {
            av[i] = all_args[i].c_str();
        }
        po::store(po::parse_command_line(all_args.size(), av, desc), vm);
        po::notify(vm);
        if (vm.count("help"))
        {
            cout << desc << "\n";
            exit(0);
        }
    }
    
    struct castSNPInfo
    {
        unsigned long coord;
        string rsID;
        char alt;
        int alt_count;
        char ref;
        int ref_count;
    };
    
    //Construct a set containing query SNPs
    struct castSNPs
    {
        map<string, map<string, vector<castSNPInfo> > > chr_dataset_SNPs;
        void construct(string &fileName);
    };
    
    void castSNPs::construct(string fileName, string dataset)
    {
        ifstream infile;
        infile.open(fileName.c_str());
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            if (currentLine[0] != '#')
            {
                split(info, currentLine, is_any_of("\t"));
                trim(info[0]);
                trim(info[1]);
                trim(info[2]);
                castSNPInfo tmp;
                tmp.coord = (info[1].c_str(), NULL, 0);
                tmp.rsID = info[2];
                tmp.ref = info[3][0];
                tmp.alt = info[4][0];
                tmp.ref_count = atoi(info[5].c_str());
                tmp.alt_count = atoi(info[6].c_str());
                chr_SNPs[info[0]][dataset].push_back(tmp);
            }
        }
    }
    
    //Construct a set containing GWAS SNPs
    
    struct castMotifInfo
    {
        unsigned long begin;//1 based
        unsigned long end;//1 based
        string name;
    };
    
    struct castMotif
    {
        vector<castMotifInfo> motifVector;
        void constructOverChr(bool exper, string &chrNum);
    };
    
    void castMotif::constructOverChr(bool exper, string &chrNum)
    {
        motifVector.clear();
        ifstream infile;
        string fileName;
        if (exper)
        {
            fileName = CHIPSEQDIRECTORY + chrNum;
        }
        else
        {
            fileName = MOTIFDIRECTORY + chrNum;
        }
        infile.open(fileName.c_str());
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            castMotifInfo tmp;
            tmp.begin = (info[2].c_str(), NULL, 0);
            if (exper)
            {
                ++ tmp.begin;
            }
            tmp.end = (info[3].c_str(), NULL, 0);
            tmp.name = info[0];
            motifVector.push_back(tmp);
        }
    }
    
    void computeOverlap(vector<castMotifInfo> &motifVector, vector<castSNPInfo> &SNPs, map<string, int> &motif_counts, ofstream &outfile, string chrNum)
    {
        int pointer = 0;
        for (int i = 0; i < motifVector.size(); ++ i)
        {
            while (1)
            {
                if (SNPs[pointer].coord >= motifVector[i].begin)
                {
                    break;
                }
                ++ pointer;
            }
            int searchPointer = pointer;
            while (1)
            {
                if (SNPs[searchPointer] > motifVector[i].end)
                {
                    break;
                }
                map<string, int>::iterator pIt = motif_counts.find(motifVector[i].name);
                if (pIt == motifVector.end())
                {
                    motif_counts[motifVector[i].name] = 1;
                }
                else
                {
                    ++ pIt -> second;
                }
                outfile << chrNum << "\t" << motifVector[i].name << "\t" << motifVector[i].begin << "\t" << motifVector[i].end << "\t" << SNPs[searchPointer].coord << "\t" << SNPs[searchPointer].rsID << "\t" << SNPs[searchPointer].ref << "\t" << SNPs[searchPointer].ref_count << "\t" << SNPs[searchPointer].alt << "\t" << SNPs[searchPointer].alt_count << endl;
                ++ searchPointer;
            }
        }
    }
}

using namespace computeOverlappingMotif;

int main_computeOverlappingMotif(const vector<string> &all_args)
{
    init(all_args);
    
    castSNPs SNPs;
    for (int i = 0; i < inputFileNames.size(); ++ i)
    {
        SNPs.construct(inputFileNames[i], inputDataSetNames[i]);
    }
    
    if (ChIPseqOn)
    {
        vector<ofstream *> outfileVector;
        for (int i = 0; i < inputFileNames.size(); ++ i)
        {
            string fileName = outputDirectory + "/" + inputDataSetNames[i] + "_" + "ChIP-seqPeakOverlapInfo";
            ofstream tmp;
            tmp.open(fileName.c_str());
            outfileVector.push_back(&tmp);
        }
        
        map<string, map<string, int> > dataset_motif_counts;
        for (int i = 1; i <= 23; i++)
        {
            string chrNum;
            if (i == 23)
            {
                chrNum = "X";
            }
            else
            {
                stringstream convert;
                convert << i;
                chrNum = convert.str();
            }
            
            castMotif motifs;
            motifs.constructOverChr(1, chrNum);
            for (int j = 0; j < inputFileNames.size(); ++ j)
            {
                computeOverlap(motifs.motifVector, SNPs.chr_dataset_SNPs[chrNum][inputDataSetNames[j]], dataset_motif_counts[inputDataSetNames[j]], outfileVector[j], chrNum);
            }
        }
        
        for (int i = 0; i < inputFileNames.size(); ++ i)
        {
            string fileName = outputDirectory + "/" + inputDataSetNames[i] + "_" + "ChIP-seqPeakOverlapHits";
            ofstream outfile;
            outfile.open(fileName.c_str());
            for (map<string, int>::iterator it = dataset_motif_counts[inputDataSetNames[i]].begin(); it != dataset_motif_counts[inputDataSetNames[i]].end(); ++ it)
            {
                outfile << it -> first << "\t" << it -> second << endl;
            }
            outfile.close();
        }
        for (int i = 0; i < inputFileNames.size(); ++ i)
        {
            outfileVector[i] -> close();
        }
    }
    
    if (motifOn)
    {
        vector<ofstream *> outfileVector;
        for (int i = 0; i < inputFileNames.size(); ++ i)
        {
            string fileName = outputDirectory + "/" + inputDataSetNames[i] + "_" + "motifOverlapInfo";
            ofstream tmp;
            tmp.open(fileName.c_str());
            outfileVector.push_back(&tmp);
        }
        
        map<string, map<string, int> > dataset_motif_counts;
        for (int i = 1; i <= 23; i++)
        {
            string chrNum;
            if (i == 23)
            {
                chrNum = "X";
            }
            else
            {
                stringstream convert;
                convert << i;
                chrNum = convert.str();
            }
            
            castMotif motifs;
            motifs.constructOverChr(0, chrNum);
            for (int j = 0; j < inputFileNames.size(); ++ j)
            {
                computeOverlap(motifs.motifVector, SNPs.chr_dataset_SNPs[chrNum][inputDataSetNames[j]], dataset_motif_counts[inputDataSetNames[j]], outfileVector[j], chrNum);
            }
        }
        
        for (int i = 0; i < inputFileNames.size(); ++ i)
        {
            string fileName = outputDirectory + "/" + inputDataSetNames[i] + "_" + "motifOverlapHits";
            ofstream outfile;
            outfile.open(fileName.c_str());
            for (map<string, int>::iterator it = dataset_motif_counts[inputDataSetNames[i]].begin(); it != dataset_motif_counts[inputDataSetNames[i]].end(); ++ it)
            {
                outfile << it -> first << "\t" << it -> second << endl;
            }
            outfile.close();
        }
        
        for (int i = 0; i < inputFileNames.size(); ++ i)
        {
            outfileVector[i] -> close();
        }
    }
    
    return 0;
}
