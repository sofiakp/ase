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

namespace computeOverlapping
{
    //User arguments
    vector<string> inputFileNames, inputDataSetNames;
    string population, ldFileDirectory, gwasLocation;
    
    void init(const vector<string> &all_args)
    {
        po::options_description desc("Compute the number of overlapping SNPs without filtering. You can specify multiple datasets, but the subjects should be of the same population. Output is directed to stdout.");
        desc.add_options()
        ("help,h", "display help message")
        ("queryInputfiles,iv", po::value< vector<string> >(&inputFileNames), "specify all the files separated by space that score the SNPs which need to be overlapped")
        ("datasets,ds", po::value< vector<string> >(&inputDataSetNames), "specify all the name of the datasets in the same order as the queryInputfiles")
        ("population,pl", po::value<string>(&population), "specify the population.");
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
        
        ldFileDirectory = "/home/yulingl/ase_diseases/hapmap_ld/"+ population + "_1";
        gwasLocation = "/home/yulingl/ase_diseases/gwasCatalog/processedGWASCatalog";
    }
    
    //Construct LD information for a particular population
    struct castLD
    {
        map<string, set<string> > tag_SNPinLD;
        void construct();
        void constructOverChr(string &chrNum);
    };
    
    void castLD::construct()
    {
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
            constructOverChr(chrNum);
        }
    }
    
    void castLD::constructOverChr(string &chrNum)
    {
        ifstream infile;
        string fileName = ldFileDirectory + "/" + chrNum;
        infile.open(fileName.c_str());
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            string ID1 = chrNum + "_" + info[0];
            string ID2 = chrNum + "_" + info[1];
            tag_SNPinLD[ID1].insert(ID2);
            tag_SNPinLD[ID1].insert(ID1);
            tag_SNPinLD[ID2].insert(ID1);
            tag_SNPinLD[ID2].insert(ID2);
        }
    }
    
    //Construct a set containing query SNPs
    struct castSNPs
    {
        set<string> SNPs;
        void construct(string &fileName);
    };
    
    void castSNPs::construct(string &fileName)
    {
        ifstream infile;
        infile.open(fileName.c_str());
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            if (currentLine[0] != '#')
            {
                string SNPID;
                split(info, currentLine, is_any_of("\t"));
                trim(info[0]);
                trim(info[2]);
                SNPID = info[0].substr(3) + "_" + info[2];
                SNPs.insert(SNPID);
            }
        }
    }
    
    //Construct a set containing GWAS SNPs
    struct castGWASCatalog
    {
        set<string> SNPs;
        map<string, set<string> > leadSNP_SNPsinLD;
        void construct();
        void constructOverChr(string &chrNum);
        void mapSNPinLD(map<string, set<string> > &tag_SNPinLD);
    };
    
    void castGWASCatalog::construct()
    {
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
            constructOverChr(chrNum);
        }
    }
    
    void castGWASCatalog::constructOverChr(string &chrNum)
    {
        ifstream infile;
        string fileName = gwasLocation + "/" + chrNum;
        infile.open(fileName.c_str());
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            SNPs.insert(info[0]);
        }
    }
    
    void castGWASCatalog::mapSNPinLD(map<string, set<string> > &tag_SNPinLD)
    {
        leadSNP_SNPsinLD.clear();
        for (set<string>::iterator it = SNPs.begin(); it != SNPs.end(); it ++)
        {
            map<string, set<string> >::iterator subIt = tag_SNPinLD.find(*it);
            if (subIt == tag_SNPinLD.end())
            {
                leadSNP_SNPsinLD[*it].insert(*it);
            }
            else
            {
                leadSNP_SNPsinLD[*it] = subIt -> second;
            }
        }
    }
    
    //Compute the overlaps
    int intersectionSize(set<string> &set1, set<string> &set2)
    {
        set<string> intersection;
        set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), inserter(intersection, intersection.end()));
        return intersection.size();
    }
    
    int computeOverlappingInLD(castSNPs &SNPs, castGWASCatalog &gwasCatalog)
    {
        int count = 0;
        for (map<string, set<string> >::iterator it1 = gwasCatalog.leadSNP_SNPsinLD.begin(); it1 != gwasCatalog.leadSNP_SNPsinLD.end(); it1 ++)
        {
            for (set<string>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++ it2)
            {
                if (SNPs.SNPs.find(*it2) != SNPs.SNPs.end())
                {
                    ++ count;
                    break;
                }
            }
        }
        return count;
    }
}

using namespace computeOverlapping;

int main_computeOverlapping(const vector<string> &all_args)
{
    init(all_args);
    
    castLD LD;
    LD.construct();
    
    castGWASCatalog gwasCatalog;
    gwasCatalog.construct();
    gwasCatalog.mapSNPinLD(LD.tag_SNPinLD);
    
    cout << "Data Set Names" << "\t" << "# of Direct Overlapping GWAS SNPs" << "\t" << "# of GWAS SNPs in LD" << endl;
    for (int i = 0; i < inputFileNames.size(); ++ i)
    {
        castSNPs SNPs;
        SNPs.construct(inputFileNames[i]);
        
        cout << inputDataSetNames[i] << "\t" << intersectionSize(gwasCatalog.SNPs, SNPs.SNPs) << "\t" << computeOverlappingInLD(SNPs, gwasCatalog) << endl;
    }
    
    return 0;
}
