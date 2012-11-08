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

#define EQTLINPUTDIRECTORY "/home/yulingl/ase_diseases/eQTL/Lymphoblastoid/mapped/"
#define GWASINPUTDIRECTORY "/home/yulingl/ase_diseases/gwasCatalog/processedGWASCatalog/"

namespace computeOverlapping
{
    //User arguments
    vector<string> inputFileNames, inputDataSetNames;
    string population;
    bool eQTLOn = 0, GWASOn = 0;
    
    void init(const vector<string> &all_args)
    {
        po::options_description desc("Compute the number of overlapping SNPs without filtering. You can specify multiple datasets, but the subjects should be of the same population. Output is directed to stdout.");
        desc.add_options()
        ("help,h", "display help message")
        ("queryInputfiles,i", po::value< vector<string> >(&inputFileNames), "specify all the files separated by space that store the SNPs which need to be overlapped")
        ("datasets,d", po::value< vector<string> >(&inputDataSetNames), "specify all the name of the datasets in the same order as the queryInputfiles")
        ("population,p", po::value<string>(&population), "specify the population.")
        ("eQTLOn,e", po::value<bool>(&eQTLOn), "if eQTL analysis is desired, specify 1.")
        ("GWASOn,g", po::value<bool>(&GWASOn), "if GWAS analysis is desired, specify 1.");
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
        string fileName = "/home/yulingl/ase_diseases/hapmap_ld/"+ population + "_1/" + chrNum;
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
    struct castAssociationSNPs
    {
        set<string> SNPs;
        map<string, set<string> > key_SNPsinLD;
        void construct(bool GWASoreQTL);
        void constructOverChr(string &chrNum, bool GWASoreQTL);
        void mapSNPinLD(map<string, set<string> > &tag_SNPinLD);
    };
    
    void castAssociationSNPs::construct(bool GWASoreQTL)
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
            constructOverChr(chrNum, GWASoreQTL);
        }
    }
    
    void castAssociationSNPs::constructOverChr(string &chrNum, bool GWASoreQTL)
    {
        ifstream infile;
        string fileName;
        if (GWASoreQTL == 1)
        {
            fileName = GWASINPUTDIRECTORY + chrNum;
        }
        else
        {
            fileName = EQTLINPUTDIRECTORY + chrNum;
        }
        infile.open(fileName.c_str());
        if (infile.fail())
        {
            cerr << "Failed to open " << fileName << endl;
            exit(1);
        }
        
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            SNPs.insert(info[0]);
        }
    }
    
    void castAssociationSNPs::mapSNPinLD(map<string, set<string> > &tag_SNPinLD)
    {
        key_SNPsinLD.clear();
        for (set<string>::iterator it = SNPs.begin(); it != SNPs.end(); it ++)
        {
            map<string, set<string> >::iterator subIt = tag_SNPinLD.find(*it);
            if (subIt == tag_SNPinLD.end())
            {
                key_SNPsinLD[*it].insert(*it);
            }
            else
            {
                key_SNPsinLD[*it] = subIt -> second;
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
    
    int computeOverlappingInLD(castSNPs &SNPs, castAssociationSNPs &associationSNPs)
    {
        int count = 0;
        for (map<string, set<string> >::iterator it1 = associationSNPs.key_SNPsinLD.begin(); it1 != associationSNPs.key_SNPsinLD.end(); it1 ++)
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
    
    if (GWASOn)
    {
        castAssociationSNPs GWASSNPs;
        GWASSNPs.construct(1);
        GWASSNPs.mapSNPinLD(LD.tag_SNPinLD);
        cout << "Data Set Names" << "\t" << "# of Direct Overlapping GWAS SNPs" << "\t" << "# of GWAS SNPs in LD" << endl;
        for (int i = 0; i < inputFileNames.size(); ++ i)
        {
            castSNPs SNPs;
            SNPs.construct(inputFileNames[i]);
            
            cout << inputDataSetNames[i] << "\t" << intersectionSize(GWASSNPs.SNPs, SNPs.SNPs) << "\t" << computeOverlappingInLD(SNPs, GWASSNPs) << endl;
        }
    }
    
    if (eQTLOn)
    {
        castAssociationSNPs eQTLs;
        eQTLs.construct(0);
        eQTLs.mapSNPinLD(LD.tag_SNPinLD);
        cout << "Data Set Names" << "\t" << "# of Direct Overlapping eQTLs" << "\t" << "# of eQTLs in LD" << endl;
        for (int i = 0; i < inputFileNames.size(); ++ i)
        {
            castSNPs SNPs;
            SNPs.construct(inputFileNames[i]);
            
            cout << inputDataSetNames[i] << "\t" << intersectionSize(eQTLs.SNPs, eQTLs.SNPs) << "\t" << computeOverlappingInLD(SNPs, eQTLs) << endl;
        }
    }
    
    return 0;
}

#undef EQTLINPUTDIRECTORY
#undef GWASINPUTDIRECTORY
