#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;

namespace processedGWASCatalog
{
    struct castSNPInfo
    {
        string allelFreq;
        string func;
        string pos;
    };
    
    struct castSNP
    {
        map<string, castSNPInfo> ucscToInfo;
        void construct(string &chrNum);
    };
    
    void castSNP::construct(string &chrNum)
    {
        ifstream infile;
        string prefix = "/home/yulingl/ase_diseases/dbSNP/dbSNP135_hg19_chr";
        string fileName;
        fileName = prefix + chrNum;
        infile.open(fileName.c_str());
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            castSNPInfo tmpSNPInfo;
            tmpSNPInfo.allelFreq = info[3];
            tmpSNPInfo.func = info[4];
            tmpSNPInfo.pos = info[2];
            string key = chrNum + "_" + info[1];
            ucscToInfo[key] = tmpSNPInfo;
        }
        infile.close();
    }
    
    struct GWASSNPInfo
    {
        double pvalue;
        vector<string> platforms;
        string function;
    };
    
    struct castGWAS
    {
        map<string, set<string> > chr_keys;
        map<string, GWASSNPInfo> key_info;
        void construct();
        void addInfoOutput(castSNP &SNPs, string &chrNum);
    };
    
    void castGWAS::construct()
    {
        ifstream infile;
        string fileName = "/home/yulingl/ase_diseases/gwasCatalog/gwascatalog.txt";
        infile.open(fileName.c_str());
        vector<string> L0Info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            split(L0Info, currentLine, is_any_of("\t"));
            if (L0Info[21] != "" && L0Info[11] != "")
            {
                if (L0Info[21].substr(0,2) == "rs")
                {
                    string chrNum;
                    if (L0Info[11] == "23")
                    {
                        chrNum = "X";
                    }
                    else
                    {
                        chrNum = L0Info[11];
                    }
                    string key = chrNum + "_" + L0Info[21];
                    GWASSNPInfo tmp;
                    tmp.pvalue = atof(L0Info[28].c_str());
                    tmp.function = L0Info[24];
                    
                    replace_first(L0Info[32], "{", "[");
                    vector<string> L1Info;
                    split(L1Info, L0Info[32], is_any_of("["));
                    vector<string> L2Info;
                    
                    replace_first(L1Info[0], "and", "&");
                    erase_all(L1Info[0], " ");
                    split(L2Info, L1Info[0], is_any_of("&"));
                    tmp.platforms = L2Info;
                    key_info[key] = tmp;
                    
                    if (chr_keys.find(chrNum) == chr_keys.end())
                    {
                        set<string> tmpSet;
                        tmpSet.insert(key);
                        chr_keys[chrNum] = tmpSet;
                    }
                    else
                    {
                        chr_keys[chrNum].insert(key);
                    }
                }
            }
        }
        infile.close();
    }
    
    void castGWAS::addInfoOutput(castSNP &SNPs, string &chrNum)
    {
        ofstream outfile;
        string fileName;
        string command = "mkdir -p /home/yulingl/ase_diseases/gwasCatalog/processedGWASCatalog";
        if (system(command.c_str()) != 0)
        {
            cerr << "Failed to execute command " << command << endl;
            exit(1);
        };
        
        fileName = "/home/yulingl/ase_diseases/gwasCatalog/processedGWASCatalog/" + chrNum;
        outfile.open(fileName.c_str());
        for (set<string>::iterator it = chr_keys[chrNum].begin(); it != chr_keys[chrNum].end(); ++ it)
        {
            map<string, castSNPInfo>::iterator searchIt = SNPs.ucscToInfo.find(*it);
            if (searchIt != SNPs.ucscToInfo.end())
            {
                map<string, GWASSNPInfo>::iterator pIt = key_info.find(*it);
                outfile << *it << "\t" << searchIt -> second.allelFreq << "\t" << searchIt -> second.pos << "\t" << pIt -> second.pvalue << "\t" << pIt -> second.function << "\t";
                for (int i = 0; i < pIt -> second.platforms.size() - 1; ++ i)
                {
                    outfile << pIt -> second.platforms[i] << ",";
                }
                outfile << pIt -> second.platforms[pIt -> second.platforms.size() - 1] << endl;
            }
        }
        outfile.close();
    }
}

using namespace processedGWASCatalog;

int main_processGWASCatalog(const vector<string> &all_args)
{
    castGWAS gwas;
    gwas.construct();
    for (int i = 1; i <= 23; i++)
    {
        castSNP SNPs;
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
        SNPs.construct(chrNum);
        gwas.addInfoOutput(SNPs, chrNum);
    }
    return 0;
}
