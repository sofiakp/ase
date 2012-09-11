#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <math.h>

using namespace std;
using namespace boost;

namespace splitLD
{
    struct castTSS
    {
        vector<unsigned long> TSS;
        void construct(string &chrNum);
        unsigned long calculateDistance(unsigned long &pos);
    };
    
    int computeBinNum(double &floatAllelFreq, unsigned long distance, bool &illumina, bool &affymetrix);
    
    void castTSS::construct(string &chrNum)
    {
        ifstream infile;
        string prefix = "/home/yulingl/ase_diseases/geneStart/geneStartinChr";
        string fileName;
        fileName = prefix + chrNum;
        infile.open(fileName.c_str());
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            unsigned long pos;
            pos = strtoul(info[2].c_str(), NULL, 0);
            TSS.push_back(pos);
        }
        sort(TSS.begin(), TSS.end());
        infile.close();
    }
    
    unsigned long castTSS::calculateDistance(unsigned long &pos)
    {
        int start = 0, end = TSS.size() - 1;
        if (pos < TSS[start])
        {
            return (TSS[start] - pos);
        }
        if (pos >= TSS[end])
        {
            return (pos - TSS[end]);
        }
        while (start < end - 1)
        {
            int middle = (start + end) / 2;
            if (pos >= TSS[start] && pos < TSS[middle])
            {
                end = middle;
            }
            else
            {
                start = middle;
            }
        }
        if (pos - TSS[start] <= TSS[end] - pos)
        {
            return (pos - TSS[start]);
        }
        else
        {
            return (TSS[end] - pos);
        }
    }
    
    int computeBinNum(double &floatAllelFreq, unsigned long distance, bool &illumina, bool &affymetrix)
    {
        int allelFreqBinNum;
        int logDistanceBinNum;
        int platform;
        
        srand(time(NULL));
        int randNum = rand() % 2;
        if (floatAllelFreq == 0.5)
        {
            allelFreqBinNum = 4;
        }
        else
        {
            allelFreqBinNum = int(floor(floatAllelFreq / 0.1));
        }
        
        double logDistance = log10(distance);
        if (logDistance <= 3.5)
        {
            logDistanceBinNum = 0;
        }
        else if (logDistance > 6)
        {
            logDistanceBinNum = 26;
        }
        else
        {
            logDistanceBinNum = int(ceil((logDistance - 3.5) / 0.1));
        }
        if (illumina && !affymetrix)
        {
            platform = 0;
        }
        else if (!illumina && affymetrix)
        {
            platform = 1;
        }
        else
        {
            platform = randNum;
        }
        return logDistanceBinNum + allelFreqBinNum * 27 + platform * 27 * 5;
    }
    
    struct castSNPInfo
    {
        string allelFreq;
        string func;
        string pos;
    };
    
    struct castSNP
    {
        set<string> ucsc;
        map<string, castSNPInfo> ucscToInfo;
        set<string> onChip;
        set<string> onAffy, onIllumina;
        set<string> onHapMap;
        void construct(std::string &chrNum, std::string &population);
    };
    
    void castSNP::construct(string &chrNum, string &population)
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
            string key = chrNum + "_" + info[1];
            ucsc.insert(key);
            castSNPInfo tmpSNPInfo;
            tmpSNPInfo.allelFreq = info[3];
            tmpSNPInfo.func = info[4];
            tmpSNPInfo.pos = info[2];
            ucscToInfo[key] = tmpSNPInfo;
        }
        infile.close();
        
        prefix = "/home/yulingl/ase_diseases/SNPonChips/allSNPonChip/SNPOnChipInChr";
        fileName = prefix + chrNum;
        infile.open(fileName.c_str());
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            string key = chrNum + "_" + info[1];
            onChip.insert(key);
        }
        infile.close();
        
        fileName = "/home/yulingl/ase_diseases/SNPonChips/Affy.chr" + chrNum + ".rsids";
        infile.open(fileName.c_str());
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            string key = chrNum + "_" + info[1];
            onAffy.insert(key);
        }
        infile.close();
        
        fileName = "/home/yulingl/ase_diseases/SNPonChips/Illumina.chr" + chrNum + ".rsids";
        infile.open(fileName.c_str());
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            string key = chrNum + "_" + info[1];
            onIllumina.insert(key);
        }
        infile.close();
        
        fileName = "/home/yulingl/ase_diseases/hapmap_ld/" + population + "_SNPSet/" + chrNum;
        infile.open(fileName.c_str());
        while (!getline(infile, currentLine).eof())
        {
            string key = chrNum + "_" + currentLine;
            onHapMap.insert(key);
        }
        infile.close();
        
    }
    
    struct castOutputFilePointersForProcessingLDToBins
    {
        map<int, ofstream *> outputFilePointers;
        void open(string &population);
        void close();
    };
    
    void castOutputFilePointersForProcessingLDToBins::open(string &population)
    {
        string command = "mkdir -p /home/yulingl/ase_diseases/hapmap_ld/" + population + "_splitByBin";
        if (system(command.c_str()) != 0)
        {
            exit(1);
        }
        
        for (int i = 0; i <= 269; ++ i)
        {
            ofstream *tmpOutfile = new ofstream;
            
            stringstream convert;
            convert << i;
            
            string fileName = "/home/yulingl/ase_diseases/hapmap_ld/" + population + "_splitByBin/" + convert.str();
            tmpOutfile -> open(fileName.c_str());
            outputFilePointers[i] = tmpOutfile;
        }
    }
    
    void castOutputFilePointersForProcessingLDToBins::close()
    {
        for (map<int, ofstream *>::iterator it = outputFilePointers.begin(); it != outputFilePointers.end(); it ++)
        {
            it -> second -> close();
        }
    }
    
    struct NULLSNPInfo
    {
        double allelFreq;
        string function;
        bool illumina, affymetrix;
        unsigned long coord;
    };
    
    int computeTagForKey(string &key, castSNP &SNPs, castTSS &TSS)
    {
        map<string, castSNPInfo>::iterator pItucscToInfo = SNPs.ucscToInfo.find(key);
        NULLSNPInfo tmpInfo;
        tmpInfo.allelFreq = atof(pItucscToInfo->second.allelFreq.c_str());
        tmpInfo.coord = strtoul(pItucscToInfo->second.pos.c_str(), NULL, 0);
        if (tmpInfo.allelFreq > 0.5)
        {
            tmpInfo.allelFreq = 1 - tmpInfo.allelFreq;
        }
        tmpInfo.function = pItucscToInfo->second.func;
        
        if (SNPs.onAffy.find(key) != SNPs.onAffy.end())
        {
            tmpInfo.affymetrix = 1;
        }
        else
        {
            tmpInfo.affymetrix = 0;
        }
        
        if (SNPs.onIllumina.find(key) != SNPs.onIllumina.end())
        {
            tmpInfo.illumina = 1;
        }
        else
        {
            tmpInfo.illumina = 0;
        }
        
        return computeBinNum(tmpInfo.allelFreq, TSS.calculateDistance(tmpInfo.coord), tmpInfo.illumina, tmpInfo.affymetrix);
    }
    
    void separateLDByBin(string &chrNum, castOutputFilePointersForProcessingLDToBins &outputFilePointersForProcessingLDToBins, string &population, castSNP &SNPs, castTSS &TSS)
    {
        ifstream infile;
        string fileName = "/home/yulingl/ase_diseases/hapmap_ld/" + population + "_0.1/" + chrNum;
        infile.open(fileName.c_str());
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            string ID1 = chrNum + "_" + info[0];
            string ID2 = chrNum + "_" + info[1];
            
            if (SNPs.ucsc.find(ID1) != SNPs.ucsc.end() && SNPs.ucsc.find(ID2) != SNPs.ucsc.end())
            {
                int binNum1 = computeTagForKey(ID1, SNPs, TSS);
                int binNum2 = computeTagForKey(ID2, SNPs, TSS);
                
                if (binNum1 > 269 or binNum1 < 0)
                {
                    cout << binNum1 << "\t" << ID1 << endl;
                    continue;
                }
                if (binNum2 > 269 or binNum2 < 0)
                {
                    cout << binNum2 << "\t" << ID2 << endl;
                    continue;
                }
                
                if (binNum1 == binNum2)
                {
                    *(outputFilePointersForProcessingLDToBins.outputFilePointers[binNum1]) << ID1 << "\t" << ID2 << "\t" << endl;
                }
                else
                {
                    *(outputFilePointersForProcessingLDToBins.outputFilePointers[binNum1]) << ID1 << "\t" << ID2 << "\t" << endl;
                    *(outputFilePointersForProcessingLDToBins.outputFilePointers[binNum2]) << ID1 << "\t" << ID2 << "\t" << endl;
                }
            }
        }
    }
}

using namespace splitLD;

int main_splitLD(const vector<string> &all_args)
{
    if (all_args.size() != 3)
    {
        cout << "Usage: ase_disease splitLD population";
        exit(0);
    }
    
    castOutputFilePointersForProcessingLDToBins outputFilePointersForProcessingLDToBins;
    string population = all_args[2];
    outputFilePointersForProcessingLDToBins.open(population);
    
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
        
        castTSS TSS;
        TSS.construct(chrNum);
        
        castSNP SNP;
        SNP.construct(chrNum, population);
        
        separateLDByBin(chrNum, outputFilePointersForProcessingLDToBins, population, SNP, TSS);
    }
    outputFilePointersForProcessingLDToBins.close();
    return 0;
}
