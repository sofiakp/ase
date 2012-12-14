#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <boost/algorithm/string.hpp>
#include <math.h>
#include "boost/program_options.hpp"
namespace po = boost::program_options;

using namespace std;
using namespace boost;

#define EQTLINPUTDIRECTORY "/home/yulingl/ase_diseases/dsQTLs/mergedResult/fairBineddsQTLs/"
#define GWASINPUTDIRECTORY "/home/yulingl/ase_diseases/gwasCatalog/fairBinedGWASCatalog/"
#define EQTLOUTPUTDIRECTORY "/home/yulingl/ase_diseases/tmpFiles/matchedNULLSetsTodsQTLByItByBin/"
#define GWASOUTPUTDIRECTORY "/home/yulingl/ase_diseases/tmpFiles/matchedNULLSetsToGWASByItByBin/"

namespace generateBinMatchNULLSet
{
    vector<string> inputFileNames, subjects;
    string population;
    bool eQTLOn = 0, GWASOn = 0;
    map<string, map<int, set<string> > > subject_bin_key;
    
    void init(const vector<string> &all_args)
    {
        po::options_description desc("Prepare the NULL set for matching. The input SNPs should be separated by chromosomes.");
        desc.add_options()
        ("help,h", "display help message")
        ("inputs,i", po::value< vector<string> >(&inputFileNames), "specify all the SNP file directories.")
        ("subjects,s", po::value< vector<string> >(&subjects), "specify all the names for the subjects with the same order as inputs.")
        ("population,p", po::value<string>(&population), "specify the population.")
        ("eQTLOn,e", po::value<bool>(&eQTLOn), "if eQTL analysis is desired, specify 1.")
        ("GWASOn,g", po::value<bool>(&GWASOn), "if GWAS analysis is desired, specify 1.");
        po::variables_map vm;
        const char* av[all_args.size()];
        for (int i = 0; i < all_args.size(); ++ i)
        {
            av[i] = all_args[i].c_str();
        }
        po::store(po::parse_command_line(all_args.size(), av, desc), vm);
        po::notify(vm);
        
        if (all_args.size() <= 2)
        {
            cout << desc << "\n";
            exit(0);
        }
        if (vm.count("help"))
        {
            cout << desc << "\n";
            exit(0);
        }
    }
}

namespace generateBinNullSet
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
        void construct(string &chrNum, string &population);
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
    
    struct NULLSNPInfo
    {
        double allelFreq;
        string function;
        bool illumina, affymetrix;
        unsigned long coord;
    };
    
    struct NULLSet
    {
        set<string> initial_keys;
        set<string> filteredKeys;
        map<string, NULLSNPInfo> key_info;
        
        void construct(string &fileDirectory, string &chrNum);
        void interSect(castSNP &SNPs);
        
        void mapInfo(castSNP &SNPs);
        void generateBin_key(castTSS &TSS, map<int, set<string> > &bin_key);
    };
    
    void NULLSet::construct(string &fileDirectory, string &chrNum)
    {
        ifstream infile;
        
        string fileName = fileDirectory + "/" + chrNum;
        infile.open(fileName.c_str());
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            if (currentLine[0] != '#')
            {
                string SNPID;
                split(info, currentLine, is_any_of("\t"));
                SNPID = info[0].substr(3) + "_" + info[2];
                initial_keys.insert(SNPID);
            }
        }
    }
    
    void NULLSet::interSect(castSNP &SNPs)
    {
        set<string> tmp;
        set_intersection(initial_keys.begin(), initial_keys.end(), SNPs.ucsc.begin(), SNPs.ucsc.end(), inserter(filteredKeys, filteredKeys.end()));
        set_intersection(filteredKeys.begin(), filteredKeys.end(), SNPs.onChip.begin(), SNPs.onChip.end(), inserter(tmp, tmp.end()));
        filteredKeys.clear();
        set_intersection(tmp.begin(), tmp.end(), SNPs.onHapMap.begin(), SNPs.onHapMap.end(), inserter(filteredKeys, filteredKeys.end()));
    }
    
    void NULLSet::mapInfo(castSNP &SNPs)
    {
        for (set<string>::iterator it = filteredKeys.begin(); it != filteredKeys.end(); ++ it)
        {
            map<string, castSNPInfo>::iterator pItucscToInfo = SNPs.ucscToInfo.find(*it);
            NULLSNPInfo tmpInfo;
            tmpInfo.allelFreq = atof(pItucscToInfo->second.allelFreq.c_str());
            tmpInfo.coord = strtoul(pItucscToInfo->second.pos.c_str(), NULL, 0);
            if (tmpInfo.allelFreq > 0.5)
            {
                tmpInfo.allelFreq = 1 - tmpInfo.allelFreq;
            }
            tmpInfo.function = pItucscToInfo->second.func;
            
            if (SNPs.onAffy.find(*it) != SNPs.onAffy.end())
            {
                tmpInfo.affymetrix = 1;
            }
            else
            {
                tmpInfo.affymetrix = 0;
            }
            
            if (SNPs.onIllumina.find(*it) != SNPs.onIllumina.end())
            {
                tmpInfo.illumina = 1;
            }
            else
            {
                tmpInfo.illumina = 0;
            }
            
            key_info[*it] = tmpInfo;
        }
    }
    
    void NULLSet::generateBin_key(castTSS &TSS, map<int, set<string> > &bin_key)
    {
        for (map<string, NULLSNPInfo>::iterator it = key_info.begin(); it != key_info.end(); ++ it)
        {
            bin_key[computeBinNum(it->second.allelFreq, TSS.calculateDistance(it->second.coord), it->second.illumina, it->second.affymetrix)].insert(it->first);
        }
    }
}

namespace generateMatchedNULLSet
{
    struct castLD
    {
        map<string, set<string> > tag_SNPinLD;
        void loadByBin(string &binNum, string &population);
    };
    
    void castLD::loadByBin(string &binNum, string &population)
    {
        ifstream infile;
        string fileName = "/home/yulingl/ase_diseases/hapmap_ld/" + population + "_splitByBin/" + binNum;
        infile.open(fileName.c_str());
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            tag_SNPinLD[info[0]].insert(info[0]);
            tag_SNPinLD[info[0]].insert(info[1]);
            tag_SNPinLD[info[1]].insert(info[0]);
            tag_SNPinLD[info[1]].insert(info[1]);
        }
        infile.close();
    }
    
    struct castTargetSNP
    {
        map<string, int> binNum_SNPNum;
        void accumulateOverChr(string &chrNum, string &subject, bool GWASoreQTL);
        void construct(string &subject, bool GWASoreQTL);
    };
    
    void castTargetSNP::construct(string &subject, bool GWASoreQTL)
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
            accumulateOverChr(chrNum, subject, GWASoreQTL);
        }
    }
    
    void castTargetSNP::accumulateOverChr(string &chrNum, string &subject, bool GWASoreQTL)
    {
        ifstream infile;
        string fileName;
        if (GWASoreQTL == 1)
        {
            fileName = GWASINPUTDIRECTORY + subject + "/" + chrNum;
        }
        else
        {
            fileName = EQTLINPUTDIRECTORY + subject + "/" + chrNum;
        }
        
        infile.open(fileName.c_str());
        if (infile.fail())
        {
            cerr << "Failed to open " << fileName << endl;
            exit(1);
        }
        
        vector<string> info;
        string currentLine;
        string SNPID;
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            map<string, int>::iterator it;
            it = binNum_SNPNum.find(info[1]);
            if (it == binNum_SNPNum.end())
            {
                binNum_SNPNum[info[1]] = 1;
            }
            else
            {
                it -> second ++;
            }
        }
    }
    
    struct castNULLSNP
    {
        vector<map<string, vector<string> > > binNum_SNPVector;
        map<string, string> SNP_binNum;
        set<string> NULLSetByBin;
        void construct(int &n, map<int, set<string> > &bin_key);
        void outputNULLSNPSetByIterationByBin(int &iteration, string &binNum, string &subject, bool GWASoreQTL);
    };
    
    void castNULLSNP::construct(int &n, map<int, set<string> > &bin_key)
    {
        map<string, vector<string> > originalBinNum_SNP;
        
        for (map<int, set<string> >::iterator it = bin_key.begin(); it != bin_key.end(); ++ it)
        {
            stringstream convert;
            convert << it -> first;
            for (set<string>::iterator subIt = it -> second.begin(); subIt != it -> second.end(); ++ subIt)
            {
                originalBinNum_SNP[convert.str()].push_back(*subIt);
            }
        }
        
        for (int i = 0; i < n; i ++)
        {
            map<string, vector<string> > swapBinNum_SNP;
            swapBinNum_SNP = originalBinNum_SNP;
            binNum_SNPVector.push_back(swapBinNum_SNP);
        }
    }
    
    void castNULLSNP::outputNULLSNPSetByIterationByBin(int &iteration, string &binNum, string &subject, bool GWASoreQTL)
    {
        stringstream convert;
        convert << iteration;
        
        string command;
        if (GWASoreQTL == 1)
        {
            command = string("mkdir -p ") + GWASOUTPUTDIRECTORY + subject + "/" + convert.str();
        }
        else
        {
            command = string("mkdir -p ") + EQTLOUTPUTDIRECTORY + subject + "/" + convert.str();
        }
        if (system(command.c_str()) != 0)
        {
            cerr << "Failed to execute command " << command << endl;
            exit(1);
        };
        
        ofstream outfile;
        string fileName;
        if (GWASoreQTL == 1)
        {
            fileName = GWASOUTPUTDIRECTORY + subject + "/" + convert.str() + "/" + binNum;
        }
        else
        {
            fileName = EQTLOUTPUTDIRECTORY + subject + "/" + convert.str() + "/" + binNum;
        }
        
        outfile.open(fileName.c_str());
        if (outfile.fail())
        {
            cerr << "Failed to open " << fileName << endl;
            exit(1);
        }
        
        for (set<string>::iterator it = NULLSetByBin.begin(); it != NULLSetByBin.end(); it ++)
        {
            outfile << *it << endl;
        }
        outfile.close();
    }
    
    struct castBinNum_eQTLNum
    {
        string binNum;
        int eQTLNum;
    };
    
    void generateNULLSet(map<string, int> &binNum_num, string &subject, map<int, set<string> > &bin_key, string &population, bool GWASoreQTL)
    {
        srand(time(NULL));
        vector<castBinNum_eQTLNum> binNum_eQTLNumVector;
        for (map<string, int>::iterator it = binNum_num.begin(); it != binNum_num.end(); it ++)
        {
            castBinNum_eQTLNum tmpBinNum_eQTLNum;
            tmpBinNum_eQTLNum.binNum = it -> first;
            tmpBinNum_eQTLNum.eQTLNum = it -> second;
            binNum_eQTLNumVector.push_back(tmpBinNum_eQTLNum);
        }
        random_shuffle(binNum_eQTLNumVector.begin(), binNum_eQTLNumVector.end());
        castNULLSNP NULLSNP;
        int n = 100;
        NULLSNP.construct(n, bin_key);
        for (vector<castBinNum_eQTLNum>::iterator it = binNum_eQTLNumVector.begin(); it != binNum_eQTLNumVector.end(); it ++)
        {
            castLD LD;
            LD.loadByBin(it -> binNum, population);
            for (int i = 0; i < n; i ++)
            {
                NULLSNP.NULLSetByBin.clear();
                map<string, vector<string> >::iterator seachItForBinNum_SNP = NULLSNP.binNum_SNPVector[i].find(it -> binNum);
                map<string, int> SNP_index;
                for (int j = 0; j < seachItForBinNum_SNP -> second.size(); j ++)
                {
                    if (seachItForBinNum_SNP -> second[j] != "")
                    {
                        SNP_index[seachItForBinNum_SNP -> second[j]] = j;
                    }
                }
                for (int j = 0; j < it -> eQTLNum; j ++)
                {
                    int indexToDelete;
                    while (true)
                    {
                        indexToDelete = rand() % seachItForBinNum_SNP -> second.size();
                        if (seachItForBinNum_SNP -> second[indexToDelete] != "")
                        {
                            break;
                        }
                    }
                    string SNPtoInsert = seachItForBinNum_SNP -> second[indexToDelete];
                    NULLSNP.NULLSetByBin.insert(SNPtoInsert);
                    map<string, set<string> >::iterator searchItForTag_SNPinLD = LD.tag_SNPinLD.find(SNPtoInsert);
                    if (searchItForTag_SNPinLD != LD.tag_SNPinLD.end())
                    {
                        for (set<string>::iterator subIt = searchItForTag_SNPinLD -> second.begin(); subIt != searchItForTag_SNPinLD -> second.end(); subIt ++)
                        {
                            map<string, string>::iterator searchItForSNP_binNum = NULLSNP.SNP_binNum.find(*subIt);
                            if (searchItForSNP_binNum != NULLSNP.SNP_binNum.end())
                            {
                                map<string, int>::iterator searchItForSNP_index = SNP_index.find(*subIt);
                                if (searchItForSNP_index != SNP_index.end())
                                {
                                    NULLSNP.binNum_SNPVector[i][searchItForSNP_binNum -> second][searchItForSNP_index -> second] = "";
                                    SNP_index.erase(searchItForSNP_index);
                                }
                            }
                        }
                    }
                    else
                    {
                        seachItForBinNum_SNP -> second[indexToDelete] = "";
                        SNP_index.erase(SNPtoInsert);
                    }
                }
                NULLSNP.outputNULLSNPSetByIterationByBin(i, it -> binNum, subject, GWASoreQTL);
            }
        }
    }
}

using namespace generateBinMatchNULLSet;
int main_generateBinFilterNULLSet(const vector<string> &all_args)
{
    init(all_args);
    
    {
        using namespace generateBinNullSet;
        
        for (int i = 1; i <= 23; ++ i)
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
            
            castSNP SNPs;
            SNPs.construct(chrNum, population);
            
            for (int j = 0; j < inputFileNames.size(); ++ j)
            {
                NULLSet nullSNP;
                nullSNP.construct(inputFileNames[j], chrNum);
                nullSNP.interSect(SNPs);
                nullSNP.mapInfo(SNPs);
                nullSNP.generateBin_key(TSS, subject_bin_key[subjects[j]]);
            }
        }
    }
    
    {
        using namespace generateMatchedNULLSet;
        
        for (int i = 0; i < subjects.size(); ++ i)
        {
            if (GWASOn)
            {
                castTargetSNP GWASSNPs;
                GWASSNPs.construct(subjects[i], 1);
            
                generateNULLSet(GWASSNPs.binNum_SNPNum, subjects[i], subject_bin_key[subjects[i]], population, 1);
            }
            
            if (eQTLOn)
            {
                castTargetSNP eQTLs;
                eQTLs.construct(subjects[i], 0);
                
                generateNULLSet(eQTLs.binNum_SNPNum, subjects[i], subject_bin_key[subjects[i]], population, 0);
            }
        }
    }
    
    return 0;
}

#undef EQTLINPUTDIRECTORY
#undef GWASINPUTDIRECTORY
#undef EQTLOUTPUTDIRECTORY
#undef GWASOUTPUTDIRECTORY
