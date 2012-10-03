#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <math.h>
#include "boost/program_options.hpp"
namespace po = boost::program_options;

using namespace std;
using namespace boost;

#define GWASINPUTLOCATION "/home/yulingl/ase_diseases/gwasCatalog/fairBinedGWASCatalog/"
#define EQTLINPUTLOCATION "/home/yulingl/ase_diseases/dsQTLs/mergedResult/fairBineddsQTLs/"
#define MAPPEDNULLSNPSTOGWASLOCATION "/home/yulingl/ase_diseases/tmpFiles/matchedNULLSetsToGWASByItByBin/"
#define MAPPEDNULLSNPSTOEQTLLOCATION "/home/yulingl/ase_diseases/tmpFiles/matchedNULLSetsTodsQTLByItByBin/"

double single_sample_t_test(double M, double Sm, double Sd, unsigned Sn)
{
    double diff = Sm - M;
    unsigned v = Sn - 1;
    double t_stat = diff / Sd;
    math::students_t dist(v);
    double q = cdf(complement(dist, fabs(t_stat)));
    return q * 2;
}

namespace prepareData
{
    vector<string> inputFileNames, inputDataSetNames;
    string subject, population, outputDirectory;
    bool eQTLOn = 0, GWASOn = 0;
    
    void init(const vector<string> &all_args)
    {
        po::options_description desc("Compute the statistics for enrichment. You can specify multiple datasets, but the datasets should be from the same subject.");
        desc.add_options()
        ("help,h", "display help message")
        ("inputs,i", po::value< vector<string> >(&inputFileNames), "specify all the dataset separated by space that store the SNPs which need to be overlapped.")
        ("datasets,d", po::value< vector<string> >(&inputDataSetNames), "specify all the name of the subjects (the name of the NULL set) in the same order as the inputs.")
        ("subject,s", po::value<string>(&subject), "specify the name for the subject.")
        ("population,p", po::value<string>(&population), "specify the population.")
        ("output,o", po::value<string>(&outputDirectory), "specify the directory for the output files.")
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
        
        string command = string("mkdir -p ") + outputDirectory;
        if (system(command.c_str()) != 0)
        {
            cerr << "Failed to execute command " << command << endl;
            exit(1);
        };
    }
    
    struct castLD
    {
        map<string, set<string> > tag_SNPinLD;
        void constructOverChr(string &chrNum, string &population);
    };
    
    void castLD::constructOverChr(string &chrNum, string &population)
    {
        ifstream infile;
        string fileName = "/home/yulingl/ase_diseases/hapmap_ld/" + population + "_1/" + chrNum;
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
    
    struct castAssociationSNP
    {
        set<string> SNPs;
        set<string> binNum;
        map<string, set<string> > associationSNP_SNPsinLD;
        void accumulateOverChr(string &chrNum, string &subject, bool GWASoreQTL);
        void mapSNPinLD(castLD &LD);
        void construct(string &subject, bool GWASoreQTL);
    };
    
    void castAssociationSNP::construct(string &subject, bool GWASoreQTL)
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
    
    void castAssociationSNP::accumulateOverChr(string &chrNum, string &subject, bool GWASoreQTL)
    {
        ifstream infile;
        string fileName;
        if (GWASoreQTL == 1)
        {
            fileName = GWASINPUTLOCATION + subject + "/" + chrNum;
        }
        else
        {
            fileName = EQTLINPUTLOCATION + subject + "/" + chrNum;
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
            SNPs.insert(info[0]);
            binNum.insert(info[1]);
        }
    }
    
    void castAssociationSNP::mapSNPinLD(castLD &LD)
    {
        associationSNP_SNPsinLD.clear();
        for (set<string>::iterator it = SNPs.begin(); it != SNPs.end(); ++ it)
        {
            map<string, set<string> >::iterator subIt = LD.tag_SNPinLD.find(*it);
            if (subIt == LD.tag_SNPinLD.end())
            {
                associationSNP_SNPsinLD[*it].insert(*it);
            }
            else
            {
                associationSNP_SNPsinLD[*it] = subIt -> second;
            }
        }
    }
    
    struct castNULLSNP
    {
        set<string> NULLSNPSet;
        map<string, set<string> > NULLSNP_SNPsinLD;
        void construct(castAssociationSNP &associationSNP, int &iteration, string subject, bool GWASoreQTL);
        void mapSNPinLD(castLD &LD);
    };
    
    void castNULLSNP::construct(castAssociationSNP &associationSNP, int &iteration, string subject, bool GWASoreQTL)
    {
        stringstream convert;
        convert << iteration;
        for (set<string>::iterator it = associationSNP.binNum.begin(); it != associationSNP.binNum.end(); it ++)
        {
            ifstream infile;
            string fileName;
            if (GWASoreQTL == 1)
            {
                fileName = MAPPEDNULLSNPSTOGWASLOCATION + subject + "/" + convert.str() + "/" + *it;
            }
            else
            {
                fileName = MAPPEDNULLSNPSTOEQTLLOCATION + subject + "/" + convert.str() + "/" + *it;
            }
            infile.open(fileName.c_str());
            if (infile.fail())
            {
                cerr << "Failed to open " << fileName << endl;
                exit(1);
            }
            
            string rsID;
            while (!getline(infile, rsID).eof())
            {
                NULLSNPSet.insert(rsID);
            }
            infile.close();
        }
    }
    
    void castNULLSNP::mapSNPinLD(castLD &LD)
    {
        NULLSNP_SNPsinLD.clear();
        for (set<string>::iterator it = NULLSNPSet.begin(); it != NULLSNPSet.end(); it ++)
        {
            map<string, set<string> >::iterator subIt = LD.tag_SNPinLD.find(*it);
            if (subIt == LD.tag_SNPinLD.end())
            {
                NULLSNP_SNPsinLD[*it].insert(*it);
            }
            else
            {
                NULLSNP_SNPsinLD[*it] = subIt -> second;
            }
        }
    }
    
    struct castAnnotatedSNPs
    {
        map<string, set<string> > assay_SNPs;
        void fillIndividualAssay(string &fileDirectory, string &datasetName);
    };
    
    void castAnnotatedSNPs::fillIndividualAssay(string &fileDirectory, string &datasetName)
    {
        ifstream infile;
        infile.open(fileDirectory.c_str());
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
                assay_SNPs[datasetName].insert(SNPID);
            }
        }
        infile.close();
    }
}

namespace computeStatWithoutLD
{
    using namespace prepareData;
    struct castStatisticalValue
    {
        double mean;
        double variance;
    };
    
    struct castStatistics
    {
        map<string, int> observed;
        map<string, set<string> > overlaps;
        map<string, vector<int> > expected;
        map<string, castStatisticalValue> expectedStatistics;
        void computeObserved(castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP);
        void computeExpected(castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP, int &n, string subject, bool GWASoreQTL);
        void computeStatisticsForExpected(int &n);
        void output(bool GWASoreQTL);
        int intersectionSize(set<string> &set1, set<string> &set2, set<string> &intersection);
        int intersectionSize(set<string> &set1, set<string> &set2);
    };
    
    int castStatistics::intersectionSize(set<string> &set1, set<string> &set2, set<string> &intersection)
    {
        set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), inserter(intersection, intersection.end()));
        return intersection.size();
    }
    
    int castStatistics::intersectionSize(set<string> &set1, set<string> &set2)
    {
        set<string> intersection;
        set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), inserter(intersection, intersection.end()));
        return intersection.size();
    }
    
    void castStatistics::computeObserved(castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP)
    {
        for (map<string, set<string> >::iterator it = annotatedSNPs.assay_SNPs.begin(); it != annotatedSNPs.assay_SNPs.end(); it ++)
        {
            observed[it -> first] = intersectionSize(it -> second, associationSNP.SNPs, overlaps[it -> first]);
        }
    }
    
    void castStatistics::computeExpected(castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP, int &n, string subject, bool GWASoreQTL)
    {
        for (map<string, set<string> >::iterator it = annotatedSNPs.assay_SNPs.begin(); it != annotatedSNPs.assay_SNPs.end(); it ++)
        {
            for (int i = 0; i < n; i ++)
            {
                castNULLSNP NULLSNP;
                if (GWASoreQTL)
                {
                    NULLSNP.construct(associationSNP, i ,subject, 1);
                }
                else
                {
                    NULLSNP.construct(associationSNP, i ,subject, 0);
                }
                expected[it -> first].push_back(intersectionSize(it -> second, NULLSNP.NULLSNPSet));
            }
        }
    }
    
    void castStatistics::computeStatisticsForExpected(int &n)
    {
        for (map<string, vector<int> >::iterator it = expected.begin(); it != expected.end(); it ++)
        {
            double EX = 0, EX2 = 0;
            for (vector<int>::iterator subIt = it -> second.begin(); subIt != it -> second.end(); subIt ++)
            {
                EX += *subIt;
                EX2 += pow(*subIt, 2.0);
            }
            EX = EX / n;
            EX2 = EX2 / n;
            double variance = EX2 - pow(EX, 2.0);
            castStatisticalValue tmpStatisticalValue;
            tmpStatisticalValue.mean = EX;
            tmpStatisticalValue.variance = variance;
            expectedStatistics[it -> first] = tmpStatisticalValue;
        }
    }
    
    void castStatistics::output(bool GWASoreQTL)
    {
        ofstream outfile;
        string fileName;
        if (GWASoreQTL == 1)
        {
            fileName = outputDirectory + "/" + subject + "_GWAS_withoutLD_statistics.txt";
        }
        else
        {
            fileName = outputDirectory + "/" + subject + "_ds_withoutLD_statistics.txt";
        }
        
        outfile.open(fileName.c_str());
        if (outfile.fail())
        {
            cerr << "Failed to open " << fileName << endl;
            exit(1);
        }
        
        outfile << "Dataset" << "\t" << "Observed" << "\t" << "Expected_mean" << "\t" << "Expected_variance" << "\t" << "Enrichment" << "\t" << "P-value" << endl;
        for (map<string, int>::iterator it = observed.begin(); it != observed.end(); it ++)
        {
            outfile << it -> first << "\t" << observed[it -> first] << "\t" << expectedStatistics[it -> first].mean << "\t" << expectedStatistics[it -> first].variance << "\t" << double(observed[it -> first])/double(expectedStatistics[it -> first].mean) << "\t" << single_sample_t_test(observed[it -> first], expectedStatistics[it -> first].mean, sqrt(expectedStatistics[it -> first].variance), 100) << endl;
        }
        
        outfile.close();
        
        if (GWASoreQTL == 1)
        {
            fileName = outputDirectory + "/" + subject + "_GWAS_overlapping_rsID.txt";
        }
        else
        {
            fileName = outputDirectory + "/" + subject + "_ds_overlapping_rsID.txt";
        }
        
        outfile.open(fileName.c_str());
        if (outfile.fail())
        {
            cerr << "Failed to open " << fileName << endl;
            exit(1);
        }
        outfile << "Dataset" << "\t" << "chr" << "\t" << "rsID" << endl;
        vector<string> tmp;
        for (map<string, set<string> >::iterator it = overlaps.begin(); it != overlaps.end(); it ++)
        {
            for (set<string>::iterator subIt = it -> second.begin(); subIt != it -> second.end(); ++ subIt)
            {
                split(tmp, *subIt, is_any_of("_"));
                outfile << it -> first << "\tchr" << tmp[0] << "\t" << tmp[1] << endl;
            }
        }
        outfile.close();
    }
}

namespace computeStatWithLD
{
    using namespace prepareData;
    struct castStatisticalValue
    {
        double mean;
        double variance;
    };
    
    struct castStatistics
    {
        map<string, int> observed;
        map<string, set<string> > overlaps;
        map<string, vector<int> > expected;
        map<string, castStatisticalValue> expectedStatistics;
        void computeObserved(castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP);
        void computeExpected(castLD &LD, castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP, int &n, string subject, bool GWASoreQTL);
        void computeStatisticsForExpected(int &n);
        void output(bool GWASoreQTL);
    };
    
    void castStatistics::computeObserved(castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP)
    {
        for (map<string, set<string> >::iterator it = annotatedSNPs.assay_SNPs.begin(); it != annotatedSNPs.assay_SNPs.end(); it ++)
        {
            int count = 0;
            for (map<string, set<string> >::iterator it1 = associationSNP.associationSNP_SNPsinLD.begin(); it1 != associationSNP.associationSNP_SNPsinLD.end(); it1 ++)
            {
                for (set<string>::iterator it2 = it1 -> second.begin(); it2 != it1 -> second.end(); it2 ++)
                {
                    set<string>::iterator it3 = it -> second.find(*it2);
                    if (it3 != it -> second.end())
                    {
                        count ++;
                        break;
                    }
                }
            }
            observed[it -> first] = count;
        }
    }
    
    void castStatistics::computeExpected(castLD &LD, castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP, int &n, string subject, bool GWASoreQTL)
    {
        for (map<string, set<string> >::iterator it = annotatedSNPs.assay_SNPs.begin(); it != annotatedSNPs.assay_SNPs.end(); it ++)
        {
            for (int i = 0; i < n; i ++)
            {
                castNULLSNP NULLSNP;
                if (GWASoreQTL)
                {
                    NULLSNP.construct(associationSNP, i ,subject, 1);
                }
                else
                {
                    NULLSNP.construct(associationSNP, i ,subject, 0);
                }
                NULLSNP.mapSNPinLD(LD);
                int count = 0;
                for (map<string, set<string> >::iterator it1 = NULLSNP.NULLSNP_SNPsinLD.begin(); it1 != NULLSNP.NULLSNP_SNPsinLD.end(); it1 ++)
                {
                    for (set<string>::iterator it2 = it1 -> second.begin(); it2 != it1 -> second.end(); it2 ++)
                    {
                        set<string>::iterator it3 = it -> second.find(*it2);
                        if (it3 != it -> second.end())
                        {
                            count ++;
                            break;
                        }
                    }
                }
                expected[it -> first].push_back(count);
            }
        }
    }
    
    void castStatistics::computeStatisticsForExpected(int &n)
    {
        for (map<string, vector<int> >::iterator it = expected.begin(); it != expected.end(); it ++)
        {
            double EX = 0, EX2 = 0;
            for (vector<int>::iterator subIt = it -> second.begin(); subIt != it -> second.end(); subIt ++)
            {
                EX += *subIt;
                EX2 += pow(*subIt, 2.0);
            }
            EX = EX / n;
            EX2 = EX2 / n;
            double variance = EX2 - pow(EX, 2.0);
            castStatisticalValue tmpStatisticalValue;
            tmpStatisticalValue.mean = EX;
            tmpStatisticalValue.variance = variance;
            expectedStatistics[it -> first] = tmpStatisticalValue;
        }
    }
    
    void castStatistics::output(bool GWASoreQTL)
    {
        ofstream outfile;
        string fileName;
        if (GWASoreQTL == 1)
        {
            fileName = outputDirectory + "/" + subject + "_GWAS_withLD_statistics.txt";
        }
        else
        {
            fileName = outputDirectory + "/" + subject + "_ds_withLD_statistics.txt";
        }
        
        outfile.open(fileName.c_str());
        if (outfile.fail())
        {
            cerr << "Failed to open " << fileName << endl;
            exit(1);
        }
        
        outfile << "Dataset" << "\t" << "Observed" << "\t" << "Expected_mean" << "\t" << "Expected_variance" << "\t" << "Enrichment" << "\t" << "P-value" << endl;
        for (map<string, int>::iterator it = observed.begin(); it != observed.end(); it ++)
        {
            outfile << it -> first << "\t" << observed[it -> first] << "\t" << expectedStatistics[it -> first].mean << "\t" << expectedStatistics[it -> first].variance << "\t" << double(observed[it -> first])/double(expectedStatistics[it -> first].mean) << "\t" << single_sample_t_test(observed[it -> first], expectedStatistics[it -> first].mean, sqrt(expectedStatistics[it -> first].variance), 100) << endl;
        }
    }
}

using namespace prepareData;
int main_computeStatistics(const vector<string> &all_args)
{
    init(all_args);
    
    castAnnotatedSNPs annotatedSNPs;
    for (int i = 0; i < inputFileNames.size(); ++ i)
    {
        annotatedSNPs.fillIndividualAssay(inputFileNames[i], inputDataSetNames[i]);
    }
    
    castAssociationSNP GWASSNPs;
    if (GWASOn)
    {
        GWASSNPs.construct(subject, 1);
    }
    
    castAssociationSNP eQTLs;
    if (eQTLOn)
    {
        eQTLs.construct(subject, 0);
    }
    
    {
        using namespace computeStatWithoutLD;
        
        if (GWASOn)
        {
            castStatistics statistics;
            statistics.computeObserved(annotatedSNPs, GWASSNPs);
            int n = 100;
            statistics.computeExpected(annotatedSNPs, GWASSNPs, n, subject, 1);
            statistics.computeStatisticsForExpected(n);
            statistics.output(1);
        }
        
        if (eQTLOn)
        {
            castStatistics statistics;
            statistics.computeObserved(annotatedSNPs, eQTLs);
            int n = 100;
            statistics.computeExpected(annotatedSNPs, eQTLs, n, subject, 0);
            statistics.computeStatisticsForExpected(n);
            statistics.output(0);
        }
    }
    
    {
        using namespace computeStatWithLD;
        castLD LD;
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
            LD.constructOverChr(chrNum, population);
        }
        
        if (GWASOn)
        {
            castStatistics statistics;
            GWASSNPs.mapSNPinLD(LD);
            statistics.computeObserved(annotatedSNPs, GWASSNPs);
            int n = 100;
            statistics.computeExpected(LD, annotatedSNPs, GWASSNPs, n, subject, 1);
            statistics.computeStatisticsForExpected(n);
            statistics.output(1);
        }
        
        if (eQTLOn)
        {
            castStatistics statistics;
            eQTLs.mapSNPinLD(LD);
            statistics.computeObserved(annotatedSNPs, eQTLs);
            int n = 100;
            statistics.computeExpected(LD, annotatedSNPs, eQTLs, n, subject, 0);
            statistics.computeStatisticsForExpected(n);
            statistics.output(0);
        }
    }
    return 0;
}

#undef GWASINPUTLOCATION
#undef EQTLINPUTLOCATION
#undef MAPPEDNULLSNPSTOGWASLOCATION
#undef MAPPEDNULLSNPSTOEQTLLOCATION
