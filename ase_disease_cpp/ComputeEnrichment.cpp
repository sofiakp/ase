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
    vector<string> inputFileNames, inputDataSetNames, subjects;
    string population;
    
    void init(const vector<string> &all_args)
    {
        po::options_description desc("Compute the statistics for enrichment. You can specify multiple datasets, but the datasets should be from the same subject.");
        desc.add_options()
        ("help,h", "display help message")
        ("inputs,i", po::value< vector<string> >(&inputFileNames), "specify all the dataset separated by space that store the SNPs which need to be overlapped")
        ("datasets,d", po::value< vector<string> >(&inputDataSetNames), "specify all the name of the datasets in the same order as the inputs")
        ("subject,s", po::value< vector<string> >(&subjects), "specify the name for the subject.")
        ("population,p", po::value<string>(&population), "specify the population.");
        po::variables_map vm;
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
        void accumulateOverChr(string &chrNum);
        void mapSNPinLD(castLD &LD);
        void construct();
    };
    
    void castAssociationSNP::construct()
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
            accumulateOverChr(chrNum);
        }
    }
    
    void castAssociationSNP::accumulateOverChr(string &chrNum)
    {
        ifstream infile;
        string fileName = "/home/yulingl/ase_diseases/gwasCatalog/fairBinedGWASCatalog/" + chrNum;
        infile.open(fileName.c_str());
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
        void construct(castAssociationSNP &associationSNP, int &iteration);
        void mapSNPinLD(castLD &LD);
    };
    
    void castNULLSNP::construct(castAssociationSNP &associationSNP, int &iteration)
    {
        stringstream convert;
        convert << iteration;
        for (set<string>::iterator it = associationSNP.binNum.begin(); it != associationSNP.binNum.end(); it ++)
        {
            ifstream infile;
            string fileName = "/home/yulingl/ase_diseases/tmpFiles/matchedNULLSetsByItByBin/" + convert.str() + "/" + *it;
            infile.open(fileName.c_str());
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
        void fillIndividualAssay(string &fileDirectory);
    };
    
    void castAnnotatedSNPs::fillIndividualAssay(string &fileDirectory)
    {
        ifstream infile;
        infile.open(fileDirectory.c_str());
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of("\t"));
            assay_SNPs[info[0]].insert(info[1]);
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
        map<string, vector<int> > expected;
        map<string, castStatisticalValue> expectedStatistics;
        void computeObserved(castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP);
        void computeExpected(castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP, int &n);
        void computeStatisticsForExpected(int &n);
        void output();
        int intersectionSize(set<string> &set1, set<string> &set2);
    };
    
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
            observed[it -> first] = intersectionSize(it -> second, associationSNP.SNPs);
        }
    }
    
    void castStatistics::computeExpected(castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP, int &n)
    {
        for (map<string, set<string> >::iterator it = annotatedSNPs.assay_SNPs.begin(); it != annotatedSNPs.assay_SNPs.end(); it ++)
        {
            for (int i = 0; i < n; i ++)
            {
                castNULLSNP NULLSNP;
                NULLSNP.construct(associationSNP, i);
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
    
    void castStatistics::output()
    {
        cout << "Dataset" << "\t" << "Observed" << "\t" << "Expected_mean" << "\t" << "Expected_variance" << "\t" << "Enrichment" << "\t" << "P-value" << endl;
        for (map<string, int>::iterator it = observed.begin(); it != observed.end(); it ++)
        {
            cout << it -> first << "\t" << observed[it -> first] << "\t" << expectedStatistics[it -> first].mean << "\t" << expectedStatistics[it -> first].variance << "\t" << double(observed[it -> first])/double(expectedStatistics[it -> first].mean) << "\t" << single_sample_t_test(observed[it -> first], expectedStatistics[it -> first].mean, sqrt(expectedStatistics[it -> first].variance), 100) << endl;
        }
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
        map<string, vector<int> > expected;
        map<string, castStatisticalValue> expectedStatistics;
        void computeObserved(castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP);
        void computeExpected(castLD &LD, castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP, int &n);
        void computeStatisticsForExpected(int &n);
        void output();
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
    
    void castStatistics::computeExpected(castLD &LD, castAnnotatedSNPs &annotatedSNPs, castAssociationSNP &associationSNP, int &n)
    {
        for (map<string, set<string> >::iterator it = annotatedSNPs.assay_SNPs.begin(); it != annotatedSNPs.assay_SNPs.end(); it ++)
        {
            for (int i = 0; i < n; i ++)
            {
                castNULLSNP NULLSNP;
                NULLSNP.construct(associationSNP, i);
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
    
    void castStatistics::output()
    {
        cout << "Dataset" << "\t" << "Observed" << "\t" << "Expected_mean" << "\t" << "Expected_variance" << "\t" << "Enrichment" << "\t" << "P-value" << endl;
        for (map<string, int>::iterator it = observed.begin(); it != observed.end(); it ++)
        {
            cout << it -> first << "\t" << observed[it -> first] << "\t" << expectedStatistics[it -> first].mean << "\t" << expectedStatistics[it -> first].variance << "\t" << double(observed[it -> first])/double(expectedStatistics[it -> first].mean) << "\t" << single_sample_t_test(observed[it -> first], expectedStatistics[it -> first].mean, sqrt(expectedStatistics[it -> first].variance), 100) << endl;
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
        annotatedSNPs.fillIndividualAssay(inputFileNames[i]);
    }
    
    castAssociationSNP associationSNP;
    associationSNP.construct();
    
    {
        using namespace computeStatWithoutLD;
        castStatistics statistics;
        statistics.computeObserved(annotatedSNPs, associationSNP);
        int n = 100;
        statistics.computeExpected(annotatedSNPs, associationSNP, n);
        statistics.computeStatisticsForExpected(n);
        statistics.output();
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
        
        castStatistics statistics;
        associationSNP.mapSNPinLD(LD);
        statistics.computeObserved(annotatedSNPs, associationSNP);
        int n = 100;
        statistics.computeExpected(LD, annotatedSNPs, associationSNP, n);
        statistics.computeStatisticsForExpected(n);
        statistics.output();
    }
    return 0;
}
