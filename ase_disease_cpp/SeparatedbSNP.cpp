#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include "boost/algorithm/string.hpp"
#include "boost/program_options.hpp"
namespace po = boost::program_options;

using namespace std;
using namespace boost;

namespace separatedbSNP
{
    string inputFileDirectory, outputFileDirectory;
    
    struct castOutputFilePointers
    {
        ofstream outputFilePointerVector[23];
        void open(string &outputFileDirectory);
        void close();
    };
    
    void init(const vector<string> &all_args)
    {
        po::options_description desc("Separate dbSNP into sub-datasets according to chromesomes to facilitate further analysis");
        desc.add_options()
        ("help,h", "display help message")
        ("inputfile-directory,i", po::value<string>(&inputFileDirectory), "specify the dbSNP file directory")
        ("output,o", po::value<string>(&outputFileDirectory), "specify the splited files directory");
        po::variables_map vm;
        const char* av[all_args.size()];
        for (int i = 0; i < all_args.size(); ++ i)
        {
            av[i] = all_args[i].c_str();
        }
        po::store(po::parse_command_line(all_args.size(), av, desc), vm);
        po::notify(vm);
        if (vm.count("help")) {
            cout << desc << "\n";
            exit(0);
        }
    }
    
    void castOutputFilePointers::open(string &outputFileDirectory)
    {
        string commandLine = "mkdir -p " + outputFileDirectory;
        system(commandLine.c_str());
        string fileName, prefix = "/dbSNP135_hg19_chr";
        for (int i = 1; i <= 22; i++)
        {
            stringstream number;
            number << i;
            fileName = outputFileDirectory + prefix + number.str();
            outputFilePointerVector[i - 1].open(fileName.c_str());
        }
        fileName = outputFileDirectory + prefix + "X";
        outputFilePointerVector[22].open(fileName.c_str());
    }
    
    void castOutputFilePointers::close()
    {
        for (int i = 0; i < 23; i++)
        {
            outputFilePointerVector[i].close();
        }
    }
    
    void processFile(ifstream &infile, castOutputFilePointers &outputFilePointers)
    {
        set<string> validChrName;
        string prefix = "chr";
        for (int i = 1; i <= 22; i++)
        {
            stringstream number;
            number << i;
            string tmp;
            tmp = prefix + number.str();
            validChrName.insert(tmp);
        }
        validChrName.insert("chrX");
        vector<string> info;
        string currentLine;
        while (!getline(infile, currentLine).eof())
        {
            if (currentLine[0] == '#')
            {
                continue;
            }
            split(info, currentLine, is_any_of("\t"));
            if (validChrName.find(info[0]) != validChrName.end())
            {
                if (info[0] == "chrX")
                {
                    outputFilePointers.outputFilePointerVector[22] << currentLine << endl;
                }
                else
                {
                    int chrNum;
                    string tmp = info[0];
                    chrNum = atoi(tmp.erase(0, 3).c_str());
                    outputFilePointers.outputFilePointerVector[chrNum -1] << currentLine << endl;
                }
            }
        }
    }
}

using namespace separatedbSNP;

int main_separatedbSNP(const vector<string> &all_args)
{
    ifstream infile;
    init(all_args);
    castOutputFilePointers outputFilePointers;
    outputFilePointers.open(outputFileDirectory);
    infile.open(inputFileDirectory.c_str());
    processFile(infile, outputFilePointers);
    infile.close();
    outputFilePointers.close();
    return 0;
}

