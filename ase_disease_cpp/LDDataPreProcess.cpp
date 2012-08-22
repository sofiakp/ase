#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;

namespace LDDataPreProcess
{
    void processFile(string &chrNum, string &population)
    {
        string command;
        command = "gunzip -c /home/yulingl/ase_diseases/hapmap_ld/ld_chr" + chrNum + "_" + population + ".txt.gz > /home/yulingl/ase_diseases/hapmap_ld/tempFile";
        system(command.c_str());
        command = "mkdir -p /home/yulingl/ase_diseases/hapmap_ld/" + population + "_1";
        system(command.c_str());
        command = "mkdir -p /home/yulingl/ase_diseases/hapmap_ld/" + population + "_0.1";
        system(command.c_str());
        command = "mkdir -p /home/yulingl/ase_diseases/hapmap_ld/" + population + "_SNPSet";
        system(command.c_str());
        
        set<string> SNPonHapMap;
        ifstream infile;
        ofstream outfileSet, outfileLD1, outfileLD01;
        
        string fileName = "/home/yulingl/ase_diseases/hapmap_ld/" + population + "_1/" + chrNum;
        outfileLD1.open(fileName.c_str());
        fileName = "/home/yulingl/ase_diseases/hapmap_ld/" + population + "_0.1/" + chrNum;
        outfileLD01.open(fileName.c_str());
        fileName = "/home/yulingl/ase_diseases/hapmap_ld/" + population + "_SNPSet/" + chrNum;
        outfileSet.open(fileName.c_str());
        
        infile.open("/home/yulingl/ase_diseases/hapmap_ld/tempFile");
        
        vector<string> info;
        string currentLine;
        double r2;
        while (!getline(infile, currentLine).eof())
        {
            split(info, currentLine, is_any_of(" "));
            SNPonHapMap.insert(info[3]);
            SNPonHapMap.insert(info[4]);
            r2 = atof(info[6].c_str());
            if (r2 >= 0.1)
            {
                outfileLD01 << info[3] << "\t" << info[4] << endl;
                if (r2 >= 1)
                {
                    outfileLD1 << info[3] << "\t" << info[4] << endl;
                }
            }
        }
        
        for (set<string>::iterator it = SNPonHapMap.begin(); it != SNPonHapMap.end(); it ++)
        {
            outfileSet << *it << endl;
        }
        infile.close();
        outfileSet.close();
        outfileLD01.close();
        outfileLD1.close();
        
        system("rm -f /home/yulingl/ase_diseases/hapmap_ld/tempFile");
    }
}

using namespace LDDataPreProcess;

int main_LDDataPreProcess(const vector<string> &all_args)
{
    if (all_args.size() != 3)
    {
        cout << "Usage: ase_disease LDDataPreProcess population" << endl;
        exit(0);
    }
    
    string population = all_args[2];
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
        processFile(chrNum, population);
    }
    return 0;
}

