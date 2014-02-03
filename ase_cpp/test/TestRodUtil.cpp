#include "swak/System.h"
#include "swak/StringUtil.h"
#include "swak/Helpers.h"
#include "swak/StlUtil.h"
#include "yaml-cpp/yaml.h"
#include "RodUtil.h"

using namespace RodUtil;

namespace TestRodUtil
{
  // User args
  string bam_file;
  string gtf_file;

  // User options
  string cleaned_gtf;

  // Program vars
  vector<string> args;
  YAML::Node options;

  void ProcessCmdLine(const vector<string> &all_args)
  {
    int num_args = 0;
    int num_prog_names = 2;

    if (!ProcessInput(all_args, num_prog_names, num_args, options, args))
    {
      string usage = Basename(all_args[0]) + " test_rod_util [options]";
      string desc = "Tests RodUtil";

      VecS opt_lines;
      //opt_lines.push_back("cleaned_gtf=FILE Output the genes after cleaning (removing dupes) to this file [null]");

      PrintUsage(usage, desc, opt_lines); 
      exit(0);
    }

    //YAML::SetVar(options, "cleaned_gtf", cleaned_gtf, false);

    //gtf_file = args[0];
  }

  void Init()
  {
  }

};


using namespace TestRodUtil;

struct Rod
{
  Rod(int64 s, int64 e) { st = s; en = e; }
  int64 start() const { return st; } 
  int64 end() const { return en; }
  int64 st, en;
};

struct AnotherRod
{
  AnotherRod(int64 s_, int64 e_) { s = s_; e = e_; }
  int64 s, e;
};

int64 astart(const AnotherRod * a)
{
  return a->s;
}

int64 aend(const AnotherRod * a)
{
  return a->e;
}


int main_test_rod_util(const vector<string> &all_args)
{
  /*
  ProcessCmdLine(all_args);
  Init();

  vector<Rod> rods;

  rods.push_back(Rod(10, 30));
  rods.push_back(Rod(20, 40));
  rods.push_back(Rod(20, 35));
  rods.push_back(Rod(23, 26));
  rods.push_back(Rod(28, 29));
  rods.push_back(Rod(30, 70));
  rods.push_back(Rod(50, 60));
  
  RodWalker<Rod, vector<Rod>::const_iterator > walker(rods.begin(), rods.end());

  LessCompEnd<Rod> comp_end; // An object to help us compare the ends of objects
  set<const Rod *, LessCompEnd<Rod> > overlaps(comp_end);
  vector<const Rod *> added, removed;

  StlFor(i, rods)
    cerr << "[" << rods[i].start() << ", " << rods[i].end() << "] ";
  cerr << endl;

  StlFor(i, rods)
  {
    for (int j = 0; j < rods[i].end(); ++j)
    {
      if (j == rods[i].start())
        cerr << j;
      if (j >= rods[i].start() && j < rods[i].end())
        cerr << "-";
      else
        cerr << " ";
    }
    cerr << rods[i].end();
    cerr << endl;
  }

  for (int64 pos = 0; pos < 85; pos += 5)
  {
    walker.JumpTo(pos);
    bool got_overlaps = walker.GetOverlaps(overlaps, added, removed);

    cout << "got?: " << got_overlaps << " pos: " << pos << "\toverlaps:\t";
    //StlForIterConst(iter, set<const Rod *>, overlaps)
    for (set<const Rod *>::const_iterator iter = overlaps.begin(); iter != overlaps.end(); ++iter)
      cout << "[" << (*iter)->start() << ", " << (*iter)->end() << "] ";

    cout << "\tadded:\t";
    StlFor(i, added)
      cout << "[" << added[i]->start() << ", " << added[i]->end() << "] ";
    cout << "\tremoved:\t";
    StlFor(i, removed)
      cout << "[" << removed[i]->start() << ", " << removed[i]->end() << "] ";
    cout << endl;
  }

  vector<AnotherRod> a_rods;
  a_rods.push_back(AnotherRod(10, 20));
  a_rods.push_back(AnotherRod(15, 20));
  a_rods.push_back(AnotherRod(30, 40));


  RodWalker<AnotherRod, vector<AnotherRod>::const_iterator, astart, aend> a_walker;

  a_walker.Init(a_rods.begin(), a_rods.end());
  a_walker.JumpTo(0);
  vector<const AnotherRod *> a_overlaps;
  a_walker.GetOverlaps(a_overlaps);

  */
  return 0;
}
