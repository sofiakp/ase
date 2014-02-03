#ifndef SWAK_HELPERS_H
#define SWAK_HELPERS_H

#include "Swak.h"
#include "System.h"
#include "yaml-cpp/yaml.h"
#include "OptionParser.h"
#include "StlUtil.h"

// Some constants to help color output
namespace Color
{
   const char reset[]   = "\x1b[0m";
   const char red[]     = "\x1b[31m";
   const char green[]   = "\x1b[32m";
   const char yellow[]  = "\x1b[33m";
   const char blue[]    = "\x1b[34m";
   const char magenta[] = "\x1b[35m";
   const char cyan[]    = "\x1b[36m";
}

template <typename T>
string Colorize(const T &item, const string &color_name)
{
  if (color_name == "reset")
    return Colorize(item, 0);
  else if (color_name == "red")
    return Colorize(item, 1);
  else if (color_name == "green")
    return Colorize(item, 2);
  else if (color_name == "yellow")
    return Colorize(item, 3);
  else if (color_name == "blue")
    return Colorize(item, 4);
  else if (color_name == "magenta")
    return Colorize(item, 5);
  else if (color_name == "cyan")
    return Colorize(item, 6);

  throw runtime_error("Illegal color name");
}

template <typename T>
string Colorize(const T &item, char color_letter)
{
  int colorint = 0;
  string color_str = Color::reset;

  if (color_letter == 'e')
    colorint = 0;
  else if (color_letter == 'r')
    colorint = 1;
  else if (color_letter == 'g')
    colorint = 2;
  else if (color_letter == 'y')
    colorint = 3;
  else if (color_letter == 'b')
    colorint = 4;
  else if (color_letter == 'm')
    colorint = 5;
  else if (color_letter == 'c')
    colorint = 6;
  else
    throw runtime_error("Illegal color letter");

  if (colorint > 0)
    color_str = "\x1b[3" + string(1, '0' + colorint) + "m";

  stringstream ss;
  ss << color_str << item << Color::reset;
  return ss.str();
}

void ConvertCArgs(int argc, char **argv, vector<string> &out_args);

vector<string> ExtractOptions(const vector<string> &args, int ignore_prefix=1, int ignore_suffix=0);
vector<string> ExtractArgs(const vector<string> &args, int args_count);
map<string,string> ParseGetoptArgs(int argc, char ** argv);

int64 NowSec();
int64 NowMicroSec();
int64 NowMilliSec();

class SwakTimer
{
public:
  SwakTimer() { reset(); }
  void reset() { start_time = -1; total_time = 0; }
  void start() { start_time = NowMicroSec(); }
  void stop() { if (start_time > 0) { total_time += NowMicroSec() - start_time; } }
  double sec_time() { return double(total_time) / 1e6; }
  int64 time() { return total_time; }

private:
  int64 start_time;
  int64 total_time;
};

inline double LogAdd(const double num1, const double num2)
{
  if (num1 > num2)
    return log1p(exp(num2 - num1)) + num1;
  else if (num1 < num2)
    return log1p(exp(num1 - num2)) + num2;

  return log(2) + num1;
}

inline ostream& operator<< (ostream &os, const YAML::Node &node)
{
  YAML::Emitter emitter;
  emitter << YAML::BeginMap;
  for(YAML::Iterator it=node.begin();it!=node.end();++it) {
     emitter << YAML::Key << it.first();
     emitter << YAML::Value;
     /*
     if(it.first() == "foo")
        emitter << "bar";
     else
     */
     emitter << it.second();
  }

  emitter << YAML::EndMap;
  os << emitter.c_str();
  return os;
}

namespace YAML
{
  void Parse(const string &str, Node &node);

  void ParseFile(const string &filename, Node &node);

  void ParseFile(istream &is, Node &node);

  bool HasKey(const Node &node, const string &key);

  template <typename T>
  void SetVar(const Node &node, const string &key, T &var, bool required=true)
  {
    if (HasKey(node, key))
    {
      try
      {
        var = node[key].to<T>();
      }
      catch (InvalidScalar e)
      {
       throw runtime_error(string("Unexpected value type for key ") + key);
      }
    }
    else if (required)
     throw runtime_error(string("Expected to find key ") + key);
  }


  void CmdLineOptionsToYaml(const vector<string> &options_vec, Node &options);
};

bool ProcessInput(const vector<string> &all_args, int num_prog_names, int num_req_args, YAML::Node &out_options, vector<string> &out_args);

// Use this version if your options contain spaces like "-t INT"
void PrintUsage(const string &usage_line, const string &desc, const VecS &lefts, const VecS &rights, int width=100, string bot_section = "Options");

// Use this version if your options are like "t=INT"
void PrintUsage(const string &usage_line, const string &desc, const VecS &lines, int width=100, string bot_section = "Options");

int ExtractDigits(int64 num, int base, vector<int> &digits, bool largest_first=false);

double Round(double d, int digits=0);

// Note that we include StlUtil which includes ToStr for vectors as well as opterator<< on vectors

#endif
