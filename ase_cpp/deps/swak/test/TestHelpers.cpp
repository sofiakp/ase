#include "swak/Helpers.h"
#include "swak/ConfigReader.h"

int main_helpers(const vector<string> &args)
{
  cerr << "* Testing Helpers..." << endl;

  string a_str1 = "default_str1";
  string b= "default_str2";
  string str3 = "default_str3";
  string config_file = "";
  bool c_flag1 = false;
  bool d= false;
  bool flag3 = false;
  int w_int1 = -10;
  int x = -20;
  int int3 = -30;

  double y_num1 = -10.1;
  double z = -20.2;
  double num3 = -30.3;

  Swak::OptionParser op("<some_arg>", "Tests helpers");

  op.AddOpt(a_str1, 'a', "str1", "STR", "desc for a/str1");
  op.AddOpt(b, 'b', "",     "STR", "desc for b xxxxxxxxxxxxxxx xxxxx xxxx xxxxxxxxxxxxxxxxxxx xxxxxxxxxxxxx xxxxxx x x xxxxxxxxxxxxx xxxxxx");
  op.AddOpt(str3,  0,  "str3", "STR", "desc for str3");

  op.AddOpt(c_flag1, 'c', "flag1", "FLAG", "desc for c/flag1");
  op.AddOpt(d, 'd', "",      "FLAG", "desc for d");
  op.AddOpt(flag3,  0,  "flag3", "FLAG", "desc for flag3");

  op.AddOpt(w_int1, 'w', "int1", "INT", "desc for int1");
  op.AddOpt(x, 'x', "",     "INT", "desc for x");
  op.AddOpt(int3,  0,  "int3", "INT", "desc for int3");

  op.AddOpt(y_num1, 'y', "num1", "NUM", "desc for y/num1");
  op.AddOpt(z, 'z', "",     "NUM", "desc for z");
  op.AddOpt(num3,  0,  "num3_long", "NUM", "desc for num3");

  op.AddOpt(config_file,  0,  "config", "FILE", "Config file containing options. (Command line opts override)");

  /*
  args.clear();
  args.push_back("--long1");
  args.push_back("--long1");
  args.push_back("--flag1");
  args.push_back("arg!");

  //args.push_back("-f");
  //args.push_back("-s");
  //args.push_back("short_val");
  */

  VecS args_only;
  if (op.GetValue(args, "config", config_file, 2))
  {
    cerr << "READING CONFIG FILE '" << config_file << "'" << endl;

    Swak::ConfigReader conf_reader;
    conf_reader.Init(config_file);
    if (!conf_reader.ReadConfig(op, args, 2))
    {
      cout << "Error when processing config file." << endl;
      op.PrintUsage();
      exit(1);
    }
    cerr << "------ Config file parsed, opts now: -------" << endl;

    pperr(a_str1);
    pperr(b);
    pperr(str3);
    pperr(c_flag1);
    pperr(d);
    pperr(flag3);
    pperr(w_int1);
    pperr(x);
    pperr(int3);
    pperr(y_num1);
    pperr(z);
    pperr(num3);
    pperr(config_file);
  }

  cerr << "-------------" << endl;
  if (!op.Parse(args, args_only, 2) || args_only.size() < 1)
  {
    op.PrintUsage(true);
    exit(1);
  }
  cerr << "-------------" << endl;

  pperr(a_str1);
  pperr(b);
  pperr(str3);
  pperr(c_flag1);
  pperr(d);
  pperr(flag3);
  pperr(w_int1);
  pperr(x);
  pperr(int3);
  pperr(y_num1);
  pperr(z);
  pperr(num3);
  pperr(config_file);

  cerr << "-------------" << endl;
  if (args_only.size() > 0)
  {
    cerr << "Non-option args:" << endl;
    cerr << ToStr(args_only) << endl;
  }


  return 0;
}
