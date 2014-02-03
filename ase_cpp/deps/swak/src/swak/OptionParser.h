#include "Swak.h"
#include "System.h"

#ifndef SWAK_OPTION_PARSER_H
#define SWAK_OPTION_PARSER_H

namespace Swak
{
  class OptionParser
  {
  private:
    vector<char> shorts;
    vector<string> longs;
    vector<string> types;
    vector<string> prefixes; // like -x, --exclude || -x INT, --exclude INT
    vector<string> opt_descs;
    vector<void *> var_ptrs;

    string argv0;
    string sub_prog_name;
    string args_str;
    string prog_desc;

    bool help_flag_used;

    int FindShort(char c);
    int FindLong(const string &s);
    int Find(const string &s);

    void CheckOpt(string &, const string &type_name);
    void CheckOpt(int &, const string &type_name);
    void CheckOpt(int64 &, const string &type_name);
    void CheckOpt(float &, const string &type_name);
    void CheckOpt(double &, const string &type_name);
    void CheckOpt(bool &, const string &type_name);

 

  public:

    string error_message;

    OptionParser(const string &args_str_, const string &prog_desc_, bool auto_help_opt=true);

    bool ParseAsStrings(const VecS &args, VecS &new_args, int start_index, map<string, string> &key_val_pairs);
    bool Parse(const VecS &args, VecS &new_args, int start_index=1);
    bool GetValue(const vector<string> &args, const string &opt_name, string &val_str, int start_index=1);

    void PrintUsage(bool compact=true, int width=100, string bot_title = "Options");

    string GetType(string opt)
    {
      int oid = Find(NormalizeOptName(opt));

      AssertMsg(oid < int(types.size()), "Something went wrong finding long opt name");

      if (oid >= 0)
        return types[oid];
      else
        return "";
    }

    string NormalizeOptName(string opt)
    {
      // Trim any leading dashes
      if (opt.size() >= 1 && opt[0] == '-')
        opt = opt.substr(1);
      if (opt.size() >= 1 && opt[0] == '-')
        opt = opt.substr(1);

      int oid = -1;
      if (opt.size() == 1)
        oid = FindShort(opt[0]);
      else
        oid = FindLong(opt);

      AssertMsg(oid < int(types.size()), "Something went wrong finding opts in NormalizeOptName");
      if (oid >= 0)
      {
        if (opt.size() == 1)
          return string(1, shorts[oid]);
        else
          return longs[oid];
      }
      else
        return "";
    }   

    bool Contains(const string &s)
    {
      return Find(s) >= 0;
    }

    template <typename T>
    void AddOpt(T &var, char short_name, const string &long_name, const string &type_name, const string &opt_desc);
  };
};

#endif
