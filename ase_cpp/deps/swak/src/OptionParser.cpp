#include "swak/OptionParser.h"
#include "swak/StlUtil.h"

namespace Swak
{
  int OptionParser::FindShort(char c)
  {
    StlFor(i, shorts)
    {
      if (shorts[i] == c)
        return i;
    }
    return -1;
  }

  int OptionParser::FindLong(const string &s)
  {
    StlFor(i, longs)
    {
      if (longs[i] == s)
        return i;
    }
    return -1;
  }

  int OptionParser::Find(const string &s)
  {
    if (s.size() == 1)
      return FindShort(s[0]);
    else if (s.size() > 1)
      return FindLong(s);
    else
      return -1;
  }

  void OptionParser::CheckOpt(string &, const string &type_name) { AssertMsg(type_name != "INT" && type_name != "NUM" && type_name != "FLAG" && type_name != "LNG" && type_name != "FLT",
      "string-like type mismatch in OptionParser::AddOpt"); }
  void OptionParser::CheckOpt(int &, const string &type_name) { AssertMsg(type_name == "INT", "int type mismatch in OptionParser::AddOpt"); }
  void OptionParser::CheckOpt(int64 &, const string &type_name) { AssertMsg(type_name == "LNG", "long long type mismatch in OptionParser::AddOpt");}
  void OptionParser::CheckOpt(float &, const string &type_name) { AssertMsg(type_name == "FLT", "float type mismatch in OptionParser::AddOpt"); }
  void OptionParser::CheckOpt(double &, const string &type_name) { AssertMsg(type_name == "NUM", "double Type mismatch in OptionParser::AddOpt"); }
  void OptionParser::CheckOpt(bool &, const string &type_name) { AssertMsg(type_name == "FLAG", "bool (FLAG) Type mismatch in OptionParser::AddOpt"); }

  OptionParser::OptionParser(const string &args_str_, const string &prog_desc_, bool auto_help_opt)
  {
    args_str = args_str_;
    prog_desc = prog_desc_;

    help_flag_used = false; // The auto-generated -h var was not set (will be set to true if -h is used)

    // TODO: Handle default values somehow

    if (auto_help_opt)
      AddOpt(help_flag_used, 'h', "help", "FLAG", "Display this help message");
  }

  bool OptionParser::ParseAsStrings(const VecS &args, VecS &new_args, int start_index, map<string, string> &key_val_pairs)
  {
    key_val_pairs.clear();

    AssertMsg(start_index >= 0 && start_index <= args.size(), "Illegal start index or args passed to Swak::OptionParser::ParseAsStrings");

    argv0 = args[0];

    if (args.size() >= 2 && start_index >= 2)
      sub_prog_name = args[1];

    new_args.clear();

    int i;

    vector<bool> seen_flag_inds(longs.size(), false); // So we don't flip a flag twice!

    for (i = start_index; i < args.size(); ++i)
    {
      const string &arg = args[i];

      //cerr << "* Reading arg: " << arg << endl;;

      if (arg[0] != '-' || arg == "-" || arg == "--")
      {
        if (arg == "--")
          ++i;
        //cerr << "* Breaking bc arg=" << arg << endl;
        break;
      }

      if (arg[1] == '-') // Long!
      {
        //cerr << "* Processing as long" << endl;
        string opt_name = arg.substr(2);
        int oid = FindLong(opt_name);

        if (oid < 0)
        {
          error_message = "Unrecognized option " + arg;
          return false;
        }

        // If no flag, then there must be an argument next
        if (types[oid] != "FLAG")
        {
          if (i >= args.size() - 1 || args[i+1] == "--")
          {
            error_message = "Option " + arg + " requires an argument";
            return false;
          }
          else
          {
            string val_str = args[i+1];
            string norm_opt_name = NormalizeOptName(opt_name);
            AssertMsg(norm_opt_name != "", "Could not normalize opt_name: '" + opt_name + "'");
            key_val_pairs[norm_opt_name] = val_str;
          }

          // So we don't think the argument is the next option name
          ++i;
        }
        else // It's a flag, flip the value
        {
          string norm_opt_name = NormalizeOptName(opt_name);
          AssertMsg(norm_opt_name != "", "Could not normalize opt_name: '" + opt_name + "'");
          key_val_pairs[norm_opt_name] = "true";
        }
      }
      else // Short!
      {
        //cerr << "* Processing as short" << endl;

        string opt_name = string(1, arg[1]);
        int oid = FindShort(opt_name[0]);

        if (oid < 0)
        {
          error_message = "Unrecognized option " + arg;
          return false;
        }

        // If no flag, then there must be an argument next
        if (types[oid] != "FLAG")
        {
          string val_str;

          if (arg.size() > 2) // Like -a10
          {
            val_str = arg.substr(2);
          }
          else if (i >= args.size() - 1 || args[i+1] == "--")
          {
            error_message = "Option " + arg + " requires an argument";
            return false;
          }
          else // Like -a 10
          {
            val_str = args[i+1];

            // So we don't think the argument is the next option name
            ++i;
          }

          string norm_opt_name = NormalizeOptName(opt_name);
          AssertMsg(norm_opt_name != "", "Could not normalize opt_name: '" + opt_name + "'");
          key_val_pairs[norm_opt_name] = val_str;
        }
        else // It's a flag, flip the value
        {
          // Skip this flag, we have seen it before
          if (seen_flag_inds[oid])
            continue;

          string norm_opt_name = NormalizeOptName(opt_name);
          AssertMsg(norm_opt_name != "", "Could not normalize opt_name: '" + opt_name + "'");
          key_val_pairs[norm_opt_name] = "true";

          seen_flag_inds[oid] = true;
        }
      }
    }

    if (help_flag_used)
      return false;

    // Copy non-option args
    AssertMsg(i <= args.size(), "Something went wrong parsing options");
    for (int j = i; j < args.size(); ++j)
      new_args.push_back(args[j]);

    return true;
  }

  bool OptionParser::Parse(const VecS &args, VecS &new_args, int start_index) 
  {
    AssertMsg(start_index >= 0 && start_index <= args.size(), "Illegal start index or args passed to Swak::OptionParser::Parse");
    map<string, string> key_val_pairs;
    if (!ParseAsStrings(args, new_args, start_index, key_val_pairs))
      return false;

    StlForMap(string, string, iter, key_val_pairs)
    {
      int oid = Find(iter->first);
      AssertMsg(oid >= 0, "Something went wrong parsing options");

      string val_str = iter->second;

      if (types[oid] == "INT")
      {
        int new_val = StrTo<int>(val_str);
        *((int *)var_ptrs[oid]) = new_val;
        //cerr << "* Set INT short=" << shorts[oid] << " longs=" << longs[oid] << " to be " << new_val << endl;
      }
      else if (types[oid] == "NUM")
      {
        double new_val = StrTo<double>(val_str);
        *((double *)var_ptrs[oid]) = new_val;
        //cerr << "* Set NUM short=" << shorts[oid] << " longs=" << longs[oid] << " to be " << new_val << endl;
      }
      else if (types[oid] == "FLT")
      {
        float new_val = StrTo<float>(val_str);
        *((float *)var_ptrs[oid]) = new_val;
        //cerr << "* Set FLT short=" << shorts[oid] << " longs=" << longs[oid] << " to be " << new_val << endl;
      }
      else if (types[oid] == "FLAG")
      {
        bool new_val = !(*((bool *)var_ptrs[oid]));
        *((bool *)var_ptrs[oid]) = new_val;
        //cerr << "* Set short=" << shorts[oid] << " longs=" << longs[oid] << " to be " << new_val << endl;
      }
      else
      {
        string new_val = val_str;

        *((string *)var_ptrs[oid]) = new_val;
        //cerr << "* Set " << types[oid] << " short=" << shorts[oid] << " longs=" << longs[oid] << " to be " << new_val << endl;
      }

    }

    return true;
  }

  bool OptionParser::GetValue(const vector<string> &args, const string &opt_name, string &val_str, int start_index)
  {
    AssertMsg(start_index >= 0 && start_index <= args.size(), "Illegal start index or args passed to Swak::OptionParser::GetValue");
    map<string, string> key_val_pairs;
    vector<string> dummy_new_args;
    if (!ParseAsStrings(args, dummy_new_args, start_index, key_val_pairs))
      return false;

    string norm_opt_name = NormalizeOptName(opt_name);

    if (Map::Contains(key_val_pairs, norm_opt_name))
    {
      val_str = key_val_pairs[norm_opt_name];
      return true;
    }
    else
      return false;
  }

  void OptionParser::PrintUsage(bool compact, int width, string bot_title)
  {
    cout << endl;
    cout << "Usage: " << Basename(argv0);
    AssertMsg(shorts.size() == longs.size(), "Something went wrong printing usage");

    if (!sub_prog_name.empty()) 
      cout << " " << sub_prog_name;
    if (shorts.size() > 0 || opt_descs.size() > 0)
      cout << " [options] ";

    cout << args_str << endl;

    if (prog_desc.size() > 0)
    {
      cout << endl;
      cout << WrapParagraph(prog_desc, width) << endl;
    }
    cout << endl;

    if (shorts.size() == 0 && longs.size() == 0)
      return;

    cout << bot_title << ":" << endl;

    // Build option names (like "-v INT, --verbosity INT")
    VecS opt_names(shorts.size());

    StlFor(i, shorts)
    {
      stringstream ss;

      if (shorts[i] != 0 && longs[i] != "")
      {
        // cerr << "x " << shorts[i] << "_" << longs[i] << endl;
        ss << "-" << shorts[i] << ", --" << longs[i];
      }
      else if (shorts[i] != 0)
      {
        // cerr << "y " << shorts[i] << endl;
        ss << "-" << shorts[i];
      }
      else if (longs[i] != "")
      {
        // cerr << "z " << longs[i] << endl;
        ss << string(4, ' ') << "--" << longs[i];
      }
      else
        throw runtime_error("An option must have either a short or long name!");

      if (types[i] != "FLAG")
        ss << " " << types[i];

      opt_names[i] = ss.str();
    }


    size_t max_left_size = 0;
    StlFor(i, opt_names)
      max_left_size = max(opt_names[i].size(), max_left_size);

    int lead_size = 2;

    StlFor(i, opt_names)
    {
      cout << string(lead_size, ' ') << ToStrL(opt_names[i], max_left_size);
      if (opt_descs[i].size() > 0)
      {
        int total_lead = 2 + lead_size + max_left_size;
        VecS lines = SplitString(WrapParagraph(opt_descs[i], width, width - total_lead), "\n");

        cout << "  " << lines[0] << endl; // First line of text after the "left" item

        // Print any additional text that may have wrapped
        for(int l = 1; l < lines.size(); ++l)
        {
          if (lines[l].size() > 0)
            cout << string(total_lead, ' ') << lines[l] << endl;
        }

        if (!compact)
          cout << endl;
      }
      else
        cout << endl;
    }

    if (error_message != "")
    {
      cout << endl << "Options error: " << error_message << endl;
    }
  }

  template <typename T>
  void OptionParser::AddOpt(T &var, char short_name, const string &long_name, const string &type_name, const string &opt_desc)
  {
    AssertMsg(short_name > 0 || long_name != "", "OptionParser::AddOpt requires at least one short or long name");
    CheckOpt(var, type_name);
    var_ptrs.push_back(&var);
    shorts.push_back(short_name);
    longs.push_back(long_name);
    types.push_back(type_name);
    opt_descs.push_back(opt_desc);
  }

  // Explicit instantiation of all the different types we expect to use with this method.
  // This allows template code to be in the cpp file and not the .h file
  template void OptionParser::AddOpt<string>(string &, char, const string &, const string &, const string &);
  template void OptionParser::AddOpt<int>(int &, char, const string &, const string &, const string &);
  template void OptionParser::AddOpt<int64>(int64 &, char, const string &, const string &, const string &);
  template void OptionParser::AddOpt<float>(float &, char, const string &, const string &, const string &);
  template void OptionParser::AddOpt<double>(double &, char, const string &, const string &, const string &);
  template void OptionParser::AddOpt<bool>(bool &, char, const string &, const string &, const string &);
};
