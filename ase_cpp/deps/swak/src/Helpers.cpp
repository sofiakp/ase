#include "swak/Helpers.h"
#include "swak/StlUtil.h"

void ConvertCArgs(int argc, char **argv, vector<string> &out_args)
{
  out_args.clear();
  out_args.assign(argv, argv + argc);
}

vector<string> ExtractOptions(const vector<string> &args, int ignore_prefix, int ignore_suffix)
{
  vector<string> options;
  options.assign(args.begin() + ignore_prefix, args.end() - ignore_suffix);
  return options;
}

vector<string> ExtractArgs(const vector<string> &all_args, int args_count)
{
  AssertMsg(all_args.size() >= args_count, string("Expected at least ") + ToStr(args_count) + " args, instead had " + ToStr(all_args.size()));
  vector<string> args(args_count);
  args.assign(all_args.begin() + (all_args.size() - args_count), all_args.end());
  return args;
}

int64 NowSeconds()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return int64(tv.tv_sec);
}

int64 NowMicroSeconds()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return int64(1000000) * int64(tv.tv_sec) + int64(tv.tv_usec);
}

int64 CurrentTimeMillis()
{
  return NowMicroSeconds() / int64(1000);
}

namespace YAML
{
  void Parse(const string &str, Node &node)
  {
    stringstream ss(str);
    Parser parser(ss);
    AssertMsg(parser.GetNextDocument(node), string("Couldn't parse YAML: ") + str);
  }

  void ParseFile(const string &filename, Node &node)
  {
    ifstream is(filename.c_str());
    AssertMsg(is.good(), string("Couldn't open YAML file: ") + filename);
    Parser parser(is);
    AssertMsg(parser.GetNextDocument(node), string("Couldn't parse YAML file: ") + filename);
  }

  void ParseFile(istream &is, Node &node)
  {
    AssertMsg(is.good(), string("Couldn't open YAML stream"));
    Parser parser(is);
    AssertMsg(parser.GetNextDocument(node), string("Couldn't parse YAML stream"));
  }

  /* ovec: command line args assumed to contain all the options for the program
   * options must be in the form of key= value or key=value NOT key = value
   * both key and value may contain spaces but they have to be quoted at the 
   * command line for this to work
   */
  void CmdLineOptionsToYaml(const vector<string> &ovec, Node &options)
  {
    VecS key_val_pairs;

    StlFor(i, ovec)
    {
      if (!StringContains(ovec[i], '='))
        throw runtime_error("Unexpected token in options: " + ovec[i]);

      string key;
      string value = "";

      if (StringEndsWith(ovec[i], '='))
      {
        key = ovec[i].substr(0, ovec[i].size() - 1);

        if (i < ovec.size() - 1 && !StringContains(ovec[i+1], '='))
        {
          value = ovec[i+1];
          ++i;
        }
      }
      else // the = must be in middle like foo=bar
      {
        size_t eq_pos = ovec[i].find('=');
        key = ovec[i].substr(0, eq_pos);
        value = ovec[i].substr(eq_pos + 1);
      }

      key_val_pairs.push_back(key + ": " + value);
    }

    string yaml_str = "{" + JoinFields(key_val_pairs, ", ") + "}";

    YAML::Parse(yaml_str, options);
  }

  bool HasKey(const Node &node, const string &key)
  {
    return (node.Type() == NodeType::Map && node.FindValue(key));
  }

};

/* Returns whether we had enough args to process input (or if -h was requested)
 * all_args: argv
 * out_args: output of required args
 * out_options: output of options
 * num_prog_names: how many names does this program have, usually 1 or 2 (ie samtools index has two names)
 * num_req_args: number of required args
 */
bool ProcessInput(const vector<string> &all_args, int num_prog_names, int num_req_args, YAML::Node &out_options, vector<string> &out_args)
{
  if (all_args.size() < num_prog_names + num_req_args || (all_args.size() > num_prog_names && all_args[num_prog_names] == "-h"))
    return false;

  out_args = ExtractArgs(all_args, num_req_args);
  vector<string> opts_vec = ExtractOptions(all_args, num_prog_names, num_req_args);
  YAML::CmdLineOptionsToYaml(opts_vec, out_options);

  return true;
}

void PrintUsage(const string &usage_line, const string &desc, const VecS &lefts, const VecS &rights, int width, string bot_section)
{
  cout << endl;
  cout << "Usage: " << usage_line << endl;
  if (desc.size() > 0)
  {
    cout << endl;
    cout << WrapParagraph(desc, width) << endl;
  }
  cout << endl;

  if (lefts.size() == 0)
    return;

  cout << bot_section << ":" << endl;

  size_t max_left_size = 0;
  StlFor(i, lefts)
    max_left_size = max(lefts[i].size(), max_left_size);

  int lead_size = 2;

  StlFor(i, lefts)
  {
    cout << string(lead_size, ' ') << ToStrL(lefts[i], max_left_size);
    if (rights[i].size() > 0)
    {
      int total_lead = 2 + lead_size + max_left_size;
      VecS lines = SplitString(WrapParagraph(rights[i], width, width - total_lead), "\n");
      
      cout << "  " << lines[0] << endl; // First line of text after the "left" item

      // Print any additional text that may have wrapped
      for(int l = 1; l < lines.size(); ++l)
      {
        if (lines[l].size() > 0)
          cout << string(total_lead, ' ') << lines[l] << endl;
      }

      cout << endl;
    }
    else
      cout << endl;
  }
}

void PrintUsage(const string &usage_line, const string &desc, const VecS &lines, int width, string bot_section)
{
  // Split lines into lefts and rights
  VecS lefts, rights;
  StlFor(i, lines)
  {
    if (StringContains(lines[i], ' '))
    {
      lefts.push_back(PrefixBefore(lines[i], ' '));
      rights.push_back(SuffixAfter(lines[i], ' '));
    }
    else
    {
      lefts.push_back(lines[i]);
      rights.push_back("");
    }
  }

  PrintUsage(usage_line, desc, lefts, rights, width, bot_section);
}

int ExtractDigits(int64 num, int base, vector<int> &digits, bool largest_first)
{
  digits.clear();

  do
  {
    int digit = num % base;
    num /= base;

    digits.push_back(digit);

  } while( num > 0);
  
  if (largest_first)
    reverse(digits.begin(), digits.end());

  return digits.size();
}

double Round(double d, int digits)
{
  long factor = pow(10, digits); 
  d *= factor;
  d = floor(d + 0.5);
  d /= factor;

  return d;
}
