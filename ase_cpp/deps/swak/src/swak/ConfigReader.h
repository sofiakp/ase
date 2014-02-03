#ifndef SWAK_CONFIG_READER_H
#define SWAK_CONFIG_READER_H

#include "Swak.h"
#include "FieldReader.h"
#include "OptionParser.h"
#include "StringUtil.h"
#include "StlUtil.h"

namespace Swak
{
  class ConfigReader
  {
    private:

      FieldReader field_reader;

      bool GetNext(string &key, string &value)
      {
        vector<string> fields;

        bool got_line = field_reader.GetNext(fields);

        if (!got_line)
          return false;

        AssertMsg(fields.size() == 2, "Config line missing a colon");

        key = fields[0];

        // Merge all stuff after the first colon to be the "value"
        vector<string> new_fields;
        for (int i = 1; i < fields.size(); ++i)
          new_fields.push_back(fields[i]);

        value = JoinFields(new_fields, "");

        return true;
      }

    public:

      ConfigReader();

      void Init(const string & fn);
      
      /* Must have called Init() prior to each call to ReadConfig
       * if cmd_line_args are sent, we assume that we want to exclude the options set on the command line
       * from consideration when we are parsing the file. cmd_line_start_index is where we will start
       * parsing options in the cmd_line_args vector.
       */
      bool ReadConfig(OptionParser &opt_parser, const VecS &cmd_line_args=VecS(), int cmd_line_start_index=1)
      {
        // If there are some command line args, then extract which options are set
        MapSS cmd_line_key_val_pairs;
        if (!cmd_line_args.empty())
        {
          VecS dummy_new_args;
          if (!opt_parser.ParseAsStrings(cmd_line_args, dummy_new_args, cmd_line_start_index, cmd_line_key_val_pairs))
            return false;
        }

        vector<string> constructed_args; // We will use this vec to set the variables through OptionParser
        string key, value;

        // Read the key value pairs
        while (GetNext(key, value))
        {
          // Construct an arg vector to send to OptionParser
          key = StripString(key);
          value = StripString(value);

          AssertMsg(key.size() >= 1, "Missing key in config file for value " + value);

          if (!opt_parser.Contains(key))
          {
            opt_parser.error_message = "Unrecognized option '" + key + "' in config file";
            return false;
          }

          // Skip any options that also show up in the command line
          if (Map::Contains(cmd_line_key_val_pairs, opt_parser.NormalizeOptName(key)))
            continue;

          if (key.size() == 1)
            key = "-" + key;
          else
            key = "--" + key;

          // If this is a flag type, then set the presence or absence of the flag  with "on/off/true/false"
          bool is_flag = opt_parser.GetType(key) == "FLAG";
          bool flag_is_on = false;

          if (is_flag)
          {
            if (value == "true" || value == "on" || value == "TRUE" || value == "ON")
              flag_is_on = true;
            else if (value == "false" || value == "off" || value == "FALSE" || value == "OFF")
              flag_is_on = false;
            else
              throw runtime_error("Unexpected value '" + value + "' for flag '" + key + "' in config file");
          }

          if (!is_flag)
          {
            constructed_args.push_back(key);
            constructed_args.push_back(value);
          }
          else if (flag_is_on)
            constructed_args.push_back(key);
        }

        // Parse as if you sent at command line
        vector<string> dummy_out_args;
        pperr(constructed_args);
        AssertMsg(opt_parser.Parse(constructed_args, dummy_out_args, 0), "There was a problem with ConfigReader: " + opt_parser.error_message);

        return true;
      }
  };

};

#endif
