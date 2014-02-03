#ifndef SWAK_STRINGUTIL_H
#define SWAK_STRINGUTIL_H


template <typename FloatType>
void FloatToString(FloatType val, string &s, int decimals);

template <typename FloatType>
string FloatToString(FloatType val, int decimals);

template <typename FloatType>
void FloatToString(FloatType val, string &s, int width, int decimals, 
		bool trailing_zeros = true, char pad = ' ');

template <typename FloatType>
string FloatToString(FloatType val, int width, int decimals,
		  bool trailing_zeros = true, char pad = ' ');

template <typename FloatType>
string FloatToSci(FloatType val, int decimals);

template <typename StringType>
StringType &ToUpper(const StringType &in, StringType &out);

template <typename StringType>
StringType ToUpper(const StringType &s);

template <typename StringType>
StringType &ToLower(const StringType &in, StringType &out);

template <typename StringType>
StringType ToLower(const StringType &s);

template <typename StringType>
StringType &StripString(const StringType &in, StringType &out);

template <typename StringType>
StringType StripString(const StringType &s);

template <typename StringType>
StringType &StripString(const StringType &in, StringType &out, char c);

template <typename StringType>
StringType StripString(const StringType &s, char c);

template <typename StringType, typename StringType2>
StringType &PrefixBefore(const StringType &in, const StringType2 &search, StringType &out);

template <typename StringType, typename StringType2>
StringType PrefixBefore(const StringType &s, const StringType2 &search);

template <typename StringType, typename StringType2>
StringType &SuffixAfter(const StringType &in, const StringType2 &search, StringType &out);

template <typename StringType, typename StringType2>
StringType SuffixAfter(const StringType &s, const StringType2 &search);

// Splits on delims, considering adjacent delims as one delim, ignores leading and trailing delims
template <typename StringType, typename StringType2>
int SplitString2(const StringType &s, const StringType2 &delim, vector<StringType> &fields);

template <typename StringType, typename StringType2>
vector<StringType> SplitString2(const StringType &s, const StringType2 &delim);

// Splits on delims, creating a field between every delim and the beginning/end of the string
template <typename StringType, typename StringType2>
int SplitString(const StringType &s, const StringType2 &delim, vector<StringType> &fields);

template <typename StringType, typename StringType2>
vector<StringType> SplitString(const StringType &s, const StringType2 &delim);

template <typename StringType, typename StringType2>
bool FastSplitString2(const StringType &s, const StringType2 &delim, StringType &field, uint32 index);

template <typename StringType, typename StringType2>
bool FastSplitString(const StringType &s, const StringType2 &delim, StringType &field, uint32 index);

// String conversion

template <typename NonStringType>
void StrTo(const string &s, NonStringType &val);

template <typename NonStringType>
NonStringType StrTo(const string &s);

template <typename NonStringType>
void ToStr(const NonStringType &val, string &s);

template <typename NonStringType>
string ToStr(const NonStringType &val);

template <typename StringType>
string ToStrL(const StringType &val, int width, char pad = ' ');

template <typename StringType>
string ToStrR(const StringType &val, int width, char pad = ' ');

//////////////////////////

template <typename NonStringType>
void StrTo(const string &s, NonStringType &val)
{
  stringstream ss(s);
  ss >> val;
}

template <typename NonStringType>
NonStringType StrTo(const string &s)
{
  NonStringType val;
  stringstream ss(s);
  ss >> val;
  return val;
}

template <typename NonStringType>
void ToStr(const NonStringType &val, string &s)
{
  stringstream ss;
  ss << val;
  s = ss.str();
}

template <typename NonStringType>
string ToStr(const NonStringType &val)
{
  stringstream ss;
  ss << val;
  return ss.str();
}

template <typename NonStringType>
string ToStrL(const NonStringType &val, int width, char pad)
{
  string s;
  ToStr(val, s);
  int s_len = min<int>(width, s.length());
  return s.substr(0, s_len) + string(width - s_len, pad);
}

template <typename NonStringType>
string ToStrR(const NonStringType &val, int width, char pad)
{
  string s;
  ToStr(val, s);
  int s_len = min<int>(width, s.length());
  return string(width - s_len, pad) + s.substr(s.length() - s_len, s_len);
}

template <typename FloatType>
void FloatToString(FloatType val, string &s, int decimals)
{
  stringstream ss;
  ss << fixed << setprecision(decimals) << val;
  s = ss.str();
}

template <typename FloatType>
string FloatToString(FloatType val, int decimals)
{
  string s;
  FloatToString(val, s, decimals);
  return s;
}

template <typename FloatType>
void FloatToString(FloatType val, string &s, int width, int decimals,
		bool trailing_zeros, char pad)
{
  stringstream ss;
  ss << fixed << setprecision(decimals) << val;
  s = ss.str();

  if ((decimals > 0) && !trailing_zeros)
  {
    while (s[s.length() - 1] == '0')
      s.resize(s.length() - 1);
    if (s[s.length() - 1] == '.')
      s += '0';
  }

  if (s.length() <= width)
    s = string(width - s.length(), pad) + s;
  else
    s = string(width, '#');
}

template <typename FloatType>
string FloatToString(FloatType val, int width, int decimals, 
		  bool trailing_zeros, char pad)
{
  string s;
  FloatToString(val, s, width, decimals, trailing_zeros, pad);
  return s;
}

template <typename FloatType>
string FloatToSci(FloatType val, int decimals)
{
  stringstream ss;
  ss << setprecision(decimals) <<  scientific << val;
  return ss.str();
}

template <typename StringType>
StringType &ToUpper(const StringType &in, StringType &out)
{
  out.resize(in.size());
  for (int i = 0; i < out.size(); ++i)
    out[i] = toupper(in[i]);
  return out;
}

template <typename StringType>
StringType ToUpper(const StringType &s)
{
  StringType result;
  ToUpper(s, result);
  return result;
}

template <typename StringType>
StringType &ToLower(const StringType &in, StringType &out)
{
  out.resize(in.size());
  for (int i = 0; i < out.size(); ++i)
    out[i] = tolower(in[i]);
  return out;
}

template <typename StringType>
StringType ToLower(const StringType &s)
{
  StringType result;
  ToLower(s, result);
  return result;
}

template <typename StringType>
StringType &StripString(const StringType &in, StringType &out)
{
  int begin = 0;
  while ((begin < in.size()) && isspace(in[begin]))
    ++begin;

  int end = in.size();
  while ((end > begin) && isspace(in[end - 1]))
    --end;

  out = in.substr(begin, end - begin);
  return out;
}

template <typename StringType>
StringType StripString(const StringType &s)
{
  StringType result;
  StripString(s, result);
  return result;
}

// Removes trailing \r
void TrimTrailingCR(string &s);

template <typename StringType>
StringType &StripString(const StringType &in, StringType &out, char c)
{
  if (in.size() == 0)
    return (out = "");

  int begin = (in[0] == c) ? 1 : 0;
  int end = (in[in.size() - 1] == c) ? in.size() - 2 : in.size() - 1;

  if (begin <= end)
    out = in.substr(begin, end - begin + 1);
  else
    out = "";

  return out;
}

template <typename StringType>
StringType StripString(const StringType &s, char c)
{
  StringType result;
  StripString(s, result, c);
  return result;
}

template <typename StringType, typename StringType2>
StringType &BeforeFirst(const StringType &in, const StringType2 &search, StringType &out)
{
  size_t pos = in.find(search);
  if (pos == StringType::npos)
    out = in;
  else
    out = in.substr(0, pos);
  return out;
}

template <typename StringType, typename StringType2>
StringType BeforeFirst(const StringType &s, const StringType2 &search)
{
  StringType result;
  BeforeFirst(s, search, result);
  return result;
}

template <typename StringType, typename StringType2>
StringType &PrefixBefore(const StringType &in, const StringType2 &search, StringType &out)
{
  size_t pos = in.find(search);
  if (pos == StringType::npos)
    out = in;
  else
    out = in.substr(0, pos);
  return out;
}

template <typename StringType, typename StringType2>
StringType PrefixBefore(const StringType &s, const StringType2 &search)
{
  StringType result;
  BeforeFirst(s, search, result);
  return result;
}

template <typename StringType, typename StringType2>
StringType &SuffixAfter(const StringType &in, const StringType2 &search, StringType &out)
{
  size_t pos = in.find(search);
  if (pos == StringType::npos)
    out.clear();
  else
    out = in.substr(pos + 1, in.length() - (pos + 1));
  return out;
}

template <typename StringType, typename StringType2>
StringType SuffixAfter(const StringType &s, const StringType2 &search)
{
  StringType result;
  SuffixAfter(s, search, result);
  return result;
}

template <typename StringType, typename StringType2>
int SplitString2(const StringType &s, const StringType2 &delim, vector<StringType> &fields)
{
  fields.clear();
  size_t begin = 0;
  while (begin < s.size())
  {
    begin = s.find_first_not_of(delim, begin);
    if (begin == StringType::npos)
      break;
    size_t end = s.find_first_of(delim, begin);
    if (end == StringType::npos)
      end = s.size();
    fields.push_back(s.substr(begin, end - begin));
    begin = end;
  }
  return fields.size();
}

template <typename StringType, typename StringType2>
vector<StringType> SplitString2(const StringType &s, const StringType2 &delim)
{
  vector<StringType> fields;
  SplitString2(s, delim, fields);
  return fields;
}

template <typename StringType, typename StringType2>
bool FastSplitString2(const StringType &s, const StringType2 &delim, StringType &field, uint32 index)
{
  uint32 i = 0; // index of current field
  size_t begin = 0;
  while (begin < s.size())
  {
    begin = s.find_first_not_of(delim, begin);
    if (begin == StringType::npos)
      break;
    size_t end = s.find_first_of(delim, begin);
    if (end == StringType::npos)
      end = s.size();
    if (i == index)
    {
      field = s.substr(begin, end - begin);
      return true;
    }
    begin = end;
    ++i;
  }
  return false;
}

template <typename StringType, typename StringType2>
int SplitString(const StringType &s, const StringType2 &delim, vector<StringType> &fields)
{
  fields.clear();
  if (s.empty())
    return 0;
  size_t begin = 0;
  while (begin <= s.size())
  {
    size_t end = s.find_first_of(delim, begin);
    if (end == StringType::npos)
      end = s.size();
    fields.push_back(s.substr(begin, end - begin));
    begin = end + 1;
  }
  return fields.size();
}

template <typename StringType, typename StringType2>
vector<StringType> SplitString(const StringType &s, const StringType2 &delim)
{
  vector<StringType> fields;
  SplitString(s, delim, fields);
  return fields;
}

template <typename StringType, typename StringType2>
bool FastSplitString(const StringType &s, const StringType2 &delim, StringType &field, uint32 index)
{
  if (s.empty())
    return false;
  uint32 i = 0; // index of current field
  size_t begin = 0;
  while (begin <= s.size())
  {
    size_t end = s.find_first_of(delim, begin);
    if (end == StringType::npos)
      end = s.size();
    if (i == index)
    {
      field = s.substr(begin, end - begin);
      return true;
    }
    begin = end + 1;
    ++i;
  }
  return false;
}

template <typename NumType1, typename NumType2>
string ToPercent(NumType1 val, NumType2 whole, int decimals=1)
{
  stringstream ss;
  double percent = double(val) / whole * 100;
  ss << fixed << setprecision(decimals) << percent << "%";
  return ss.str();
}

template <typename StringType, typename StringType2>
string JoinFields(const vector<StringType> &fields, const StringType2 &delim)
{
  stringstream ss;
  if (fields.size() == 0)
    return "";

  ss << fields[0];
  for (int i = 1; i < fields.size(); ++i)
  {
    ss << delim << fields[i];
  }

  return ss.str();
}

void FindAndReplace(string &str, const string &string_to_find, const string &replace_with);

template <typename StringType, typename StringType2>
bool StringContains(const StringType &str, const StringType2 &target)
{
  size_t pos = str.find(target);
  return (pos != StringType::npos);
}

bool StringEndsWith(const string &s, char c);
bool StringStartsWith(const string &s, char c);

string WrapParagraph(const string &s, int width, int first_line_width = -1);




#endif
