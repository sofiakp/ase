#include "swak/System.h"
#include "swak/StringUtil.h"

void FindAndReplace(string &str, const string &string_to_find, const string &replace_with)
{
  size_t pos = 0;
  size_t find_len = string_to_find.length();
  size_t uReplaceLen = replace_with.length();

  if( find_len == 0 )
    return;

  for( ;(pos = str.find( string_to_find, pos )) != string::npos; )
  {
    str.replace( pos, find_len, replace_with );
    pos += uReplaceLen;
  }
}

bool StringEndsWith(const string &s, char c)
{
  return s[s.size() - 1] == c;
}

bool StringStartsWith(const string &s, char c)
{
  return s[0] == c;
}

string WrapParagraph(const string &in_str, int width, int first_line_width)
{
  if (first_line_width < 0)
    first_line_width = width;

  string str = in_str;
  size_t curWidth = first_line_width;

  while( curWidth < str.length() ) {
    std::string::size_type spacePos = str.rfind( ' ', curWidth );
    if( spacePos == std::string::npos )
      spacePos = str.find( ' ', curWidth );
    if( spacePos != std::string::npos ) {
      str[ spacePos ] = '\n';
      curWidth = spacePos + width + 1;
    }
  }

  return str;
}
