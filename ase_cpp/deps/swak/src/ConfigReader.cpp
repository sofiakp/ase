#include "swak/ConfigReader.h"

using namespace Swak;

ConfigReader::ConfigReader()
{
}

void ConfigReader::Init(const string & fn)
{
  field_reader.Init(fn);
  Assert(field_reader.IsOpen());
  field_reader.SetDelim(":");
}
