#include "swak/System.h"
#include "swak/Helpers.h"

int main_system(const vector<string> &)
{
  cerr << "* Testing System...";
  string dir;
  string file;
  SplitPath("./dir/file", dir, file);
  Assert(dir == "./dir");
  Assert(file == "file");
  cout << "OK" << endl;
  return 0;
}
