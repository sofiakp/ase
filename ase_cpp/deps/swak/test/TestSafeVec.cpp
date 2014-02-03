#include "swak/Swak.h"
#include "swak/SafeVec.h"

int main_safe_vec(const vector<string> &)
{
  cerr << "* Testing SafeVec..." << endl;

  SafeVec<int> v;

  v.push_back(1);

  cout << v[0] << endl;

  return 0;
}
