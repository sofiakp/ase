#include "swak/Swak.h"
#include "swak/Helpers.h"

int main_extract_digits(const vector<string> &args)
{
  cerr << "* Testing ExtractDigits..." << endl;

  Swak::OptionParser op("<base> <num>", "Tests ExtractDigits");

  vector<string> args_only;

  if (!op.Parse(args, args_only, 2) || args_only.size() < 2)
  {
    op.PrintUsage();
    exit(1);
  }

  int base = StrTo<int>(args_only[0]);
  int num = StrTo<int>(args_only[1]);

  vector<int> digits;

  ExtractDigits(num, base, digits, true);

  cout << digits << endl;

  return 0;
}
