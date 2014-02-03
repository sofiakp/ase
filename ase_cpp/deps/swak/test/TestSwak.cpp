#include "swak/Swak.h"
#include "swak/PolyExe.h"
#include "swak/FoldedOutStream.h"
#include "swak/System.h"
#include "swak/StringUtil.h"

void Usage(const string &exe)
{

  cout << "usage: " << exe << "<test>" << endl;
  cout << "===============" << endl;
  cout << "Available tests" << endl;
  cout << "===============" << endl;
  for (int i = 0; i < PolyExe::NumMains(); ++i)
  {
    cout << ToStrL(PolyExe::MainNames(i), 20) << " " << PolyExe::MainDescs(i) << endl;
  }
}

int main_fold(const vector<string> &args)
{ 
  if (args.size() < 3 || args[2] == "-h")
  { 
    cout << "usage: " << Basename(args[0]) << " " << args[1] << " <line_width>" << endl;
    cout << "  Reads from stdin writes folded output to stdout" << endl;
    return 0;
  }

  int line_width = StrTo<int>(args[2]);
  FoldedOutStream out(cout, line_width);

  string line;
  while (getline(cin, line))
  { 
    out.Write(line + '\n');
  }

  return 0;
}


int main(int argc, char ** argv)
{
  DeclarePolyExe(main, helpers, "Tests Helpers.h");
  DeclarePolyExe(main, system, "Tests System.h");
  DeclarePolyExe(main, fold, "Tests FoldedOutStream.h");
  DeclarePolyExe(main, safe_vec, "Tests SafeVec.h");
  DeclarePolyExe(main, extract_digits, "Tests ExtractDigits in Helpers.h");
  DeclarePolyExe(main, binary, "Tests Binary IO");

  if (argc <= 1 || string(argv[1]) == string("-h") || !PolyExe::HaveMain(argv[1]))
  {
    Usage(argv[0]);
    return 0;
  }

  return PolyExe::CallMain(argc, argv);
}
