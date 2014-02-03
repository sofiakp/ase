#include "swak/System.h"
#include "swak/Helpers.h"
#include "swak/BinaryIO.h"

using namespace Swak;

int main_binary(const vector<string> &)
{
  cerr << "* Testing BinaryIO..." << endl;
  {
    cerr << "* Testing writing to out.bin..." << endl;
    int i = 30;
    long l = 40;
    string s = "Hello";

    BinWriter writer;
    writer.Init("out.bin", "TEST");
    writer.Write(i);
    writer.Write(l);
    writer.Write(s);
  }

  {
    cerr << "* Testing reading from out.bin..." << endl;
    int i = 0;
    long l = 0;
    string s = "";

    BinReader reader;
    reader.Init("out.bin", "TEST");
    reader.Read(i);
    reader.Read(l);
    reader.Read(s);

    pperr(i);
    pperr(l);
    pperr(s);
  }
  
  cerr << "* Done" << endl;
  return 0;
}
