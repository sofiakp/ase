#include "swak/Swak.h"
#include "swak/PolyExe.h"
#include "Gff.h"
#include "GffReader.h"
#include "swak/StringUtil.h"
#include "swak/FoldedOutStream.h"
#include "swak/StringUtil.h"
#include "swak/Helpers.h"
#include "GeneReader.h"

void Usage(const string &exe)
{
  cout << "usage: " << exe << "<test>" << endl;
  cout << "===============" << endl;
  cout << "Available tests" << endl;
  cout << "===============" << endl;
  for (int i = 0; i < PolyExe::NumMains(); ++i)
    cout << PolyExe::MainNames(i) << endl;
}

int main_gff(const vector<string> &args)
{
  if (args.size() < 3 || args[2] == "-h")
  {
    cout << "usage: " << args[0] << " " << args[1] << " <file.gff>" << endl;
    cout << "  Reads and prints out the file" << endl;
    return 0;
  }

  GffReader reader(args[2]);
  GffEntry g;
  while (reader.GetNext(g))
  {
    cout << g << endl;

    cout << "Attribs:" << endl;
    MapSS m = g.GetAttributes();
    for (MapSS::const_iterator iter = m.begin(); iter != m.end(); ++iter)
    {
      cout << "\t" << iter->first << ": " << iter->second << endl;
    }
    cout << endl;
  }

  return 0;
}

int main_gene_reader(const vector<string> &args)
{
  if (args.size() < 3 || args[2] == "-h")
  {
    cout << "usage: " << args[0] << " " << args[1] << " <file.gtf>" << endl;
    cout << "  Reads and prints out the file" << endl;
    return 0;
  }

  GeneReader reader(args[2]);

  map<string, vector<Gene> > genes_by_chrom;
  bool print_warnings = true;
  reader.GetGenesByChrom(genes_by_chrom, print_warnings);
  Assert(GeneReader::GenesAreValid(genes_by_chrom));
  cout << GeneReader::GetNumWarnings(genes_by_chrom, print_warnings) << " warnings" << endl;

  for (map<string, vector<Gene> >::iterator iter = genes_by_chrom.begin();
      iter != genes_by_chrom.end(); ++iter)
  {
    cout << "==============" << endl;
    cout << "chrom: " << iter->first << endl;
    cout << "==============" << endl;
    StlFor(g, iter->second)
      cout << iter->second[g];
  }

  return 0;
}

int main_splices(const vector<string> &args)
{
  if (args.size() < 3 || args[2] == "-h")
  {
    cout << "usage: " << args[0] << " " << args[1] << " <file.gtf>" << endl;
    cout << "  Prints out splice junctions" << endl;
    return 0;
  }

  GeneReader reader(args[2]);

  map<string, vector<Gene> > genes_by_chrom;
  bool print_warnings = true;
  reader.GetGenesByChrom(genes_by_chrom, print_warnings);
  Assert(GeneReader::GenesAreValid(genes_by_chrom));
  cout << GeneReader::GetNumWarnings(genes_by_chrom, print_warnings) << " warnings" << endl;

  StlForMap(string, vector<Gene>, iter, genes_by_chrom)
  {
    StlFor(g, iter->second)
    {
      const Gene &gene = iter->second[g];
      cout << Color::green << gene << Color::reset;
      vector<Splice> splices = gene.GetSplices(10);
      StlFor(s, splices)
        cout << Color::cyan << splices[s].ToGff() << Color::reset << endl;
    }
  }

  return 0;
}

int main(int argc, char ** argv)
{
  DeclarePolyExe(main, gff, "Tests GFF class");
  DeclarePolyExe(main, gene_reader, "Tests reading GTF genes");
  DeclarePolyExe(main, splices, "Tests splice class");

  if (argc <= 1 || string(argv[1]) == string("-h") || !PolyExe::HaveMain(argv[1]))
  {
    Usage(argv[0]);
    return 0;
  }

  return PolyExe::CallMain(argc, argv);
}
