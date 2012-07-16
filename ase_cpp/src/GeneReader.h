#ifndef GENE_READER_H
#define GENE_READER_H

#include "GffReader.h"
#include "Gene.h"
#include "swak/StlUtil.h"
#include "swak/Helpers.h"

class GeneReader
{
public:
  GeneReader();
  GeneReader(const string &filename);

  void Init(const string &filename);

  // chr_gene_map: chrom name to list of genes on that chrom
  int GetGenesByChrom(map<string, vector<Gene> > &chr_gene_map, bool print_warnings=false);
  int GetGenes(vector<Gene> &out_genes, bool print_warnings=false);
  int GetGenesById(map<string, Gene> &out_genes, bool print_warnings=false);
  int GetTranscripts(vector<Transcript> &out_tx, bool print_warnings=false);
  int GetTranscriptsById(map<string, Transcript> &out_tx, bool print_warnings=false);

  static bool GenesAreValid(const map<string, vector<Gene> > &genes_by_chrom);

  static int GetNumWarnings(const map<string, vector<Gene> > &genes_by_chrom, bool print_warnings);

private:
  GffReader gff_reader;
  GffEntry gff;
};

#endif
