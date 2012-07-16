#include "GeneReader.h"

GeneReader::GeneReader()
{
}

GeneReader::GeneReader(const string &filename)
{
  Init(filename);
}

void GeneReader::Init(const string &filename)
{
  gff_reader.Init(filename);
}

// chr_gene_map: chrom name to list of genes on that chrom
// Returns number of warnings encountered
int GeneReader::GetGenesByChrom(map<string, vector<Gene> > &chr_gene_map, bool print_warnings)
{
  int num_warning_messages = 0;
  map<string, Gene> gene_map; // Name of gene to gene object
  map<string, Transcript * > tx_map; // Name of tx to tx object

  //=============================
  // Read genes into memory
  //=============================
  while (gff_reader.GetNext(gff))
  {
    if (gff.feature != "CDS" && gff.feature != "exon")
      continue;

    MapSS attribs = gff.GetAttributes();
    string &gene_id = attribs["gene_id"];
    string &tx_id = attribs["transcript_id"];
    string &gene_name = attribs["gene_name"];

    if (gene_id.empty())
    {
      cerr << gff << endl;
      throw runtime_error("Gff error: Empty gene_id");
    }
    if (tx_id.empty())
      throw runtime_error("Gff error: Empty transcript_id");

    // New gene, set name and strand
    if (!Map::Contains(gene_map, gene_id))
    {
      gene_map[gene_id].chrom = gff.chrom;
      gene_map[gene_id].id = gene_id;
      gene_map[gene_id].name = gene_name;
      gene_map[gene_id].strand = gff.strand;
    }

    Gene &gene = gene_map[gene_id];

    if (gene.chrom != gff.chrom)
    {
      if (print_warnings)
        cerr << "Warning: gene " << gene_id << " exists on two chromosomes: " << gene.chrom << " and " << gff.chrom << ". Skipping entries on " << gff.chrom << endl;
      ++num_warning_messages;
      continue;
    }

    if (!Map::Contains(tx_map, tx_id))
    {
      // Add a transcript to the the gene's transcript vector and get a pointer to it
      gene.transcripts.resize(gene.transcripts.size() + 1);

      // Reassign all pointers in case transcripts.resize() caused vector to be moved
      StlFor(i, gene.transcripts)
        tx_map[gene.transcripts[i].id] = &gene.transcripts[i];

      Transcript * tx_ptr = &(gene.transcripts[gene.transcripts.size() - 1]);
      tx_map[tx_id] = tx_ptr;

      // Setup this new transcript info
      tx_ptr->id = tx_id;
      tx_ptr->chrom = gff.chrom;
      tx_ptr->parent_gene_id = gene.id;
      tx_ptr->strand = gff.strand;
    }

    Transcript * tx_ptr = tx_map[tx_id];
    Exon exon(gff.start, gff.stop + 1);
    if (gff.feature == "exon")
      tx_ptr->exons.push_back(exon);
    else if (gff.feature == "CDS")
      tx_ptr->cds_regions.push_back(exon);
  }

  //=========================================================================
  // Prepare output
  //=========================================================================
  // 1. Sort exons/CDS inside each transcript
  StlForMap(string, Transcript *, iter, tx_map)
  {
    Transcript * tx_ptr = iter->second;

    sort(tx_ptr->exons.begin(), tx_ptr->exons.end());
    sort(tx_ptr->cds_regions.begin(), tx_ptr->cds_regions.end());
  }

  // 2. Sort transcripts within genes, add genes to vector per chrom
  chr_gene_map.clear();

  StlForMap(string, Gene, iter, gene_map)
  {
    Gene &gene = iter->second;
    sort(gene.transcripts.begin(), gene.transcripts.end());
    chr_gene_map[gene.chrom].push_back(gene);
  }

  // Sort genes within the chromosome
  for (map<string, vector<Gene> >::iterator iter = chr_gene_map.begin();
      iter != chr_gene_map.end(); ++iter)
  {
    vector<Gene> &gene_vec = iter->second;
    sort(gene_vec.begin(), gene_vec.end());
  }

  return num_warning_messages;
}

int GeneReader::GetGenes(vector<Gene> &out_genes, bool print_warnings)
{
  out_genes.clear();
  map<string, vector<Gene> > genes_by_chrom;
  int num_warnings = GetGenesByChrom(genes_by_chrom, print_warnings);

  StlForMap(string, vector<Gene>, iter, genes_by_chrom)
  {
    StlFor(g, iter->second)
      out_genes.push_back(iter->second[g]);
  }
  return num_warnings;
}

int GeneReader::GetGenesById(map<string, Gene> &out_genes, bool print_warnings)
{
  out_genes.clear();

  vector<Gene> gene_vec;
  int num_warnings = GetGenes(gene_vec, print_warnings);

  StlFor(i, gene_vec)
    out_genes[gene_vec[i].id] = gene_vec[i];

  return num_warnings;
}

int GeneReader::GetTranscripts(vector<Transcript> &out_tx, bool print_warnings)
{
  out_tx.clear();

  vector<Gene> gene_vec;
  int num_warnings = GetGenes(gene_vec, print_warnings);

  StlFor(g, gene_vec)
  {
    const Gene &gene = gene_vec[g];
    StlFor(t, gene.transcripts)
    {
      const Transcript &tx = gene.transcripts[t];
      out_tx.push_back(tx);
    }
  }

  return num_warnings;
}

int GeneReader::GetTranscriptsById(map<string, Transcript> &out_tx, bool print_warnings)
{
  out_tx.clear();

  vector<Transcript> tx_vec;
  int num_warnings = GetTranscripts(tx_vec, print_warnings);

  StlFor(t, tx_vec)
    out_tx[tx_vec[t].id] = tx_vec[t];

  return num_warnings;
}


bool GeneReader::GenesAreValid(const map<string, vector<Gene> > &genes_by_chrom)
{
  map<string, int> counts_of_tx_names;

  StlForMapConst(string, vector<Gene>, iter, genes_by_chrom)
  {
    const string &chrom_name = iter->first;
    const vector<Gene> &genes = iter->second;

    StlFor(g, genes)
    {
      const Gene &gene = genes[g];

      Assert(gene.chrom == chrom_name);

      StlFor(t, gene.transcripts)
      {
        const Transcript &tx = gene.transcripts[t];

        counts_of_tx_names[tx.id] += 1;

        Assert(tx.parent_gene_id == gene.id);
        Assert(tx.exons.size() > 0); // Note we must have exons but we may not have CDS if it's a non-coding RNA gene

        // Make sure all exons/CDS regions are in order and not overlapping (< operator ensures non-overlapping)
        for (int i = 1; i < tx.exons.size(); ++i)
          Assert(tx.exons[i - 1] < tx.exons[i]);
        for (int i = 1; i < tx.cds_regions.size(); ++i)
          Assert(tx.cds_regions[i - 1] < tx.cds_regions[i]);
      }
    }
  }

  StlForMap(string, int, iter, counts_of_tx_names)
    AssertMsg(iter->second == 1, "Transcript ids must be unique, but the transcript " + iter->first + " appears in multiple genes");

  return true;
}

int GeneReader::GetNumWarnings(const map<string, vector<Gene> > &genes_by_chrom, bool print_warnings)
{
  int num_warning_messages = 0;

  StlForMapConst(string, vector<Gene>, iter, genes_by_chrom)
  {
    const vector<Gene> &genes = iter->second;

    StlFor(g, genes)
    {
      const Gene &gene = genes[g];

      StlFor(t, gene.transcripts)
      {
        const Transcript &tx = gene.transcripts[t];

        if (tx.strand != gene.strand) 
        {
          ++num_warning_messages;

          if (print_warnings)
            cerr << "--- Warning: Transcript doesn't have same strand as parent gene.\nGene:\n" << gene << "Transcript:\n" << tx << endl;
        }
      }
    }
  }
  return num_warning_messages;
}
