include_directories(. ../deps/swak/src ../deps/swak/deps/yaml-cpp/include ${bamtools_INCLUDE_DIR} ../deps/vcftools)

add_library(gene
    Gene.cpp
    Gff.cpp
    GeneReader.cpp
)
target_link_libraries(gene swak)

# Combine my bam code with the bamtools library
add_library(bam
            BamUtil.cpp)
# Use the first one for source and the second for static lib
target_link_libraries(bam BamTools)
#target_link_libraries(bam ${bamtools_LIBRARY})

# Combine my vcf code with the vcftools lib
add_library(vcf
            VcfUtil.cpp)
target_link_libraries(vcf vcftools z)
 
add_library(rod
            RodUtil.cpp)

add_executable(Ase 
  Ase.cpp
  Reconcile.cpp
  AseQuantMultiBAM.cpp
  AseQuantMultiRG.cpp
  AseRegion.cpp
)

add_executable(AseDisease
  AseDisease.cpp
  SeparatedbSNP.cpp
  LDDataPreProcess.cpp
  SplitLD.cpp
  ProcessGWASCatalog.cpp
  GenerateFairBinGWASCatalog.cpp
  ComputeOverlapping.cpp
  GenerateBinFilterNullSet.cpp
  ComputeEnrichment.cpp
)

target_link_libraries(Ase yaml-cpp swak gene bam vcf rod)
target_link_libraries(AseDisease yaml-cpp swak gene bam vcf rod boost_program_options)
