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
  AseQuantMultiRG.cpp
  AseQuantMultiBAM.cpp
  AseRegion.cpp
)

target_link_libraries(Ase yaml-cpp swak gene bam vcf rod)
